#!/usr/bin/env python

import flywheel_gear_toolkit
import flywheel
import logging
import pip
import smtplib
import sys
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime, date, timedelta
import time
import pandas as pd
import re
from flywheel import (
    ProjectOutput,
    AcquisitionListOutput,
    FileListOutput,
    Gear
)
from heudiconv.heuristics.reproin import parse_series_spec

pip.main(["install", "--upgrade", "git+https://github.com/poldracklab/wbhi-utils.git"])
from wbhiutils import parse_dicom_hdr
from wbhiutils.constants import (
    ADMIN_EMAIL
)

log = logging.getLogger(__name__)

WAIT_TIMEOUT = 3600 * 2

DATAVIEW_COLUMNS = [
    'subject.id',
    'subject.label',
    'session.id',
    'session.tags',
    'session.info.BIDS',
    'acquisition.label',
    'acquisition.id',
    'file.classification.Intent',
    'file.classification.Measurement',
    'file.classification.Features',
    'file.modality',
    'file.file_id',
    'file.tags',
    'file.type',
    'file.created',
    'file.name',
    'acquisition.timestamp'
]

def create_view_df(container, columns: list, filter=None) -> pd.DataFrame:
    """Get unique labels for all acquisitions in the container.

    This is done using a single Data View which is more efficient than iterating through
    all acquisitions, sessions, and subjects. This prevents time-out errors in large projects.
    """

    builder = flywheel.ViewBuilder(
        container='acquisition',
        filename="*.*",
        match='all',
        filter=filter,
        process_files=False,
        include_ids=False,
        include_labels=False
    )
    for c in columns:
        builder.column(src=c)
   
    view = builder.build()
    return client.read_view_dataframe(view, container.id)

def pydeface_filter(df: pd.DataFrame) -> bool:
    """Returns True if any acquisitions containing a Structural T1 or T2 dicom
    DO NOT also contain a nifti with a 'pydeface' label. Otherwise, returns False."""
    
    # Check for any Structural T1 or T2 dicoms
    dcm_df = df.copy()
    try:
        dcm_df = dcm_df[dcm_df['file.type'] == 'dicom']
        dcm_df = dcm_df[~dcm_df['file.classification.Intent'].isna()]
        dcm_df = dcm_df[~dcm_df['file.classification.Measurement'].isna()]
        mask = dcm_df['file.classification.Intent'].apply(lambda x: 'Structural' in x)
        dcm_df = dcm_df[mask]
        mask = dcm_df['file.classification.Measurement'].apply(lambda x: 'T1' in x or 'T2' in x)
        dcm_df = dcm_df[mask]
    except KeyError:
        return False
    if dcm_df.empty:
        return False

    # If structurals present, check for niftis without 'pydeface' tag
    try:
        nii_df = df.copy()
        nii_df = nii_df[nii_df['file.type'] == 'nifti']
        mask = nii_df['file.tags'].apply(lambda x: 'pydeface' not in x).any()
    except KeyError:
        return False
    return mask

def get_subjects(project: ProjectOutput) -> pd.DataFrame:
    """Returns a list of subject IDs of subjects that have not yet been bidsified
    and are ready to be."""

    file_df = create_view_df(project, DATAVIEW_COLUMNS) 
    file_df = file_df[file_df['file.type'].isin(('dicom', 'nifti'))]
    
    # Skip subject if all sessions have a 'bidsified' tag
    file_df = file_df.groupby('subject.id').filter(
            lambda x: x['session.tags'].apply(lambda x: 'bidsified' not in x).any()
    )
    nonbids_df = file_df.copy()
    
    # Skip subjects if not all structural T1 and T2 niftis have 'pydeface' tag
    nii_df = file_df.copy()
    nii_df = nii_df.groupby('acquisition.id').filter(pydeface_filter)
    file_df = file_df[~file_df['subject.id'].isin(nii_df['subject.id'])]

    # Add column of filename without the suffix to find matching dicoms and niftis
    file_df = file_df[~file_df['file.name'].isnull()]
    file_df['file.name'] = file_df['file.name'].str.replace('_', '.')
    file_df['file_root'] = file_df['file.name'].apply(
        lambda x: x.split('.dicom')[0].split('.dcm')[0]
    )
    
    # Raise warning if some non-bidsified sessions are filtered out
    filt_out_set = set(nonbids_df['subject.id']) - set(file_df['subject.id'])
    if filt_out_set:
        log.warning("The following subjects are not yet bidsified and filtered out: "
                    f"\n{filt_out_set}"
        )

    if file_df.empty:
        log.info("No subjects were found that need to be bidsified.")
        sys.exit(0)

    return file_df


def classify(file_df: pd.DataFrame) -> pd.DataFrame:
    """Determines which acquisitions should be bidsified and converts the acquisition labels
    to repronim standard. For acquisitions that won't be bidsified, adds '_ignore-BIDS' to
    the end of the acquisition label, signaling the bids-curate gear to skip."""

    def reproin_filter(row: pd.Series) -> str:
        measurement = row['file.classification.Measurement']
        intent = row['file.classification.Intent']
        label = row['acquisition.label']
        label = re.sub(r'_\d$', '', label)
        
        validator = parse_series_spec(label)
        reproin = None

        if 'bids-gear-ignore'in row['file.tags']:
            # Above tag can be used to mark files that can't be niftified and should be ignored
            return None
        elif '_ignore-BIDS' in label:
            return None
        elif label.startswith('GOBRAIN_'):
            return row['acquisition.label'] + '_ignore-BIDS'
        elif intent and 'Localizer' in intent:
            return row['acquisition.label'] + '_ignore-BIDS'
        elif measurement and 'Perfusion' in measurement:
            reproin = 'asl'
        elif measurement and intent and 'Structural' in intent:
            if 'T1' in measurement:
                reproin = 'anat-T1w'
            elif 'T2' in measurement:
                reproin = 'anat-T2w'   
            elif 'Diffusion' in measurement:
                reproin = 'dwi'
        elif (intent and 'Functional' in intent and
            label and 'rest' in label.lower()):
            reproin = 'func_task-rest'
        elif intent and 'Fieldmap' in intent:
            reproin = 'fmap'
        
        if reproin:
            if validator:
                return label
            return reproin
        
        # Use the pre-regex acq-label to avoid duplicate filenames since runs won't be added
        return row['acquisition.label'] + '_ignore-BIDS'

    classified_df = file_df.copy()
    breakpoint()
    classified_df = classified_df[classified_df['session.tags'].apply(lambda x: 'bidsified' not in x)]

    reproin_series = classified_df.apply(reproin_filter, axis=1)
    if reproin_series.empty:
        log.info("No acquisitions were found that need to be bidsified.")
        sys.exit(0)

    classified_df['reproin'] = reproin_series

    # Remove subjects if any included acqs don't have corresponding niftis
    filenames_joined = ','.join(classified_df[
        classified_df['file.type'] == 'nifti']['file.name'].tolist())
    bidsify_df = classified_df.copy()
    bidsify_df = bidsify_df[~bidsify_df['reproin'].isna()]
    bidsify_df = bidsify_df[bidsify_df['file.modality'] != 'SR']
    bidsify_df = bidsify_df[~bidsify_df['reproin'].str.contains('_ignore-BIDS')]
    bidsify_df = bidsify_df[bidsify_df['file.type'] == 'dicom'].groupby('subject.id').filter(
        lambda x: x['file_root'].apply(
            lambda y: y in filenames_joined
        ).all()
    )

    filt_out_set = set(classified_df['subject.id']) - set(bidsify_df['subject.id'])
    if filt_out_set:
        log.warning("The following subjects are not yet bidsified and filtered out: "
                    f"\n{filt_out_set}"
        )
    classified_df = classified_df[classified_df['subject.id'].isin(bidsify_df['subject.id'])]

    return classified_df

def add_repronim(classified_df: pd.DataFrame) -> pd.DataFrame:
    repronim_df = classified_df.copy()
    repronim_df = repronim_df[[
        'subject.id',
        'session.id',
        'acquisition.id',
        'acquisition.label',
        'acquisition.timestamp',
        'reproin'
    ]]
    # Produce a df with one row per acq
    repronim_df = repronim_df.drop_duplicates()

    def add_run(df: pd.DataFrame()):
        if df.shape[0] > 1:
            df = df.sort_values(by='acquisition.timestamp')
            n_digits = max(2, len(str(df.shape[0])))
            padded_list = [str(n).zfill(n_digits) for n in range(1, df.shape[0] + 1)]
            df.insert(df.shape[1], 'run', padded_list)
            df['reproin'] = df['reproin'] + '_run-' + df['run']
        return df

    repronim_df = repronim_df.groupby(['session.id', 'reproin']).apply(add_run).reset_index(drop=True)

    for name, ses_df in repronim_df.groupby('session.id'):
        for i, row in ses_df.iterrows():
            if row['reproin'] != row['acquisition.label']:
                acq_match_df = ses_df[ses_df['acquisition.label'] == row['reproin']]
                if not acq_match_df.empty:
                    acq_tmp = client.get_acquisition(acq_match_df['acquisition.id'].iloc[0])
                    acq_tmp.update({'label': row['reproin'] + '_tmp'})
                acq = client.get_acquisition(row['acquisition.id'])
                acq.update({'label': row['reproin']})
            
    return repronim_df

def rename_sessions(repronim_df: pd.DataFrame()) -> None:
    """Renames session to '01', or the next available number if existing sessions."""
    for sub_id in repronim_df['subject.id'].unique():
        subject = client.get_subject(sub_id)
        sub_sessions = subject.sessions()

        if len(sub_sessions) == 1:
            new_session_label = '01'
            session = sub_sessions[0]
            session.update({'label': new_session_label})
        else:
            sub_sessions_sorted = sorted(sub_sessions, key=lambda d: d.timestamp)
            num_digits = len(str(len(sub_sessions_sorted)))
            zero_pad = max(num_digits, 2)
            for i, session in enumerate(sub_sessions_sorted, 1):
                new_session_label = str(i).zfill(zero_pad)
                session.update({'label': new_session_label})

def submit_bids_jobs(repronim_df: pd.DataFrame()):
    curate_bids_gear = client.lookup('gears/curate-bids')
    for sub_id in repronim_df['subject.id'].unique():
        subject = client.get_subject(sub_id)
        run_gear(curate_bids_gear, {}, {}, subject)

def wait_for_jobs(project: ProjectOutput) -> None:
    query = f'parents.project={project.id},state=running,gear_info.name=curate-bids'
    start_time = time.time()
    while client.jobs.find(query):
        if time.time() - start_time > WAIT_TIMEOUT:
            log.error("Wait timeout for copy to complete")
            sys.exit(1)
        time.sleep(5)

def tag_and_email(project: ProjectOutput) -> None:
    cutoff = (datetime.now() - timedelta(hours=6)).isoformat()

    failed_q = f'parents.project={project.id},state=failed,gear_info.name=curate-bids,created>{cutoff}'
    failed_jobs = client.jobs.find(failed_q)
    failed_job_subjects = [job.parents.subject for job in failed_jobs]
    for sub_id in failed_job_subjects:
        sub = client.get_subject(sub_id)
        for ses in sub.sessions(): 
            if 'bids-failed' not in ses.tags:
                ses.add_tag('bids-failed')

    complete_q = f'parents.project={project.id},state=complete,gear_info.name=curate-bids,created>{cutoff}'
    complete_jobs = client.jobs.find(complete_q)
    complete_job_subjects = [job.parents.subject for job in complete_jobs]
    for sub_id in complete_job_subjects:
        sub = client.get_subject(sub_id)
        for ses in sub.sessions():
            if 'bids-failed' in ses.tags:
                ses.delete_tag('bids-failed')
            if 'bidsified' not in ses.tags:
                ses.add_tag('bids-failed')
    
    send_email(
        'failed bids-curate jobs',
        ', '.join(failed_job_subjects),
        config["gmail_address"],
        ADMIN_EMAIL,
        config["gmail_password"]
    )

def send_email(subject, html_content, sender, recipients, password):
    msg = MIMEMultipart()
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = ', '.join(recipients)
    msg.attach(MIMEText(html_content, 'html'))

    with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp_server:
       smtp_server.login(sender, password)
       smtp_server.sendmail(sender, recipients, msg.as_string())

    log.info(f"Email sent to {recipients}")

def run_gear(
    gear: Gear,
    inputs: dict,
    config: dict,
    dest,
    tags=None) -> str:
    """Submits a job with specified gear and inputs. dest can be any type of container
    that is compatible with the gear (project, subject, session, acquisition)"""
   
    for i in range(0,3):
        try:
            # Run the gear on the inputs provided, stored output in dest constainer and returns job ID
            gear_job_id = gear.run(inputs=inputs, config=config, destination=dest, tags=tags)
            log.debug('Submitted job %s', gear_job_id)
            return gear_job_id
        except flywheel.rest.ApiException:
            #log.exception('An exception was raised when attempting to submit a job for %s', gear.name)
            time.sleep(1)
def main():
    gtk_context.init_logging()
    gtk_context.log_config()

    deid_project = client.lookup("joe_test/deid_joe")
    #deid_project = client.lookup("wbhi/deid")

    try:
        file_df = get_subjects(deid_project)
    except (KeyError, NameError):
        log.info("All sessions have already been bidsified. Exiting")
        sys.exit(0)

    file_df = classify(file_df)
    file_df = add_repronim(file_df)
    rename_sessions(file_df)
    submit_bids_jobs(file_df)
    wait_for_jobs(deid_project)
    tag_and_email(deid_project)
    breakpoint()

if __name__ == "__main__":
    with flywheel_gear_toolkit.GearToolkitContext() as gtk_context:
        config = gtk_context.config
        client = gtk_context.client
        
        main()

