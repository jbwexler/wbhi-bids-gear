#!/usr/bin/env python

import flywheel_gear_toolkit
import logging
import pip
import sys
from datetime import datetime, timedelta
import time
import pandas as pd
from flywheel import ProjectOutput
from flywheel.client import Client
from utils.reproin import (
    pydeface_filter,
    reproin_filter,
    add_nifti_filenames,
    add_sbref,
    add_fmap,
    add_run,
)
from utils.flywheel import (
    create_view_df,
    send_email,
    run_gear
)
from utils.constants import WAIT_TIMEOUT, DATAVIEW_COLUMNS

pip.main(["install", "--upgrade", "git+https://github.com/poldracklab/wbhi-utils.git"])
from wbhiutils.constants import ( # noqa: E402
    ADMIN_EMAIL
)

log = logging.getLogger(__name__)



def get_subjects(project: ProjectOutput, client: Client) -> pd.DataFrame:
    """Returns a list of subject IDs of subjects that have not yet been bidsified
    and are ready to be."""

    file_df = create_view_df(project, DATAVIEW_COLUMNS, client) 
    #file_df = file_df[file_df['file.type'].isin(('dicom', 'nifti'))]
    file_df['acquisition.timestamp'] = pd.to_datetime(file_df['acquisition.timestamp'])
    file_df['sbref'] = None
    
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

    # Skip subjects that have any acquisitions without niftis
    has_nii_df = file_df.copy()

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

def split_multiple_dicoms(df: pd.DataFrame, client):
    """Finds any acquisitions with multiple dicoms and splits them into
    multiple acquisitions"""
    mult_df = df.groupby('acquisition.id').filter(
        lambda x: len(x[x['file.type'] == 'dicom']) > 1
    )
    for acq_name,acq_group in mult_df.groupby('acquisition.id'):
        for i,dcm in acq_group[acq_group['file.type'] == 'dicom'].iloc[1:].iterrows():
            new_acq_df = acq_group[acq_group['file.name'].apply(
                lambda x, dcm=dcm: dcm['file_root'] in x
            )]

def classify(file_df: pd.DataFrame) -> pd.DataFrame:
    """Determines which acquisitions should be bidsified and converts the acquisition labels
    to reproin standard. For acquisitions that won't be bidsified, adds '_ignore-BIDS' to
    the end of the acquisition label, signaling the bids-curate gear to skip."""
    classified_df = file_df.copy()
    classified_df = classified_df[classified_df['session.tags'].apply(lambda x: 'bidsified' not in x)]

    reproin_series = classified_df.apply(reproin_filter, axis=1)
    if reproin_series.empty:
        log.info("No acquisitions were found that need to be bidsified.")
        sys.exit(0)

    classified_df['reproin'] = reproin_series

    # Remove subjects if any included acqs don't have corresponding niftis
    bidsify_df = classified_df.copy()
    bidsify_df = bidsify_df[~bidsify_df['reproin'].isna()]
    bidsify_df = bidsify_df[~bidsify_df['reproin'].str.contains('_ignore-BIDS')]
    bidsify_df = bidsify_df.groupby('subject.id').apply(add_nifti_filenames).reset_index(drop=True)
    breakpoint()

    filt_out_set = set(classified_df['subject.label']) - set(bidsify_df['subject.label'])
    if filt_out_set:
        log.warning("The following subjects were not bidsified: "
          f"\n{filt_out_set}"
        )

    return bidsify_df

def add_reproin(classified_df: pd.DataFrame) -> pd.DataFrame:
    reproin_df = classified_df.copy()

    print(reproin_df.shape)
    reproin_df = reproin_df.groupby('session.id').apply(add_fmap).reset_index(drop=True)
    reproin_df = reproin_df.groupby(['session.id', 'reproin']).apply(add_run).reset_index(drop=True)
    reproin_df = reproin_df.groupby('session.id').apply(add_sbref).reset_index(drop=True)

    #####################
    columns = [
        'subject.label',
        'reproin',
        'file.info.header.dicom.SeriesDescription',
        'file.info.header.dicom_array.ImageType.0',
        'file.classification.Intent',
        'file.classification.Measurement',
        'file.classification.Features',
        'acquisition.timestamp',
        'file.modality',
        'file.type'
    ]
    output_path = gtk_context.output_dir.joinpath('reproin_dryrun.csv')
    reproin_df[columns].to_csv(output_path)
    breakpoint()
    #####################

    for name, ses_df in reproin_df.groupby('session.id'):
        for i, row in ses_df.iterrows():
            if row['reproin'] != row['acquisition.label']:
                acq_match_df = ses_df[ses_df['acquisition.label'] == row['reproin']]
                if not acq_match_df.empty:
                    acq_tmp = client.get_acquisition(acq_match_df['acquisition.id'].iloc[0])
                    acq_tmp.update({'label': row['reproin'] + '_tmp'})
                acq = client.get_acquisition(row['acquisition.id'])
                acq.update({'label': row['reproin']})
    return reproin_df

def rename_sessions(reproin_df: pd.DataFrame()) -> None:
    """Renames session to '01', or the next available number if existing sessions."""
    for sub_id in reproin_df['subject.id'].unique():
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

def submit_bids_jobs(reproin_df: pd.DataFrame()):
    curate_bids_gear = client.lookup('gears/curate-bids')
    for sub_id in reproin_df['subject.id'].unique():
        subject = client.get_subject(sub_id)
        run_gear(curate_bids_gear, {}, {'reset': True}, subject)

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


def main():
    gtk_context.init_logging()
    gtk_context.log_config()

    deid_project = client.lookup("joe_test/deid_joe")
    #deid_project = client.lookup("wbhi/deid")

    try:
        file_df = get_subjects(deid_project, client)
    except (KeyError, NameError):
        log.info("All sessions have already been bidsified. Exiting")
        sys.exit(0)

    #file_df = split_multiple_dicoms(file_df)
    file_df = classify(file_df)
    file_df = add_reproin(file_df)
    submit_bids_jobs(file_df)
    wait_for_jobs(deid_project)
    tag_and_email(deid_project)

if __name__ == "__main__":
    with flywheel_gear_toolkit.GearToolkitContext() as gtk_context:
        config = gtk_context.config
        client = gtk_context.client
        
        main()

