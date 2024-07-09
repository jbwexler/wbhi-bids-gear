#!/usr/bin/env python

import flywheel_gear_toolkit
import flywheel
import logging
import pip
import smtplib
import sys
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime
import pandas as pd
from flywheel import (
    ProjectOutput,
    AcquisitionListOutput,
    FileListOutput,
    Gear
)

pip.main(["install", "--upgrade", "git+https://github.com/poldracklab/wbhi-utils.git"])
from wbhiutils import parse_dicom_hdr
from wbhiutils.constants import (
    REDCAP_API_URL,
    DATE_FORMAT_FW
)

log = logging.getLogger(__name__)

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
    'file.file_id',
    'file.tags',
    'file.type',
    'file.created',
    'file.name'
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
        mask = nii_df['file.tags'].apply(lambda x: 'dcm2niix' not in x).any()
    except KeyError:
        return False
    return mask

def get_subjects(project: ProjectOutput) -> list:
    """Returns a list of subject IDs of subjects that have not yet been bidsified
    and are ready to be."""

    file_df = create_view_df(project, DATAVIEW_COLUMNS) 
    file_df = file_df[file_df['file.type'].isin(('dicom', 'nifti'))]

    # Skip subject if all sessions have been bidsified
    file_df = file_df.groupby('subject.id').filter(
        lambda x: x['session.info.BIDS'].isna().any()
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
        lambda x: x.split('.dicom')[0]
    )

    # Skips subjects if not all dicoms have corresponding niftis
    filenames_joined = ','.join(file_df[file_df['file.type'] == 'nifti']['file.name'].tolist())
    dcm_df = file_df.copy()
    dcm_df = dcm_df[dcm_df['file.type'] == 'dicom'].groupby('subject.id').filter(
        lambda x: x['file_root'].apply(
            lambda y: y in filenames_joined
        ).all()
    )
    file_df = file_df[file_df['subject.id'].isin(dcm_df['subject.id'])]

    # Raise warning if some non-bidsified sessions are filtered out
    filt_out_set = set(nonbids_df['subject.id']) - set(file_df['subject.id'])
    if filt_out_set:
        log.warning("The following subjects are not yet bidsified and filtered out: "
                    f"\n{filt_out_set}"
        )

    return file_df['subject.id'].unique().tolist()


def get_acq_path(acq: AcquisitionListOutput) -> str:
    """Takes an acquisition and returns its path:
    project/subject/session/acquisition"""
    project_label = client.get_project(acq.parents.project).label
    sub_label = client.get_subject(acq.parents.subject).label
    ses_label = client.get_session(acq.parents.session).label

    return f"{project_label}/{sub_label}/{ses_label}/{acq.label}"

def get_hdr_fields(dicom: FileListOutput, site: str) -> dict:
    """Get relevant fields from dicom header of an acquisition"""
    if "file-classifier" not in dicom.tags or "header" not in dicom.info:
        log.error(f"File-classifier gear has not been run on {get_acq_path(acq)}")
        return {"error": "FILE_CLASSIFIER_NOT_RUN"}

    dcm_hdr = dicom.reload().info["header"]["dicom"]
    return {
        "error": None,
        "site": site,
        "date": datetime.strptime(dcm_hdr["StudyDate"], DATE_FORMAT_FW),
        "am_pm": "am" if float(dcm_hdr["StudyTime"]) < 120000 else "pm",
        "pi_id": parse_dicom_hdr.parse_pi(dcm_hdr, site).casefold(),
        "sub_id": parse_dicom_hdr.parse_sub(dcm_hdr, site).casefold()
    }

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

    try:
        # Run the gear on the inputs provided, stored output in dest constainer and returns job ID
        gear_job_id = gear.run(inputs=inputs, config=config, destination=dest, tags=tags)
        log.debug('Submitted job %s', gear_job_id)
        return gear_job_id
    except flywheel.rest.ApiException:
        log.exception('An exception was raised when attempting to submit a job for %s', gear.name)

def main():
    gtk_context.init_logging()
    gtk_context.log_config()

    deid_project = client.lookup("wbhi/deid")
    try:
        sub_list = get_subjects(deid_project)
    except (KeyError, NameError):
        log.info("All sessions have already been bidsified. Exiting")
        sys.exit(0)

if __name__ == "__main__":
    with flywheel_gear_toolkit.GearToolkitContext() as gtk_context:
        config = gtk_context.config
        client = gtk_context.client
        
        main()

