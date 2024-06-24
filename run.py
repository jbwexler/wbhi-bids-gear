#!/usr/bin/env python

import csv
import flywheel_gear_toolkit
import flywheel
import logging
import pip
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime, timedelta
import pandas as pd
from redcap import Project
from flywheel import (
    ProjectOutput,
    SessionListOutput,
    AcquisitionListOutput,
    FileListOutput,
    Gear
)

pip.main(["install", "--upgrade", "git+https://github.com/poldracklab/wbhi-utils.git"])
from wbhiutils import parse_dicom_hdr
from wbhiutils.constants import *

log = logging.getLogger(__name__)

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

def create_file_view_df(container) -> pd.DataFrame:
    columns = [
        'subject.id',
        'subject.label',
        'session.id',
        'session.tags',
        'session.info.BIDS',
        'acquisition.label',
        'file.file_id',
        'file.tags',
        'file.type',
        'file.created'
    ]
    return create_view_df(container, columns)

def create_session_view_df(container) -> pd.DataFrame:
    columns = [
        'subject.label',
        'session.id',
        'session.tags',
        'acquisition.label',
        'acquisition.timestamp'
    ]
    breakpoint()
    return create_view_df(container, columns, filter)

def get_deid_sessions(file_df: pd.DataFrame) -> list:
    dcm_tags = {
            'file-metadata-importer-deid',
            'file-classifier-deid',
            'dcm2niix'
    }
    nii_tag = 'pydeface'

    # Skip subject if all sessions have been bidsified
    nonbids_df = file_df.groupby('subject.id').filter(
        lambda x: x['session.info.BIDS'].isna().any()
    )
    
    # Skip sessions if any dicoms or niftis missing proper tags
    nii_df = nonbids_df[nonbids_df['file.type'] == 'nifti']
    nii_df = nii_df[nii_df['file.tags'].apply(lambda x: nii_tag in x)]
    filt_df = nonbids_df[nonbids_df['session.id'].isin(nii_df['session.id'])]
    dcm_df = nii_df[nii_df['file.type'] == 'dicom']
    dcm_df = dcm_df[dcm_df['file.tags'].apply(
        lambda x, tags=dcm_tags: dcm_tags.issubset(x)
    )]
    if not dcm_df.empty:
        filt_df = filt_df[filt_df['session.id'].isin(dcm_df['session.id'])]
    else:
        filt_df = filt_df[0:0]

    # Raise warning if some non-bidsified sessions are filtered out
    nonbids_ses_set = set(nonbids_df['subject.label'])
    filt_ses_set = set(filt_df['subject.label'])
    filt_out_set = nonbids_ses_set - filt_ses_set
    if filt_out_set:
        log.warning("The following subjects are not yet bidsified but "
                    f"contain files that are missing necessary tags: \n{filt_out_set}"
        )
    return filt_df['session.id'].tolist()

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

def main():
    gtk_context.init_logging()
    gtk_context.log_config()

    redcap_api_key = config["redcap_api_key"]
    redcap_project = Project(REDCAP_API_URL, redcap_api_key)
    
    deid_project = client.lookup("wbhi/deid")

if __name__ == "__main__":
    with flywheel_gear_toolkit.GearToolkitContext() as gtk_context:
        config = gtk_context.config
        client = gtk_context.client
        
        main()

