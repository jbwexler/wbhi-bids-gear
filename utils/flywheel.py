#!/usr/bin/env python

import flywheel
from flywheel import Gear
import pandas as pd
import time
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import smtplib
import logging

log = logging.getLogger(__name__)


def create_view_df(container, columns: list, client, filter=None) -> pd.DataFrame:
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



