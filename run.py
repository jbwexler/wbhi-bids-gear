#!/usr/bin/env python

import flywheel_gear_toolkit
import logging
import os
import pip
import sys
import time
import pandas as pd
from flywheel import ProjectOutput
from utils.reproin import (
    check_pydeface,
    check_nifti_presence,
    check_asl_and_spectroscopy,
    check_anat_presence,
    reproin_filter,
    add_sbref,
    add_fmap,
    add_run,
    add_rec,
    build_acq_df,
    check_duplicate_reproin,
    rm_acq_labels,
)
from utils.flywheel import create_view_df, send_email, run_gear, mv_session
from utils.constants import WAIT_TIMEOUT, DATAVIEW_COLUMNS

pip.main(["install", "--upgrade", "git+https://github.com/poldracklab/wbhi-utils.git"])
from wbhiutils.constants import (  # noqa: E402
    ADMIN_EMAIL,
)

log = logging.getLogger(__name__)

PKL_PATH = "output/file_df.pkl"


def get_subjects(project: ProjectOutput) -> pd.DataFrame:
    """Returns a df of files from subjects that are ready to be bidsified."""
    if config.get("test_mode") and os.path.exists(PKL_PATH):
        log.info(f"Reusing df from {PKL_PATH}.")
        file_df = pd.read_pickle(PKL_PATH)
    else:
        file_df = create_view_df(project, DATAVIEW_COLUMNS.keys(), client)
        if config.get("test_mode"):
            log.info(f"Saving df as {PKL_PATH}.")
            file_df.to_pickle(PKL_PATH)

    if file_df.empty:
        log.info("No subjects were found in deid project.")
        sys.exit(0)

    bidsified_sessions = file_df[
        file_df["session.tags"].apply(lambda x: "bidsified" in x)
    ]["session.id"].unique()
    for session_id in bidsified_sessions:
        client.get_session(session_id).delete_tag("bidsified")

    file_df["sbref"] = None
    file_df["error"] = ""

    if file_df.empty:
        log.info("No subjects were found that need to be bidsified.")
        sys.exit(0)

    return file_df


def get_acq_datetime(df: pd.DataFrame()) -> pd.Series():
    """Returns a Series containing the acquisition datetime."""
    df = df.copy()

    df["acq_datetime"] = df["file.info.header.dicom.AcquisitionDateTime"]

    acq_date_plus_time = df["file.info.header.dicom.AcquisitionDate"].fillna("") + df[
        "file.info.header.dicom.AcquisitionTime"
    ].fillna("")
    df["acq_datetime"] = (
        df["acq_datetime"].replace("", pd.NA).fillna(acq_date_plus_time)
    )

    df["acq_datetime"] = pd.to_datetime(df["acq_datetime"], format="%Y%m%d%H%M%S.%f")
    return df["acq_datetime"].replace("", pd.NA)


def rename_sessions(df: pd.DataFrame()) -> None:
    """Renames session to '01', or the next available number if existing sessions."""
    df_relabel = df[df["session.label"] != "01"]

    for sub_id in df_relabel["subject.id"].unique():
        subject = client.get_subject(sub_id)
        sub_sessions = subject.sessions()

        if len(sub_sessions) == 1:
            session = sub_sessions[0]

            log.debug("Renaming session %s to 01" % (session.id,))
            session.update({"label": "01"})
        else:
            sub_sessions_sorted = sorted(sub_sessions, key=lambda d: d.timestamp)
            num_digits = len(str(len(sub_sessions_sorted)))
            zero_pad = max(num_digits, 2)

            for i, session in enumerate(sub_sessions_sorted, 1):
                new_session_label = str(i).zfill(zero_pad)

                if session.label != new_session_label:
                    log.debug(
                        "Renaming session %s to %s" % (session.id, new_session_label)
                    )
                    session.update({"label": new_session_label})


def skip_multiple_dicoms(file_df: pd.DataFrame) -> pd.DataFrame:
    """Removes any subjects containing any acquisitions with multiple dicoms."""
    file_df = file_df.copy()

    mult_df = file_df[
        ~file_df["file.classification.Intent"].apply(
            lambda x: type(x) is list and "Localizer" in x
        )
    ]
    mult_df = mult_df.groupby("acquisition.id").filter(
        lambda x: len(x[x["file.type"] == "dicom"]) > 1
    )
    mult_df_subs = mult_df["subject.label"].unique()

    if mult_df_subs.size > 0:
        log.info(
            """Skipping the following subjects because they contain files with
            multiple dicoms: %s"""
            % mult_df_subs
        )

    return file_df[~file_df["subject.label"].isin(mult_df_subs)]


def partition_ignored(
    acq_df: pd.DataFrame, ignored_df: pd.DataFrame = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Moves acquisitions whose reproin contains '_ignore-BIDS' from acq_df into
    ignored_df. Returns (acq_df, ignored_df)."""
    is_ignored = acq_df["reproin"].str.contains("_ignore-BIDS", na=False)
    ignored_df = pd.concat([ignored_df, acq_df[is_ignored]])
    return acq_df[~is_ignored], ignored_df


def classify(file_df: pd.DataFrame) -> pd.DataFrame:
    """Determines which acquisitions should be bidsified and adds columns containing
    appropriate reproin labels, as well as fmap, rec, run and sbref info. For acquisitions
    that won't be bidsified, adds '_ignore-BIDS' to the end of reproin label."""
    file_df = file_df.copy()

    file_df = file_df[~file_df["acquisition.label"].str.endswith("_ignore-BIDS")]
    file_df = skip_multiple_dicoms(file_df)

    if file_df.empty:
        log.info("No acquisitions were found that need to be bidsified.")
        sys.exit(0)

    acq_df = build_acq_df(file_df)
    acq_df["acq_datetime"] = get_acq_datetime(acq_df)

    acq_df["reproin"] = acq_df.apply(reproin_filter, axis=1)
    acq_df = (
        acq_df.groupby("session.id")[acq_df.columns]
        .apply(rm_acq_labels)
        .reset_index(drop=True)
    )

    acq_df = check_nifti_presence(acq_df)
    acq_df = check_pydeface(acq_df)
    acq_df = check_asl_and_spectroscopy(acq_df)
    acq_df = check_anat_presence(acq_df)

    acq_df, ignored_df = partition_ignored(acq_df)

    acq_df = (
        acq_df.groupby("session.id")[acq_df.columns]
        .apply(add_fmap)
        .reset_index(drop=True)
    )

    acq_df, ignored_df = partition_ignored(acq_df, ignored_df)

    acq_df = (
        acq_df.groupby("session.id")[acq_df.columns]
        .apply(add_rec)
        .reset_index(drop=True)
    )

    acq_df, ignored_df = partition_ignored(acq_df, ignored_df)

    acq_df = (
        acq_df.groupby(["session.id", "reproin"])[acq_df.columns]
        .apply(add_run)
        .reset_index(drop=True)
    )
    acq_df = (
        acq_df.groupby("session.id")[acq_df.columns]
        .apply(add_sbref)
        .reset_index(drop=True)
    )

    return pd.concat([acq_df, ignored_df], ignore_index=True)


def add_reproin(file_df: pd.DataFrame) -> pd.DataFrame:
    file_df = file_df.copy()
    breakpoint()

    file_df = (
        file_df.groupby("session.id")[file_df.columns]
        .apply(check_duplicate_reproin)
        .reset_index(drop=True)
    )
    file_df = file_df.groupby("subject.id").filter(
        lambda x: x["error"].apply(lambda x: not x).all()
    )

    for name, ses_df in file_df.groupby("session.id"):
        rename_df = ses_df[ses_df["reproin"] != ses_df["acquisition.label"]]
        acqs = [
            client.get_acquisition(acq_id) for acq_id in rename_df["acquisition.id"]
        ]

        # Rename to temporary labels first so no acquisition is renamed to a label
        # that another acquisition still holds.
        for acq, reproin in zip(acqs, rename_df["reproin"]):
            acq.update({"label": reproin + "_tmp"})
        for acq, reproin in zip(acqs, rename_df["reproin"]):
            acq.update({"label": reproin})

    return file_df


def submit_bids_jobs(file_df: pd.DataFrame()) -> list[str]:
    """Submits a curate-bids job for each subject. Returns the job ids."""
    curate_bids_gear = client.lookup("gears/curate-bids")
    job_ids = []
    for sub_id in file_df["subject.id"].unique():
        subject = client.get_subject(sub_id)
        job_id = run_gear(curate_bids_gear, {}, {"reset": True}, subject)
        if job_id:
            job_ids.append(job_id)
        else:
            log.error("Failed to submit curate-bids job for subject %s" % subject.label)
    return job_ids


def wait_for_jobs(job_ids: list[str]) -> list:
    """Waits until all jobs have finished. Returns the finished jobs."""
    pending = set(job_ids)
    finished = []
    start_time = time.time()
    while pending:
        for job_id in list(pending):
            job = client.get_job(job_id)
            if job.state in ("complete", "failed", "cancelled"):
                finished.append(job)
                pending.remove(job_id)
        if pending:
            if time.time() - start_time > WAIT_TIMEOUT:
                log.error("Timed out waiting for curate-bids jobs to complete")
                sys.exit(1)
            time.sleep(30)
    return finished


def tag_and_email(jobs: list) -> None:
    """Tags each job's subject's sessions as bidsified or bids-failed, and emails
    the labels of subjects whose job didn't complete."""
    failed_subjects = []
    for job in jobs:
        sub = client.get_subject(job.parents.subject)
        if job.state == "complete":
            add_tag, rm_tag = "bidsified", "bids-failed"
        else:
            add_tag, rm_tag = "bids-failed", "bidsified"
            failed_subjects.append(sub.label)

        for ses in sub.sessions():
            if add_tag not in ses.tags:
                ses.add_tag(add_tag)
            if rm_tag in ses.tags:
                ses.delete_tag(rm_tag)

    if failed_subjects:
        send_email(
            "failed bids-curate jobs",
            ", ".join(failed_subjects),
            config["gmail_address"],
            ADMIN_EMAIL,
            os.environ["GMAIL_PASSWORD"],
        )


def mv_bidsified(src_project: ProjectOutput, dst_project: ProjectOutput) -> None:
    columns = ["session.id", "session.tags"]
    df = create_view_df(src_project, columns, client, container_type="session")
    bidsified_df = df[
        df["session.tags"].apply(lambda x: "bidsified" in x and "bids-failed" not in x)
    ]

    for i, ses_id in bidsified_df["session.id"].items():
        session = client.get_session(ses_id)
        mv_session(session, dst_project, client)


def main():
    gtk_context.init_logging()
    gtk_context.log_config()

    if config["test_mode"]:
        group = "joe_test"
    else:
        group = "wbhi"

    deid_project = client.lookup(f"{group}/deid")
    file_df = get_subjects(deid_project)
    rename_sessions(file_df)

    file_df = classify(file_df)
    file_df = add_reproin(file_df)

    job_ids = submit_bids_jobs(file_df)
    jobs = wait_for_jobs(job_ids)
    tag_and_email(jobs)

    staging_project = client.lookup(f"{group}/staging")
    mv_bidsified(deid_project, staging_project)


if __name__ == "__main__":
    with flywheel_gear_toolkit.GearToolkitContext() as gtk_context:
        config = gtk_context.config
        client = gtk_context.client

        main()
