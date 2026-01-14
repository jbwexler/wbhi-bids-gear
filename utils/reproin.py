#!/usr/bin/env python

import re
import pandas as pd
from heudiconv.heuristics.reproin import parse_series_spec
import logging
from datetime import timedelta
from utils.constants import (
    ALLOWED_DATATYPES,
    SBREF_DELTA,
    FMAP_DELTA,
    REC_DELTA,
    FMAP_SUFFIXES,
    IGNORE_SCANS,
    DWI_SUFFIXES_IGNORE,
)

log = logging.getLogger(__name__)

REC_LIST = []


def pydeface_filter(df: pd.DataFrame) -> bool:
    """Returns True if any acquisitions containing a Structural T1 or T2 dicom
    DO NOT also contain a nifti with a 'pydeface' label. Otherwise, returns False."""

    # Check for any Structural T1 or T2 dicoms
    dcm_df = df.copy()
    try:
        dcm_df = dcm_df[dcm_df["file.type"] == "dicom"]
        dcm_df = dcm_df[~dcm_df["file.classification.Intent"].isna()]
        dcm_df = dcm_df[~dcm_df["file.classification.Measurement"].isna()]
        mask = dcm_df["file.classification.Intent"].apply(lambda x: "Structural" in x)
        dcm_df = dcm_df[mask]
        mask = dcm_df["file.classification.Measurement"].apply(
            lambda x: "T1" in x or "T2" in x
        )
        dcm_df = dcm_df[mask]
    except KeyError:
        return False
    if dcm_df.empty:
        return False

    # If structurals present, check for niftis without 'pydeface' tag
    try:
        nii_df = df.copy()
        nii_df = nii_df[nii_df["file.type"] == "nifti"]
        mask = nii_df["file.tags"].apply(lambda x: "pydeface" not in x).any()
    except KeyError:
        return False
    return mask


def check_pydeface(file_df: pd.DataFrame) -> pd.DataFrame:
    """Returns a df of files from subjects for which all structural T1 and T2 niftis
    have 'pydeface' tag"""
    file_df = file_df.copy()

    if file_df.empty:
        return file_df

    no_pydeface_df = file_df.groupby("acquisition.id").filter(pydeface_filter)
    file_df = file_df[
        ~file_df["subject.label"].isin(no_pydeface_df["subject.label"].unique())
    ]

    if not no_pydeface_df.empty:
        log.warning(
            "The following subjects contain anatomical nifti(s) without the 'pydeface' tag: "
            "\n%s" % no_pydeface_df['subject.label'].unique()
        )

    return file_df


def check_nifti_presence(file_df: pd.DataFrame) -> pd.DataFrame:
    """Removes subjects that have any acquisitions without niftis"""
    file_df = file_df.copy()
    no_nii_df = file_df.copy()

    no_nii_df = no_nii_df[
        ~no_nii_df["reproin"].str.endswith(("_ignore-BIDS"), na=False)
    ]
    no_nii_df = no_nii_df[no_nii_df["file.modality"] == "MR"]
    no_nii_df = (
        no_nii_df.groupby("acquisition.id")
        .filter(lambda x: not x["file.type"].str.contains("nifti").any())
        .reset_index(drop=True)
    )

    file_df = file_df[~file_df["subject.label"].isin(no_nii_df["subject.label"])]
    if not no_nii_df.empty:
        pd.options.display.max_rows = 5000
        log.warning(
            "The following subjects contain acquisition(s) without nifti(s): "
            "\n%s" % no_nii_df[['subject.label', 'acquisition.label']]
        )

    return file_df


def reproin_filter(row: pd.Series) -> str:
    label = row["acquisition.label"]
    if row["file.type"] != "dicom":
        return label

    # Utilize acq label if already in reproin format. Otherwise, create a reproin label
    label_re = re.sub(r"_\d$", "", label)
    label_ignore = label + "_ignore-BIDS"
    if label_re in IGNORE_SCANS:
        return label_ignore
    image_type = row["file.info.header.dicom.ImageType"]
    # to-do: add exception for mp2rage
    if image_type and "derived" in image_type:
        return label_ignore
    try:
        validator = parse_series_spec(label)
    except IndexError:
        validator = {}
    if (
        validator
        and validator["datatype"] == "anat"
        and "datatype_suffix" not in validator
    ):
        validator = None

    if validator:
        if validator["datatype"] not in ALLOWED_DATATYPES:
            return label_ignore
        elif validator["datatype"] == "func":
            if "task" in validator and (
                "rest" in validator["task"] or validator["task"] == "rs"
            ):
                return re.sub("task-.*?(_|$)", r"task-rest\1", label_re)
            else:
                return label_ignore
        return label_re
    else:
        measurement = row["file.classification.Measurement"]
        intent = row["file.classification.Intent"]
        features = row["file.classification.Features"]

        if label.startswith("GOBRAIN_"):
            return label_ignore
        elif label == "t2_tse_tra_hi-res_hippocampus":
            return "anat-T2w_acq-hippocampus"
        elif label_re.endswith(("PhysioLog", "setter", "TENSOR")):
            return label_ignore
        elif intent and ("Localizer" in intent or "Spectroscopy" in intent):
            return label_ignore
        elif measurement and "Perfusion" in measurement:
            return "asl"
        elif measurement and intent and "Structural" in intent:
            if "T1" in measurement:
                return "anat-T1w"
            elif "T2" in measurement:
                if features and "FLAIR" in features:
                    return "anat-FLAIR"
                else:
                    return "anat-T2w"
            elif "Diffusion" in measurement:
                if label_re.endswith(DWI_SUFFIXES_IGNORE):
                    return label_ignore
                else:
                    return "dwi"
            else:
                return label_ignore
        elif (
            intent
            and "Functional" in intent
            and label
            and ("rest" in label.lower() or "rsfmri" in label.lower())
        ):
            return "func_task-rest"
        elif intent and "Fieldmap" in intent:
            return "fmap"
        else:
            return label_ignore


def classify_fmap_acq(image_type_list: list) -> str:
    """Determines the type of fieldmap file (phase, magnitude or diffusion) based on
    the presence of corresponding strings in the ImageType dict."""
    if not image_type_list:
        return None
    elif "M" in image_type_list:
        return "magnitude"
    elif "P" in image_type_list:
        return "phase"
    elif "DIFFUSION" in image_type_list:
        return "diffusion"
    elif "FMRI" in image_type_list:
        return "fmri"


def get_runs(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    run_df = df.sort_values(by="acquisition.timestamp")
    n_digits = max(2, len(str(run_df.shape[0])))
    padded_list = [str(n).zfill(n_digits) for n in range(1, run_df.shape[0] + 1)]
    run_df.insert(run_df.shape[1], "run", padded_list)
    return run_df


def add_fmap_runs(df: pd.DataFrame()) -> pd.DataFrame():
    df = df.copy()
    group_df = df[["acquisition.timestamp", "group"]].drop_duplicates(
        subset="group", keep="first"
    )
    if len(group_df) > 1:
        group_mapping = get_runs(group_df).set_index("group")
        df["run"] = df["group"].map(group_mapping["run"])
        df["reproin"] = df["reproin"] + "_run-" + df["run"]

    return df


def timestamp_group(df: pd.DataFrame, delta: timedelta) -> pd.DataFrame:
    """Groups rows in a dataframe that fall within a time delta. Returns a df with
    added 'time_diff' and 'group' columns"""
    df = df.copy()
    df = df.sort_values(by="acquisition.timestamp")
    df["time_diff"] = df["acquisition.timestamp"].diff().fillna(pd.Timedelta(seconds=0))
    df["group"] = (df["time_diff"] > delta).cumsum()
    return df


def get_case_four(group: pd.DataFrame(), fmap_file_types: dict, subject: str):
    group = group.copy()

    if not fmap_file_types["diffusion"].empty:
        non_phase = fmap_file_types["diffusion"]

        if (
            group[group["file.type"] == "bval"].empty
            or group[group["file.type"] == "bvec"].empty
        ):
            err_msg = "Subject %s has diffusion fmaps without bval/bvecs." % subject
            log.error(err_msg)
            group.loc[:, "error"] = err_msg
            return group

        group["fmap_acq"] = "dwi"
    elif not fmap_file_types["magnitude"].empty:
        non_phase = fmap_file_types["magnitude"]
        group["fmap_acq"] = "bold"
    elif not fmap_file_types["fmri"].empty:
        non_phase = fmap_file_types["fmri"]
        group["fmap_acq"] = "bold"

    non_phase["direction"] = non_phase["acquisition.label"].str.extract(
        r".*(LR|RL|AP|PA)(?=[^a-zA-Z0-9]|$)"
    )[0]

    if non_phase["direction"].isnull().any():
        err_msg = "Subject %s has acquisition with 2 epis but no direction." % subject
        log.error(err_msg)
        group.loc[:, "error"] = err_msg
        return group
    elif len(non_phase) == 2 and len(non_phase["acquisition.id"].unique()) == 1:
        err_msg = "Subject %s has 2 epis in one acquisition." % subject
        log.error(err_msg)
        group.loc[:, "error"] = err_msg
        return group

    reproin_mapping = non_phase.groupby("acquisition.id")[non_phase.columns].apply(
        lambda x: "fmap-epi_dir-" + x["direction"].iloc[0]
    )
    group.loc[:, "reproin"] = group["acquisition.id"].map(reproin_mapping)

    return group


def add_fmap(df: pd.DataFrame()):
    df = df.copy()
    df = df[~df["reproin"].str.endswith("_ignore-BIDS", na=False)]
    if df.empty:
        return df

    subject = df["subject.label"].iloc[0]
    fmap_df = df[df["file.modality"] == "MR"]
    fmap_df = (
        fmap_df.groupby("acquisition.id")
        .filter(
            lambda x: (
                not x[x["file.type"] == "dicom"].empty
                and not x[x["file.type"] == "dicom"]["file.classification.Intent"]
                .isnull()
                .values.any()
                and "Fieldmap"
                in x[x["file.type"] == "dicom"]["file.classification.Intent"].iloc[0]
            )
        )
        .reset_index(drop=True)
    )
    fmap_df["reproin_dict"] = fmap_df["reproin"].apply(lambda x: parse_series_spec(x))

    if fmap_df.empty:
        return df

    # Determine fieldmap file types (magnitude, phase, etc.)
    type_mapping = fmap_df.groupby("acquisition.id")[fmap_df.columns].apply(
        lambda x: classify_fmap_acq(
            x[x["file.type"] == "dicom"]["file.info.header.dicom.ImageType"].iloc[0]
        )
    )
    if type_mapping.empty:
        err_msg = (
            "Subject %s doesn't contain any files with proper ImageType." % subject
        )
        log.error(err_msg)
        df.loc[:, "error"] = err_msg
        return df
    fmap_df["fmap_file_type"] = fmap_df["acquisition.id"].map(type_mapping)
    fmap_df = fmap_df[
        fmap_df["file.type"].apply(lambda x: x in ("nifti", "bval", "bvec"))
    ]
    fmap_df["magnitude_type"] = fmap_df["file.name"].apply(
        lambda x: 1
        if ("_e1.nii" in x or "_magnitude1.nii" in x)
        else (2 if ("_e2.nii" in x or "_magnitude2.nii" in x) else None)
    )

    # Group acquisitions by timestamp
    fmap_df = timestamp_group(fmap_df, FMAP_DELTA)
    fmap_df["fmap_acq"] = None
    fmap_df = fmap_df.set_index("file.file_id")

    # Determine bids fmap case and assign reproin labels
    for name, group in fmap_df.groupby("group"):
        if (
            group.iloc[0]["acquisition.timestamp"]
            - group.iloc[-1]["acquisition.timestamp"]
        ) > FMAP_DELTA:
            err_msg = "Fieldmap group in subject %s spans more than %s" % (
                    fmap_df.iloc[0]['subject.label'],
                    FMAP_DELTA,
                )
            log.error(err_msg)
            group.loc[:, "error"] = err_msg
            return group
        elif (
            group["reproin_dict"]
            .apply(
                lambda x: "datatype_suffix" in x
                and x["datatype_suffix"] in FMAP_SUFFIXES
            )
            .all()
        ):
            continue

        niftis = group[group["file.type"] == "nifti"]
        magnitude = niftis[niftis["fmap_file_type"] == "magnitude"]
        phase = niftis[niftis["fmap_file_type"] == "phase"]
        diffusion = niftis[niftis["fmap_file_type"] == "diffusion"]
        fmri = niftis[niftis["fmap_file_type"] == "fmri"]
        fmap_file_types = {
            "magnitude": magnitude,
            "phase": phase,
            "diffusion": diffusion,
            "fmri": fmri,
        }

        if magnitude.empty and diffusion.empty and fmri.empty:
            err_msg = "Subject %s has fieldmaps with no magnitude, diffusion or fmri files." % subject
            log.error(err_msg)
            group.loc[:, "error"] = err_msg
        elif sum((magnitude.empty, diffusion.empty, fmri.empty)) < 2:
            err_msg = "Subject %s has fieldmap(s) with more than one type of files." % subject
            log.error(err_msg)
            group.loc[:, "error"] = err_msg
        elif not magnitude.empty:
            if len(magnitude) > 2:
                err_msg = "Subject %s has fieldmap(s) with > 2 magnitude files." % subject
                log.error(err_msg)
                group.loc[:, "error"] = err_msg
            elif len(phase) > 2:
                err_msg = "Subject %s has fieldmap(s) with > 2 phase files." % subject
                log.error(err_msg)
                group.loc[:, "error"] = err_msg
            elif len(phase) in (1, 2):
                pattern = re.compile(
                    r"(?i)(spinecho|(^|[^a-zA-Z0-9])se($|[^a-zA-Z0-9])|(?-i:AP|PA))"
                )
                magnitude.loc[:, "spinecho"] = magnitude["acquisition.label"].str.match(
                    pattern
                )
                num_spinecho = len(magnitude[magnitude["spinecho"]])

                if num_spinecho not in (0, len(magnitude)):
                    err_msg = "Subject %s has magnitude(s) that are both spinecho and not." % subject
                    log.error(err_msg)
                    group.loc[:, "error"] = err_msg
                elif num_spinecho == len(magnitude) and (
                    len(magnitude[magnitude["acquisition.label"].str.contains("AP")])
                    == len(magnitude[magnitude["acquisition.label"].str.contains("PA")])
                ):
                    # Case 4
                    group = get_case_four(group, fmap_file_types, subject)
                elif magnitude["magnitude_type"].isna().any():
                    err_msg = "Subject %s has magnitude(s) not matching expected filenames." % subject
                    log.error(err_msg)
                    group.loc[:, "error"] = err_msg
                elif len(phase) == 1:
                    # Case 1
                    group["fmap_acq"] = "phasediff"
                else:
                    # Case 2
                    group["fmap_acq"] = "twophase"
                    group["reproin"] = "fmap-gre"
            elif len(phase) == 0:
                # Case 4
                group = get_case_four(group, fmap_file_types, subject)
        else:
            # Case 4
            group = get_case_four(group, fmap_file_types, subject)

        fmap_df.update(group)

    if fmap_df["fmap_acq"].nunique() > 1:
        fmap_df["reproin"] += "_acq-" + fmap_df["fmap_acq"]
    fmap_df["reproin"] = (
        fmap_df.groupby("fmap_acq")[fmap_df.columns]
        .apply(add_fmap_runs)
        .reset_index(level=0, drop=True)["reproin"]
    )

    fmap_df = fmap_df.set_index("acquisition.id")
    df = df.set_index("acquisition.id")
    fmap_mapping = fmap_df.groupby(fmap_df.index).first()
    fmap_mapping = fmap_mapping["reproin"]
    fmap_mapping = fmap_mapping.dropna()
    df["reproin"] = df["reproin"].where(
        ~df.index.isin(fmap_mapping.index), df.index.map(fmap_mapping)
    )

    return df.reset_index()


def add_rec(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    rec_df = df[df["reproin"].str.startswith(("anat", "asl"))]

    if rec_df.empty:
        return df

    rec_df = (
        rec_df.groupby("reproin")[rec_df.columns]
        .apply(lambda x: timestamp_group(x, REC_DELTA))
        .reset_index(drop=True)
    )

    subject = df["subject.label"].iloc[0]

    for i, group in rec_df.groupby(["reproin", "group"]):
        if len(group) == 1:
            continue

        REC_LIST.append(group)

        if len(group) > 2:
            err_msg = "Subject %s has more than 2 acquisitions in the same rec group." % subject
            log.error(err_msg)
            df.loc[df["acquisition.id"].isin(group["acquisition.id"]), "error"] = err_msg
        elif len(group) == 2:
            group["norm"] = group["file.info.header.dicom.ImageType"].apply(
                lambda x: x is not None and "NORM" in x
            )

            if len(group[group["norm"]]) != 1:
                err_msg = "Subject %s has a rec group with 0 or 2 NORM acquisitions" % subject
                log.error(err_msg)
                df.loc[df["acquisition.id"].isin(group["acquisition.id"]), "error"] = err_msg
            else:
                acq_id_update = group[group["norm"]]["acquisition.id"].iloc[0]
                df.loc[df["acquisition.id"] == acq_id_update, "reproin"] += "_rec-norm"

    return df


def add_run(df: pd.DataFrame) -> pd.DataFrame:
    skip_df = df[df["acquisition.label"].str.contains("sbref", case=False)]
    run_df = df[~df["file.file_id"].isin(skip_df["file.file_id"])]

    # TODO: Handle recs properly (make sure rec pairs get same run number)
    if run_df.shape[0] > 1:
        run_df = get_runs(run_df)
        if run_df["reproin"].iloc[[0]].str.startswith("fmap").item():
            run_df["reproin"] = run_df["reproin"] + "_" + run_df["run"]
        else:
            run_df["reproin"] = run_df["reproin"] + "_run-" + run_df["run"]
    return pd.concat([run_df, skip_df])


def add_sbref(df: pd.DataFrame()):
    df = df.copy()
    sbref_df = df[df["acquisition.label"].str.contains("sbref", case=False)]

    if not sbref_df.empty:
        func_dwi = df[df["reproin"].str.startswith(("func", "dwi"))]
        func_dwi = func_dwi[
            ~func_dwi["acquisition.label"].str.contains("sbref", case=False)
        ]

        for i, sbref in sbref_df.iterrows():
            func_dwi["timedelta"] = (
                sbref["acquisition.timestamp"] - func_dwi["acquisition.timestamp"]
            )
            match = func_dwi[
                (func_dwi["timedelta"] <= SBREF_DELTA)
                & (func_dwi["timedelta"] >= timedelta(seconds=-3))
            ]

            if len(match) == 1:
                df.loc[i, "sbref"] = match["acquisition.id"].iloc[0]
                df.loc[i, "reproin"] = (
                    match["reproin"].iloc[0].removesuffix("SBRef") + "_SBRef"
                )
            elif len(match) > 1:
                err_msg = (
                    "%s in subject %s should match exactly 1 bold image, but it matched 2 %s" % (
                        sbref["acquisition.label"],
                        sbref["subject.label"],
                        match["acquisition.label"],
                    )
                )
                log.error(err_msg)
                df.loc[i, "error"] = err_msg 
            else:
                err_msg = "%s in subject %s didn't match any bold images." % (
                        sbref["acquisition.label"],
                        sbref["subject.label"],
                    )
                log.error(err_msg)
                df.loc[i, "error"] = err_msg 

    return df
