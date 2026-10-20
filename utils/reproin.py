#!/usr/bin/env python

import re
import pandas as pd
from heudiconv.heuristics.reproin import parse_series_spec
import logging
from datetime import timedelta
from utils.constants import (
    ALLOWED_DATATYPES,
    ANAT_SUFFIX_CASE,
    ACQ_PARAMS,
    OPTIONAL_ACQ_PARAMS,
    ACQ_PARAMS_FALLBACK,
    SBREF_DELTA,
    FMAP_DELTA,
    REC_DELTA,
    FMAP_SUFFIXES,
    IGNORE_SCANS,
    DWI_SUFFIXES_IGNORE,
)

log = logging.getLogger(__name__)


_NIFTI_DERIV_RE = re.compile(r"(?:_ph|_e[2-9]|_imaginary|_real|_magnitude[2-9])\.nii")

# Shape of each dict in acq_df["nifti_files"]; keys = file_df columns, values = dict keys
_NIFTI_FILE_FIELDS = {
    "file.file_id": "file_id",
    "file.name": "name",
    "file.tags": "tags",
    "file.info.SeriesDescription": "SeriesDescription",
    "ImageType": "ImageType",
    "file.info.EchoTime": "EchoTime",
    "magnitude_type": "magnitude_type",
}


def _magnitude_type(name: str):
    if "_e1.nii" in name or "_magnitude1.nii" in name:
        return 1
    if "_e2.nii" in name or "_magnitude2.nii" in name:
        return 2
    return None


def _combine_image_type(image_type, image_type_text):
    """Merges ImageType and ImageTypeText into one list, dropping duplicates.
    Returns image_type unchanged if neither is a list."""
    combined = []
    for val in (image_type, image_type_text):
        if isinstance(val, list):
            combined += [x for x in val if x not in combined]
    return combined or image_type


def build_acq_df(file_df: pd.DataFrame) -> pd.DataFrame:
    """Reduces a per-file dataframe to one row per acquisition.

    Uses the canonical nifti as the base row (first non-derivative nifti, or
    first alphabetically if all are derivatives). Falls back to the first
    dicom row when an acquisition has no niftis. Attaches a nifti_files
    list-of-dicts column and has_bval / has_bvec presence flags so downstream
    stages never need to fall back to the raw file_df.
    """
    # Per-acquisition lookups
    types_by_acq = file_df.groupby("acquisition.id")["file.type"].apply(set)
    dicom_rows = (
        file_df[file_df["file.type"] == "dicom"].groupby("acquisition.id").first()
    )

    niftis = file_df[file_df["file.type"] == "nifti"].copy()
    niftis["ImageType"] = [
        _combine_image_type(*vals)
        for vals in zip(
            niftis["file.info.ImageType"], niftis["file.info.ImageTypeText"]
        )
    ]
    niftis["magnitude_type"] = niftis["file.name"].apply(_magnitude_type)
    niftis["_is_deriv"] = niftis["file.name"].str.contains(
        _NIFTI_DERIV_RE, regex=True, na=False
    )

    # Build nifti_files list-of-dicts in a single pass
    records = (
        niftis[list(_NIFTI_FILE_FIELDS)]
        .rename(columns=_NIFTI_FILE_FIELDS)
        .to_dict(orient="records")
    )
    nifti_files: dict = {}
    for acq_id, rec in zip(niftis["acquisition.id"], records):
        nifti_files.setdefault(acq_id, []).append(rec)

    # Pick canonical row per acquisition (non-deriv preferred, then alphabetical)
    canonical_niftis = (
        niftis.sort_values(["_is_deriv", "file.name"])
        .drop_duplicates("acquisition.id", keep="first")
        .set_index("acquisition.id")
    )
    no_nifti_rows = dicom_rows.loc[
        dicom_rows.index.difference(canonical_niftis.index)
    ].copy()
    no_nifti_rows["file.info.SeriesDescription"] = no_nifti_rows[
        "file.info.header.dicom.SeriesDescription"
    ]
    no_nifti_rows["ImageType"] = no_nifti_rows["file.info.header.dicom.ImageType"]
    acq_df = pd.concat([canonical_niftis, no_nifti_rows])

    # Backfill fields that only exist on dicom
    for field in ("AcquisitionTime", "AcquisitionDate", "AcquisitionDateTime"):
        col = f"file.info.header.dicom.{field}"
        acq_df[col] = acq_df.index.map(dicom_rows[col])

    # Attach derived columns
    acq_df["nifti_files"] = acq_df.index.map(lambda acq_id: nifti_files.get(acq_id, []))
    acq_df["has_bval"] = types_by_acq.apply(lambda s: "bval" in s).reindex(
        acq_df.index, fill_value=False
    )
    acq_df["has_bvec"] = types_by_acq.apply(lambda s: "bvec" in s).reindex(
        acq_df.index, fill_value=False
    )

    # Drop file-level + internal columns
    drop_cols = (
        "file.file_id",
        "file.type",
        "file.name",
        "file.tags",
        "file.created",
        "file.info.ImageType",
        "file.info.ImageTypeText",
        "_is_deriv",
        "magnitude_type",
    )
    acq_df = acq_df.drop(columns=[c for c in drop_cols if c in acq_df.columns])

    return acq_df.reset_index()


def check_pydeface(acq_df: pd.DataFrame) -> pd.DataFrame:
    """Removes subjects with any Structural T1 or T2 acquisitions containing a
    nifti without the 'pydeface' tag, unless the acquisition's reproin ends with
    _ignore-BIDS."""
    acq_df = acq_df.copy()

    is_struct = acq_df["file.classification.Intent"].apply(
        lambda x: type(x) is list and "Structural" in x
    ) & acq_df["file.classification.Measurement"].apply(
        lambda x: type(x) is list and ("T1" in x or "T2" in x)
    )
    no_pydeface_df = acq_df[
        is_struct
        & acq_df["nifti_files"].apply(
            lambda files: any("pydeface" not in f["tags"] for f in files)
        )
        & ~acq_df["reproin"].str.endswith("_ignore-BIDS", na=False)
    ]

    acq_df = acq_df[~acq_df["subject.label"].isin(no_pydeface_df["subject.label"])]
    if not no_pydeface_df.empty:
        log.warning(
            "The following subjects contain T1/T2 acquisition(s) without the 'pydeface' tag "
            "and will be skipped: \n%s"
            % no_pydeface_df[["subject.label", "acquisition.label"]].to_string()
        )

    return acq_df


def check_nifti_presence(acq_df: pd.DataFrame) -> pd.DataFrame:
    """Removes subjects with any MR acquisitions missing niftis, unless the
    acquisition's reproin ends with _ignore-BIDS."""
    acq_df = acq_df.copy()

    no_nii_df = acq_df[
        (acq_df["file.modality"] == "MR")
        & (acq_df["nifti_files"].apply(len) == 0)
        & ~acq_df["reproin"].str.endswith("_ignore-BIDS", na=False)
    ]

    acq_df = acq_df[~acq_df["subject.label"].isin(no_nii_df["subject.label"])]
    if not no_nii_df.empty:
        log.warning(
            "The following subjects contain acquisition(s) without nifti(s): "
            "\n%s" % no_nii_df[["subject.label", "acquisition.label"]].to_string()
        )

    return acq_df


def check_asl_and_spectroscopy(file_df: pd.DataFrame) -> pd.DataFrame:
    """Removes subjects that have any ASL/perfusion or spectroscopy acquisitions"""
    file_df = file_df.copy()

    excluded_df = file_df.groupby("subject.label").filter(
        lambda x: (
            x["reproin"].str.startswith("asl", na=False).any()
            or x["file.classification.Intent"]
            .apply(lambda i: isinstance(i, list) and "Spectroscopy" in i)
            .any()
            or x["acquisition.label"].str.startswith("csf_flow", na=False).any()
        )
    )

    file_df = file_df[
        ~file_df["subject.label"].isin(excluded_df["subject.label"].unique())
    ]

    if not excluded_df.empty:
        log.warning(
            "The following subjects have ASL/perfusion or spectroscopy acquisitions and will be skipped: "
            "\n%s" % excluded_df["subject.label"].unique()
        )

    return file_df


def check_anat_presence(file_df: pd.DataFrame) -> pd.DataFrame:
    """Removes sessions that have no anatomical acquisitions"""
    file_df = file_df.copy()

    no_anat_df = file_df.groupby("session.id").filter(
        lambda x: not x["reproin"].str.startswith("anat", na=False).any()
    )

    file_df = file_df[~file_df["session.id"].isin(no_anat_df["session.id"].unique())]

    if not no_anat_df.empty:
        log.warning(
            "The following sessions have no anatomical acquisitions and will be skipped: "
            "\n%s" % no_anat_df["session.id"].unique()
        )

    return file_df


def reproin_filter(row: pd.Series) -> str:
    label = row["file.info.SeriesDescription"]
    if not isinstance(label, str):
        label = row["acquisition.label"]
    label_lower = label.lower()

    # Utilize acq label if already in reproin format. Otherwise, create a reproin label
    label_re = re.sub(r"_\d$", "", label_lower)
    label_ignore = label + "_ignore-BIDS"
    if label_re in IGNORE_SCANS:
        return label_ignore

    image_type = row["ImageType"]
    # to-do: add exception for mp2rage
    # if image_type and "DERIVED" in image_type:
    #    return label_ignore
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
        if validator.get("datatype_suffix") == "scout":
            return label_ignore
        if validator["datatype"] == "func":
            task = validator.get("task", "")
            if task and (
                "rest" in task.lower()
                or "rsfmri" in task.lower()
                or task.lower() == "rs"
            ):
                validator["task"] = "rest"
            else:
                return label_ignore

        new_label = validator["datatype"]
        if validator.get("datatype_suffix") and validator["datatype_suffix"] != "None":
            suffix = validator["datatype_suffix"]
            suffix = ANAT_SUFFIX_CASE.get(suffix, suffix)
            new_label += "-" + suffix
        if validator.get("task") and validator["task"] != "None":
            new_label += "_task-" + validator["task"]
        if validator.get("acq") and validator["acq"] != "None":
            new_label += "_acq-" + validator["acq"]
        return new_label
    else:
        measurement = row["file.classification.Measurement"]
        intent = row["file.classification.Intent"]
        features = row["file.classification.Features"]

        if label_lower.startswith("gobrain_"):
            return label_ignore
        elif label_re.endswith(("physiolog", "setter", "tensor")):
            return label_ignore
        elif intent and ("Localizer" in intent or "Spectroscopy" in intent):
            return label_ignore
        elif label_lower.startswith("se_distortion"):
            return "fmap"
        elif measurement and "Perfusion" in measurement:
            return "asl"
        elif intent and "Structural" in intent:
            if features and "FLAIR" in features:
                return "anat-FLAIR"
            elif measurement:
                if "T1" in measurement:
                    return "anat-T1w"
                elif "T2" in measurement:
                    if "hippocampus" in label_lower:
                        return "anat-T2w_acq-hippocampus"
                    else:
                        return "anat-T2w"
                elif "Diffusion" in measurement:
                    if not label_re.endswith(DWI_SUFFIXES_IGNORE):
                        return "dwi"
        elif (
            intent
            and "Functional" in intent
            and label_lower
            and ("rest" in label_lower or "rsfmri" in label_lower)
        ):
            return "func_task-rest"
        elif intent and "Fieldmap" in intent:
            return "fmap"

    return label_ignore


def rm_acq_labels(df: pd.DataFrame()) -> pd.DataFrame():
    """Removes certain blacklisted acq labels from reproin columns when appropriate."""
    df = df.copy()

    df["reproin_no_acq"] = df["reproin"].apply(
        lambda x: re.sub("_acq-.*?(_|$)", r"\1", x)
    )
    if df["reproin"].equals(df["reproin_no_acq"]):
        return df

    df["validator"] = df["reproin"].apply(lambda x: parse_series_spec(x.lower()))
    df["acq"] = df["validator"].apply(lambda x: x.get("acq"))
    df["datatype_suffix"] = df["validator"].apply(lambda x: x.get("datatype_suffix"))
    group_df = df.groupby("reproin_no_acq")[df.columns]

    for _, group in group_df:
        mm_matches = group[group["acq"].str.match(r"^\d+[.p]?\d+mm$", na=False)]
        if not mm_matches.empty:
            if (len(mm_matches) == len(group)) and (mm_matches["acq"].nunique() == 1):
                df.loc[mm_matches.index, "reproin"] = mm_matches["reproin_no_acq"]
            continue

        mprage_matches = group[
            (group["acq"] == "mprage") & (group["datatype_suffix"] == "t1w")
        ]
        space_matches = group[
            (group["acq"] == "space") & (group["datatype_suffix"] == "t2w")
        ]

        for matches in (mprage_matches, space_matches):
            if not matches.empty:
                if len(matches) == len(group):
                    df.loc[matches.index, "reproin"] = matches["reproin_no_acq"]

    return df


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
    run_df = df.sort_values(by="acq_datetime")
    n_digits = max(2, len(str(run_df.shape[0])))
    padded_list = [str(n).zfill(n_digits) for n in range(1, run_df.shape[0] + 1)]
    run_df.insert(run_df.shape[1], "run", padded_list)
    return run_df


def add_fmap_runs(df: pd.DataFrame()) -> pd.DataFrame():
    df = df.copy()
    group_df = df[["acq_datetime", "group"]].drop_duplicates(
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
    df = df.sort_values(by="acq_datetime")
    df["time_diff"] = df["acq_datetime"].diff().fillna(pd.Timedelta(seconds=0))
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

    non_phase["direction"] = non_phase[
        "file.info.header.dicom.SeriesDescription"
    ].str.extract(r".*(LR|RL|AP|PA)(?=[^a-zA-Z0-9]|$)")[0]

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
    if df.empty:
        return df

    has_fmap = df["file.classification.Intent"].apply(
        lambda x: isinstance(x, list) and "Fieldmap" in x
    ) & (df["file.modality"] == "MR")

    if not has_fmap.any():
        return df.reset_index(drop=True)

    subject = df["subject.label"].iloc[0]
    err_msg = "Subject %s has fieldmap(s) and will be skipped." % subject
    log.warning(err_msg)
    df.loc[:, "error"] = err_msg
    return df.reset_index(drop=True)

    # ------------------------------------------------------------------
    # NEEDS PORT TO acq_df: everything below this early return is
    # preserved as a reference for the file-level fmap classification
    # logic. It still assumes the old data model (file.type / file.name /
    # bval / bvec rows) and will fail at runtime under acq_df. Rewrite
    # before removing the early return above.
    # ------------------------------------------------------------------
    fmap_df = df.copy()
    fmap_related_df = df[df["reproin"].str.startswith(("asl", "func", "dwi"), na=False)]
    if fmap_related_df.empty:
        fmap_df["reproin"] += "_ignore-BIDS"
        df.update(fmap_df)
        return df.reset_index(drop=True)

    fmap_df["reproin_dict"] = fmap_df["reproin"].apply(lambda x: parse_series_spec(x))

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
        return df.reset_index(drop=True)
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
            group.iloc[0]["acq_datetime"] - group.iloc[-1]["acq_datetime"]
        ) > FMAP_DELTA:
            err_msg = "Fieldmap group in subject %s spans more than %s" % (
                fmap_df.iloc[0]["subject.label"],
                FMAP_DELTA,
            )
            log.error(err_msg)
            group.loc[:, "error"] = err_msg
            return group.reset_index(drop=True)
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
            err_msg = (
                "Subject %s has fieldmaps with no magnitude, diffusion or fmri files."
                % subject
            )
            log.error(err_msg)
            group.loc[:, "error"] = err_msg
        elif sum((magnitude.empty, diffusion.empty, fmri.empty)) < 2:
            err_msg = (
                "Subject %s has fieldmap(s) with more than one type of files." % subject
            )
            log.error(err_msg)
            group.loc[:, "error"] = err_msg
        elif not magnitude.empty:
            if len(magnitude) > 2:
                err_msg = (
                    "Subject %s has fieldmap(s) with > 2 magnitude files." % subject
                )
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
                magnitude.loc[:, "spinecho"] = magnitude[
                    "file.info.header.dicom.SeriesDescription"
                ].str.match(pattern)
                num_spinecho = len(magnitude[magnitude["spinecho"]])

                if num_spinecho not in (0, len(magnitude)):
                    err_msg = (
                        "Subject %s has magnitude(s) that are both spinecho and not."
                        % subject
                    )
                    log.error(err_msg)
                    group.loc[:, "error"] = err_msg
                elif num_spinecho == len(magnitude) and (
                    len(
                        magnitude[
                            magnitude[
                                "file.info.header.dicom.SeriesDescription"
                            ].str.contains("AP")
                        ]
                    )
                    == len(
                        magnitude[
                            magnitude[
                                "file.info.header.dicom.SeriesDescription"
                            ].str.contains("PA")
                        ]
                    )
                ):
                    # Case 4
                    group = get_case_four(group, fmap_file_types, subject)
                elif magnitude["magnitude_type"].isna().any():
                    err_msg = (
                        "Subject %s has magnitude(s) not matching expected filenames."
                        % subject
                    )
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

    # Note: need to handle ASL/Angio (eg sub-G5GRGM)
    for i, group in rec_df.groupby(["reproin", "group"]):
        if len(group) == 1:
            continue

        reproin = group["reproin"].iloc[0]

        rec_type_dict = {
            "derived": {"DERIVED"},
            "norm": {"NORM"},
            "dis": {"DIS3D", "DIS2D"},
        }
        for label, keywords in rec_type_dict.items():
            group[label] = group["ImageType"].apply(
                lambda x: isinstance(x, list) and bool(keywords & set(x))
            )
        memprage_type_dict = {
            "memprage": "memprage|mpr",
            "rms": "rms",
        }
        for label, regex in memprage_type_dict.items():
            group[label] = group["acquisition.label"].str.contains(regex, case=False)

        ignore = pd.Series(False, index=group.index)

        if group["derived"].any():
            if not reproin.startswith("anat-T1w"):
                err_msg = "Subject %s's %s group contains DERIVED but is not T1w" % (
                    subject,
                    reproin,
                )
                log.error(err_msg)
                df.loc[df["acquisition.id"].isin(group["acquisition.id"]), "error"] = (
                    err_msg
                )

            ignore |= group["derived"]
        if group["norm"].any():
            ignore |= ~group["norm"]
        if group["dis"].any():
            if (~ignore).sum() > 2:
                err_msg = (
                    "Subject %s has more than 2 acquisitions in the same DIS rec group."
                    % subject
                )
                log.error(err_msg)
                df.loc[df["acquisition.id"].isin(group["acquisition.id"]), "error"] = (
                    err_msg
                )
            ignore |= ~group["dis"]
        if (group["memprage"] & group["rms"]).any():
            ignore |= group["memprage"] & ~group["rms"]

        ignore_ids = group.loc[ignore, "acquisition.id"]
        df.loc[df["acquisition.id"].isin(ignore_ids), "reproin"] += "_ignore-BIDS"

        if (~ignore).sum() != 1:
            err_msg = (
                "Subject %s does not have exactly 1 included acq in its %s rec group: %s"
                % (subject, reproin, group["acquisition.label"].tolist())
            )
            log.error(err_msg)
            df.loc[df["acquisition.id"].isin(group["acquisition.id"]), "error"] = (
                err_msg
            )

    return df


def _is_empty(val) -> bool:
    if isinstance(val, list):
        return len(val) == 0
    return pd.isna(val)


def _describe_missing_params(
    labels: pd.Series, empty_mask: pd.DataFrame, missing_cols: list
) -> str:
    """Lists each missing parameter with the acquisitions missing it, e.g.
    'EchoTime (rest_1, rest_2)'."""
    return "; ".join(
        "%s (%s)"
        % (col.removeprefix("file.info."), ", ".join(labels[empty_mask[col].values]))
        for col in missing_cols
    )


def _group_mismatched_runs(labels: pd.Series, params_df: pd.DataFrame) -> str:
    """Groups runs with identical parameters, e.g. '[T1w_1] vs [T1w_2, T1w_3]'."""
    groups = {}
    for key, label in zip(params_df.astype(str).itertuples(index=False), labels):
        groups.setdefault(tuple(key), []).append(label)
    return " vs ".join("[%s]" % ", ".join(lbls) for lbls in groups.values())


def get_acq_params(row: pd.Series) -> pd.Series:
    """Returns ACQ_PARAMS values for a row, falling back through any
    ACQ_PARAMS_FALLBACK paths in order when the top-level field is empty.
    Single-element lists are unwrapped so list-typed fallbacks (e.g.
    InversionTimes) compare equal to scalar top-level fields."""
    result = {}
    for param in ACQ_PARAMS:
        val = row.get(param)
        if _is_empty(val):
            for fallback in ACQ_PARAMS_FALLBACK.get(param, ()):
                val = row.get(fallback)
                if not _is_empty(val):
                    break
        if isinstance(val, list) and len(val) == 1:
            val = val[0]
        result[param] = val
    return pd.Series(result)


def add_run(df: pd.DataFrame) -> pd.DataFrame:
    skip_df = df[df["file.info.SeriesDescription"].str.contains("sbref", case=False)]
    run_df = df[~df["acquisition.id"].isin(skip_df["acquisition.id"])]

    # TODO: Handle recs properly (make sure rec pairs get same run number)
    if run_df.shape[0] > 1:
        run_df = get_runs(run_df)

        if run_df["reproin"].iloc[[0]].str.startswith("fmap").item():
            run_df["reproin"] = run_df["reproin"] + "_" + run_df["run"]
        else:
            unique_df = run_df.apply(get_acq_params, axis=1)

            subject = df["subject.label"].iloc[0]
            reproin = df["reproin"].iloc[0]

            unique_req_df = unique_df[
                [c for c in unique_df.columns if c not in OPTIONAL_ACQ_PARAMS]
            ]
            empty_mask = unique_req_df.map(_is_empty)
            missing_cols = empty_mask.columns[empty_mask.any()].tolist()
            if missing_cols:
                err_msg = "Subject %s is missing acquisition parameters for %s: %s" % (
                    subject,
                    reproin,
                    _describe_missing_params(
                        run_df["acquisition.label"], empty_mask, missing_cols
                    ),
                )
                log.error(err_msg)
                df.loc[:, "error"] = err_msg
                return df

            unique_df = unique_df.apply(
                lambda col: col.map(lambda x: tuple(x) if isinstance(x, list) else x)
            )
            if len(unique_df.drop_duplicates()) > 1:
                err_msg = (
                    "Subject %s has runs with mismatched acquisition parameters for %s: %s"
                    % (
                        subject,
                        reproin,
                        _group_mismatched_runs(run_df["acquisition.label"], unique_df),
                    )
                )
                log.error(err_msg)
                df.loc[:, "error"] = err_msg
                return df

            run_df["reproin"] = run_df["reproin"] + "_run-" + run_df["run"]

    return pd.concat([run_df, skip_df])


def add_sbref(df: pd.DataFrame()):
    df = df.copy()
    sbref_df = df[df["file.info.SeriesDescription"].str.contains("sbref", case=False)]

    if not sbref_df.empty:
        func_dwi = df[df["reproin"].str.startswith(("func", "dwi"))]
        func_dwi = func_dwi[
            ~func_dwi["file.info.SeriesDescription"].str.contains("sbref", case=False)
        ]

        if (~func_dwi["acq_datetime"].isna()).all() and (
            ~sbref_df["acq_datetime"].isna()
        ).all():
            use_series = False
        elif (~func_dwi["file.info.SeriesNumber"].isna()).all() and (
            ~sbref_df["file.info.SeriesNumber"].isna()
        ).all():
            use_series = True
        else:
            err_msg = (
                "Subject %s contains files relevant to sbref matching that don't contain "
                "acq_datetime or SeriesNumber."
            ) % df.iloc[0]["subject.label"]
            log.error(err_msg)
            df.loc[:, "error"] = err_msg
            return df

        for i, sbref in sbref_df.iterrows():
            match = func_dwi[
                sbref["file.info.ProtocolName"] == func_dwi["file.info.ProtocolName"]
            ].copy()

            if use_series:
                match = match[
                    (
                        match["file.info.SeriesNumber"]
                        - sbref["file.info.SeriesNumber"]
                        == 1
                    )
                ]
            else:
                match["timedelta"] = sbref["acq_datetime"] - match["acq_datetime"]
                match = match[
                    (match["timedelta"] <= SBREF_DELTA)
                    & (match["timedelta"] >= -SBREF_DELTA)
                ]

            if len(match) == 1:
                df.loc[i, "sbref"] = match["acquisition.id"].iloc[0]
                df.loc[i, "reproin"] = (
                    match["reproin"].iloc[0].removesuffix("SBRef") + "_SBRef"
                )
            elif len(match) > 1:
                err_msg = (
                    "%s in subject %s should match exactly 1 bold image, but it matched multiple %s"
                    % (
                        sbref["acquisition.label"],
                        sbref["subject.label"],
                        match["acquisition.label"].tolist(),
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


def check_duplicate_reproin(df: pd.DataFrame()) -> pd.DataFrame():
    """Returns df if there are not duplicates in the 'reproin' column"""
    df = df.copy()
    duplicated = df[df["reproin"].duplicated()]

    if not duplicated.empty:
        subject_label = df["subject.label"].iloc[0]
        err_msg = "Subject %s contains a session with duplicate reproins: %s" % (
            subject_label,
            duplicated,
        )
        log.error(err_msg)
        return None

    return df
