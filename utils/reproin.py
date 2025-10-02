#!/usr/bin/env python

import re
import pandas as pd
from heudiconv.heuristics.reproin import parse_series_spec
import sys
import logging
from datetime import timedelta
from utils.constants import (
    ALLOWED_DATATYPES,
    SBREF_DELTA,
    FMAP_DELTA,
    REC_DELTA,
    FMAP_SUFFIXES,
    IGNORE_SCANS,
    DWI_SUFFIXES_IGNORE

)

log = logging.getLogger(__name__)

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

def check_pydeface(file_df: pd.DataFrame) -> pd.DataFrame:
    """Removes subjects if not all structural T1 and T2 niftis have 'pydeface' tag"""
    no_pydeface_df = file_df.copy()
    no_pydeface_df = no_pydeface_df.groupby('acquisition.id').filter(pydeface_filter)
    file_df = file_df[~file_df['subject.label'].isin(no_pydeface_df['subject.label'])]

    if not no_pydeface_df.empty:
        log.warning("The following subjects contain anatomical nifti(s) without the 'pydeface' tag: "
                    f"\n{no_pydeface_df['subject.label'].unique()}"
        )
    
    return file_df

def check_nifti_presence(file_df: pd.DataFrame) -> pd.DataFrame:
    """Removes subjects that have any acquisitions without niftis"""
    no_nii_df = file_df.copy()
    no_nii_df = no_nii_df.groupby('acquisition.id').filter(
        lambda x: x[x['reproin'].str.endswith('_ignore-BIDS')].empty
    ).reset_index(drop=True)
    no_nii_df = no_nii_df[no_nii_df['file.modality'] == 'MR']
    no_nii_df = no_nii_df.groupby('acquisition.id').filter(
            lambda x: not x['file.type'].str.contains('nifti').any()
    ).reset_index(drop=True)

    file_df = file_df[~file_df['subject.label'].isin(no_nii_df['subject.label'])]
    if not no_nii_df.empty:
        log.warning("The following subjects contain acquisition(s) without nifti(s): "
                    f"\n{no_nii_df['subject.label'].unique()}"
        )
    
    return file_df


def reproin_filter(row: pd.Series) -> str:
    label = row['acquisition.label']
    if row['file.type'] != 'dicom':
        return label
    elif '_ignore-BIDS' in label:
        return label
    
    # Utilize acq label if already in reproin format. Otherwise, create a reproin label
    label_re = re.sub(r'_\d$', '', label)
    label_ignore = label + '_ignore-BIDS'
    if label_re in IGNORE_SCANS:
        return label_ignore
    image_type = row['file.info.header.dicom_array.ImageType.0']
    # to-do: add exception for mp2rage
    if image_type and 'derived' in image_type:
        return label_ignore
    validator = parse_series_spec(label)
    if validator and validator['datatype'] == 'anat' and 'datatype_suffix' not in validator:
        validator = None

    if validator:
        if validator['datatype'] not in ALLOWED_DATATYPES:
            return label_ignore
        elif validator['datatype'] == 'func':
            if 'task' in validator and ('rest' in validator['task'] or validator['task'] == 'rs'):
                return re.sub('task-.*?(_|$)', r'task-rest\1', label_re)
            else:
                return label_ignore
        return label_re
    else:
        measurement = row['file.classification.Measurement']
        intent = row['file.classification.Intent']
        features = row['file.classification.Features']

        if label.startswith('GOBRAIN_'):
            return label_ignore
        elif label == 't2_tse_tra_hi-res_hippocampus':
            return 'anat-T2w_acq-hippocampus'
        elif label_re.endswith('PhysioLog'):
            return label_ignore
        elif intent and ('Localizer' in intent or 'Spectroscopy' in intent):
            return label_ignore
        elif measurement and 'Perfusion' in measurement:
            return 'asl'
        elif measurement and intent and 'Structural' in intent:
            if 'T1' in measurement:
                return 'anat-T1w'
            elif 'T2' in measurement:
                if features and 'FLAIR' in features:
                    return 'anat-FLAIR'
                else:
                    return 'anat-T2w'   
            elif 'Diffusion' in measurement:
                if label_re.endswith(DWI_SUFFIXES_IGNORE):
                    return label_ignore
                else:
                    return 'dwi'
        elif (intent and 'Functional' in intent and label
            and ('rest' in label.lower() or 'rsfmri' in label.lower())):
            return 'func_task-rest'
        elif intent and 'Fieldmap' in intent:
            return 'fmap'
        else:
            return label_ignore


def classify_fmap_acq(image_type_list: list) -> str:
    """Determines the type of fieldmap file (phase, magnitude or diffusion) based on
    the presence of corresponding strings in the ImageType dict."""

    if 'M' in image_type_list:
        return 'magnitude'
    elif 'P' in image_type_list:
        return 'phase'
    elif 'DIFFUSION' in image_type_list:
        return 'diffusion'

def get_runs(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    run_df = df.sort_values(by='acquisition.timestamp')
    n_digits = max(2, len(str(run_df.shape[0])))
    padded_list = [str(n).zfill(n_digits) for n in range(1, run_df.shape[0] + 1)]
    run_df.insert(run_df.shape[1], 'run', padded_list)
    return run_df

def add_fmap_runs(df: pd.DataFrame()) -> pd.DataFrame():
    df = df.copy()
    group_df = df[['acquisition.timestamp', 'group']].drop_duplicates(subset='group', keep='first')
    if len(group_df) > 1:
        group_mapping = get_runs(group_df).set_index('group')
        df['run'] = df['group'].map(group_mapping['run'])
        df['reproin'] = df['reproin'] + '_run-' + df['run']

    return df

def timestamp_group(df: pd.DataFrame, delta: timedelta) -> pd.DataFrame:
    """Groups rows in a dataframe that fall within a time delta. Returns a df with
    added 'time_diff' and 'group' columns"""
    df = df.copy()
    df = df.sort_values(by='acquisition.timestamp')
    df['time_diff'] = df['acquisition.timestamp'].diff().fillna(
         pd.Timedelta(seconds=0)
    )
    df['group'] = (df['time_diff'] > delta).cumsum()
    return df
    

def add_fmap(df: pd.DataFrame()):
    df = df.copy()
    subject = df['subject.label'].iloc[0]
    fmap_df = df[df['file.modality'] == 'MR']
    fmap_df = fmap_df.groupby('acquisition.id').filter(
        lambda x: (
            not x[x['file.type'] == 'dicom']['file.classification.Intent'].isnull().values.any()
            and 'Fieldmap' in x[x['file.type'] == 'dicom']['file.classification.Intent'].iloc[0]
        )
    ).reset_index(drop=True)
    fmap_df['reproin_dict'] = fmap_df['reproin'].apply(
        lambda x: parse_series_spec(x)
    )

    if fmap_df.empty:
        return df

    # Determine fieldmap type (magnitude, phase, etc.)
    type_mapping = fmap_df.groupby('acquisition.id')[fmap_df.columns].apply(
        lambda x: classify_fmap_acq(
            x[x['file.type'] == 'dicom']['file.info.header.dicom_array.ImageType.0'].iloc[0]
        )
    )
    fmap_df['fmap_file_type'] = fmap_df['acquisition.id'].map(type_mapping)
    fmap_df = fmap_df[fmap_df['file.type'].apply(lambda x: x in ('nifti', 'bval', 'bvec'))]
    fmap_df['magnitude_type'] = fmap_df['file.name'].apply(
        lambda x: 1 if ('_e1.nii' in x or '_magnitude1.nii' in x)
            else (2 if ('_e2.nii' in x or '_magnitude2.nii' in x) else None)
    )
    
    # Group acquisitions by timestamp
    fmap_df = timestamp_group(fmap_df, FMAP_DELTA)
    fmap_df['fmap_acq'] = None
    fmap_df = fmap_df.set_index('file.file_id')
    for name, group in fmap_df.groupby('group'):
        if ((group.iloc[0]['acquisition.timestamp'] - group.iloc[-1]['acquisition.timestamp'])
             > FMAP_DELTA):
             log.error(f"Fieldmap group in subject {fmap_df.iloc[0]['subject.label']} spans more than {FMAP_DELTA}")
             return None
        elif group['reproin_dict'].apply(
                lambda x:'datatype_suffix' in x and x['datatype_suffix'] in FMAP_SUFFIXES
        ).all(): 
            continue

        niftis = group[group['file.type'] == 'nifti']
        magnitude = niftis[niftis['fmap_file_type'] == 'magnitude']
        phase = niftis[niftis['fmap_file_type'] == 'phase']
        diffusion = niftis[niftis['fmap_file_type'] == 'diffusion']


        if magnitude.empty and diffusion.empty:
            log.error(f"Subject {subject} has fieldmap(s) with no magnitude or diffusion files.")
            return None
        elif not magnitude.empty and not diffusion.empty:
            log.error(f"Subject {subject} has fieldmap(s) with both magnitude and diffusion files.")
            return None
        elif len(magnitude) > 2:
            log.error(f"Subject {subject} has fieldmap(s) with > 2 magnitude files.")
            return None
        elif len(phase) > 2:
            log.error(f"Subject {subject} has fieldmap(s) with > 2 phase files.")
            return None
        elif len(phase) in (1, 2):
            if magnitude['magnitude_type'].isna().any():
                log.error(f"Subject {subject} has magnitude(s) not matching expected filenames.")
                return None
            if len(phase) == 1:
                group['fmap_acq'] = 'phasediff'
            else:
                group['fmap_acq'] = 'twophase'
            group['reproin'] = 'fmap-gre'
        elif len(phase) == 0:
            if not diffusion.empty:
                mag_dif = diffusion
                if (group[group['file.type'] == 'bval'].empty 
                    or group[group['file.type'] == 'bvec'].empty):
                    log.error(f"Subject {subject} has diffusion fmaps without bval/bvecs.")
                    return None
                group['fmap_acq'] = 'dwi'
            else:
                mag_dif = magnitude
                group['fmap_acq'] = 'bold'

            mag_dif['direction'] = mag_dif['acquisition.label'].str.extract(
                r'.*(LR|RL|AP|PA)(?=[^a-zA-Z0-9]|$)'
            )[0]

            if len(mag_dif) == 1:
                if not mag_dif['direction'].isnull().all():
                    log.error(f"Subject {subject} has acquisition with 1 epis but direction.")
                    return None
                group['reproin'] = 'fmap-epi'
            elif len(mag_dif == 2):
                if mag_dif['direction'].isnull().any():
                    log.error(f"Subject {subject} has acquisition with 2 epis but no direction.")
                    return None
                elif len(mag_dif['acquisition.id'].unique()) == 1:
                    log.error(f"Subject {subject} has acquisition with 2 epis in one acquisition.")
                    return None

                reproin_mapping = mag_dif.groupby('acquisition.id')[mag_dif.columns].apply(
                    lambda x: 'fmap-epi_dir-' + x['direction'].iloc[0]
                )
                group['reproin'] = group['acquisition.id'].map(reproin_mapping)

        fmap_df.update(group)
    
    if fmap_df['fmap_acq'].nunique() > 1:
        fmap_df['reproin'] += '_acq-' + fmap_df['fmap_acq']
    fmap_df['reproin'] = fmap_df.groupby('fmap_acq')[fmap_df.columns].apply(
        add_fmap_runs).reset_index(level=0, drop=True)['reproin']
    
    fmap_df = fmap_df.set_index('acquisition.id')
    df = df.set_index('acquisition.id')
    fmap_mapping = fmap_df.groupby(fmap_df.index).first()
    fmap_mapping = fmap_mapping['reproin']
    fmap_mapping = fmap_mapping.dropna()
    df['reproin'] = df['reproin'].where(
        ~df.index.isin(fmap_mapping.index), df.index.map(fmap_mapping))

    return df.reset_index()

def add_rec(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    subject = df['subject.label'].iloc[0]
    rec_df = df[df['reproin'].str.startswith(('anat', 'asl'))]
    if rec_df.empty:
        return df
    rec_df = rec_df.groupby('reproin')[rec_df.columns].apply(
        lambda x: timestamp_group(x, REC_DELTA)
    ).reset_index(drop=True)
    
    for i, groupby in rec_df.groupby(['reproin', 'group']):
        if len(groupby) > 2:
            log.error(f"Subject {subject} has more than 2 acquisitions in the same rec group.")
            return None
        elif len(groupby) == 2:
            groupby['norm'] = groupby['file.info.header.dicom_array.ImageType.0'].apply(
                lambda x: 'NORM' in x)
            if len(groupby[groupby['norm']]) != 1:
                log.error(f"Subject {subject} has a rec group with 0 or 2 NORM acquisitions")
                return None
            acq_id_update = groupby[groupby['norm']]['acquisition.id'].iloc[0]
            df.loc[df['acquisition.id'] == acq_id_update, 'reproin'] += '_rec-norm' 
    
    return df

def add_run(df: pd.DataFrame) -> pd.DataFrame:
    skip_df = df[df['acquisition.label'].str.contains('sbref', case=False)]
    run_df = df[~df['file.file_id'].isin(skip_df['file.file_id'])]
    
    # TODO: Handle recs properly (make sure rec pairs get same run number)
    if run_df.shape[0] > 1:
        run_df = get_runs(run_df)
        if run_df['reproin'].iloc[[0]].str.startswith('fmap').item():
            run_df['reproin'] = run_df['reproin'] + '_' + run_df['run']
        else:
            run_df['reproin'] = run_df['reproin'] + '_run-' + run_df['run']
    return pd.concat([run_df, skip_df])

def add_sbref(df: pd.DataFrame()):
    df = df.copy()
    sbref_df = df[df['acquisition.label'].str.contains('sbref', case=False)]
    sbref_df = sbref_df[~sbref_df['acquisition.label'].str.endswith('_ignore-BIDS')]
    if not sbref_df.empty:
        func_dwi = df[df['reproin'].str.startswith(('func', 'dwi'))]
        func_dwi = func_dwi[~func_dwi['acquisition.label'].str.contains('sbref', case=False)]
        for i,sbref in sbref_df.iterrows():
            func_dwi['timedelta'] = func_dwi['acquisition.timestamp'] - sbref['acquisition.timestamp']
            match = func_dwi[
                (func_dwi['timedelta'] <= SBREF_DELTA) & (func_dwi['timedelta'] >= timedelta(seconds=0))
            ]

            if len(match) == 1:
                df.loc[i, 'sbref'] = match['acquisition.id'].iloc[0]
                df.loc[i, 'reproin'] = match['reproin'].iloc[0].removesuffix('SBRef') + '_SBRef'
            elif len(match) > 1:
                log.error(
                    "%s in subject %s should match exactly 1 bold image, but it matched 2: %s" % ({sbref['acquisition.label']}, {sbref['subject.label']}, {len(match)}, {match['acquisition.label']},)
                )
            else:
                log.error(
                    "%s in subject %s didn't match any bold images" % ({sbref['acquisition.label']}, {sbref['subject.label']})
                )

                return None
    return df




