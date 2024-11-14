#!/usr/bin/env python

import re
import pandas as pd
from heudiconv.heuristics.reproin import parse_series_spec
import sys
import logging
from datetime import timedelta
from utils.constants import ALLOWED_DATATYPES, SBREF_DELTA, FMAP_DELTA, FMAP_SUFFIXES

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

def reproin_filter(row: pd.Series) -> str:
    label = row['acquisition.label']
    if row['file.type'] != 'dicom':
        return None
    elif '_ignore-BIDS' in label:
        return None
    
    # Utilize acq label if already in reproin format. Otherwise, create a reproin label.
    label_re = re.sub(r'_\d$', '', label)
    label_ignore = label + '_ignore-BIDS'
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

        if label.startswith('GOBRAIN_'):
            return label_ignore
        elif label == 't2_tse_tra_hi-res_hippocampus':
            return 'anat-T2w_acq-hippocampus'
        elif intent and 'Localizer' in intent:
            return label_ignore
        elif measurement and 'Perfusion' in measurement:
            return 'asl'
        elif measurement and intent and 'Structural' in intent:
            if 'T1' in measurement:
                return 'anat-T1w'
            elif 'T2' in measurement:
                return 'anat-T2w'   
            elif 'Diffusion' in measurement:
                return 'dwi'
        elif (intent and 'Functional' in intent and label
            and ('rest' in label.lower() or 'rsfmri' in label.lower())):
            return 'func_task-rest'
        elif intent and 'Fieldmap' in intent:
            return 'fmap'
        else:
            return label_ignore

def add_nifti_filenames(sub_df: pd.DataFrame) -> pd.DataFrame:
    sub_df['niftis'] = sub_df.apply(
        lambda x: [y.split('.nii')[0] for y in filenames_list if x['file_root'] in y],
        axis=1
    )
    print(sub_df[sub_df['file.modality'] != 'SR'][['subject.label','niftis']])
    if sub_df[sub_df['file.modality'] != 'SR']['niftis'].apply(lambda x: x == []).any():
        return None
    else:
        return sub_df

def add_sbref(df: pd.DataFrame()):
    sbref_df = df[df['file.classification.Features'].apply(
        lambda x: x is not None and 'SBRef' in x
    )]
    if not sbref_df.empty:
        bold = df[df['reproin'].str.startswith('func')]
        bold = bold[bold['file.classification.Features'].apply(
            lambda x: x is None or 'SBRef' not in x
        )]
        for i,sbref in sbref_df.iterrows():
            bold['timedelta'] = bold['acquisition.timestamp'] - sbref['acquisition.timestamp']
            match = bold[
                (bold['timedelta'] <= SBREF_DELTA) & (bold['timedelta'] >= timedelta(seconds=0))
            ]

            if len(match) == 1:
                df.loc[i, 'sbref'] = match['acquisition.id'].iloc[0]
                df.loc[i, 'reproin'] = match['reproin'].iloc[0] + '_SBRef'
            else:
                log.error(f"{sbref['acquisition.label']} in subject {sbref['subject.id']} "
                f"should match exactly 1 bold image, but it matched {len(match)}: "
                f"{match['acquisition.label']}")
                sys.exit(1)
    return df

def classify_fmap_file(image_type_dict: dict) -> str:
    """Determines the type of fieldmap file (phase, magnitude or diffusion) based on
    the presence of corresponding strings in the ImageType dict."""
    if 'M' in image_type_dict:
        return 'magnitude'
    elif 'P' in image_type_dict:
        return 'phase'
    elif 'DIFFUSION' in image_type_dict:
        return 'diffusion'

def add_fmap(df: pd.DataFrame()):
    fmap_df = df[df['file.classification.Intent'].apply(
        lambda x: x is not None and 'Fieldmap' in x
    )]
    fmap_df['reproin_dict'] = fmap_df['reproin'].apply(
        lambda x: parse_series_spec(x)
    )

    print(df.iloc[0]['session.id'])
    if fmap_df.empty:
        return df

    fmap_df = fmap_df.sort_values(by='acquisition.timestamp')
    fmap_df['time_diff'] = fmap_df['acquisition.timestamp'].diff().fillna(
         pd.Timedelta(seconds=0)
    )
    fmap_df['group'] = (fmap_df['time_diff'] > FMAP_DELTA).cumsum()
    fmap_df['fmap_type'] = fmap_df['file.info.header.dicom_array.ImageType.0'].apply(
        classify_fmap_file
    )
    for name,group in fmap_df.groupby('group'):
        if (
             group.iloc[0]['acquisition.timestamp'] - 
             group.iloc[-1]['acquisition.timestamp']
        ) > FMAP_DELTA:
             log.error(f"Fieldmap group in subject {fmap_df.iloc[o]['subject.label']} spans more than {FMAP_DELTA}")
        elif group['reproin_dict'].apply(
                lambda x:'datatype_suffix' in x and x['datatype_suffix'] in FMAP_SUFFIXES
        ).all(): 
            continue

        
        magnitude = [x for xs in group[group['fmap_type'] == 'magnitude']['niftis'] for x in xs]
        phase = [x for xs in group[group['fmap_type'] == 'phase']['niftis'] for x in xs]
        diffusion = [x for xs in group[group['fmap_type'] == 'diffusion']['niftis'] for x in xs]
         
        if not magnitude:
            session = df.iloc[0]['session.label']
            log.error(f"Session {session} has fieldmap(s) with no magnitude files.")
            return None
        elif len(magnitude) > 2:
            log.error(f"Session {session} has fieldmap(s) with > 2 magnitude files.")
        elif len(phase):
            pass

        breakpoint()
    return df

def add_run(df: pd.DataFrame()):
    skip_df = df[df.apply(
        lambda x: 
            (x['file.classification.Features'] is not None
                and 'SBRef' in x['file.classification.Features'])
            or (x['file.classification.Intent'] is not None
                and'Fieldmap' in x['file.classification.Intent']),
        axis=1
    )]
    if not skip_df.empty:
        print(skip_df)
    run_df = df[~df.isin(skip_df).all(axis=1)]

    if run_df.shape[0] > 1:
        run_df = run_df.sort_values(by='acquisition.timestamp')
        n_digits = max(2, len(str(run_df.shape[0])))
        padded_list = [str(n).zfill(n_digits) for n in range(1, run_df.shape[0] + 1)]
        run_df.insert(run_df.shape[1], 'run', padded_list)
        if run_df['reproin'].iloc[[0]].str.startswith('fmap').item():
            run_df['reproin'] = run_df['reproin'] + '_' + run_df['run']
        else:
            run_df['reproin'] = run_df['reproin'] + '_run-' + run_df['run']
    return pd.concat([run_df, skip_df])




