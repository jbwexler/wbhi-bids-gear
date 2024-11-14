#!/usr/bin/env python

from datetime import timedelta

WAIT_TIMEOUT = 3600 * 2
SBREF_DELTA = timedelta(seconds=30)
FMAP_DELTA = timedelta(minutes=1)
DATAVIEW_COLUMNS = (
    'subject.id',
    'subject.label',
    'session.id',
    'session.tags',
    'session.info.BIDS',
    'acquisition.label',
    'acquisition.id',
    'file.info.header.dicom.SeriesDescription',
    'file.info.header.dicom_array.ImageType.0',
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
)
ALLOWED_DATATYPES = ('func', 'anat', 'dwi', 'fmap', 'asl')
FMAP_SUFFIXES = ('gre', 'epi', 'phasediff', 'phase1', 'phase2', 'magnitude1', 'magnitude2',
                 'magnitude', 'fieldmap', 'm0scan')
