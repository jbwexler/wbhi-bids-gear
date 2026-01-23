#!/usr/bin/env python

from datetime import timedelta

WAIT_TIMEOUT = 3600 * 2
SBREF_DELTA = timedelta(seconds=40)
FMAP_DELTA = timedelta(minutes=1)
REC_DELTA = timedelta(minutes=1)
DATAVIEW_COLUMNS = {
    "subject.id": str,
    "subject.label": str,
    "session.id": str,
    "session.label": str,
    "session.tags": list,
    "session.info.BIDS": dict,
    "acquisition.label": str,
    "acquisition.id": str,
    "file.info.header.dicom.SeriesDescription": str,
    "file.info.header.dicom.ImageType": list,
    "file.info.header.dicom_array.ImageType.0": list,
    "file.classification.Intent": list,
    "file.classification.Measurement": list,
    "file.classification.Features": list,
    "file.modality": str,
    "file.file_id": str,
    "file.tags": list,
    "file.type": str,
    "file.created": str,
    "file.name": str,
    "acquisition.timestamp": str,
}
ALLOWED_DATATYPES = ("func", "anat", "dwi", "fmap", "asl")
FMAP_SUFFIXES = (
    "gre",
    "epi",
    "phasediff",
    "phase1",
    "phase2",
    "magnitude1",
    "magnitude2",
    "magnitude",
    "fieldmap",
    "m0scan",
)
IGNORE_SCANS = ("BIAS_BC", "BIAS_32CH", "BIAS_BodyCoil", "BIAS_HeadCoil")
DWI_SUFFIXES_IGNORE = ("ADC", "TRACEW", "FA", "COLFA")
