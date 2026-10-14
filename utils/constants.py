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
    "file.info.SeriesDescription": str,
    "file.info.ProtocolName": str,
    "file.info.header.dicom.AcquisitionTime": str,
    "file.info.header.dicom.AcquisitionDate": str,
    "file.info.header.dicom.AcquisitionDateTime": str,
    "file.info.ImageType": list,
    "file.info.ManufacturerModelName": str,
    "file.info.SoftwareVersions": str,
    "file.info.ScanningSequence": str,
    "file.info.SequenceVariant": str,
    "file.info.SequenceName": str,
    "file.info.PulseSequenceName": str,
    "file.info.SliceThickness": float,
    "file.info.EchoTime": float,
    "file.info.RepetitionTime": float,
    "file.info.FlipAngle": float,
    "file.info.SeriesNumber": int,
    "file.info.InversionTime": float,
    "file.info.BaseResolution": int,
    "file.info.ParallelReductionFactorInPlane": float,
    "file.info.header.nifti.dim.0": int,
    "file.info.header.nifti.dim.1": int,
    "file.info.header.nifti.dim.2": int,
    "file.info.header.nifti.dim.3": int,
    "file.info.header.dicom.SeriesDescription": str,
    "file.info.header.dicom.ImageType": list,
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
ACQ_PARAMS = (
    "file.info.ImageType",
    "file.info.ManufacturerModelName",
    "file.info.SoftwareVersions",
    "file.info.ScanningSequence",
    "file.info.SequenceVariant",
    "file.info.SequenceName",
    "file.info.SliceThickness",
    "file.info.EchoTime",
    "file.info.RepetitionTime",
    "file.info.FlipAngle",
    "file.info.InversionTime",
    "file.info.BaseResolution",
    "file.info.ParallelReductionFactorInPlane",
)
# Fallback paths to try in order when the top-level ACQ_PARAMS field is empty.
ACQ_PARAMS_FALLBACK = {
    "file.info.SequenceName": (
        "file.info.PulseSequenceName",
    ),
}
ALLOWED_DATATYPES = ("func", "anat", "dwi", "fmap", "asl")
# Map lowercase anat datatype_suffix from parse_series_spec to proper BIDS case.
ANAT_SUFFIX_CASE = {
    "t1w": "T1w",
    "t2w": "T2w",
    "flair": "FLAIR",
}
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
IGNORE_SCANS = ("bias_bc", "bias_32ch", "bias_bodycoil", "bias_headcoil")
DWI_SUFFIXES_IGNORE = ("adc", "tracew", "fa", "colfa")
