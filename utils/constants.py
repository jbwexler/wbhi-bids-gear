#!/usr/bin/env python

from datetime import timedelta

WAIT_TIMEOUT = 3600 * 2
SBREF_DELTA = timedelta(seconds=40)
FMAP_DELTA = timedelta(minutes=1)
REC_DELTA = timedelta(minutes=1)
_SHARED = "file.info.header.dicom.SharedFunctionalGroupsSequence.0"
_PERFRAME = "file.info.header.dicom.PerFrameFunctionalGroupsSequence"
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
    "file.info.header.dicom.AcquisitionTime": str,
    "file.info.header.dicom.AcquisitionDate": str,
    "file.info.header.dicom.AcquisitionDateTime": str,
    "file.info.header.dicom.ImageType": list,
    "file.info.header.dicom.ManufacturerModelName": str,
    "file.info.header.dicom.SoftwareVersions": str,
    "file.info.header.dicom.ScanningSequence": str,
    "file.info.header.dicom.SequenceVariant": str,
    "file.info.header.dicom.SequenceName": str,
    "file.info.header.dicom.SliceThickness": float,
    "file.info.header.dicom.EchoTime": float,
    "file.info.header.dicom.RepetitionTime": float,
    "file.info.header.dicom.InversionTime": float,
    "file.info.header.dicom.FlipAngle": float,
    "file.info.header.dicom.BaseResolution": int,
    "file.info.header.dicom.ParallelReductionFactorInPlane": float,
    f"{_SHARED}.MRImageFrameTypeSequence.0.FrameType": list,
    f"{_SHARED}.PixelMeasuresSequence.0.SliceThickness": float,
    f"{_SHARED}.MREchoSequence.0.EffectiveEchoTime": float,
    f"{_SHARED}.MRTimingAndRelatedParametersSequence.0.RepetitionTime": float,
    f"{_SHARED}.MRTimingAndRelatedParametersSequence.0.FlipAngle": float,
    f"{_SHARED}.MRModifierSequence.0.InversionTimes": list,
    f"{_SHARED}.MRModifierSequence.0.ParallelReductionFactorInPlane": float,
    f"{_SHARED}.MRFOVGeometrySequence.0.MRAcquisitionFrequencyEncodingSteps": int,
    f"{_PERFRAME}.PixelMeasuresSequence.0.SliceThickness": float,
    f"{_PERFRAME}.MREchoSequence.0.EffectiveEchoTime": float,
    "file.info.header.dicom.PulseSequenceName": str,
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
ACQ_PARAMS = (
    "file.info.header.dicom.ImageType",
    "file.info.header.dicom.ManufacturerModelName",
    "file.info.header.dicom.SoftwareVersions",
    "file.info.header.dicom.ScanningSequence",
    "file.info.header.dicom.SequenceVariant",
    "file.info.header.dicom.SequenceName",
    "file.info.header.dicom.SliceThickness",
    "file.info.header.dicom.EchoTime",
    "file.info.header.dicom.RepetitionTime",
    "file.info.header.dicom.InversionTime",
    "file.info.header.dicom.FlipAngle",
    "file.info.header.dicom.BaseResolution",
    "file.info.header.dicom.ParallelReductionFactorInPlane",
)
# ACQ_PARAMS that are allowed to be empty/NaN without raising a missing-params error.
OPTIONAL_ACQ_PARAMS = (
    "file.info.header.dicom.ScanningSequence",
    "file.info.header.dicom.SequenceVariant",
)
# Fallback paths to try in order when the top-level ACQ_PARAMS field is empty.
ACQ_PARAMS_FALLBACK = {
    "file.info.header.dicom.ImageType": (
        f"{_SHARED}.MRImageFrameTypeSequence.0.FrameType",
    ),
    "file.info.header.dicom.SliceThickness": (
        f"{_SHARED}.PixelMeasuresSequence.0.SliceThickness",
        f"{_PERFRAME}.PixelMeasuresSequence.0.SliceThickness",
    ),
    "file.info.header.dicom.EchoTime": (
        f"{_SHARED}.MREchoSequence.0.EffectiveEchoTime",
        f"{_PERFRAME}.MREchoSequence.0.EffectiveEchoTime",
    ),
    "file.info.header.dicom.RepetitionTime": (
        f"{_SHARED}.MRTimingAndRelatedParametersSequence.0.RepetitionTime",
    ),
    "file.info.header.dicom.FlipAngle": (
        f"{_SHARED}.MRTimingAndRelatedParametersSequence.0.FlipAngle",
    ),
    "file.info.header.dicom.InversionTime": (
        f"{_SHARED}.MRModifierSequence.0.InversionTimes",
    ),
    "file.info.header.dicom.BaseResolution": (
        f"{_SHARED}.MRFOVGeometrySequence.0.MRAcquisitionFrequencyEncodingSteps",
    ),
    "file.info.header.dicom.ParallelReductionFactorInPlane": (
        f"{_SHARED}.MRModifierSequence.0.ParallelReductionFactorInPlane",
    ),
    "file.info.header.dicom.SequenceName": (
        "file.info.header.dicom.PulseSequenceName",
    ),
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
