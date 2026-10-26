"""BIDS sidecar decoding: `dataset_description.json`, JSON metadata sidecars
and the tabular `.tsv` files.

The Brain Imaging Data Structure stores the science in ordinary neuroimaging
and electrophysiology files and everything about it in plain sidecars beside
them. A NIfTI volume does not record its repetition time, a BrainVision
triplet does not record the task, and no raw file records who the subject was:
those live in the sidecars, so a pipeline that reads only the data files sees
a session it cannot describe.

`.json` and `.tsv` are maximally generic extensions, so both decoders decline
anything that is not BIDS-shaped: a JSON sidecar must carry at least one key
the standard defines, and a TSV must have one of the standard's column
signatures (`onset`+`duration` for events, `participant_id` for the subject
table, and so on). Anything else costs one bounded read and returns nothing.

BIDS key names are kept exactly as the standard spells them rather than
recased to this library's usual convention: these are the names the standard,
its validator and every downstream tool use, and renaming them would make a
record impossible to match against the specification.

Spec: https://bids-specification.readthedocs.io/
"""

import json
from typing import Dict, List, Optional, Tuple

from .decoders import DecodedHeader, Read, _finite, register_decoder

# A participants table or an events file runs to thousands of rows; what is
# read is reported as a sample, never as a file total.
BIDS_READ_BYTES = 256 * 1024

_KEY_CAP = 128
_LIST_CAP = 16
_COLUMN_CAP = 64
_LEVEL_CAP = 32
_TEXT_CAP = 512

# Keys the standard defines that no other JSON is likely to carry. One of
# these must be present for a `.json` to be claimed as a BIDS sidecar.
_SIDECAR_KEYS = frozenset({
    # MRI acquisition
    "RepetitionTime", "EchoTime", "EchoTime1", "EchoTime2", "FlipAngle", "InversionTime",
    "MagneticFieldStrength", "ScanningSequence", "SequenceVariant", "SequenceName",
    "PulseSequenceType", "SliceTiming", "SliceEncodingDirection", "PhaseEncodingDirection",
    "EffectiveEchoSpacing", "TotalReadoutTime", "DwellTime", "MultibandAccelerationFactor",
    "ParallelReductionFactorInPlane", "ParallelAcquisitionTechnique", "PartialFourier",
    "RepetitionTimeExcitation", "RepetitionTimePreparation", "MTState", "NonlinearGradientCorrection",
    "B0FieldIdentifier", "B0FieldSource", "IntendedFor", "AnatomicalLandmarkCoordinates",
    # Task and electrophysiology
    "TaskName", "SamplingFrequency", "PowerLineFrequency", "RecordingDuration", "RecordingType",
    "EpochLength", "SoftwareFilters", "HardwareFilters", "EEGReference", "EEGChannelCount",
    "EEGGround", "EEGPlacementScheme", "iEEGReference", "iEEGGround", "iEEGPlacementScheme",
    "iEEGElectrodeGroups", "ECOGChannelCount", "SEEGChannelCount", "MEGChannelCount",
    "MEGREFChannelCount", "DewarPosition", "DigitizedLandmarks", "DigitizedHeadPoints",
    "ContinuousHeadLocalization", "TriggerChannelCount", "MiscChannelCount", "EMGChannelCount",
    "NIRSChannelCount", "NIRSSourceOptodeCount", "NIRSDetectorOptodeCount",
    # Coordinate systems
    "EEGCoordinateSystem", "iEEGCoordinateSystem", "MEGCoordinateSystem",
    "NIRSCoordinateSystem", "FiducialsCoordinateSystem",
    # PET
    "TracerName", "TracerRadionuclide", "InjectedRadioactivity", "ModeOfAdministration",
    "FrameTimesStart", "FrameDuration", "ReconMethodName", "AcquisitionMode", "ScanStart",
    # Microscopy
    "PixelSize", "PixelSizeUnits", "ChunkTransformationMatrix", "SampleEnvironment",
    "SampleStaining", "ImageAcquisitionProtocol",
    # Continuous recordings beside an imaging run (physio, stim, eye tracking)
    "Columns", "StartTime", "PhysioType",
    # Derivatives
    "SpatialReference", "SkullStripped", "RawSources",
})

# A coordinate-system sidecar names its space in a key suffixed this way
# (`EEGCoordinateSystem`, `EMGCoordinateUnits`, ...). Matching the suffix
# rather than listing every modality keeps a new one from being missed.
_COORDINATE_SUFFIXES = ("CoordinateSystem", "CoordinateUnits", "CoordinateSystemDescription")

# Carried verbatim into the record when present, in this order.
_SIDECAR_REPORTED = (
    "TaskName", "TaskDescription", "Modality", "Manufacturer", "ManufacturersModelName",
    "InstitutionName", "DeviceSerialNumber", "StationName", "SoftwareVersions",
    "MagneticFieldStrength", "RepetitionTime", "EchoTime", "FlipAngle", "InversionTime",
    "ScanningSequence", "SequenceVariant", "SequenceName", "PulseSequenceType",
    "SliceEncodingDirection", "PhaseEncodingDirection", "EffectiveEchoSpacing",
    "TotalReadoutTime", "MultibandAccelerationFactor", "ParallelReductionFactorInPlane",
    "PartialFourier", "BodyPartExamined", "MRAcquisitionType",
    "SamplingFrequency", "PowerLineFrequency", "RecordingDuration", "RecordingType",
    "EpochLength", "EEGReference", "EEGGround", "EEGPlacementScheme",
    "iEEGReference", "iEEGGround", "iEEGPlacementScheme", "iEEGElectrodeGroups",
    "EEGChannelCount", "ECOGChannelCount", "SEEGChannelCount", "EOGChannelCount",
    "ECGChannelCount", "EMGChannelCount", "MiscChannelCount", "TriggerChannelCount",
    "MEGChannelCount", "MEGREFChannelCount", "DewarPosition", "ContinuousHeadLocalization",
    "EEGCoordinateSystem", "iEEGCoordinateSystem", "MEGCoordinateSystem",
    "NIRSCoordinateSystem", "FiducialsCoordinateSystem",
    "TracerName", "TracerRadionuclide", "InjectedRadioactivity", "ModeOfAdministration",
    "ReconMethodName", "AcquisitionMode", "ScanStart", "Units",
    "SpatialReference", "SkullStripped", "RawSources",
    "Columns", "StartTime", "PhysioType",
    "PixelSize", "PixelSizeUnits", "BodyPart", "SampleEnvironment", "SampleStaining",
)

# Which part of the standard a sidecar belongs to, decided by the first
# matching group. Order matters: an iEEG sidecar also carries EEG keys.
_MODALITIES: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("ieeg", ("iEEGReference", "iEEGGround", "iEEGPlacementScheme", "iEEGElectrodeGroups",
              "iEEGCoordinateSystem", "ECOGChannelCount", "SEEGChannelCount")),
    ("meg", ("MEGChannelCount", "MEGREFChannelCount", "DewarPosition", "MEGCoordinateSystem",
             "ContinuousHeadLocalization")),
    ("nirs", ("NIRSChannelCount", "NIRSSourceOptodeCount", "NIRSDetectorOptodeCount",
              "NIRSCoordinateSystem")),
    ("eeg", ("EEGReference", "EEGGround", "EEGChannelCount", "EEGPlacementScheme",
             "EEGCoordinateSystem")),
    ("pet", ("TracerName", "TracerRadionuclide", "InjectedRadioactivity", "FrameTimesStart",
             "ReconMethodName", "ModeOfAdministration")),
    ("micr", ("PixelSize", "PixelSizeUnits", "ChunkTransformationMatrix", "SampleStaining",
              "SampleEnvironment")),
    ("mri", ("RepetitionTime", "EchoTime", "FlipAngle", "MagneticFieldStrength",
             "ScanningSequence", "PulseSequenceType", "SliceTiming", "PhaseEncodingDirection",
             "EffectiveEchoSpacing", "InversionTime", "MRAcquisitionType")),
)

# Column signatures, in order; the first whose required columns are all
# present names the file.
_TABLES: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("events", ("onset", "duration")),
    ("participants", ("participant_id",)),
    ("sessions", ("session_id",)),
    ("samples", ("sample_id",)),
    ("scans", ("filename",)),
    ("electrodes", ("name", "x", "y", "z")),
    ("optodes", ("name", "type", "x", "y", "z")),
    ("channels", ("name", "type", "units")),
)


def _scalar(value: object) -> bool:
    return isinstance(value, (str, int, float, bool))


def _carry(fields: Dict[str, object], key: str, value: object) -> None:
    """Copy one JSON value into the record: scalars verbatim (long strings
    truncated), short scalar lists verbatim, anything longer as a count."""
    if isinstance(value, str):
        fields[key] = value[:_TEXT_CAP]
    elif isinstance(value, bool) or isinstance(value, int):
        fields[key] = value
    elif isinstance(value, float):
        finite = _finite(value)
        if finite is not None:
            fields[key] = finite
    elif isinstance(value, list):
        if len(value) <= _LIST_CAP and all(_scalar(v) for v in value):
            fields[key] = value
        else:
            fields[key + "Count"] = len(value)


def _dataset_description(obj: dict) -> DecodedHeader:
    fields: Dict[str, object] = {"bidsFileType": "datasetDescription"}
    for key in ("Name", "BIDSVersion", "HEDVersion", "DatasetType", "License", "DatasetDOI"):
        if key in obj:
            _carry(fields, key, obj[key])
    for key, name in (
        ("Authors", "authors"), ("Funding", "funding"),
        ("ReferencesAndLinks", "referencesAndLinks"), ("EthicsApprovals", "ethicsApprovals"),
    ):
        value = obj.get(key)
        if isinstance(value, list):
            fields[name + "Count"] = len(value)
            entries = [v for v in value if isinstance(v, str)]
            if entries:
                fields[name] = [e[:_TEXT_CAP] for e in entries[:_LIST_CAP]]
    generated = obj.get("GeneratedBy")
    if isinstance(generated, list):
        names = []
        for entry in generated:
            if not isinstance(entry, dict):
                continue
            label = entry.get("Name")
            version = entry.get("Version")
            if isinstance(label, str):
                names.append(f"{label} {version}" if isinstance(version, str) else label)
        if names:
            fields["generatedBy"] = names[:_LIST_CAP]
    sources = obj.get("SourceDatasets")
    if isinstance(sources, list):
        fields["sourceDatasetCount"] = len(sources)
    for key in ("Acknowledgements", "HowToAcknowledge"):
        if isinstance(obj.get(key), str):
            fields[key] = obj[key][:_TEXT_CAP]
    fields["keyCount"] = len(obj)

    name = fields.get("Name")
    parts = ["BIDS dataset"]
    if isinstance(name, str) and name:
        parts.append(name)
    version = fields.get("BIDSVersion")
    if isinstance(version, str):
        parts.append(f"BIDS {version}")
    kind = fields.get("DatasetType")
    if isinstance(kind, str):
        parts.append(kind)
    return DecodedHeader(format="bids-dataset", summary=", ".join(parts), fields=fields)


def _data_dictionary(obj: dict) -> Optional[DecodedHeader]:
    """A sidecar that describes the columns of a `.tsv` beside it: every value
    is an object carrying `Description`, `LongName`, `Units` or `Levels`."""
    columns = [
        k for k, v in obj.items()
        if isinstance(v, dict) and ({"Description", "LongName", "Units", "Levels", "HED"} & set(v))
    ]
    if not columns or len(columns) * 2 < len(obj):
        return None
    fields: Dict[str, object] = {
        "bidsFileType": "dataDictionary",
        "columnCount": len(columns),
        "columns": columns[:_COLUMN_CAP],
    }
    described: Dict[str, object] = {}
    for name in columns[:_COLUMN_CAP]:
        entry = obj[name]
        out: Dict[str, object] = {}
        for key in ("LongName", "Description", "Units"):
            if isinstance(entry.get(key), str):
                out[key] = entry[key][:_TEXT_CAP]
        levels = entry.get("Levels")
        if isinstance(levels, dict):
            out["levelCount"] = len(levels)
            out["levels"] = sorted(levels)[:_LEVEL_CAP]
        if out:
            described[name] = out
    if described:
        fields["columnDescriptions"] = described
    return DecodedHeader(
        format="bids-dictionary",
        summary=f"BIDS data dictionary, {len(columns)} columns",
        fields=fields,
    )


def _coordinate_keys(obj: dict) -> List[str]:
    return [k for k in obj if k.endswith(_COORDINATE_SUFFIXES)]


def _sidecar(obj: dict) -> Optional[DecodedHeader]:
    coordinates = _coordinate_keys(obj)
    if not (_SIDECAR_KEYS & set(obj)) and not coordinates:
        return None
    fields: Dict[str, object] = {"bidsFileType": "sidecar"}
    if coordinates:
        fields["bidsModality"] = "coordsystem"
    else:
        for modality, keys in _MODALITIES:
            if any(k in obj for k in keys):
                fields["bidsModality"] = modality
                break
    for key in _SIDECAR_REPORTED:
        if key in obj:
            _carry(fields, key, obj[key])
    for key in coordinates:
        if key not in fields:
            _carry(fields, key, obj[key])
    timing = obj.get("SliceTiming")
    if isinstance(timing, list):
        fields["sliceCount"] = len(timing)
    frames = obj.get("FrameTimesStart")
    if isinstance(frames, list):
        fields["frameCount"] = len(frames)
    intended = obj.get("IntendedFor")
    if isinstance(intended, list):
        fields["intendedForCount"] = len(intended)
    elif isinstance(intended, str):
        fields["intendedForCount"] = 1
    # A physio or eye-tracking sidecar is a hybrid: recording metadata plus a
    # description object per column of the `.tsv.gz` beside it.
    described = [k for k, v in obj.items() if isinstance(v, dict) and ({"Description", "Units", "Levels"} & set(v))]
    if described:
        fields["columnDescriptionCount"] = len(described)
    fields["keyCount"] = len(obj)
    fields["keys"] = sorted(obj)[:_KEY_CAP]

    parts = ["BIDS sidecar"]
    modality = fields.get("bidsModality")
    if isinstance(modality, str):
        parts.append(modality)
    task = fields.get("TaskName")
    if isinstance(task, str) and task:
        parts.append(f"task {task}")
    for key, unit in (("RepetitionTime", "s TR"), ("SamplingFrequency", "Hz")):
        value = fields.get(key)
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            parts.append(f"{_number(float(value))} {unit}")
    return DecodedHeader(format="bids-sidecar", summary=", ".join(parts), fields=fields)


def decode_bids_json(data: bytes) -> Optional[DecodedHeader]:
    try:
        # A BOM is common in sidecars written on Windows and is not valid JSON.
        obj = json.loads(data.decode("utf-8", "replace").lstrip("\ufeff"))
    except ValueError:
        return None
    if not isinstance(obj, dict):
        return None
    if "BIDSVersion" in obj:
        return _dataset_description(obj)
    sidecar = _sidecar(obj)
    if sidecar is not None:
        return sidecar
    return _data_dictionary(obj)


def _rows(data: bytes) -> Tuple[List[List[str]], bool]:
    """Tab-separated rows from the leading read, and whether the read holds the
    whole file. A row cut in half by the read boundary is dropped."""
    text = data.decode("utf-8", "replace").lstrip("\ufeff")
    complete = len(data) < BIDS_READ_BYTES
    lines = text.split("\n")
    if not complete and lines:
        lines = lines[:-1]
    rows = [line.rstrip("\r").split("\t") for line in lines if line.strip()]
    return rows, complete


def _floats(rows: List[List[str]], index: int) -> List[float]:
    out = []
    for row in rows:
        if index >= len(row):
            continue
        try:
            value = _finite(float(row[index]))
        except ValueError:
            continue
        if value is not None:
            out.append(value)
    return out


def decode_bids_tsv(data: bytes) -> Optional[DecodedHeader]:
    rows, complete = _rows(data)
    if not rows:
        return None
    columns = [c.strip() for c in rows[0]]
    # A trailing tab on the header line leaves an empty name; that is a
    # formatting slip, not a different file, so it is dropped.
    while columns and not columns[-1]:
        columns.pop()
    if not columns or any(not c for c in columns):
        return None
    lookup = {c: i for i, c in enumerate(columns)}
    kind = next((name for name, required in _TABLES if all(c in lookup for c in required)), None)
    if kind is None:
        return None
    body = [r for r in rows[1:] if len(r) > 1]

    fields: Dict[str, object] = {
        "bidsFileType": kind,
        "columnCount": len(columns),
        "columns": columns[:_COLUMN_CAP],
        "rowsRead": len(body),
        "rowsComplete": complete,
    }
    if kind == "events":
        onsets = _floats(body, lookup["onset"])
        if onsets:
            fields["onsetFirstSeconds"] = onsets[0]
            fields["onsetLastSeconds"] = onsets[-1]
        durations = _floats(body, lookup["duration"])
        if durations:
            fields["durationMinSeconds"] = min(durations)
            fields["durationMaxSeconds"] = max(durations)
        if "trial_type" in lookup:
            fields.update(_counts(body, lookup["trial_type"], "trialTypes"))
    elif kind == "channels":
        fields.update(_counts(body, lookup["type"], "channelTypes"))
        if "sampling_frequency" in lookup:
            rates = sorted({v for v in _floats(body, lookup["sampling_frequency"])})
            if len(rates) == 1:
                fields["samplingRateHz"] = rates[0]
            elif rates:
                fields["samplingRatesHz"] = rates[:_LIST_CAP]
    elif kind in ("participants", "samples", "sessions", "electrodes", "optodes"):
        fields[{"participants": "participantCount", "samples": "sampleCount",
                "sessions": "sessionCount", "electrodes": "electrodeCount",
                "optodes": "optodeCount"}[kind]] = len(body)

    noun = {"events": "event", "participants": "participant", "samples": "sample",
            "sessions": "session", "scans": "scan", "channels": "channel",
            "electrodes": "electrode", "optodes": "optode"}[kind]
    if len(body) != 1:
        noun += "s"
    suffix = "" if complete else " read"
    return DecodedHeader(
        format="bids-table",
        summary=f"BIDS {kind} table, {len(body)} {noun}{suffix}, {len(columns)} columns",
        fields=fields,
    )


def _counts(rows: List[List[str]], index: int, name: str) -> Dict[str, object]:
    """Value counts for one column, most common first."""
    counts: Dict[str, int] = {}
    for row in rows:
        if index >= len(row):
            continue
        value = row[index].strip()
        if not value or value == "n/a":
            continue
        counts[value] = counts.get(value, 0) + 1
    if not counts:
        return {}
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return {name: dict(ordered[:_LEVEL_CAP]), name + "Distinct": len(counts)}


def _number(value: float) -> str:
    return str(int(value)) if value == int(value) else f"{value:g}"


register_decoder(["json"], decode_bids_json, read=Read(leading=BIDS_READ_BYTES))
register_decoder(["tsv"], decode_bids_tsv, read=Read(leading=BIDS_READ_BYTES))
