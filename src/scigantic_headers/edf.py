"""EDF, EDF+, BDF and BDF+ header decoding.

The European Data Format is what clinical and research electrophysiology rigs
write: polysomnography, EEG, EMG, ECG. BDF is BioSemi's 24-bit variant of the
same layout. Both open with a 256-byte ASCII record describing the recording,
followed by one 256-byte block per signal holding that signal's label, units,
calibration and sample count. Everything needed to describe a session is in
those bytes; the samples that follow are never read here.

The signal header is column-major, not one block per signal: all labels, then
all transducers, then all physical dimensions, and so on. A short read is
therefore partial in every column at once, which `headerComplete` reports.

EDF+ keeps structured subfields in the two free-text fields (patient and
recording) and carries an unabbreviated start date there, which is the only
way to tell 1989 from 2089 in a file whose header date is `24.04.89`.

Specs: https://www.edfplus.info/specs/edf.html and
https://www.biosemi.com/faq/file_format.htm
"""

import re
from typing import Dict, List, Optional, Tuple

from .decoders import DecodedHeader, Read, _finite, register_decoder

# 256 bytes fixed plus 256 per signal. 256 KiB covers 1023 signals; a clinical
# montage is tens of channels and a dense BDF cap is 256, so this is headroom,
# and a file with more signals still decodes as far as the read goes.
EDF_READ_BYTES = 256 * 1024

_FIXED_BYTES = 256
_SIGNAL_BYTES = 256
# Signals described in full. A long montage is listed by label instead.
_SIGNAL_DETAIL = 32
_LABEL_CAP = 256

_BDF_MAGIC = b"\xffBIOSEMI"
_EDF_MAGIC = b"0       "

_DATE = re.compile(r"^(\d{2})\.(\d{2})\.(\d{2})$")
_TIME = re.compile(r"^(\d{2})\.(\d{2})\.(\d{2})$")
# EDF+ recording field: "Startdate 24-APR-1989 hospital technician equipment".
_STARTDATE = re.compile(r"^Startdate (\d{2})-([A-Za-z]{3})-(\d{4})\b")
_MONTHS = {
    "JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5, "JUN": 6,
    "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10, "NOV": 11, "DEC": 12,
}
# The reserved label of the EDF+/BDF+ annotations channel. A channel is
# recognised by this label alone, including in a legacy file whose reserved
# field does not say EDF+: the label is reserved either way, and a channel
# carrying it holds timestamped text, not samples, so it is kept out of the
# recording's sampling rate.
_ANNOTATIONS = ("EDF Annotations", "BDF Annotations")


def _text(data: bytes, offset: int, width: int) -> str:
    """One fixed-width ASCII field, trailing padding removed. EDF is ASCII but
    writers do emit Latin-1 accents in names, which never raises here."""
    return data[offset:offset + width].decode("latin-1").strip()


def _int(text: str) -> Optional[int]:
    try:
        return int(text.strip())
    except ValueError:
        return None


def _float(text: str) -> Optional[float]:
    try:
        return _finite(float(text.strip()))
    except ValueError:
        return None


def _start_datetime(date_text: str, time_text: str, recording: str) -> Optional[str]:
    """ISO-8601 local start, or None when either field is not a real date.

    The header's two-digit year is clipped by the spec: 85-99 is 19xx, 00-84 is
    20xx. An EDF+ recording field overrides it with the four-digit year, which
    is what files written before 1985 or after 2084 rely on."""
    d, t = _DATE.match(date_text), _TIME.match(time_text)
    if not d or not t:
        return None
    day, month, yy = int(d.group(1)), int(d.group(2)), int(d.group(3))
    year = 1900 + yy if yy >= 85 else 2000 + yy
    m = _STARTDATE.match(recording)
    if m and m.group(2).upper() in _MONTHS:
        day, month, year = int(m.group(1)), _MONTHS[m.group(2).upper()], int(m.group(3))
    hour, minute, second = int(t.group(1)), int(t.group(2)), int(t.group(3))
    if not (1 <= month <= 12 and 1 <= day <= 31 and hour < 24 and minute < 60 and second < 61):
        return None
    return f"{year:04d}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:02d}"


def _patient_subfields(patient: str) -> Dict[str, object]:
    """EDF+ splits the patient field into code, sex, birthdate and name, space
    separated, `X` for unknown. A plain EDF file keeps free text there, which
    has no subfields and is left alone."""
    parts = patient.split(" ")
    if len(parts) < 4:
        return {}
    code, sex, birth, name = parts[0], parts[1], parts[2], " ".join(parts[3:])
    if sex not in ("M", "F", "X"):
        return {}
    out: Dict[str, object] = {}
    if code != "X":
        out["subjectCode"] = code
    if sex != "X":
        out["subjectSex"] = sex
    if birth != "X":
        out["subjectBirthdate"] = birth
    if name and name != "X":
        out["subjectName"] = name
    return out


def _recording_subfields(recording: str) -> Dict[str, object]:
    """EDF+ recording field: `Startdate <date> <hospital> <technician>
    <equipment>`, again `X` for unknown."""
    if not _STARTDATE.match(recording):
        return {}
    parts = recording.split(" ")
    out: Dict[str, object] = {}
    for key, idx in (("hospitalCode", 2), ("technician", 3), ("equipment", 4)):
        if len(parts) > idx and parts[idx] != "X":
            out[key] = parts[idx]
    return out


def _column(data: bytes, offset: int, width: int, count: int) -> Tuple[List[str], int]:
    """One column of the signal header: `count` fields of `width` bytes from
    `offset`. Returns the fields the read actually holds and the offset of the
    next column, which is computed from the full count so a short read does not
    shift the columns that follow."""
    out = []
    for i in range(count):
        start = offset + i * width
        if start + width > len(data):
            break
        out.append(data[start:start + width].decode("latin-1").strip())
    return out, offset + count * width


def _flavour(kind: str, reserved: str) -> str:
    """`EDF`, `EDF+C`, `EDF+D` or the BDF equivalents. The reserved field names
    the continuity: C for one uninterrupted recording, D for interrupted."""
    tag = reserved.split(" ")[0].upper() if reserved else ""
    if tag in (f"{kind}+C", f"{kind}+D"):
        return tag
    return kind


def decode_edf(data: bytes) -> Optional[DecodedHeader]:
    if len(data) < _FIXED_BYTES:
        return None
    version = data[:8]
    if version == _EDF_MAGIC:
        kind, bit_depth = "EDF", 16
    elif version == _BDF_MAGIC:
        kind, bit_depth = "BDF", 24
    else:
        return None

    header_bytes = _int(_text(data, 184, 8))
    signal_count = _int(_text(data, 252, 4))
    if signal_count is None or signal_count < 0 or signal_count > 65536:
        return None
    # The two counts are redundant by construction, which is what makes them a
    # cheap validity check on bytes that are otherwise all free text.
    if header_bytes != _FIXED_BYTES * (signal_count + 1):
        return None

    patient = _text(data, 8, 80)
    recording = _text(data, 88, 80)
    reserved = _text(data, 192, 44)
    records = _int(_text(data, 236, 8))
    record_seconds = _float(_text(data, 244, 8))

    fields: Dict[str, object] = {
        "edfType": _flavour(kind, reserved),
        "bitDepth": bit_depth,
        "headerBytes": header_bytes,
        "signalCount": signal_count,
    }
    start = _start_datetime(_text(data, 168, 8), _text(data, 176, 8), recording)
    if start:
        fields["startDateTime"] = start
    if patient:
        fields["patient"] = patient
    fields.update(_patient_subfields(patient))
    if recording:
        fields["recording"] = recording
    fields.update(_recording_subfields(recording))
    if reserved:
        fields["reserved"] = reserved
    # -1 is the spec's "still being written, count unknown".
    if records is not None and records >= 0:
        fields["numDataRecords"] = records
    if record_seconds is not None:
        fields["dataRecordDurationSeconds"] = record_seconds

    labels, off = _column(data, _FIXED_BYTES, 16, signal_count)
    transducers, off = _column(data, off, 80, signal_count)
    dimensions, off = _column(data, off, 8, signal_count)
    phys_min, off = _column(data, off, 8, signal_count)
    phys_max, off = _column(data, off, 8, signal_count)
    dig_min, off = _column(data, off, 8, signal_count)
    dig_max, off = _column(data, off, 8, signal_count)
    prefilters, off = _column(data, off, 80, signal_count)
    samples, off = _column(data, off, 8, signal_count)

    fields["headerComplete"] = len(data) >= header_bytes
    if labels:
        fields["signalLabels"] = labels[:_LABEL_CAP]

    signals: List[Dict[str, object]] = []
    rates = set()
    annotation_count = 0
    for i in range(min(len(labels), signal_count)):
        label = labels[i]
        entry: Dict[str, object] = {"label": label}
        if i < len(dimensions) and dimensions[i]:
            entry["physicalDimension"] = dimensions[i]
        for key, column in (
            ("physicalMin", phys_min), ("physicalMax", phys_max),
        ):
            if i < len(column):
                value = _float(column[i])
                if value is not None:
                    entry[key] = value
        for key, column in (("digitalMin", dig_min), ("digitalMax", dig_max)):
            if i < len(column):
                value = _int(column[i])
                if value is not None:
                    entry[key] = value
        per_record = _int(samples[i]) if i < len(samples) else None
        if per_record is not None:
            entry["samplesPerRecord"] = per_record
            if record_seconds:
                rate = _finite(per_record / record_seconds)
                if rate is not None:
                    entry["samplingRateHz"] = rate
                    if label not in _ANNOTATIONS:
                        rates.add(rate)
        if i < len(transducers) and transducers[i]:
            entry["transducer"] = transducers[i]
        if i < len(prefilters) and prefilters[i]:
            entry["prefiltering"] = prefilters[i]
        if label in _ANNOTATIONS:
            annotation_count += 1
            entry["annotations"] = True
        signals.append(entry)

    if signals:
        fields["signals"] = signals[:_SIGNAL_DETAIL]
    if annotation_count:
        fields["annotationSignalCount"] = annotation_count
    if len(rates) == 1:
        fields["samplingRateHz"] = rates.pop()
    duration = None
    if records is not None and records >= 0 and record_seconds:
        duration = _finite(records * record_seconds)
        if duration is not None:
            fields["durationSeconds"] = duration

    parts = [f"{fields['edfType']} recording", f"{signal_count} signals"]
    if duration is not None:
        parts.append(_duration(duration))
    rate = fields.get("samplingRateHz")
    if isinstance(rate, float):
        parts.append(f"{_number(rate)} Hz")
    return DecodedHeader(format="edf", summary=", ".join(parts), fields=fields)


def _number(value: float) -> str:
    """A rate or count without a trailing `.0`, so 100.0 Hz reads as 100 Hz."""
    return str(int(value)) if value == int(value) else f"{value:g}"


def _duration(seconds: float) -> str:
    """`h:mm:ss` for a recording that runs an hour or more, else `m:ss`."""
    total = int(seconds)
    hours, rest = divmod(total, 3600)
    minutes, secs = divmod(rest, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{secs:02d}"
    return f"{minutes}:{secs:02d}"


register_decoder(["edf", "bdf"], decode_edf, read=Read(leading=EDF_READ_BYTES))
