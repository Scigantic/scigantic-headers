"""Neuralynx Cheetah and Pegasus file header decoding.

Every Neuralynx acquisition file, whatever it holds, opens with the same
16 KiB ASCII header followed by fixed-size binary records: continuously
sampled channels (`.ncs`), events (`.nev`), single-electrode and tetrode
spikes (`.nse`, `.ntt`) and video tracking (`.nvt`). The header is the whole
description, and it is the same 16384 bytes in every file the rig writes.

Two header dialects exist and both are read here. Cheetah before 5.6.4 wrote
the file name and open/close times as `##` comment lines and the application
version as `-CheetahRev`; 5.7 onwards and Pegasus write them as ordinary
`-Key Value` properties (`-OriginalFileName`, `-TimeCreated`, `-TimeClosed`,
`-ApplicationName Cheetah "5.7.4"`).

Checked against the header grammar in python-neo's `NlxHeader`, which is the
reference reader for this format.
"""

import re
from typing import Dict, List, Optional

from .decoders import DecodedHeader, Read, _finite, register_decoder

# The header is exactly this long in every Neuralynx file, NUL padded.
NEURALYNX_HEADER_BYTES = 16 * 1024

_MAGIC = b"########"
# `-Key Value`, where the value runs to end of line and may hold spaces.
_PROPERTY = re.compile(r"^-(\S+)[ \t]*(.*)$")
# `## File Name: C:\...\CSC1.ncs` and `## File Name C:\...`, both seen.
_FILENAME = re.compile(r"^## File Name:?\s+(.*)$")
# `## Time Opened (m/d/y): 11/28/2016  (h:m:s.ms) 21:50:33.322` and the
# `Date Opened: (mm/dd/yyy): ... At Time: ...` spelling neuraview writes.
_TIMESTAMP = re.compile(
    r"^## (?:Time|Date) (Opened|Closed):*\s*\((?:m/d/y|mm/dd/yyy)\):\s*(\S+)\s+"
    r"(?:\(h:m:s\.ms\)|At Time:)\s*(\S+)"
)
# `-ApplicationName Cheetah "5.7.4"`.
_APPLICATION = re.compile(r'^(\S+)\s+"([^"]*)"$')

_FIELD_CAP = 128

_STRINGS = (
    ("FileType", "fileType"), ("FileVersion", "fileVersion"),
    ("AcqEntName", "channelName"), ("ProbeName", "probeName"),
    ("HardwareSubSystemName", "hardwareSubSystemName"),
    ("HardwareSubSystemType", "hardwareSubSystemType"),
    ("AcquisitionSystem", "acquisitionSystem"), ("ReferenceChannel", "referenceChannel"),
    ("SessionUUID", "sessionUuid"), ("FileUUID", "fileUuid"),
    ("DspLowCutFilterType", "dspLowCutFilterType"),
    ("DspHighCutFilterType", "dspHighCutFilterType"),
    ("NLX_Base_Class_Type", "baseClassType"),
)
_INTS = (
    ("RecordSize", "recordBytes"), ("ADMaxValue", "adMaxValue"),
    ("NumADChannels", "adChannelCount"), ("ADChannel", "adChannel"),
    ("WaveformLength", "waveformSamples"), ("AlignmentPt", "alignmentPoint"),
    ("DspLowCutNumTaps", "dspLowCutNumTaps"), ("DspHighCutNumTaps", "dspHighCutNumTaps"),
)
_FLOATS = (
    ("SamplingFrequency", "samplingRateHz"), ("ADBitVolts", "adBitVolts"),
    ("DspLowCutFrequency", "dspLowCutFrequencyHz"),
    ("DspHighCutFrequency", "dspHighCutFrequencyHz"),
)
_BOOLS = (
    ("DSPLowCutFilterEnabled", "dspLowCutFilterEnabled"),
    ("DSPHighCutFilterEnabled", "dspHighCutFilterEnabled"),
    ("InputInverted", "inputInverted"),
    ("DspDelayCompensation", "dspDelayCompensation"),
)
# File type to the kind of record that follows the header.
_RECORD_KINDS = {
    "CSC": "continuous", "Event": "event", "Spike": "spike", "Video": "video",
}


def _first(values: Optional[List[str]]) -> Optional[str]:
    """Several channels in one file repeat a key once per channel; the first is
    the one a single-value field can carry. Cheetah quotes any value that holds
    spaces, and the quotes are not part of it."""
    if not values:
        return None
    value = values[0]
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        return value[1:-1]
    return value


def _number(text: str) -> Optional[float]:
    """A value that may be a bare number or a space-separated per-channel list;
    the first entry is taken."""
    token = text.split(" ")[0].strip()
    try:
        return _finite(float(token))
    except ValueError:
        return None


def _iso(date: str, time: str) -> Optional[str]:
    """`11/28/2016` plus `21:50:33.322` to ISO-8601, or `2017/02/16` plus
    `17:56:04` from the property form."""
    if "/" in date:
        parts = date.split("/")
        if len(parts) != 3:
            return None
        if len(parts[0]) == 4:
            year, month, day = parts
        else:
            month, day, year = parts
    elif "-" in date:
        split = date.split("-")
        if len(split) != 3:
            return None
        year, month, day = split
    else:
        return None
    try:
        y, m, d = int(year), int(month), int(day)
    except ValueError:
        return None
    if not (1 <= m <= 12 and 1 <= d <= 31):
        return None
    clock = time.split(":")
    if len(clock) != 3:
        return None
    try:
        hh, mm = int(clock[0]), int(clock[1])
        ss = float(clock[2])
    except ValueError:
        return None
    if hh > 23 or mm > 59 or ss >= 61:
        return None
    seconds = f"{ss:06.3f}".rstrip("0").rstrip(".") if ss % 1 else f"{int(ss):02d}"
    return f"{y:04d}-{m:02d}-{d:02d}T{hh:02d}:{mm:02d}:{seconds}"


def decode_neuralynx(data: bytes) -> Optional[DecodedHeader]:
    if not data.startswith(_MAGIC):
        return None
    header = data[:NEURALYNX_HEADER_BYTES].split(b"\x00", 1)[0].decode("latin-1")
    if "Neuralynx" not in header[:256]:
        return None

    properties: Dict[str, List[str]] = {}
    fields: Dict[str, object] = {}
    for line in header.splitlines():
        line = line.rstrip()
        if not line:
            continue
        name = _FILENAME.match(line)
        if name:
            fields.setdefault("originalFileName", name.group(1).strip())
            continue
        stamp = _TIMESTAMP.match(line)
        if stamp:
            when = _iso(stamp.group(2), stamp.group(3))
            if when:
                fields["timeOpened" if stamp.group(1) == "Opened" else "timeClosed"] = when
            continue
        prop = _PROPERTY.match(line)
        if prop:
            properties.setdefault(prop.group(1), []).append(prop.group(2).strip())

    for key, name in _STRINGS:
        value = _first(properties.get(key))
        if value:
            fields[name] = value
    for key, name in _INTS:
        value = _first(properties.get(key))
        if value is None:
            continue
        number = _number(value)
        if number is not None and number == int(number):
            fields[name] = int(number)
    for key, name in _FLOATS:
        value = _first(properties.get(key))
        if value is None:
            continue
        number = _number(value)
        if number is not None:
            fields[name] = number
    for key, name in _BOOLS:
        value = _first(properties.get(key))
        if value in ("True", "False"):
            fields[name] = value == "True"
    for key, name in (("InputRange", "inputRange"), ("OriginalFileName", "originalFileName")):
        value = _first(properties.get(key))
        if value:
            fields[name] = value

    # Application and version: `-CheetahRev 5.5.1` on older files,
    # `-ApplicationName Cheetah "5.7.4"` on newer ones.
    revision = _first(properties.get("CheetahRev"))
    application = _first(properties.get("ApplicationName"))
    if revision:
        fields["applicationName"] = "Cheetah"
        fields["applicationVersion"] = revision
    elif application:
        match = _APPLICATION.match(application)
        if match:
            fields["applicationName"], fields["applicationVersion"] = match.group(1), match.group(2)
        else:
            fields["applicationName"] = application
    elif "baseClassType" in fields:
        fields["applicationName"] = "BML"

    for key, name in (("TimeCreated", "timeOpened"), ("TimeClosed", "timeClosed")):
        value = _first(properties.get(key))
        if not value:
            continue
        parts = value.split()
        if len(parts) == 2:
            when = _iso(parts[0], parts[1])
            if when:
                fields[name] = when

    bit_volts = fields.get("adBitVolts")
    if isinstance(bit_volts, float) and bit_volts > 0:
        microvolts = _finite(bit_volts * 1e6)
        if microvolts is not None:
            fields["microvoltsPerBit"] = microvolts
    file_type = fields.get("fileType")
    if isinstance(file_type, str) and file_type in _RECORD_KINDS:
        fields["recordKind"] = _RECORD_KINDS[file_type]
    fields["headerBytes"] = NEURALYNX_HEADER_BYTES
    fields["propertyCount"] = len(properties)
    if properties:
        fields["propertyNames"] = sorted(properties)[:_FIELD_CAP]

    parts = [f"Neuralynx {file_type or 'file'}"]
    channel = fields.get("channelName")
    if channel:
        parts.append(str(channel))
    rate = fields.get("samplingRateHz")
    if isinstance(rate, float):
        parts.append(f"{str(int(rate)) if rate == int(rate) else f'{rate:g}'} Hz")
    application_name = fields.get("applicationName")
    if application_name:
        version = fields.get("applicationVersion")
        parts.append(f"{application_name} {version}" if version else str(application_name))
    return DecodedHeader(format="neuralynx", summary=", ".join(parts), fields=fields)


register_decoder(
    ["ncs", "nev", "nse", "ntt", "nvt"],
    decode_neuralynx,
    read=Read(leading=NEURALYNX_HEADER_BYTES),
)
