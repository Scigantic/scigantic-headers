"""BrainVision header (`.vhdr`) and marker (`.vmrk`) decoding.

Brain Products' format is a triplet: an INI-style `.vhdr` describing the
recording, a `.vmrk` holding the event markers, and a raw `.eeg` holding the
samples. The two text files are the whole description, which is why BIDS
names BrainVision one of its three allowed EEG and iEEG container formats.

Both files are read as text and both decline anything whose first line is not
their own magic sentence, so a stray INI file under either extension costs one
bounded read and nothing more. Sample data is never touched.

Spec: BrainVision Core Data Format 1.0, Brain Products GmbH.
"""

import re
from typing import Dict, List, Optional, Tuple

from .decoders import DecodedHeader, Read, _finite, register_decoder

# A header with 256 channels and a coordinates section is tens of KB; markers
# run much longer, and a marker file past this is reported as sampled.
BRAINVISION_READ_BYTES = 256 * 1024

_HEADER_MAGIC = "Brain Vision Data Exchange Header File"
_MARKER_MAGIC = "Brain Vision Data Exchange Marker File"
_VERSION = re.compile(r"Version\s+([0-9.]+)")
_CHANNEL = re.compile(r"^Ch(\d+)$", re.IGNORECASE)
_MARKER = re.compile(r"^Mk(\d+)$", re.IGNORECASE)
# The "New Segment" marker carries the wall-clock start as YYYYMMDDhhmmssffffff.
_SEGMENT_DATE = re.compile(r"^(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})(\d{6})?$")

# Channel names and marker descriptions escape a comma as \1.
_COMMA = "\\1"

_CHANNEL_DETAIL = 64
_NAME_CAP = 512
_MARKER_DETAIL = 16
_MARKER_TYPE_CAP = 32

# Sample width in bytes, for turning a file size into a sample count.
_BINARY_FORMATS = {
    "INT_16": 2, "UINT_16": 2, "INT_32": 4, "UINT_32": 4,
    "IEEE_FLOAT_32": 4, "IEEE_FLOAT_64": 8,
}


def _sections(text: str) -> Tuple[Dict[str, List[Tuple[str, str]]], Optional[str]]:
    """The file as `{section: [(key, value), ...]}`, plus the free text of the
    `[Comment]` section. Lines starting with `;` are comments; a line with no
    `=` outside `[Comment]` is ignored, which is how the spec's own commentary
    lines (`Data orientation: MULTIPLEXED=ch1,pt1 ...`) are skipped."""
    out: Dict[str, List[Tuple[str, str]]] = {}
    comment: List[str] = []
    section = ""
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip().lower()
            out.setdefault(section, [])
            continue
        if section == "comment":
            comment.append(raw)
            continue
        if not line or line.startswith(";"):
            continue
        key, sep, value = line.partition("=")
        if not sep:
            continue
        out.setdefault(section, []).append((key.strip(), value.strip()))
    body = "\n".join(comment).strip()
    return out, body or None


def _lookup(sections: Dict[str, List[Tuple[str, str]]], section: str, key: str) -> Optional[str]:
    for k, v in sections.get(section, ()):
        if k.lower() == key.lower():
            return v
    return None


def _int(text: Optional[str]) -> Optional[int]:
    if text is None:
        return None
    try:
        return int(text.strip())
    except ValueError:
        return None


def _float(text: Optional[str]) -> Optional[float]:
    if text is None:
        return None
    try:
        return _finite(float(text.strip().replace(",", ".")))
    except ValueError:
        return None


def _unescape(text: str) -> str:
    return text.replace(_COMMA, ",")


def _split(value: str, count: int) -> List[str]:
    """A comma-separated entry padded to `count` fields; trailing fields are
    optional in the spec and simply absent in most files."""
    parts = [p.strip() for p in value.split(",")]
    return parts + [""] * (count - len(parts))


def decode_brainvision_vhdr(data: bytes) -> Optional[DecodedHeader]:
    text = data.decode("utf-8", "replace")
    if not text.lstrip("﻿").startswith(_HEADER_MAGIC):
        return None
    sections, comment = _sections(text)

    fields: Dict[str, object] = {}
    version = _VERSION.search(text.splitlines()[0])
    if version:
        fields["brainVisionVersion"] = version.group(1)

    for key, name in (
        ("DataFile", "dataFile"), ("MarkerFile", "markerFile"),
        ("DataFormat", "dataFormat"), ("DataOrientation", "dataOrientation"),
        ("DataType", "dataType"), ("Codepage", "codepage"),
    ):
        value = _lookup(sections, "common infos", key)
        if value:
            fields[name] = value

    channel_count = _int(_lookup(sections, "common infos", "NumberOfChannels"))
    if channel_count is not None:
        fields["channelCount"] = channel_count
    interval = _float(_lookup(sections, "common infos", "SamplingInterval"))
    rate = None
    if interval:
        fields["samplingIntervalMicroseconds"] = interval
        rate = _finite(1_000_000.0 / interval)
        if rate is not None:
            fields["samplingRateHz"] = rate
    points = _int(_lookup(sections, "common infos", "DataPoints"))
    if points is not None:
        fields["dataPoints"] = points
        if rate:
            duration = _finite(points / rate)
            if duration is not None:
                fields["durationSeconds"] = duration
    for key, name in (("SegmentDataPoints", "segmentDataPoints"), ("AveragedSegments", "averagedSegments")):
        value = _int(_lookup(sections, "common infos", key))
        if value is not None:
            fields[name] = value
    averaged = _lookup(sections, "common infos", "Averaged")
    if averaged:
        fields["averaged"] = averaged.upper() == "YES"

    binary_format = _lookup(sections, "binary infos", "BinaryFormat")
    if binary_format:
        fields["binaryFormat"] = binary_format
        width = _BINARY_FORMATS.get(binary_format.upper())
        if width:
            fields["bytesPerSample"] = width
    big_endian = _lookup(sections, "binary infos", "UseBigEndianOrder")
    if big_endian:
        fields["byteOrder"] = "big" if big_endian.upper() == "YES" else "little"
    elif binary_format:
        # The spec's default, and what every writer in the wild emits.
        fields["byteOrder"] = "little"
    for key, name in (
        ("DecimalSymbol", "decimalSymbol"), ("SkipLines", "skipLines"), ("SkipColumns", "skipColumns"),
    ):
        value = _lookup(sections, "ascii infos", key)
        if value:
            fields[name] = value

    names: List[str] = []
    channels: List[Dict[str, object]] = []
    resolutions = set()
    for key, value in sections.get("channel infos", ()):
        if not _CHANNEL.match(key):
            continue
        name, reference, resolution, unit = _split(value, 4)[:4]
        name = _unescape(name)
        names.append(name)
        entry: Dict[str, object] = {"name": name}
        if reference:
            entry["reference"] = _unescape(reference)
        value_resolution = _float(resolution)
        if value_resolution is not None:
            entry["resolution"] = value_resolution
            resolutions.add(value_resolution)
        if unit:
            entry["unit"] = _unescape(unit)
        channels.append(entry)
    if names:
        fields["channelNames"] = names[:_NAME_CAP]
    if channels:
        fields["channels"] = channels[:_CHANNEL_DETAIL]
        fields["channelsRead"] = len(channels)
    if len(resolutions) == 1:
        # No unit column means microvolts, which is the format's default.
        fields["resolution"] = resolutions.pop()
    if any(_CHANNEL.match(k) for k, _ in sections.get("coordinates", ())):
        fields["hasCoordinates"] = True
    if comment:
        fields["comment"] = comment[:1024]

    parts = ["BrainVision"]
    if channel_count is not None:
        parts.append(f"{channel_count} channels")
    elif channels:
        parts.append(f"{len(channels)} channels")
    if rate:
        parts.append(f"{_number(rate)} Hz")
    tail = [str(fields[k]) for k in ("binaryFormat", "dataOrientation") if k in fields]
    summary = ", ".join(parts)
    if tail:
        summary += ", " + " ".join(tail).lower()
    return DecodedHeader(format="brainvision", summary=summary, fields=fields)


def decode_brainvision_vmrk(data: bytes) -> Optional[DecodedHeader]:
    text = data.decode("utf-8", "replace")
    if not text.lstrip("﻿").startswith(_MARKER_MAGIC):
        return None
    sections, _ = _sections(text)

    fields: Dict[str, object] = {}
    version = _VERSION.search(text.splitlines()[0])
    if version:
        fields["brainVisionVersion"] = version.group(1)
    for key, name in (("DataFile", "dataFile"), ("Codepage", "codepage")):
        value = _lookup(sections, "common infos", key)
        if value:
            fields[name] = value

    markers: List[Dict[str, object]] = []
    types: Dict[str, int] = {}
    read = 0
    for key, value in sections.get("marker infos", ()):
        if not _MARKER.match(key):
            continue
        read += 1
        kind, description, position, size, channel, date = _split(value, 6)[:6]
        kind, description = _unescape(kind), _unescape(description)
        types[kind] = types.get(kind, 0) + 1
        entry: Dict[str, object] = {"type": kind}
        if description:
            entry["description"] = description
        for name, raw in (("position", position), ("size", size), ("channel", channel)):
            parsed = _int(raw)
            if parsed is not None:
                entry[name] = parsed
        if date and "segmentStartTime" not in fields:
            stamp = _SEGMENT_DATE.match(date)
            if stamp:
                y, mo, d, h, mi, s, micro = stamp.groups()
                iso = f"{y}-{mo}-{d}T{h}:{mi}:{s}"
                if micro:
                    iso += f".{micro}"
                fields["segmentStartTime"] = iso
        markers.append(entry)

    # The read is bounded, so this is what the leading bytes hold, not a file
    # total: a long marker file is cut off and says so.
    fields["markersRead"] = read
    fields["markersComplete"] = len(data) < BRAINVISION_READ_BYTES
    if markers:
        fields["markers"] = markers[:_MARKER_DETAIL]
    if types:
        ordered = sorted(types.items(), key=lambda kv: (-kv[1], kv[0]))
        fields["markerTypes"] = dict(ordered[:_MARKER_TYPE_CAP])

    noun = "marker" if read == 1 else "markers"
    summary = f"BrainVision markers, {read} {noun}"
    if not fields["markersComplete"]:
        summary += " read"
    if types:
        summary += ", " + " ".join(f"{k} x{v}" for k, v in ordered[:3])
    return DecodedHeader(format="brainvision-markers", summary=summary, fields=fields)


def _number(value: float) -> str:
    return str(int(value)) if value == int(value) else f"{value:g}"


register_decoder(["vhdr"], decode_brainvision_vhdr, read=Read(leading=BRAINVISION_READ_BYTES))
register_decoder(["vmrk"], decode_brainvision_vmrk, read=Read(leading=BRAINVISION_READ_BYTES))
