"""SpikeGLX run metadata (`.meta`) decoding.

SpikeGLX is the acquisition program for Neuropixels probes. Every binary it
writes has a `.meta` sidecar beside it: plain `key=value` lines, no sections,
holding the sample rate, the channel layout, the probe identity and the run's
size and duration. The binary itself is a headerless int16 matrix, so this
sidecar is the only description of it that exists.

`.meta` is a generic extension, so this declines anything that does not carry
SpikeGLX's own required keys (`nSavedChans` plus one of the sample-rate keys).

Keys prefixed with `~` hold a parenthesised table, one entry per channel
(`imroTbl`, `snsChanMap`, `snsGeomMap`). Those run to tens of KB and say
nothing a record needs, so only their entry counts are reported.

Reference: https://billkarsh.github.io/SpikeGLX/Sgl_help/Metadata_30.html
"""

import re
from typing import Dict, List, Optional

from .decoders import DecodedHeader, Read, _finite, register_decoder

# The per-channel tables dominate the file: a 384-channel snsGeomMap is ~20 KB,
# and a 4-shank run with every table present is larger again.
SPIKEGLX_READ_BYTES = 256 * 1024

_SAMPLE_RATE_KEYS = ("imSampRate", "niSampRate", "obSampRate")
# `<run>_g<gate>_t<trigger>`, which is how SpikeGLX names every file it writes.
_RUN = re.compile(r"^(?P<run>.+)_g(?P<gate>\d+)_t(?P<trigger>\d+|cat)$")
_TABLE_ENTRY = re.compile(r"\)\(")

# imDatPrb_type values documented in the metadata reference. Only the three
# that are unambiguous are named; any other type is reported as its number
# and its part number, which the probe table resolves exactly.
_PROBE_MODELS = {
    0: "Neuropixels 1.0",
    21: "Neuropixels 2.0 (1 shank)",
    24: "Neuropixels 2.0 (4 shank)",
}

# Verbatim keys worth carrying into a record, with the name they take there.
_STRINGS = (
    ("typeThis", "streamType"), ("appVersion", "appVersion"),
    ("fileName", "originalFileName"), ("fileCreateTime", "fileCreateTime"),
    ("gateMode", "gateMode"), ("trigMode", "trigMode"),
    ("imDatPrb_pn", "probePartNumber"), ("imDatPrb_sn", "probeSerial"),
    ("imDatApi", "imecApiVersion"), ("imDatBsc_pn", "basestationPartNumber"),
    ("imDatHs_pn", "headstagePartNumber"), ("snsSaveChanSubset", "savedChannelSubset"),
    ("imRoFile", "imroFile"), ("imCalibrated", "calibrated"),
)
_INTS = (
    ("nSavedChans", "channelCount"), ("fileSizeBytes", "fileSizeBytes"),
    ("firstSample", "firstSample"), ("imMaxInt", "maxInt"),
    ("imDatPrb_type", "probeType"), ("imDatPrb_slot", "probeSlot"),
    ("imDatPrb_port", "probePort"), ("imDatPrb_dock", "probeDock"),
)
_FLOATS = (
    ("fileTimeSecs", "durationSeconds"), ("imAiRangeMax", "analogRangeMaxVolts"),
    ("imAiRangeMin", "analogRangeMinVolts"), ("niAiRangeMax", "analogRangeMaxVolts"),
)


def _parse(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for line in text.splitlines():
        key, sep, value = line.partition("=")
        if not sep:
            continue
        out[key.strip()] = value.strip()
    return out


def _int(meta: Dict[str, str], key: str) -> Optional[int]:
    raw = meta.get(key)
    if raw is None:
        return None
    try:
        return int(raw)
    except ValueError:
        return None


def _float(meta: Dict[str, str], key: str) -> Optional[float]:
    raw = meta.get(key)
    if raw is None:
        return None
    try:
        return _finite(float(raw))
    except ValueError:
        return None


def _counts(meta: Dict[str, str], key: str, names: List[str]) -> Dict[str, object]:
    """A comma-separated channel-count triple (`snsApLfSy` = AP, LF, sync) into
    named fields."""
    raw = meta.get(key)
    if not raw:
        return {}
    parts = raw.split(",")
    if len(parts) != len(names):
        return {}
    out: Dict[str, object] = {}
    for name, part in zip(names, parts):
        try:
            out[name] = int(part)
        except ValueError:
            return {}
    return out


def _table_entries(value: str) -> Optional[int]:
    """Entries in a `~key=(header)(e1)(e2)...` table, not counting the leading
    header element."""
    if not value.startswith("(") or not value.endswith(")"):
        return None
    return len(_TABLE_ENTRY.findall(value))


def _stream_kind(name: str) -> Dict[str, object]:
    """Run name, gate and trigger index, device and stream, read back out of
    the original file name SpikeGLX recorded."""
    stem = name.replace("\\", "/").rsplit("/", 1)[-1]
    if stem.lower().endswith(".bin"):
        stem = stem[:-4]
    parts = stem.split(".")
    out: Dict[str, object] = {}
    if len(parts) > 1 and parts[-1] in ("ap", "lf"):
        out["stream"] = parts[-1]
        parts = parts[:-1]
    if len(parts) > 1 and (parts[-1].startswith("imec") or parts[-1].startswith("obx") or parts[-1] == "nidq"):
        out["device"] = parts[-1]
        parts = parts[:-1]
    run = _RUN.match(parts[0])
    if run:
        out["runName"] = run.group("run")
        out["gateIndex"] = int(run.group("gate"))
        if run.group("trigger") != "cat":
            out["triggerIndex"] = int(run.group("trigger"))
        else:
            out["concatenated"] = True
    elif parts[0]:
        out["runName"] = parts[0]
    return out


def _microvolts_per_bit(meta: Dict[str, str], fields: Dict[str, object]) -> Optional[float]:
    """Volts per ADC count scaled to microvolts, for the two cases the metadata
    reference defines without consulting the probe table.

    Neuropixels 1.0 divides the full-scale range by 512 counts and by the AP
    gain in the first `imroTbl` entry. Neuropixels 2.0 has a fixed gain and
    divides by `imMaxInt` instead. Any other probe is left alone rather than
    guessed at."""
    full_scale = fields.get("analogRangeMaxVolts")
    if not isinstance(full_scale, float) or full_scale <= 0:
        return None
    probe_type = fields.get("probeType")
    if probe_type == 0:
        gain = fields.get("apGain")
        if not isinstance(gain, (int, float)) or gain <= 0:
            return None
        return _finite(full_scale * 1e6 / 512.0 / float(gain))
    max_int = fields.get("maxInt")
    if probe_type in (21, 24) and isinstance(max_int, int) and max_int > 0:
        gain = _float(meta, "imChan0apGain") or 80.0
        if gain <= 0:
            return None
        return _finite(full_scale * 1e6 / max_int / gain)
    return None


def decode_spikeglx_meta(data: bytes) -> Optional[DecodedHeader]:
    meta = _parse(data.decode("latin-1"))
    if "nSavedChans" not in meta:
        return None
    rate_key = next((k for k in _SAMPLE_RATE_KEYS if k in meta), None)
    if rate_key is None:
        return None

    fields: Dict[str, object] = {}
    for key, name in _STRINGS:
        if meta.get(key):
            fields[name] = meta[key]
    for key, name in _INTS:
        value = _int(meta, key)
        if value is not None:
            fields[name] = value
    for key, name in _FLOATS:
        if name in fields:
            continue
        value = _float(meta, key)
        if value is not None:
            fields[name] = value

    rate = _float(meta, rate_key)
    if rate:
        fields["samplingRateHz"] = rate
    fields.update(_counts(meta, "snsApLfSy", ["apChannels", "lfChannels", "syncChannels"]))
    fields.update(_counts(meta, "snsXaDwSy", ["analogChannels", "digitalChannels", "syncChannels"]))
    fields.update(_counts(meta, "snsMnMaXaDw", ["multiplexedChannels", "auxAnalogChannels", "analogChannels", "digitalChannels"]))
    fields.update(_counts(meta, "acqApLfSy", ["acquiredApChannels", "acquiredLfChannels", "acquiredSyncChannels"]))

    channels = fields.get("channelCount")
    size = fields.get("fileSizeBytes")
    if isinstance(channels, int) and channels > 0 and isinstance(size, int) and size > 0:
        # The binary is int16, one sample per channel per tick, no header.
        fields["sampleCount"] = size // 2 // channels

    if isinstance(fields.get("probeType"), int):
        model = _PROBE_MODELS.get(fields["probeType"])
        if model:
            fields["probeModel"] = model

    for key in ("imroTbl", "snsChanMap", "snsGeomMap", "snsShankMap"):
        entries = _table_entries(meta.get("~" + key, ""))
        if entries:
            fields[key + "Entries"] = entries
    # The first imroTbl row is `chan bank refid apgain lfgain apfilt` on a 1.0
    # probe; a 2.0 row is shorter and carries no gain, which this skips.
    imro = meta.get("~imroTbl", "")
    if imro.startswith("("):
        rows = imro[1:-1].split(")(") if imro.endswith(")") else []
        if len(rows) > 1:
            row = rows[1].split(" ")
            if len(row) >= 6:
                try:
                    fields["apGain"] = int(row[3])
                    fields["lfGain"] = int(row[4])
                except ValueError:
                    pass

    uv = _microvolts_per_bit(meta, fields)
    if uv is not None:
        fields["microvoltsPerBit"] = uv
    if meta.get("fileName"):
        fields.update(_stream_kind(meta["fileName"]))
    fields["keyCount"] = len(meta)

    kind = fields.get("streamType") or "spikeglx"
    stream = fields.get("stream")
    label = f"SpikeGLX {kind}{'.' + str(stream) if stream else ''}"
    parts = [label]
    if isinstance(channels, int):
        parts.append(f"{channels} channels")
    if rate:
        parts.append(f"{_number(rate)} Hz")
    duration = fields.get("durationSeconds")
    if isinstance(duration, float):
        parts.append(f"{_number(duration)} s")
    model = fields.get("probeModel") or fields.get("probePartNumber")
    if model:
        parts.append(str(model))
    return DecodedHeader(format="spikeglx", summary=", ".join(parts), fields=fields)


def _number(value: float) -> str:
    return str(int(value)) if value == int(value) else f"{value:g}"


register_decoder(["meta"], decode_spikeglx_meta, read=Read(leading=SPIKEGLX_READ_BYTES))
