"""Open Ephys binary-format recording manifest (`structure.oebin`) decoding.

The Open Ephys GUI's binary format writes each recording as a directory: raw
`continuous.dat` matrices, event and spike directories, and one JSON manifest
named `structure.oebin` that describes all of them. The matrices carry no
header of their own, so the manifest is the only description of the sample
rate, the channel count and the per-channel scaling that exists.

One manifest covers every stream a record node wrote, so a Neuropixels session
lists AP and LF side by side with the event lines. Per-channel entries run to
hundreds of KB; this reports what is uniform across a stream (scaling, units)
and the counts, not the per-channel table.

Format: https://open-ephys.github.io/gui-docs/User-Manual/Recording-data/Binary-format.html
"""

import json
from typing import Dict, List, Optional

from .decoders import DecodedHeader, Read, _finite, register_decoder

# A 384-channel stream is ~100 KB of per-channel entries; a record node with
# several probes is several times that. A manifest past this cannot be parsed
# as JSON at all, so the decoder declines rather than describe half of one.
OPENEPHYS_READ_BYTES = 2 * 1024 * 1024

_STREAM_DETAIL = 16
_EVENT_DETAIL = 16


def _number(value: object) -> Optional[float]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return _finite(float(value))


def _uniform(entries: List[dict], key: str) -> Optional[object]:
    """The one value every channel in a stream agrees on, or None when they
    differ. Reporting a mixed scaling as if it were a single number would be
    worse than reporting none."""
    seen = {e.get(key) for e in entries if isinstance(e, dict) and e.get(key) not in (None, "")}
    return seen.pop() if len(seen) == 1 else None


def _continuous(entry: dict) -> Dict[str, object]:
    out: Dict[str, object] = {}
    for key, name in (
        ("folder_name", "folderName"), ("stream_name", "streamName"),
        ("source_processor_name", "sourceProcessorName"),
        ("recorded_processor", "recordedProcessor"),
    ):
        value = entry.get(key)
        if isinstance(value, str) and value:
            out[name] = value
    for key, name in (
        ("source_processor_id", "sourceProcessorId"),
        ("recorded_processor_id", "recordedProcessorId"),
        ("num_channels", "channelCount"),
    ):
        value = entry.get(key)
        if isinstance(value, int) and not isinstance(value, bool):
            out[name] = value
    rate = _number(entry.get("sample_rate"))
    if rate is not None:
        out["samplingRateHz"] = rate
    channels = [c for c in entry.get("channels", []) if isinstance(c, dict)]
    if channels:
        out["channelsListed"] = len(channels)
        bit_volts = _number(_uniform(channels, "bit_volts"))
        if bit_volts is not None:
            out["bitVolts"] = bit_volts
        units = _uniform(channels, "units")
        if isinstance(units, str):
            out["units"] = units
        first = channels[0].get("channel_name")
        if isinstance(first, str) and first:
            out["firstChannelName"] = first
    return out


def _event(entry: dict) -> Dict[str, object]:
    out: Dict[str, object] = {}
    for key, name in (
        ("folder_name", "folderName"), ("channel_name", "channelName"),
        ("type", "type"), ("stream_name", "streamName"),
        ("source_processor", "sourceProcessor"), ("description", "description"),
    ):
        value = entry.get(key)
        if isinstance(value, str) and value:
            out[name] = value
    for key, name in (("num_channels", "channelCount"),):
        value = entry.get(key)
        if isinstance(value, int) and not isinstance(value, bool):
            out[name] = value
    rate = _number(entry.get("sample_rate"))
    if rate is not None:
        out["samplingRateHz"] = rate
    return out


def decode_openephys_oebin(data: bytes) -> Optional[DecodedHeader]:
    try:
        manifest = json.loads(data.decode("utf-8", "replace"))
    except ValueError:
        return None
    if not isinstance(manifest, dict):
        return None
    # A manifest always declares at least one of the three record classes, and
    # nothing else writes this shape, so requiring one keeps a stray JSON out.
    if not any(isinstance(manifest.get(k), list) for k in ("continuous", "events", "spikes")):
        return None

    fields: Dict[str, object] = {}
    version = manifest.get("GUI version")
    if isinstance(version, str) and version:
        fields["guiVersion"] = version

    continuous = [e for e in manifest.get("continuous", []) if isinstance(e, dict)]
    events = [e for e in manifest.get("events", []) if isinstance(e, dict)]
    spikes = [e for e in manifest.get("spikes", []) if isinstance(e, dict)]
    if not (continuous or events or spikes) and "GUI version" not in manifest:
        return None

    fields["continuousStreams"] = len(continuous)
    fields["eventStreams"] = len(events)
    fields["spikeStreams"] = len(spikes)

    described = [_continuous(e) for e in continuous]
    if described:
        fields["continuous"] = described[:_STREAM_DETAIL]
    total = sum(int(d["channelCount"]) for d in described if isinstance(d.get("channelCount"), int))
    if total:
        fields["totalContinuousChannels"] = total
    rates = {d["samplingRateHz"] for d in described if isinstance(d.get("samplingRateHz"), float)}
    if len(rates) == 1:
        fields["samplingRateHz"] = rates.pop()
    elif rates:
        fields["samplingRatesHz"] = sorted(rates)
    names = [d["streamName"] for d in described if isinstance(d.get("streamName"), str)]
    if names:
        fields["streamNames"] = names[:_STREAM_DETAIL]
    processors = sorted({
        d["sourceProcessorName"] for d in described if isinstance(d.get("sourceProcessorName"), str)
    })
    if processors:
        fields["sourceProcessors"] = processors[:_STREAM_DETAIL]
    if events:
        fields["events"] = [_event(e) for e in events][:_EVENT_DETAIL]
    if spikes:
        fields["spikes"] = [_event(e) for e in spikes][:_EVENT_DETAIL]

    parts = ["Open Ephys recording"]
    if continuous:
        plural = "" if len(continuous) == 1 else "s"
        parts.append(f"{len(continuous)} continuous stream{plural}")
    if total:
        parts.append(f"{total} channels")
    rate = fields.get("samplingRateHz")
    if isinstance(rate, float):
        parts.append(f"{str(int(rate)) if rate == int(rate) else f'{rate:g}'} Hz")
    if events:
        plural = "" if len(events) == 1 else "s"
        parts.append(f"{len(events)} event stream{plural}")
    return DecodedHeader(format="openephys-oebin", summary=", ".join(parts), fields=fields)


register_decoder(["oebin"], decode_openephys_oebin, read=Read(leading=OPENEPHYS_READ_BYTES))
