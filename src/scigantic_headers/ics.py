"""Decode an Image Cytometry Standard (ICS) header.

ICS (Dean et al. 1990; libics) is a text header, one keyword per line with
tab-separated fields, that describes a raw binary image. Version 1.0 keeps the
pixels in a sibling `.ids` file; version 2.0 appends them to the `.ics` after
the `end` line. The header is small (a few KB), so a bounded leading read
holds all of it.

    <tab><newline>                              separator declaration (first line)
    ics_version    1.0 | 2.0
    filename       <stem>
    layout         parameters   N              number of dimensions incl. "bits"
    layout         order        bits x y z ... dimension names, bits first
    layout         sizes        16 512 512 30  one per order entry
    layout         significant_bits  16
    representation format       integer | real | complex
    representation sign         signed | unsigned
    representation byte_order   1 2            1 2 = little-endian
    representation compression  uncompressed | gzip
    parameter      scale        1.0 0.2 0.2 0.5  sampling, one per order entry
    parameter      units        bits micrometers micrometers micrometers
    parameter      labels       intensity x y z
    sensor         type / model / s_params <name> <value>
    history        <key> <value>               software, date, free text
    end

The dimension named `bits` is the sample width, not an image axis; sizeX..T
come from the x, y, z, ch (or c), t entries. Units are converted to
micrometres for voxelSizeUm when they name a length (micrometers, nm, mm, m).
Huygens (SVI) is the common writer; its sensor block carries the objective NA
and wavelengths, which are kept under `sensorParams`.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional

from .decoders import DecodedHeader, Read, register_decoder

_TO_UM = {
    "micrometers": 1.0, "micrometer": 1.0, "microns": 1.0, "micron": 1.0, "um": 1.0, "µm": 1.0,
    "nanometers": 1e-3, "nanometer": 1e-3, "nm": 1e-3,
    "millimeters": 1e3, "millimeter": 1e3, "mm": 1e3,
    "meters": 1e6, "meter": 1e6, "m": 1e6,
}
_AXIS = {"x": "X", "y": "Y", "z": "Z", "ch": "C", "c": "C", "t": "T", "time": "T"}
_MAX_HISTORY = 40
_MAX_SENSOR = 64


def _finite(text) -> Optional[float]:
    try:
        v = float(text)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _int(text) -> Optional[int]:
    try:
        v = int(text)
    except (TypeError, ValueError):
        return None
    return v if 0 <= v < 1 << 31 else None


def decode_ics(data: bytes) -> Optional[DecodedHeader]:
    head = data[:4096]
    if b"ics_version" not in head and b"layout" not in head:
        return None
    text = data.decode("latin-1")
    # A v2 file carries binary pixels after "end": stop there.
    lines = text.split("\n")
    layout: Dict[str, List[str]] = {}
    repr_: Dict[str, List[str]] = {}
    param: Dict[str, List[str]] = {}
    history: Dict[str, str] = {}
    sensor: Dict[str, str] = {}
    sensor_params: Dict[str, str] = {}
    history_count = 0
    version = filename = None
    saw_end = False
    for raw in lines[:4000]:
        line = raw.rstrip("\r")
        if line == "end":
            saw_end = True
            break
        parts = line.split("\t")
        key = parts[0].strip().lower()
        if not key:
            continue
        rest = [p.strip() for p in parts[1:]]
        if key == "ics_version" and rest:
            version = rest[0]
        elif key == "filename" and rest:
            filename = rest[0]
        elif key == "layout" and len(rest) >= 2:
            layout[rest[0].lower()] = rest[1:]
        elif key == "representation" and len(rest) >= 2:
            repr_[rest[0].lower()] = rest[1:]
        elif key == "parameter" and len(rest) >= 2:
            param[rest[0].lower()] = rest[1:]
        elif key == "history" and rest:
            history_count += 1
            if len(history) < _MAX_HISTORY:
                k = rest[0] if len(rest) > 1 else "note%d" % history_count
                v = "\t".join(rest[1:]) if len(rest) > 1 else rest[0]
                if k in history:
                    k = "%s#%d" % (k, history_count)
                history[k] = v
        elif key == "sensor" and len(rest) >= 2:
            if rest[0].lower() in ("s_params",) and len(rest) >= 3:
                if len(sensor_params) < _MAX_SENSOR:
                    sensor_params[rest[1]] = "\t".join(rest[2:])
            elif rest[0].lower() in ("type", "model"):
                sensor[rest[0].lower()] = " ".join(rest[1:])
    if version is None and not layout:
        return None

    order = [o.lower() for o in layout.get("order", [])]
    sizes = [_int(s) for s in layout.get("sizes", [])]
    scale = [_finite(s) for s in param.get("scale", [])]
    units = param.get("units", [])
    labels = param.get("labels", [])

    fields: Dict[str, object] = {
        "icsVersion": version, "filename": filename,
        "sizeX": None, "sizeY": None, "sizeZ": None, "sizeC": None, "sizeT": None,
    }
    bits = None
    voxel: List[Optional[float]] = [None, None, None]
    axes = []
    for i, name in enumerate(order):
        size = sizes[i] if i < len(sizes) else None
        if name == "bits":
            bits = size
            continue
        axes.append(name)
        ax = _AXIS.get(name)
        if ax is not None:
            fields["size" + ax] = size
        if name in ("x", "y", "z"):
            sc = scale[i] if i < len(scale) else None
            unit = units[i].lower() if i < len(units) else ""
            factor = _TO_UM.get(unit)
            if sc is not None and factor is not None:
                voxel["xyz".index(name)] = round(sc * factor, 6)

    fmt = (repr_.get("format") or [""])[0].lower()
    sign = (repr_.get("sign") or [""])[0].lower()
    dtype = None
    if bits:
        if fmt == "real":
            dtype = "float%d" % bits
        elif fmt == "complex":
            dtype = "complex%d" % bits
        elif fmt == "integer" or fmt == "":
            dtype = ("int%d" if sign == "signed" else "uint%d") % bits
    byte_order = repr_.get("byte_order")
    if byte_order:
        byte_order = "little" if byte_order[0] == "1" else "big" if byte_order[0].isdigit() else None
    else:
        byte_order = None

    fields.update({
        "dtype": dtype,
        "bits": bits,
        "significantBits": _int((layout.get("significant_bits") or [None])[0]),
        "axes": axes,
        "sizes": [s for n, s in zip(order, sizes) if n != "bits"],
        "byteOrder": byte_order,
        "compression": (repr_.get("compression") or [None])[0],
        "voxelSizeUm": voxel,
        "units": units,
        "labels": labels,
        "sensorType": sensor.get("type"),
        "sensorModel": sensor.get("model"),
        "sensorParams": sensor_params,
        "historyCount": history_count,
        "history": history,
        "hasInlineData": saw_end and (version or "").startswith("2"),
        "software": history.get("software"),
    })
    na = _finite(sensor_params.get("NumAperture"))
    if na is not None:
        fields["numericalAperture"] = na

    dims = [fields["sizeX"], fields["sizeY"]]
    if fields["sizeZ"] and fields["sizeZ"] > 1:
        dims.append(fields["sizeZ"])
    summary = "ICS %s image" % (version or "?")
    if all(d is not None for d in dims):
        summary += " %s" % "x".join(str(d) for d in dims)
    if fields["sizeC"] and fields["sizeC"] > 1:
        summary += ", %d channels" % fields["sizeC"]
    if fields["sizeT"] and fields["sizeT"] > 1:
        summary += ", %d timepoints" % fields["sizeT"]
    if dtype:
        summary += ", " + dtype
    if all(v is not None for v in voxel[:2]):
        summary += ", %s um" % "x".join("%g" % v for v in voxel if v is not None)
    return DecodedHeader(format="ics", summary=summary, fields=fields)


register_decoder("ics", decode_ics, read=Read(leading=256 * 1024))
