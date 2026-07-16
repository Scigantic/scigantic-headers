"""Read acquisition optics from a CryoSPARC .cs metadata file.

A .cs file is an .npy with a structured dtype, one record per particle, with
fields like 'blob/psize_A' (pixel size), 'ctf/accel_kv' (voltage), 'ctf/cs_mm'
(spherical aberration). The header gives the schema; the values are in the data,
so this reads the first record and extracts the optics fields. Same job as the
RELION STAR reader (star.py), different container.

Records are a flat C-packed buffer, so a field's value is at a computed offset.
That holds as long as no field is an object dtype ('|O'), object fields make
.npy store a pickle instead of a flat buffer, so this bails (returns {}) in that
case rather than reading garbage. CryoSPARC stores paths as fixed-width strings,
not objects, so real .cs files parse.

Pure stdlib. Verified against numpy-written structured arrays (numpy is the
reference for the .npy format).
"""

from __future__ import annotations

import math
import struct
from typing import Dict, Optional

from .decoders import npy_preamble

# Header + first record for typical CryoSPARC schemas (many fields, some
# fixed-string paths). Generous; the first record is right after the header.
CS_SCAN_BYTES = 64 * 1024

# CryoSPARC field name -> (output key, rounding). All are float32 in .cs.
_CS_OPTICS = {
    "blob/psize_A": ("pixelSizeA", 4),
    "ctf/accel_kv": ("voltageKv", 1),
    "ctf/cs_mm": ("sphericalAberrationMm", 3),
}

_NUM_STRUCT = {
    ("f", 2): "e", ("f", 4): "f", ("f", 8): "d",
    ("i", 1): "b", ("i", 2): "h", ("i", 4): "i", ("i", 8): "q",
    ("u", 1): "B", ("u", 2): "H", ("u", 4): "I", ("u", 8): "Q",
    ("b", 1): "?",
}


def _base(ts: str) -> str:
    return ts[1:] if ts[:1] in "<>|=" else ts


def _itemsize(ts) -> Optional[int]:
    """Bytes one element of typestr `ts` occupies, or None if it is an object
    dtype or otherwise not a flat fixed-size element."""
    if not isinstance(ts, str):
        return None
    s = _base(ts)
    if len(s) < 2 or not s[1:].isdigit():
        return None
    kind, size = s[0], int(s[1:])
    if kind == "O":
        return None                 # object -> pickled, not a flat buffer
    if kind == "U":
        return size * 4             # unicode: 4 bytes/char
    if kind in ("S", "a", "V", "f", "i", "u", "c", "b"):
        return size
    return None


def _num_fmt(ts: str) -> Optional[str]:
    s = _base(ts)
    ch = _NUM_STRUCT.get((s[0], int(s[1:]))) if len(s) >= 2 and s[1:].isdigit() else None
    if ch is None:
        return None
    return (">" if ts[:1] == ">" else "<") + ch


def _field_layout(descr):
    """{name: (offset, typestr)} plus the record size, or None if any field is
    object/unknown dtype (so the array is not a flat readable buffer)."""
    offset = 0
    layout = {}
    for spec in descr:
        try:
            name, ts = spec[0], spec[1]
            shape = spec[2] if len(spec) > 2 else ()
        except (TypeError, IndexError):
            return None
        base = _itemsize(ts)
        if base is None:
            return None
        count = 1
        for d in (shape if isinstance(shape, (tuple, list)) else (shape,)):
            count *= int(d)
        layout[name] = (offset, ts)
        offset += base * count
    return layout, offset


def parse_cryosparc_optics(data: bytes) -> Dict[str, float]:
    """Extract {pixelSizeA, voltageKv, sphericalAberrationMm} (those present)
    from the first record of a .cs buffer. {} if not a flat structured array or
    the fields aren't in the bytes we have."""
    pre = npy_preamble(data)
    if pre is None:
        return {}
    meta, data_off = pre
    descr = meta.get("descr")
    if not isinstance(descr, list):
        return {}
    lay = _field_layout(descr)
    if lay is None:
        return {}
    layout, _rec_size = lay

    out: Dict[str, float] = {}
    for cs_name, (key, ndigits) in _CS_OPTICS.items():
        if cs_name not in layout:
            continue
        foff, ts = layout[cs_name]
        fmt = _num_fmt(ts)
        size = _itemsize(ts)
        if fmt is None or size is None:
            continue
        pos = data_off + foff
        if pos + size > len(data):
            continue                 # field lies past the bytes we read
        try:
            (val,) = struct.unpack_from(fmt, data, pos)
        except struct.error:
            continue
        if isinstance(val, float) and not math.isfinite(val):
            continue
        out[key] = round(float(val), ndigits)
    return out


def read_cryosparc_optics(path: str, max_bytes: int = CS_SCAN_BYTES) -> Dict[str, float]:
    """Read a .cs file's leading bytes (header + first record) and parse its
    optics. {} if unreadable or not a flat structured array."""
    try:
        with open(path, "rb") as fh:
            data = fh.read(max_bytes)
    except OSError:
        return {}
    return parse_cryosparc_optics(data)
