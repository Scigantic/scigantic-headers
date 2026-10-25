"""NRRD (Nearly Raw Raster Data) header decoding.

NRRD is the volume format of teem, ITK, 3D Slicer and the Allen Institute's
mouse brain atlases (the Common Coordinate Framework annotation and template
volumes are `.nrrd`). The file opens with a text header: a magic line
`NRRD000X` (X is the format version, 1 to 5), then one `field: value` line per
field and any number of `key:=value` pairs, ending at the first blank line. For
an attached file the raw, gzip, bzip2, ascii or hex encoded data follows that
blank line; a detached header (`.nhdr`) names its data file(s) instead and may
end without one.

This reads the header only. The data is never touched, so the encoding is
reported, not inflated. Spec: https://teem.sourceforge.net/nrrd/format.html
"""

import math
import re
from typing import Dict, List, Optional

from .decoders import DecodedHeader, Read, _finite, register_decoder

# 64 KiB: a plain volume header is a few hundred bytes, but 3D Slicer writes
# segmentation metadata as key/value pairs that run to tens of KB.
NRRD_READ_BYTES = 64 * 1024

_MAGIC = re.compile(rb"^NRRD000([1-5])\r?\n")

# The spec's type names (and their C-ish aliases) to a NumPy-style dtype name.
_TYPES = {
    "signed char": "int8", "int8": "int8", "int8_t": "int8",
    "uchar": "uint8", "unsigned char": "uint8", "uint8": "uint8", "uint8_t": "uint8",
    "short": "int16", "short int": "int16", "signed short": "int16", "signed short int": "int16",
    "int16": "int16", "int16_t": "int16",
    "ushort": "uint16", "unsigned short": "uint16", "unsigned short int": "uint16",
    "uint16": "uint16", "uint16_t": "uint16",
    "int": "int32", "signed int": "int32", "int32": "int32", "int32_t": "int32",
    "uint": "uint32", "unsigned int": "uint32", "uint32": "uint32", "uint32_t": "uint32",
    "longlong": "int64", "long long": "int64", "long long int": "int64", "signed long long": "int64",
    "signed long long int": "int64", "int64": "int64", "int64_t": "int64",
    "ulonglong": "uint64", "unsigned long long": "uint64", "unsigned long long int": "uint64",
    "uint64": "uint64", "uint64_t": "uint64",
    "float": "float32", "double": "float64", "block": "block",
}
_ENCODINGS = {
    "raw": "raw", "txt": "ascii", "text": "ascii", "ascii": "ascii", "hex": "hex",
    "gz": "gzip", "gzip": "gzip", "bz2": "bzip2", "bzip2": "bzip2",
}
# Fields whose value is one token per axis.
_PER_AXIS_INT = ("sizes",)
_PER_AXIS_FLOAT = ("spacings", "thicknesses", "axis mins", "axis maxs")
_PER_AXIS_WORD = ("kinds", "centers", "centerings")
_PER_AXIS_QUOTED = ("labels", "units", "space units")
_INT_FIELDS = ("dimension", "space dimension", "line skip", "lineskip", "byte skip", "byteskip")
_VECTOR_LIST = ("space directions", "measurement frame")
_VECTOR = ("space origin",)
_CANONICAL = {
    "lineskip": "line skip", "byteskip": "byte skip", "datafile": "data file",
    "centers": "centerings", "old min": "old min", "oldmin": "old min", "old max": "old max", "oldmax": "old max",
}


def _header_end(data: bytes) -> Optional[int]:
    """Offset just past the blank line that ends the header, or None when the
    read holds no blank line (a header longer than the read, or a detached
    header file that simply stops)."""
    best = None
    for sep in (b"\n\n", b"\r\n\r\n", b"\n\r\n"):
        i = data.find(sep)
        if i >= 0 and (best is None or i < best[0]):
            best = (i, len(sep))
    return best[0] + best[1] if best else None


def _vector(text: str) -> Optional[List[Optional[float]]]:
    """`(1,0,0)` to [1.0, 0.0, 0.0]; `none` to None."""
    t = text.strip()
    if t.lower() == "none":
        return None
    if not (t.startswith("(") and t.endswith(")")):
        return None
    out = []
    for part in t[1:-1].split(","):
        try:
            out.append(_finite(float(part)))
        except ValueError:
            out.append(None)
    return out


def _vector_list(text: str) -> List[Optional[List[Optional[float]]]]:
    toks = re.findall(r"\([^)]*\)|none", text, flags=re.IGNORECASE)
    return [_vector(t) for t in toks]


def _quoted(text: str) -> List[str]:
    return [m.group(1) if m.group(1) is not None else m.group(2)
            for m in re.finditer(r'"([^"]*)"|(\S+)', text)]


def _floats(text: str) -> List[Optional[float]]:
    out = []
    for tok in text.split():
        try:
            out.append(_finite(float(tok)))
        except ValueError:
            out.append(None)
    return out


def decode_nrrd(data: bytes) -> Optional[DecodedHeader]:
    m = _MAGIC.match(data)
    if not m:
        return None
    version = int(m.group(1))
    end = _header_end(data)
    complete = end is not None
    header = data[:end] if complete else data
    text = header.decode("utf-8", "replace")

    fields: Dict[str, object] = {}
    keyvals: Dict[str, str] = {}
    for line in text.splitlines()[1:]:
        line = line.rstrip("\r")
        if not line.strip() or line.startswith("#"):
            continue
        if ":=" in line:
            k, _, v = line.partition(":=")
            keyvals[k.strip()] = v.strip()
            continue
        if ": " not in line and not line.endswith(":"):
            continue
        k, _, v = line.partition(":")
        k = re.sub(r"\s+", " ", k.strip().lower())
        k = _CANONICAL.get(k, k)
        fields[k] = v.strip()

    dimension = None
    if "dimension" in fields:
        try:
            dimension = int(fields["dimension"])
        except ValueError:
            dimension = None
    sizes: List[int] = []
    if "sizes" in fields:
        try:
            sizes = [int(t) for t in fields["sizes"].split()]
        except ValueError:
            sizes = []
    if dimension is not None and not 1 <= dimension <= 16:
        return None
    if dimension is None and sizes:
        dimension = len(sizes)

    type_name = fields.get("type")
    dtype = _TYPES.get(re.sub(r"\s+", " ", type_name.lower())) if type_name else None
    if type_name and dtype is None:
        dtype = type_name
    encoding_name = fields.get("encoding")
    encoding = _ENCODINGS.get(encoding_name.lower(), encoding_name) if encoding_name else None

    directions = _vector_list(fields["space directions"]) if "space directions" in fields else None
    spacings = _floats(fields["spacings"]) if "spacings" in fields else None
    kinds = fields["kinds"].split() if "kinds" in fields else None
    # Per-axis spacing: the length of each space-direction vector when given
    # (axes with `none` are non-spatial, e.g. an RGB or vector axis), else the
    # `spacings` field. NRRD attaches no unit unless `space units` says so.
    spacing: Optional[List[Optional[float]]] = None
    if directions:
        spacing = []
        for vec in directions:
            if vec is None or any(c is None for c in vec):
                spacing.append(None)
            else:
                spacing.append(_finite(round(math.sqrt(sum(c * c for c in vec)), 6)))
    elif spacings:
        spacing = spacings

    data_files = None
    if "data file" in fields:
        data_files = fields["data file"]
    detached = data_files is not None

    n = None
    if sizes:
        n = 1
        for d in sizes:
            n *= d

    f: Dict[str, object] = {
        "nrrdVersion": version,
        "dimension": dimension,
        "shape": sizes or None,
        "dtype": dtype,
        "encoding": encoding,
        "endian": fields.get("endian"),
        "numElements": n,
        "space": fields.get("space"),
        "spaceDimension": int(fields["space dimension"]) if fields.get("space dimension", "").isdigit() else None,
        "spaceDirections": directions,
        "spaceOrigin": _vector(fields["space origin"]) if "space origin" in fields else None,
        "spaceUnits": _quoted(fields["space units"]) if "space units" in fields else None,
        "spacing": spacing,
        "kinds": kinds,
        "labels": _quoted(fields["labels"]) if "labels" in fields else None,
        "units": _quoted(fields["units"]) if "units" in fields else None,
        "centerings": fields["centerings"].split() if "centerings" in fields else None,
        "measurementFrame": _vector_list(fields["measurement frame"]) if "measurement frame" in fields else None,
        "content": fields.get("content"),
        "dataFile": data_files,
        "detached": detached,
        "byteSkip": int(fields["byte skip"]) if re.fullmatch(r"-?\d+", fields.get("byte skip", "")) else None,
        "lineSkip": int(fields["line skip"]) if fields.get("line skip", "").isdigit() else None,
        "keyValueCount": len(keyvals),
        "keyValues": dict(list(keyvals.items())[:64]),
        "headerComplete": complete,
        "headerBytes": end,
    }

    parts = ["NRRD"]
    if dimension:
        parts[0] += " %d-D" % dimension
    if sizes:
        parts.append("x".join(str(s) for s in sizes))
    if dtype:
        parts.append(dtype)
    if encoding:
        parts.append(encoding)
    if spacing and all(v is not None for v in spacing):
        unit = ""
        su = f["spaceUnits"]
        if su and len(set(su)) == 1:
            unit = " " + su[0]
        parts.append("spacing " + "x".join("%g" % v for v in spacing) + unit)
    if detached:
        parts.append("detached header")
    if not complete:
        parts.append("header truncated")
    return DecodedHeader(format="nrrd", summary=", ".join(parts), fields=f)


register_decoder(["nrrd", "nhdr"], decode_nrrd, read=Read(leading=NRRD_READ_BYTES))
