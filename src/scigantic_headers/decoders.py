"""Pure header decoders, bytes in, typed fields out. No I/O, no dependencies.

Given the leading bytes of a scientific file, decode a self-contained binary
header into typed fields plus a one-line summary an agent can read instead of
re-opening the file. This module never touches the filesystem or the network;
reading bytes is the job of `sources.py`. Keeping the decode pure makes it
trivially testable, keeps the hot path allocation-light, and lets the
TypeScript twin (backend headerDecoders.ts) mirror it exactly.

Zero runtime dependencies is deliberate: the dtype table is a plain dict rather
than numpy, so importing this costs nothing and it installs on an air-gapped
box. (The older scigantic_empiar.parse_mrc_header imported numpy solely to size
a dtype, this drops that.)

Add a format by writing a pure `bytes -> DecodedHeader | None` function and
calling `register_decoder`.
"""

from __future__ import annotations

import ast
import math
import struct
from dataclasses import asdict, dataclass, field
from typing import Callable, Dict, Optional


def _finite(x: Optional[float]) -> Optional[float]:
    """A float that is safe to put in a record, or None. Garbage bytes can
    decode into NaN/inf; those are not valid JSON (json emits a bare `NaN`
    token that strict parsers reject) and are meaningless as a measurement, so
    they become None."""
    return x if x is not None and math.isfinite(x) else None

# Leading bytes any registered decoder may inspect. Readers should fetch exactly
# this many and no more, for MRC everything needed sits below byte 212, and no
# other format we support needs more than 1 KiB.
HEADER_BYTES = 1024


@dataclass(frozen=True)
class DecodedHeader:
    """Typed result of decoding a binary header. `fields` is format-specific."""

    format: str
    summary: str
    fields: Dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


# A decoder inspects the leading bytes and returns typed fields, or None if the
# bytes do not validate for its format.
Decoder = Callable[[bytes], Optional[DecodedHeader]]

_DECODERS_BY_EXT: Dict[str, Decoder] = {}


def register_decoder(extensions, decoder: Decoder) -> None:
    """Register `decoder` for one or more lower-case extensions (no dot).

    Extending to a new format is a pure function plus one call to this, nothing
    in the dispatch, the readers, or the batch path is format-specific.
    """
    if isinstance(extensions, str):
        extensions = [extensions]
    for ext in extensions:
        _DECODERS_BY_EXT[ext.lower().lstrip(".")] = decoder


def extension_of(key: str) -> str:
    """Lower-cased final extension of a path/key, without the dot. '' if none."""
    base = key.rsplit("/", 1)[-1]
    dot = base.rfind(".")
    return base[dot + 1 :].lower() if dot > 0 else ""


def has_decoder_for(key: str) -> bool:
    """True if a decoder is registered for this key's extension. Cheap: no I/O,
    so callers filter on this before opening a file."""
    return extension_of(key) in _DECODERS_BY_EXT


def decode_bytes(key: str, data: bytes) -> Optional[DecodedHeader]:
    """Decode `data` (the leading bytes of `key`) if a decoder is registered and
    the bytes validate. Returns None when there is no decoder for the extension
    or the header fails validation (wrong magic, insane values, truncated).
    Callers treat None as "no context", never as an error."""
    decoder = _DECODERS_BY_EXT.get(extension_of(key))
    if decoder is None:
        return None
    try:
        return decoder(data)
    except Exception:
        return None


# ── MRC / MRCS (cryo-EM maps, micrographs, tilt series) ────────────────────
#
# MRC2014 layout, little-endian (the 'MAP ' stamp at byte 208 and MACHST at 212
# record endianness; effectively all cryo-EM data is little-endian). Offsets and
# the mode->dtype table match backend headerDecoders.ts and the reference
# parse_mrc_header in scigantic_empiar; all three are held to the golden file
# tests/fixtures/mrc-cases.json so they cannot silently diverge.
#
#   0   NX,NY,NZ    int32   columns, rows, sections
#   12  MODE        int32   pixel data type
#   28  MX          int32   grid sampling (X)
#   40  CELLA_X     float32  cell dimension X, angstrom
#   92  NSYMBT      int32   extended-header bytes after the 1024 base
#   208 'MAP '      char[4]  format stamp (validity check)
#
# Everything read here sits below byte 212, so 512 leading bytes already suffice.

_MRC_MODE_DTYPE = {
    0: "int8",
    1: "int16",
    2: "float32",
    3: "complex_int16",
    4: "complex_float32",
    6: "uint16",
    12: "float16",
}
_MRC_DTYPE_BYTES = {
    "int8": 1, "int16": 2, "float32": 4, "complex_int16": 4,
    "complex_float32": 8, "uint16": 2, "float16": 2,
}
_MRC_MAGIC_OFFSET = 208
_MRC_STAMP_MIN_BYTES = _MRC_MAGIC_OFFSET + 4  # need this many to read the stamp
_MRC_FIELDS_MIN_BYTES = 96  # last numeric field (NSYMBT) ends at byte 96

# One reusable compiled struct for the leading NX/NY/NZ/MODE block, compiling
# the format string once is marginal here (the read dominates), but it keeps the
# hot path free of per-call format parsing.
_MRC_DIMS = struct.Struct("<4i")


def decode_mrc_header(data: bytes, *, strict: bool = True) -> Optional[DecodedHeader]:
    """Decode an MRC2014 header.

    strict=True (default, used by the registry/dispatch): require the 'MAP '
    stamp, so a random binary file with an .mrc extension is rejected rather
    than decoded into garbage. strict=False: skip the stamp check for a caller
    that already knows the bytes are MRC, matches the historical permissive
    behavior of scigantic_empiar.parse_mrc_header, which read pre-2014 files
    that predate the stamp.
    """
    min_bytes = _MRC_STAMP_MIN_BYTES if strict else _MRC_FIELDS_MIN_BYTES
    if len(data) < min_bytes:
        return None
    # Cheapest reject first: without the stamp a random binary header decodes
    # into plausible-looking garbage. Skipped when the caller vouches for the
    # format (strict=False).
    if strict and data[_MRC_MAGIC_OFFSET:_MRC_MAGIC_OFFSET + 4] != b"MAP ":
        return None

    nx, ny, nz, mode = _MRC_DIMS.unpack_from(data, 0)
    (mx,) = struct.unpack_from("<i", data, 28)
    (cella_x,) = struct.unpack_from("<f", data, 40)
    (nsymbt,) = struct.unpack_from("<i", data, 92)

    # Guard against a stamp match on non-image bytes. Only under strict dispatch,
    # a vouched caller (strict=False) gets the historical no-validation behavior.
    if strict and not (0 < nx <= 1_000_000 and 0 < ny <= 1_000_000 and 0 < nz <= 10_000_000):
        return None

    dtype = _MRC_MODE_DTYPE.get(mode, "float32")
    bytes_per_voxel = _MRC_DTYPE_BYTES.get(dtype, 4)
    # Pixel size is derivable only when both grid sampling and cell size are set.
    # Raw movie stacks often ship CELLA = 0 (pixel size then lives in the
    # acquisition metadata / STAR, not the header), report unknown, not 0.
    # _finite guards against a garbage CELLA decoding to inf/NaN.
    pixel_size_a = _finite(round(cella_x / mx, 3)) if mx > 0 and cella_x > 0 else None
    is_stack = nz > 1

    dims = f"{nx}x{ny}x{nz}" if is_stack else f"{nx}x{ny}"
    apix = f", {pixel_size_a} A/px" if pixel_size_a is not None else ""
    summary = f"MRC {'stack' if is_stack else 'image'} {dims}, {dtype}{apix}"

    return DecodedHeader(
        format="mrc",
        summary=summary,
        fields={
            "nx": nx, "ny": ny, "nz": nz, "mode": mode, "dtype": dtype,
            "pixelSizeA": pixel_size_a,
            "frameBytes": nx * ny * bytes_per_voxel,
            "dataOffset": 1024 + nsymbt,
            "isStack": is_stack,
        },
    )


register_decoder(["mrc", "mrcs", "st", "ali", "rec", "map"], decode_mrc_header)


# ── NPY (NumPy array, ML, genomics, materials; not cryo-EM) ───────────────
#
# .npy layout (numpy format 1.0/2.0/3.0):
#   0   \x93NUMPY        6-byte magic
#   6   major, minor     two bytes
#   8   HEADER_LEN       uint16 (v1) or uint32 (v2/3), little-endian
#   10/12  header dict   ASCII/UTF-8 Python literal, e.g.
#                        {'descr': '<f8', 'fortran_order': False, 'shape': (100, 50), }
# The dict is padded to align the whole preamble to 64 bytes, so it always fits
# in the leading read. Parsed with ast.literal_eval (safe, literals only).

_NPY_MAGIC = b"\x93NUMPY"
_NPY_KIND = {"f": "float", "i": "int", "u": "uint", "c": "complex", "b": "bool"}


def _npy_dtype_name(descr) -> str:
    # descr like '<f8', '>i4', '|u1', '|b1'. Byte order is irrelevant to the name.
    if not isinstance(descr, str) or len(descr) < 2:
        return "structured" if isinstance(descr, list) else str(descr)
    s = descr[1:] if descr[0] in "<>|=" else descr
    kind, size = s[0], s[1:]
    if kind in _NPY_KIND and size.isdigit():
        return "bool" if kind == "b" else f"{_NPY_KIND[kind]}{int(size) * 8}"
    return descr


def npy_preamble(data: bytes):
    """Parse an .npy preamble into (meta_dict, data_offset), or None. meta has
    'descr'/'shape'/'fortran_order'; data_offset is where the array bytes begin
    (past the padded header). Shared by the NPY and CryoSPARC .cs paths, a .cs
    file is an .npy with a structured descr."""
    if len(data) < 10 or data[:6] != _NPY_MAGIC:
        return None
    major = data[6]
    if major == 1:
        (hlen,) = struct.unpack_from("<H", data, 8)
        start = 10
    else:  # 2.0 / 3.0 use a 4-byte length
        if len(data) < 12:
            return None
        (hlen,) = struct.unpack_from("<I", data, 8)
        start = 12
    if start + hlen > len(data):  # header longer than the bounded read (never in practice)
        return None
    try:
        meta = ast.literal_eval(data[start:start + hlen].decode("latin1"))
    except Exception:
        return None
    if not isinstance(meta, dict) or "shape" not in meta or "descr" not in meta:
        return None
    return meta, start + hlen


def decode_npy_header(data: bytes) -> Optional[DecodedHeader]:
    pre = npy_preamble(data)
    if pre is None:
        return None
    meta, _data_offset = pre

    # numpy always writes shape as a tuple of ints; anything else is not a real
    # .npy header, so reject rather than raise on the conversion below.
    raw_shape = meta["shape"]
    if not isinstance(raw_shape, tuple):
        return None
    try:
        shape = [int(d) for d in raw_shape]
    except (TypeError, ValueError):
        return None
    if any(d < 0 for d in shape):
        return None

    dtype = _npy_dtype_name(meta["descr"])
    fortran = bool(meta.get("fortran_order", False))
    n = 1
    for d in shape:
        n *= d
    dims = "x".join(str(d) for d in shape) if shape else "scalar"
    order = ", Fortran-order" if fortran else ""
    return DecodedHeader(
        format="npy",
        summary=f"NPY array {dims}, {dtype}{order}",
        fields={
            "shape": shape, "dtype": dtype, "fortranOrder": fortran,
            "ndim": len(shape), "numElements": n,
        },
    )


register_decoder("npy", decode_npy_header)


# ── NIfTI-1 (neuroimaging volumes, MRI/fMRI; not cryo-EM) ──────────────────
#
# 348-byte fixed binary header. Endianness is detected by reading sizeof_hdr:
# it must equal 348 in the file's byte order. Validated by the magic at byte 344
# ('n+1\0' single-file, 'ni1\0' header/image pair). Only raw .nii is handled;
# .nii.gz is gzip-wrapped and, like BAM, needs a decompress step first.
#
#   0    sizeof_hdr   int32   == 348
#   40   dim[8]       int16   dim[0]=ndim, dim[1..]=sizes
#   70   datatype     int16   NIfTI datatype code
#   72   bitpix       int16   bits per voxel
#   76   pixdim[8]    float32 pixdim[1..3]=voxel sizes in mm
#   344  magic        char[4] 'n+1\0' | 'ni1\0'

_NIFTI_DTYPE = {
    2: "uint8", 4: "int16", 8: "int32", 16: "float32", 32: "complex64",
    64: "float64", 256: "int8", 512: "uint16", 768: "uint32",
    1024: "int64", 1280: "uint64", 1792: "complex128",
}
_NIFTI_MAGIC = (b"n+1\x00", b"ni1\x00")


def decode_nifti_header(data: bytes) -> Optional[DecodedHeader]:
    if len(data) < 348:
        return None
    # Detect byte order: sizeof_hdr reads as 348 only in the correct one.
    order = None
    for bo in ("<", ">"):
        if struct.unpack_from(bo + "i", data, 0)[0] == 348:
            order = bo
            break
    if order is None or data[344:348] not in _NIFTI_MAGIC:
        return None

    dim = struct.unpack_from(order + "8h", data, 40)
    (datatype,) = struct.unpack_from(order + "h", data, 70)
    pixdim = struct.unpack_from(order + "8f", data, 76)
    ndim = dim[0]
    if not 1 <= ndim <= 7:
        return None

    shape = list(dim[1:1 + ndim])
    # _finite guards against uninitialized / garbage pixdim decoding to NaN/inf.
    voxel = [_finite(round(pixdim[i], 4)) for i in range(1, min(1 + ndim, 4))]
    dtype = _NIFTI_DTYPE.get(datatype, f"code{datatype}")

    dims = "x".join(str(s) for s in shape)
    # Only show the voxel spacing in the summary when it is fully known.
    mm = "x".join(str(v) for v in voxel) if all(v is not None for v in voxel) else None
    tail = f", {mm} mm" if mm else ""
    return DecodedHeader(
        format="nifti",
        summary=f"NIfTI-1 volume {dims}, {dtype}{tail}",
        fields={
            "ndim": ndim, "shape": shape, "dtype": dtype,
            "voxelSizesMm": voxel, "datatypeCode": datatype,
        },
    )


register_decoder("nii", decode_nifti_header)


# ── CryoSPARC .cs (metadata dataset, cryo-EM) ─────────────────────────────
#
# A .cs file is an .npy with a *structured* dtype: one record per particle /
# micrograph, fields like 'blob/psize_A' (pixel size), 'ctf/accel_kv' (voltage),
# 'ctf/cs_mm'. The header alone gives the record count and the field schema,
# enough to identify the file and what it carries. Reading the actual optics
# *values* means reading a data row; that lives in cryosparc.py, next to the
# STAR reader, since it is the same "recover the pixel size" job.

def decode_cryosparc_header(data: bytes) -> Optional[DecodedHeader]:
    pre = npy_preamble(data)
    if pre is None:
        return None
    meta, _ = pre
    descr = meta["descr"]
    # A CryoSPARC dataset is always a structured array (list of field specs).
    if not isinstance(descr, list) or not descr:
        return None
    try:
        fields = [str(spec[0]) for spec in descr]
    except (TypeError, IndexError):
        return None
    shape = meta.get("shape")
    num_records = int(shape[0]) if isinstance(shape, tuple) and shape else None

    n = f"{num_records} records" if num_records is not None else "records"
    return DecodedHeader(
        format="cryosparc",
        summary=f"CryoSPARC dataset, {n}, {len(fields)} fields",
        fields={
            "numRecords": num_records,
            "numFields": len(fields),
            # A preview of the schema, enough to see the optics fields exist.
            "fields": fields[:40],
        },
    )


register_decoder("cs", decode_cryosparc_header)
