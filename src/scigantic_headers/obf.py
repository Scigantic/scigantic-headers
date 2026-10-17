"""Decode an Imspector OBF / MSR (Abberior, LaVision BioTec ImSpector) header.

OBF is a linked list of stacks behind a small file header; an MSR measurement
file embeds an OBF (the OBF header is at byte 0 and Imspector's own window and
hardware records live between the header and the first stack). Little-endian,
packed structs. Verified against the Imspector file-format documentation and
the obffile reference reader on real .msr measurements.

File header (26 bytes):
    0   char[10]  "OMAS_BF\\n" + 0xff 0xff
    10  uint32    format_version
    14  uint64    first_stack_pos
    22  uint32    descr_len            then descr_len bytes of UTF-8 (XML)
    (version >= 2) uint64 meta_data_position, right after the description

Stack header (368 bytes) at first_stack_pos, then at each next_stack_pos:
    0    char[16]    "OMAS_BF_STACK\\n" + 0xff 0xff
    16   uint32      format_version
    20   uint32      rank                  valid entries in the 15-slot arrays
    24   uint32[15]  res                   samples per axis
    84   double[15]  len                   physical length per axis (metres for
                                           spatial axes, see the footer SI units)
    204  double[15]  off                   physical offset per axis
    324  uint32      dt                    data type flags (below)
    328  uint32      compression_type      0 none, 1 zlib
    332  uint32      compression_level
    336  uint32      name_len
    340  uint32      descr_len
    344  uint64      reserved
    352  uint64      data_len_disk
    360  uint64      next_stack_pos        0 after the last stack
    368  name, description, then data_len_disk bytes of data, then the footer

Stack footer (format_version >= 1) at 368 + name_len + descr_len + data_len_disk:
    0    uint32      size                  bytes to skip to reach the labels
    4    uint32[15]  has_col_positions
    64   uint32[15]  has_col_labels
    124  uint32      metadata_length
    128  SI unit     si_value              (v2+) 9 x (int32 num, int32 den) + double
    208  SI unit[15] si_dimensions         (v2+) one per axis; metres exponent
                                           1/1 marks a spatial axis
    size rank x (uint32 n, n bytes)        axis labels ("ExpControl X", ...)

dt flags: 0x1 uint8, 0x2 int8, 0x4 uint16, 0x8 int16, 0x10 uint32, 0x20 int32,
0x40 float32, 0x80 float64, 0x400 RGB (3 x uint8), 0x800 RGBA, 0x1000 uint64,
0x2000 int64, 0x10000 bool; 0x40000000 OR'd in marks complex data.

The voxel pitch is len / res per axis, reported in micrometres for the axes the
footer marks as metres (or, with no SI footer, for the first three axes whose
length is positive). Only the first stack is detailed; the stack list is walked
to a bound so the count is known.
"""

from __future__ import annotations

import math
import struct
from typing import Dict, List, Optional

from .decoders import DecodedHeader, Read, Source, register_decoder

OBF_READ_BUDGET = 1024 * 1024
_FILE_MAGIC = b"OMAS_BF\n\xff\xff"
_STACK_MAGIC = b"OMAS_BF_STACK\n\xff\xff"
_FILE_HEADER = struct.Struct("<10sIQI")
_STACK_HEADER = struct.Struct("<16sII15I15d15d5I3Q")  # 368 bytes
_MAX_DIMS = 15
_MAX_STACKS = 64
_MAX_STACKS_LISTED = 16
_MAX_NAME = 4096
_MAX_DESCR = 64 * 1024
_SI_UNIT = 9 * 8 + 8  # 9 fractions + double scale factor
_FOOTER_SI_DIMS_OFFSET = 128 + _SI_UNIT
_FOOTER_WITH_SI = _FOOTER_SI_DIMS_OFFSET + _MAX_DIMS * _SI_UNIT  # 1408 bytes

_DTYPE = {
    0x1: "uint8", 0x2: "int8", 0x4: "uint16", 0x8: "int16", 0x10: "uint32",
    0x20: "int32", 0x40: "float32", 0x80: "float64", 0x400: "rgb-uint8",
    0x800: "rgba-uint8", 0x1000: "uint64", 0x2000: "int64", 0x10000: "bool",
}
_COMPLEX = 0x40000000
_COMPRESSION = {0: "none", 1: "zlib"}


def _finite(x: float) -> Optional[float]:
    return x if math.isfinite(x) else None


def _dtype_name(dt: int) -> Optional[str]:
    base = _DTYPE.get(dt & ~_COMPLEX)
    if base is None:
        return None
    if dt & _COMPLEX:
        return {"float32": "complex64", "float64": "complex128"}.get(base, "complex-" + base)
    return base


def _text(raw: bytes) -> Optional[str]:
    s = raw.decode("utf-8", "replace").strip("\x00 \r\n\t")
    return s or None


def _read_stack(src: Source, pos: int):
    """(header dict, data_pos) for the stack header at `pos`, or None."""
    hdr = src.read_at(pos, _STACK_HEADER.size)
    if len(hdr) < _STACK_HEADER.size or hdr[:16] != _STACK_MAGIC:
        return None
    v = _STACK_HEADER.unpack(hdr)
    rank = v[2]
    if not 0 < rank <= _MAX_DIMS:
        return None
    res = list(v[3:3 + _MAX_DIMS])
    lengths = list(v[18:18 + _MAX_DIMS])
    offsets = list(v[33:33 + _MAX_DIMS])
    dt, comp_type, comp_level, name_len, descr_len = v[48:53]
    _reserved, data_len, next_pos = v[53:56]
    if name_len > _MAX_NAME or descr_len > _MAX_DESCR:
        return None
    name = _text(src.read_at(pos + _STACK_HEADER.size, name_len)) if name_len else None
    data_pos = pos + _STACK_HEADER.size + name_len + descr_len
    return {
        "version": v[1], "rank": rank, "res": res[:rank],
        "len": [_finite(x) for x in lengths[:rank]],
        "off": [_finite(x) for x in offsets[:rank]],
        "dt": dt, "compression": comp_type, "compressionLevel": comp_level,
        "name": name, "descrLen": descr_len, "dataLen": data_len,
        "next": next_pos, "dataPos": data_pos, "pos": pos,
    }


def _read_footer(src: Source, stack) -> Dict[str, object]:
    """Axis labels and which axes are metres, from the stack footer. Best
    effort: an unreadable footer leaves both unknown."""
    out: Dict[str, object] = {"labels": None, "metres": None}
    if stack["version"] < 1:
        return out
    pos = stack["dataPos"] + stack["dataLen"]
    rank = stack["rank"]
    head = src.read_at(pos, 4)
    if len(head) < 4:
        return out
    (size,) = struct.unpack_from("<I", head, 0)
    if size < 124 or size > 1 << 20:
        return out
    body = src.read_at(pos, size + 4 + rank * (4 + 256))
    if len(body) < size:
        return out
    if size >= _FOOTER_WITH_SI and len(body) >= _FOOTER_WITH_SI:
        metres = []
        for i in range(rank):
            o = _FOOTER_SI_DIMS_OFFSET + i * _SI_UNIT
            num, den = struct.unpack_from("<ii", body, o)
            metres.append(den != 0 and num / den == 1.0)
        out["metres"] = metres
    labels: List[Optional[str]] = []
    o = size
    for _ in range(rank):
        if o + 4 > len(body):
            break
        (n,) = struct.unpack_from("<I", body, o)
        o += 4
        if n > 256 or o + n > len(body):
            break
        labels.append(_text(body[o:o + n]))
        o += n
    if labels:
        out["labels"] = labels + [None] * (rank - len(labels))
    return out


def _voxel_um(stack, metres: Optional[List[bool]]) -> List[Optional[float]]:
    out: List[Optional[float]] = []
    for i in range(min(3, stack["rank"])):
        n = stack["res"][i]
        length = stack["len"][i]
        spatial = metres[i] if metres is not None and i < len(metres) else (length is not None and length > 0 and n > 1)
        if spatial and n and length:
            out.append(round(abs(length) / n * 1e6, 6))
        else:
            out.append(None)
    while len(out) < 3:
        out.append(None)
    return out


def decode_obf(src: Source) -> Optional[DecodedHeader]:
    head = src.read_at(0, _FILE_HEADER.size)
    if len(head) < _FILE_HEADER.size or head[:10] != _FILE_MAGIC:
        return None
    _magic, version, first_pos, descr_len = _FILE_HEADER.unpack(head)
    if descr_len > _MAX_DESCR:
        return None
    description = _text(src.read_at(_FILE_HEADER.size, descr_len)) if descr_len else None
    meta_pos = None
    if version >= 2:
        raw = src.read_at(_FILE_HEADER.size + descr_len, 8)
        if len(raw) == 8:
            meta_pos = struct.unpack("<Q", raw)[0]

    stacks = []
    pos = first_pos
    seen = set()
    truncated = False
    while pos and pos not in seen and len(stacks) < _MAX_STACKS:
        seen.add(pos)
        st = _read_stack(src, pos)
        if st is None:
            break
        stacks.append(st)
        pos = st["next"]
    if pos and len(stacks) >= _MAX_STACKS:
        truncated = True
    if not stacks:
        return None

    first = stacks[0]
    footer = _read_footer(src, first)
    dtype = _dtype_name(first["dt"])
    sizes = first["res"]
    fields: Dict[str, object] = {
        "formatVersion": version,
        "description": description[:500] if description else None,
        "hasMetadata": meta_pos is not None and meta_pos > 0,
        "stackCount": len(stacks),
        "stackCountTruncated": truncated,
        "stacks": [{"name": s["name"], "sizes": s["res"], "dtype": _dtype_name(s["dt"])}
                   for s in stacks[:_MAX_STACKS_LISTED]],
        "name": first["name"],
        "stackVersion": first["version"],
        "rank": first["rank"],
        "sizes": sizes,
        "axisLabels": footer["labels"],
        "sizeX": sizes[0] if len(sizes) > 0 else None,
        "sizeY": sizes[1] if len(sizes) > 1 else None,
        "sizeZ": sizes[2] if len(sizes) > 2 else None,
        "sizeT": sizes[3] if len(sizes) > 3 else None,
        "dtype": dtype,
        "dataTypeFlags": first["dt"],
        "compression": _COMPRESSION.get(first["compression"], "unknown-%d" % first["compression"]),
        "voxelSizeUm": _voxel_um(first, footer["metres"]),
        "lengths": first["len"],
        "offsets": first["off"],
        "dataBytesOnDisk": first["dataLen"],
    }
    parts = ["Imspector OBF"]
    if sizes:
        parts[0] += " " + "x".join(str(n) for n in sizes[:3] if n) if any(sizes[:3]) else ""
    if dtype:
        parts.append(dtype)
    vs = [v for v in fields["voxelSizeUm"] if v is not None]
    if len(vs) >= 2:
        parts.append("x".join("%g" % v for v in vs) + " um")
    if len(stacks) > 1:
        parts.append("%d stacks%s" % (len(stacks), "+" if truncated else ""))
    return DecodedHeader(format="obf", summary=", ".join(parts), fields=fields)


register_decoder(["msr", "obf"], decode_obf, read=Read(seek=OBF_READ_BUDGET))
