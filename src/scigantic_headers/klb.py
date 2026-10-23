"""Decode a Keller Lab Block (KLB) header.

KLB is the blocked, compressed volume format from the Keller lab's light-sheet
pipelines (SiMView, Fiji/KLB plugins). The header is fixed, little-endian,
319 bytes, and the whole description of the volume is in it:

    0    headerVersion    uint8      1 or 2
    1    xyzct[5]         uint32     image size per dimension
    21   pixelSize[5]     float32    sampling per dimension (see units below)
    41   dataType         uint8      see _KLB_DTYPE
    42   compressionType  uint8      0 none, 1 bzip2, 2 zlib
    43   metadata[256]    char       free text, NUL-terminated
    299  blockSize[5]     uint32     block (chunk) size per dimension
    319  blockOffset[N]   uint64     end offset of each compressed block, one per
                                    block, N = product of ceil(size / blockSize)
    319 + 8N             the compressed blocks

KLB does not record a unit for pixelSize; the writers (and the Fiji plugin)
treat it as micrometres, which is what `voxelSizeUm` assumes. A header is
rejected unless the version, type and compression codes are in range and
every size is a positive, plausible count.
"""

from __future__ import annotations

import math
import struct
from typing import Optional

from .decoders import DecodedHeader, register_decoder

KLB_HEADER_BYTES = 319

_KLB_DTYPE = {
    0: "uint8", 1: "uint16", 2: "uint32", 3: "uint64",
    4: "int8", 5: "int16", 6: "int32", 7: "int64",
    8: "float32", 9: "float64",
}
_KLB_COMPRESSION = {0: "none", 1: "bzip2", 2: "zlib"}
_HEADER = struct.Struct("<B5I5fBB256s5I")  # exactly 319 bytes


def decode_klb(data: bytes) -> Optional[DecodedHeader]:
    if len(data) < KLB_HEADER_BYTES:
        return None
    version, x, y, z, c, t, px, py, pz, pc, pt, dtype_code, comp_code, meta, bx, by, bz, bc, bt = \
        _HEADER.unpack_from(data, 0)
    if version not in (1, 2) or dtype_code not in _KLB_DTYPE or comp_code not in _KLB_COMPRESSION:
        return None
    sizes = (x, y, z, c, t)
    if any(s == 0 or s > (1 << 31) for s in sizes):
        return None

    def um(v: float) -> Optional[float]:
        return round(v, 6) if math.isfinite(v) and v > 0 else None

    blocks = [bx, by, bz, bc, bt]
    num_blocks = None
    if all(b > 0 for b in blocks):
        num_blocks = math.prod(-(-s // b) for s, b in zip(sizes, blocks))

    dtype = _KLB_DTYPE[dtype_code]
    metadata = meta.split(b"\x00", 1)[0].decode("utf-8", "replace").strip() or None
    voxel = [um(px), um(py), um(pz)]
    fields = {
        "headerVersion": version,
        "sizeX": x, "sizeY": y, "sizeZ": z, "sizeC": c, "sizeT": t,
        "dtype": dtype,
        "compression": _KLB_COMPRESSION[comp_code],
        "voxelSizeUm": voxel,
        "samplingCT": [um(pc), um(pt)],
        "blockSize": blocks,
        "metadata": metadata,
        # The blocks start after the offset table, not right after the header.
        "dataOffset": KLB_HEADER_BYTES + 8 * num_blocks if num_blocks is not None else None,
    }
    dims = "%dx%dx%d" % (x, y, z)
    extra = ""
    if c > 1:
        extra += ", %d channels" % c
    if t > 1:
        extra += ", %d timepoints" % t
    vs = ""
    if all(v is not None for v in voxel):
        vs = ", %s um" % "x".join("%g" % v for v in voxel)
    summary = "KLB volume %s%s, %s, %s%s" % (dims, extra, dtype, _KLB_COMPRESSION[comp_code], vs)
    return DecodedHeader(format="klb", summary=summary, fields=fields)


register_decoder("klb", decode_klb)
