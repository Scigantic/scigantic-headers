"""Decode a TIFF header: classic and BigTIFF, and the microscopy flavors that
ride on it (OME-TIFF, ImageJ hyperstacks, Micro-Manager, ImSpector/LaVision).

TIFF keeps its metadata behind pointers, so this is a seeking decoder: it is
handed a `Source` and reads where the pointers lead, within the budget.

    Classic:  0  'II' | 'MM'   byte order
              2  42            uint16 magic
              4  IFD0 offset   uint32
    BigTIFF:  2  43            uint16 magic
              4  8             uint16 offset size
              6  0             uint16
              8  IFD0 offset   uint64

An IFD is an entry count (uint16 / uint64), entries of 12 / 20 bytes (tag u16,
type u16, count u32/u64, value-or-offset u32/u64; the value is inline when it
fits 4 / 8 bytes), then the next-IFD pointer (0 ends the chain). libtiff and
Bio-Formats usually write the IFD after the image data, so IFD0 is read where
the header says, not assumed at byte 8.

The first IFD is read in full (dimensions, sample layout, resolution,
Software, DateTime, ImageDescription, Micro-Manager tag). Later IFDs are only
counted, by following the next pointers up to a cap, so a stack's page count
is exact when the chain is short and reported as truncated otherwise. Where
the writer's own metadata states the plane count (OME SizeZ*SizeC*SizeT,
ImageJ `images=`, the Micro-Manager index map) that is reported as planeCount.

Micro-Manager files add a private header in bytes 8-39 of the TIFF header
(offset-table pointers and the summary metadata JSON length, with the JSON at
byte 40), which gives the acquisition's channels/slices/frames/pixel size
without touching a single IFD beyond the first.
"""

from __future__ import annotations

import json
import math
import re
import struct
from typing import Dict, List, Optional, Tuple

from .decoders import DecodedHeader, Read, Source, register_decoder
from .ome_xml import is_ome_xml, ome_summary, parse_ome_xml

# Byte budget per decode: header + first IFD + an OME-XML description that may
# run to a few hundred KB + a few dozen later IFD blocks.
TIFF_READ_BUDGET = 4 * 1024 * 1024

# Follow at most this many IFDs. Each later IFD costs one bounded read, and the
# Source caps reads at SEEK_MAX_READS, so this keeps the whole decode within it.
IFD_CHAIN_CAP = 64
# One read per IFD. The first IFD carries the tags of interest and is read in a
# block that covers ~170 classic entries; a later IFD is only read for its count
# and next pointer, so a smaller block (40 entries) is tried first and extended
# only when the IFD is longer.
_IFD_BLOCK = 2048
_IFD_BLOCK_LATER = 512
_MAX_ENTRIES = 512
_MAX_DESCRIPTION = 1024 * 1024
_MAX_SMALL_VALUE = 64 * 1024

_TYPE_SIZE = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8, 16: 8, 17: 8, 18: 8}
_TYPE_FMT = {1: "B", 3: "H", 4: "I", 6: "b", 8: "h", 9: "i", 11: "f", 12: "d", 16: "Q", 17: "q", 18: "Q"}

T_WIDTH, T_LENGTH, T_BITS, T_COMPRESSION, T_PHOTOMETRIC = 256, 257, 258, 259, 262
T_DESCRIPTION, T_SPP, T_XRES, T_YRES, T_RESUNIT = 270, 277, 282, 283, 296
T_SOFTWARE, T_DATETIME, T_TILEWIDTH, T_SAMPLEFORMAT = 305, 306, 322, 339
T_MM_METADATA = 51123
_WANTED = {T_WIDTH, T_LENGTH, T_BITS, T_COMPRESSION, T_PHOTOMETRIC, T_DESCRIPTION, T_SPP,
           T_XRES, T_YRES, T_RESUNIT, T_SOFTWARE, T_DATETIME, T_TILEWIDTH, T_SAMPLEFORMAT,
           T_MM_METADATA}

_COMPRESSION = {
    1: "none", 2: "ccitt-rle", 3: "ccitt-t4", 4: "ccitt-t6", 5: "lzw", 6: "jpeg-old", 7: "jpeg",
    8: "deflate", 32773: "packbits", 32946: "deflate", 33003: "aperio-jp2k", 33005: "aperio-jp2k",
    34712: "jpeg2000", 34925: "lzma", 50000: "zstd", 50001: "webp", 50002: "jpegxl",
}
_PHOTOMETRIC = {
    0: "white-is-zero", 1: "black-is-zero", 2: "rgb", 3: "palette", 4: "mask",
    5: "separated", 6: "ycbcr", 8: "cielab", 32844: "logl", 32845: "logluv",
}

# Micro-Manager private header words (micro-manager.org, Multipage TIFF format).
_MM_INDEX_MAP_HEADER = 54773648
_MM_INDEX_MAP_INNER = 3453623
_MM_SUMMARY_HEADER = 2355492

# Length units ImageJ writes in `unit=`, as a factor to micrometres.
_IJ_UNIT_UM = {
    "micron": 1.0, "microns": 1.0, "um": 1.0, "µm": 1.0, "μm": 1.0, "micrometer": 1.0,
    "nm": 1e-3, "mm": 1e3, "cm": 1e4, "m": 1e6, "inch": 25400.0, "pixel": None, "pixels": None,
}


def _finite(x):
    return x if isinstance(x, (int, float)) and math.isfinite(x) else None


class _Ifd:
    """One parsed IFD: {tag: (type, count, value_or_offset, inline_bytes)} and
    the next-IFD pointer. None when the bytes do not form an IFD."""

    def __init__(self, bo: str, big: bool):
        self.bo, self.big = bo, big
        self.entries: Dict[int, Tuple[int, int, int, bytes]] = {}
        self.next = 0
        self.count = 0

    @property
    def entry_size(self):
        return 20 if self.big else 12


def _read_ifd(src: Source, offset: int, bo: str, big: bool, block_size: int = _IFD_BLOCK) -> Optional[_Ifd]:
    if offset <= 0:
        return None
    block = src.read_at(offset, block_size)
    cnt_size = 8 if big else 2
    if len(block) < cnt_size:
        return None
    count = struct.unpack_from(bo + ("Q" if big else "H"), block, 0)[0]
    if count == 0 or count > _MAX_ENTRIES:
        return None
    ifd = _Ifd(bo, big)
    ifd.count = count
    need = cnt_size + count * ifd.entry_size + (8 if big else 4)
    if len(block) < need:
        more = src.read_at(offset + len(block), need - len(block))
        block += more
        if len(block) < need:
            return None
    ent = bo + ("HHQQ" if big else "HHII")
    pos = cnt_size
    for _ in range(count):
        tag, typ, n, _raw = struct.unpack_from(ent, block, pos)
        inline = block[pos + (12 if big else 8):pos + ifd.entry_size]
        if tag in _WANTED:
            ifd.entries[tag] = (typ, n, _raw, inline)
        pos += ifd.entry_size
    ifd.next = struct.unpack_from(bo + ("Q" if big else "I"), block, pos)[0]
    return ifd


def _value_bytes(src: Source, ifd: _Ifd, tag: int, cap: int) -> Optional[bytes]:
    e = ifd.entries.get(tag)
    if e is None:
        return None
    typ, n, raw, inline = e
    size = _TYPE_SIZE.get(typ)
    if size is None or n < 0:
        return None
    total = size * n
    if total > cap:
        total = cap - cap % size if size else cap
    if size * n <= (8 if ifd.big else 4):
        return inline[:total]
    data = src.read_at(raw, total)
    return data if len(data) == total else None


def _numbers(src: Source, ifd: _Ifd, tag: int) -> Optional[List[float]]:
    e = ifd.entries.get(tag)
    if e is None:
        return None
    typ = e[0]
    data = _value_bytes(src, ifd, tag, _MAX_SMALL_VALUE)
    if data is None:
        return None
    try:
        if typ in (5, 10):  # RATIONAL / SRATIONAL: pairs
            f = ifd.bo + ("II" if typ == 5 else "ii")
            out = []
            for i in range(0, len(data) - 7, 8):
                num, den = struct.unpack_from(f, data, i)
                out.append(num / den if den else 0.0)
            return out
        fmt = _TYPE_FMT.get(typ)
        if fmt is None:
            return None
        size = _TYPE_SIZE[typ]
        return [float(v) if fmt in "fd" else v for v in struct.unpack_from(ifd.bo + fmt * (len(data) // size), data, 0)]
    except struct.error:
        return None


def _first_int(src, ifd, tag) -> Optional[int]:
    vals = _numbers(src, ifd, tag)
    if not vals:
        return None
    v = vals[0]
    return int(v) if isinstance(v, (int, float)) and math.isfinite(v) else None


def _text(src: Source, ifd: _Ifd, tag: int, cap: int) -> Optional[str]:
    data = _value_bytes(src, ifd, tag, cap)
    if not data:
        return None
    return data.split(b"\x00", 1)[0].decode("utf-8", "replace").strip() or None


def _dtype(bits: Optional[int], fmt: Optional[int]) -> Optional[str]:
    if not bits:
        return None
    kind = {2: "int", 3: "float", 5: "complex", 6: "complex"}.get(fmt or 1, "uint")
    if kind == "float" and bits not in (16, 32, 64):
        return None
    if bits == 1:
        return "bit"
    return "%s%d" % (kind, bits)


def _parse_imagej(desc: str) -> Dict[str, str]:
    out = {}
    for line in desc.split("\n")[:64]:
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def _ij_int(ij: Dict[str, str], key: str) -> Optional[int]:
    try:
        return int(ij[key])
    except (KeyError, ValueError):
        return None


def _ij_float(ij: Dict[str, str], key: str) -> Optional[float]:
    try:
        return _finite(float(ij[key]))
    except (KeyError, ValueError):
        return None


def _resolution_um(res: Optional[List[float]], unit: Optional[int], ij_unit: Optional[str]) -> Optional[float]:
    """Pixel pitch in micrometres from a TIFF resolution (pixels per unit).
    ResolutionUnit 2 = inch, 3 = cm; 1/None = no absolute unit unless ImageJ's
    `unit=` line names one (ImageJ writes micron pitches with ResolutionUnit 1)."""
    if not res or not res[0] or res[0] <= 0 or not math.isfinite(res[0]):
        return None
    per_unit = res[0]
    if ij_unit is not None:
        factor = _IJ_UNIT_UM.get(ij_unit.lower())
        if factor is None:
            return None
        return _finite(round(factor / per_unit, 6))
    if unit == 2:
        return _finite(round(25400.0 / per_unit, 6))
    if unit == 3:
        return _finite(round(10000.0 / per_unit, 6))
    return None


def _micromanager_header(src: Source, head: bytes, bo: str):
    """Micro-Manager's private words in the classic TIFF header: summary JSON
    (and length) and the index map's entry count. (None, None) if absent."""
    if len(head) < 40:
        return None, None
    try:
        words = struct.unpack_from(bo + "8I", head, 8)
    except struct.error:
        return None, None
    summary = None
    planes = None
    if words[6] == _MM_SUMMARY_HEADER and 0 < words[7] <= _MAX_DESCRIPTION:
        raw = src.read_at(40, words[7])
        if len(raw) == words[7]:
            try:
                summary = json.loads(raw.decode("utf-8", "replace"))
                if isinstance(summary, dict) and "Summary" in summary:
                    summary = summary["Summary"]
                if not isinstance(summary, dict):
                    summary = None
            except ValueError:
                summary = None
    if words[0] == _MM_INDEX_MAP_HEADER and words[1] > 0:
        hdr = src.read_at(words[1], 8)
        if len(hdr) == 8:
            inner, n = struct.unpack_from(bo + "II", hdr, 0)
            if inner == _MM_INDEX_MAP_INNER:
                planes = n
    return summary, planes


def _mm_int(d: dict, *keys) -> Optional[int]:
    for k in keys:
        v = d.get(k)
        if isinstance(v, bool):
            continue
        if isinstance(v, int):
            return v
        if isinstance(v, str) and v.isdigit():
            return int(v)
    return None


def _mm_float(d: dict, *keys) -> Optional[float]:
    for k in keys:
        v = d.get(k)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            return _finite(float(v))
        if isinstance(v, str):
            try:
                return _finite(float(v))
            except ValueError:
                pass
    return None


_MM_PIXEL_TYPE = {"GRAY8": "uint8", "GRAY16": "uint16", "GRAY32": "float32", "RGB32": "uint8", "RGB64": "uint16"}


def decode_tiff(src: Source) -> Optional[DecodedHeader]:
    head = src.read_at(0, 40)
    if len(head) < 8:
        return None
    bo = {b"II": "<", b"MM": ">"}.get(head[:2])
    if bo is None:
        return None
    magic = struct.unpack_from(bo + "H", head, 2)[0]
    if magic == 42:
        big = False
        ifd0 = struct.unpack_from(bo + "I", head, 4)[0]
    elif magic == 43:
        big = True
        if len(head) < 16 or struct.unpack_from(bo + "HH", head, 4) != (8, 0):
            return None
        ifd0 = struct.unpack_from(bo + "Q", head, 8)[0]
    else:
        return None
    if ifd0 < 8 or (src.size is not None and ifd0 >= src.size):
        return None

    ifd = _read_ifd(src, ifd0, bo, big)
    if ifd is None:
        return None

    width = _first_int(src, ifd, T_WIDTH)
    height = _first_int(src, ifd, T_LENGTH)
    if not width or not height or width > 1 << 31 or height > 1 << 31:
        return None
    bits = _first_int(src, ifd, T_BITS)
    spp = _first_int(src, ifd, T_SPP) or 1
    fmt = _first_int(src, ifd, T_SAMPLEFORMAT)
    compression = _first_int(src, ifd, T_COMPRESSION)
    photometric = _first_int(src, ifd, T_PHOTOMETRIC)
    tiled = T_TILEWIDTH in ifd.entries
    software = _text(src, ifd, T_SOFTWARE, _MAX_SMALL_VALUE)
    datetime = _text(src, ifd, T_DATETIME, 256)
    xres = _numbers(src, ifd, T_XRES)
    yres = _numbers(src, ifd, T_YRES)
    resunit = _first_int(src, ifd, T_RESUNIT)
    desc_raw = _value_bytes(src, ifd, T_DESCRIPTION, _MAX_DESCRIPTION)
    desc = desc_raw.split(b"\x00", 1)[0] if desc_raw else b""

    mm_summary, mm_planes = (None, None)
    if not big:
        mm_summary, mm_planes = _micromanager_header(src, head, bo)
    mm_tag = None
    if T_MM_METADATA in ifd.entries:
        raw = _value_bytes(src, ifd, T_MM_METADATA, _MAX_SMALL_VALUE)
        if raw:
            try:
                mm_tag = json.loads(raw.split(b"\x00", 1)[0].decode("utf-8", "replace"))
                if not isinstance(mm_tag, dict):
                    mm_tag = None
            except ValueError:
                mm_tag = None

    # Classify the description.
    ome = parse_ome_xml(desc) if is_ome_xml(desc) else None
    ij = None
    mm_desc = None
    desc_kind = None
    if ome is not None:
        desc_kind = "ome-xml"
    elif desc.startswith(b"ImageJ="):
        ij = _parse_imagej(desc.decode("utf-8", "replace"))
        desc_kind = "imagej"
    elif desc.lstrip().startswith(b"{"):
        try:
            j = json.loads(desc.decode("utf-8", "replace"))
            if isinstance(j, dict):
                mm_desc = j.get("Summary") if isinstance(j.get("Summary"), dict) else None
                desc_kind = "micromanager" if mm_desc is not None else "text"
        except ValueError:
            desc_kind = "text"
    elif desc:
        desc_kind = "text"
    if mm_summary is None and mm_desc is not None:
        mm_summary = mm_desc

    lowered = ((software or "") + " " + desc[:4096].decode("latin-1")).lower()
    imspector = "imspector" in lowered or "lavision" in lowered or "miltenyi" in lowered

    # Count pages by following the chain, bounded.
    pages = 1
    nxt = ifd.next
    seen = {ifd0}
    truncated = False
    while nxt:
        if nxt in seen or pages >= IFD_CHAIN_CAP:
            truncated = pages >= IFD_CHAIN_CAP
            break
        seen.add(nxt)
        later = _read_ifd(src, nxt, bo, big, _IFD_BLOCK_LATER)
        if later is None:
            truncated = True
            break
        pages += 1
        nxt = later.next

    dtype = _dtype(bits, fmt)
    fields: Dict[str, object] = {
        "bigTiff": big,
        "byteOrder": "little" if bo == "<" else "big",
        "width": width,
        "height": height,
        "bitsPerSample": bits,
        "samplesPerPixel": spp,
        "sampleFormat": {1: "uint", 2: "int", 3: "float", 4: "void", 5: "complex-int", 6: "complex-float"}.get(fmt, None) if fmt else "uint",
        "dtype": dtype,
        "compression": _COMPRESSION.get(compression, str(compression) if compression is not None else None),
        "photometric": _PHOTOMETRIC.get(photometric, str(photometric) if photometric is not None else None),
        "tiled": tiled,
        "pageCount": pages,
        "pageCountTruncated": truncated,
        "planeCount": None,
        "imageDescriptionKind": desc_kind,
        "software": software,
        "dateTime": datetime,
    }

    if ome is not None:
        fmt_name = "ome-tiff"
        fields.update(ome)
        z, c, t = ome.get("sizeZ") or 1, ome.get("sizeC") or 1, ome.get("sizeT") or 1
        fields["planeCount"] = z * c * t if ome.get("sizeX") else None
        if ome.get("dtype") is None:
            fields["dtype"] = dtype
        fields["flavor"] = "imspector" if imspector else "ome"
        summary = ome_summary("OME-TIFF", fields)
    elif mm_summary is not None or mm_tag is not None:
        fmt_name = "micromanager-tiff"
        s = mm_summary or {}
        intended = s.get("IntendedDimensions") if isinstance(s.get("IntendedDimensions"), dict) else {}
        ch = _mm_int(s, "Channels") or _mm_int(intended, "channel")
        zz = _mm_int(s, "Slices") or _mm_int(intended, "z")
        tt = _mm_int(s, "Frames") or _mm_int(intended, "time")
        pp = _mm_int(s, "Positions") or _mm_int(intended, "position")
        px = _mm_float(s, "PixelSize_um") or (_mm_float(mm_tag, "PixelSizeUm") if mm_tag else None)
        zstep = _mm_float(s, "z-step_um")
        names = s.get("ChNames") if isinstance(s.get("ChNames"), list) else []
        fields.update({
            "sizeX": _mm_int(s, "Width") or width,
            "sizeY": _mm_int(s, "Height") or height,
            "sizeZ": zz, "sizeC": ch, "sizeT": tt, "positions": pp,
            "channels": [str(n) for n in names[:32]],
            "channelCount": ch or (len(names) or None),
            "voxelSizeUm": [px, px, zstep],
            "pixelType": s.get("PixelType") if isinstance(s.get("PixelType"), str) else None,
            "micromanagerVersion": s.get("MicroManagerVersion") if isinstance(s.get("MicroManagerVersion"), str) else None,
            "acquisitionDate": s.get("StartTime") if isinstance(s.get("StartTime"), str) else (s.get("Time") if isinstance(s.get("Time"), str) else None),
            "camera": (mm_tag or {}).get("Camera") if isinstance((mm_tag or {}).get("Camera"), str) else None,
            "binning": _mm_int(mm_tag or {}, "Binning"),
            "exposureMs": _mm_float(mm_tag or {}, "Exposure-ms"),
            "indexMapEntries": mm_planes,
        })
        if fields["dtype"] is None and fields["pixelType"]:
            fields["dtype"] = _MM_PIXEL_TYPE.get(fields["pixelType"])
        fields["planeCount"] = mm_planes if mm_planes is not None else (
            (zz or 1) * (ch or 1) * (tt or 1) if (zz or ch or tt) else None)
        fields["flavor"] = "micromanager"
        summary = ome_summary("Micro-Manager TIFF", fields)
    elif ij is not None:
        fmt_name = "imagej-tiff"
        images = _ij_int(ij, "images")
        slices = _ij_int(ij, "slices")
        channels = _ij_int(ij, "channels")
        frames = _ij_int(ij, "frames")
        unit = ij.get("unit")
        spacing = _ij_float(ij, "spacing")
        zf = _IJ_UNIT_UM.get(unit.lower()) if unit else None
        fields.update({
            "sizeX": width, "sizeY": height,
            "sizeZ": slices, "sizeC": channels, "sizeT": frames,
            "channelCount": channels,
            "voxelSizeUm": [
                _resolution_um(xres, resunit, unit),
                _resolution_um(yres, resunit, unit),
                _finite(round(spacing * zf, 6)) if (spacing and zf) else None,
            ],
            "imagejVersion": ij.get("ImageJ"),
            "hyperstack": ij.get("hyperstack", "").lower() == "true" if "hyperstack" in ij else None,
            "unit": unit,
        })
        fields["planeCount"] = images
        fields["flavor"] = "imagej"
        summary = ome_summary("ImageJ TIFF", fields)
        if images and images > 1:
            summary += ", %d planes" % images
    else:
        fmt_name = "imspector-tiff" if imspector else "tiff"
        fields.update({
            "sizeX": width, "sizeY": height,
            "sizeZ": None if truncated else pages, "sizeC": None, "sizeT": None,
            "voxelSizeUm": [_resolution_um(xres, resunit, None), _resolution_um(yres, resunit, None), None],
        })
        if not truncated:
            fields["planeCount"] = pages
        fields["flavor"] = "imspector" if imspector else "plain"
        kind = "BigTIFF" if big else "TIFF"
        summary = "%s %dx%d" % (kind, width, height)
        if pages > 1:
            summary += ", %s%d pages" % (">=" if truncated else "", pages)
        if spp > 1:
            summary += ", %d samples/px" % spp
        if dtype:
            summary += ", " + dtype
        vs = fields["voxelSizeUm"]
        if vs[0] is not None and vs[1] is not None:
            summary += ", %gx%g um" % (vs[0], vs[1])
        if fields["compression"] not in (None, "none"):
            summary += ", " + str(fields["compression"])
        if imspector:
            summary = "ImSpector " + summary

    if imspector and fmt_name == "ome-tiff":
        summary = "ImSpector " + summary
    return DecodedHeader(format=fmt_name, summary=summary, fields=fields)


register_decoder(["tif", "tiff", "btf", "tf2", "tf8"], decode_tiff, read=Read(seek=TIFF_READ_BUDGET))
