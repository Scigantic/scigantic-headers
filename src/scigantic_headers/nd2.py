"""Decode a Nikon NIS-Elements ND2 header.

A modern ND2 (file versions 2.x and 3.x) is a sequence of named chunks with a
chunk map at the very end, so the decode is: read the last 40 bytes, jump to
the chunk map, then fetch the handful of metadata chunks by name. Image frames
are never read. Verified against the `nd2` Python package on real files of
both generations.

Chunk:
    0   uint32  magic 0x0ABECEDA
    4   uint32  nameLength
    8   uint64  dataLength
    16  name    nameLength bytes, NUL padded, e.g. b"ImageAttributesLV!"
    ..  data    dataLength bytes

File start: a chunk named b"ND2 FILE SIGNATURE CHUNK NAME01!" whose 64-byte
data begins "Ver2.1" / "Ver3.0". File end (last 40 bytes):
    char[32]  b"ND2 CHUNK MAP SIGNATURE 0000001!"
    uint64    position of the chunk-map chunk (named b"ND2 FILEMAP SIGNATURE NAME 0001!")
Chunk-map data: repeated (name up to and including '!', uint64 offset, uint64
size), terminated by the chunk-map signature as a name.

Metadata chunks, two encodings by file generation:
  * Ver2.x: b"ImageAttributes!", b"ImageCalibration|0!", b"ImageTextInfo!",
    b"ImageMetadata!", b"ImageMetadataSeq|0!" hold CLxVariant XML
    (`<variant><no_name runtype="CLxListVariant"><uiWidth runtype="lx_uint32"
    value="696"/>...`).
  * Ver3.x: the same names with an "LV" suffix (b"ImageAttributesLV!" ...) hold
    Nikon's LiteVariant binary: records of uint8 type, uint8 nameLength (UTF-16
    code units), UTF-16-LE name, then a value by type: 1 bool (u8), 2 int32,
    3 uint32, 4 int64, 5 uint64, 6 double, 7 pointer (u64), 8 string (UTF-16,
    NUL terminated), 9 byte array (u64 length), 11 level (uint32 count, uint64
    byte length of the record measured from its start, `count` child records,
    then count x uint64 offset table). Type 76 ('L') is a zlib-compressed level:
    skip 10 bytes and inflate the rest.
Both decode to the same tree of Nikon's prefixed names (uiWidth, dCalibration,
sDescription, ...), so one extraction serves both generations.

What is read:
  ImageAttributes      uiWidth, uiHeight, uiComp (components per frame),
                       uiBpcInMemory, uiBpcSignificant, uiSequenceCount (frames),
                       eCompression / uiCompression
  ImageCalibration|0   dCalibration (micrometres per pixel)
  ImageMetadataSeq|0   sPicturePlanes (channel count and names) from the first
                       32 KiB; wsObjectiveName, dObjectiveMag, dObjectiveNA scanned
                       out of the last 8 KiB (the chunk can be hundreds of KB of
                       camera settings in between)
  ImageMetadata        SLxExperiment: the acquisition loop tree. eType 1 TimeLoop,
                       2 XYPosLoop, 4 ZStackLoop (uiCount, dZStep), 6 SpectLoop,
                       7 CustomLoop, 8 NETimeLoop (pPeriod entries, pPeriodValid);
                       nested loops under ppNextLevelEx.
  ImageTextInfo        TextInfoItem_5 description, _6 capturing, _9 date, _13 optics

A legacy ND2 (NIS-Elements 2.x) is a JPEG 2000 container (starts with the
jP signature box); it is reported as format "nd2-legacy" with no fields.
"""

from __future__ import annotations

import math
import re
import struct
import xml.etree.ElementTree as ET
import zlib
from typing import Dict, List, Optional, Tuple

from .decoders import DecodedHeader, Read, Source, register_decoder

ND2_READ_BUDGET = 4 * 1024 * 1024
_MAGIC = 0x0ABECEDA
_FILE_SIG = b"ND2 FILE SIGNATURE CHUNK NAME01!"
_MAP_SIG = b"ND2 CHUNK MAP SIGNATURE 0000001!"
_FILEMAP_SIG = b"ND2 FILEMAP SIGNATURE NAME 0001!"
_JP2_SIG = b"\x00\x00\x00\x0cjP  \r\n\x87\n"
_CHUNK_HDR = struct.Struct("<IIQ")
_MAX_MAP_BYTES = 1024 * 1024
_MAX_CHUNK_BYTES = 1024 * 1024
_MAX_MAP_ENTRIES = 100000
_MAX_LV_RECORDS = 200000
_MAX_DEPTH = 24
_MAX_CHANNELS = 32

_WANTED = ("ImageAttributes", "ImageCalibration|0", "ImageMetadataSeq|0", "ImageMetadata", "ImageTextInfo")

_LOOP_NAMES = {1: "TimeLoop", 2: "XYPosLoop", 4: "ZStackLoop", 6: "SpectLoop", 7: "CustomLoop", 8: "NETimeLoop"}


def _finite(x) -> Optional[float]:
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


# ── LiteVariant (Ver3 "LV" chunks) ────────────────────────────────────────────

class _LV:
    """Bounded reader for a LiteVariant buffer -> nested dict. Repeated empty
    names (Nikon's lists) become Python lists. Never raises past `parse`."""

    def __init__(self, buf: bytes):
        self.b = buf
        self.records = 0

    def parse(self, start: int, end: int, count: Optional[int], depth: int) -> Dict[str, object]:
        out: Dict[str, object] = {}
        i = start
        n = 0
        while i < end and (count is None or n < count):
            if self.records > _MAX_LV_RECORDS or depth > _MAX_DEPTH:
                break
            self.records += 1
            if i + 2 > end:
                break
            kind, name_len = self.b[i], self.b[i + 1]
            if kind == 76:  # zlib-compressed level
                try:
                    inner = zlib.decompress(self.b[i + 12:end])
                except zlib.error:
                    break
                sub = _LV(inner)
                sub.records = self.records
                out.update(sub.parse(0, len(inner), None, depth + 1))
                self.records = sub.records
                break
            rec_start = i
            i += 2
            name = self.b[i:i + name_len * 2].decode("utf-16-le", "replace").rstrip("\x00")
            i += name_len * 2
            value: object
            if kind == 1:
                value = bool(self.b[i]) if i < end else None
                i += 1
            elif kind == 2:
                value = struct.unpack_from("<i", self.b, i)[0] if i + 4 <= end else None
                i += 4
            elif kind == 3:
                value = struct.unpack_from("<I", self.b, i)[0] if i + 4 <= end else None
                i += 4
            elif kind in (4, 5, 7):
                fmt = "<q" if kind == 4 else "<Q"
                value = struct.unpack_from(fmt, self.b, i)[0] if i + 8 <= end else None
                i += 8
            elif kind == 6:
                value = _finite(struct.unpack_from("<d", self.b, i)[0]) if i + 8 <= end else None
                i += 8
            elif kind == 8:
                j = i
                while j + 2 <= end and self.b[j:j + 2] != b"\x00\x00":
                    j += 2
                value = self.b[i:j].decode("utf-16-le", "replace")
                i = j + 2
            elif kind == 9:
                if i + 8 > end:
                    break
                (size,) = struct.unpack_from("<Q", self.b, i)
                i += 8
                value = None  # byte arrays are not needed here
                i += size
            elif kind == 11:
                if i + 12 > end:
                    break
                child_count, length = struct.unpack_from("<IQ", self.b, i)
                i += 12
                level_end = min(rec_start + length, end)
                value = self.parse(i, level_end, child_count, depth + 1)
                i = level_end + child_count * 8
            else:
                break
            if name == "" and name in out:
                cur = out[name]
                if not isinstance(cur, list):
                    out[name] = [cur]
                out[name].append(value)  # type: ignore[union-attr]
            else:
                out[name] = value
            n += 1
        return out


def _from_lv(data: bytes) -> Dict[str, object]:
    try:
        return _LV(data).parse(0, len(data), None, 0)
    except Exception:
        return {}


# ── CLxVariant XML (Ver2 chunks) ───────────────────────────────────────────────

def _xml_value(el):
    rt = el.get("runtype") or ""
    children = list(el)
    if children or rt.startswith("CLxListVariant") or rt.startswith("RLx") or rt.startswith("CLx") and el.get("value") is None and rt not in ("CLxStringW", "CLxByteArray"):
        out: Dict[str, object] = {}
        for c in children:
            v = _xml_value(c)
            key = "" if c.tag == "no_name" else c.tag
            if key in out:
                cur = out[key]
                if not isinstance(cur, list):
                    out[key] = [cur]
                out[key].append(v)  # type: ignore[union-attr]
            else:
                out[key] = v
        return out
    v = el.get("value")
    if rt == "bool":
        return (v or "").lower() in ("true", "1")
    if rt in ("double", "float"):
        return _finite(v)
    if rt.startswith("lx_"):
        try:
            return int(v)
        except (TypeError, ValueError):
            return None
    return v


def _from_xml(data: bytes) -> Dict[str, object]:
    try:
        root = ET.fromstring(data.split(b"?>", 1)[-1] if data.startswith(b"<?xml") else data)
    except ET.ParseError:
        return {}
    if root.tag != "variant":
        return {}
    out: Dict[str, object] = {}
    for c in root:
        out["" if c.tag == "no_name" else c.tag] = _xml_value(c)
    return out


def _decode_chunk(data: bytes) -> Dict[str, object]:
    head = data[:64].lstrip()
    if head.startswith(b"<?xml") or head.startswith(b"<variant"):
        return _from_xml(data)
    return _from_lv(data)


# ── tree helpers ───────────────────────────────────────────────────────────────

def _find(tree, key: str, depth: int = 0):
    """First value under `key` anywhere in a nested dict/list, bounded depth."""
    if depth > _MAX_DEPTH:
        return None
    if isinstance(tree, dict):
        if key in tree:
            return tree[key]
        for v in tree.values():
            r = _find(v, key, depth + 1)
            if r is not None:
                return r
    elif isinstance(tree, list):
        for v in tree:
            r = _find(v, key, depth + 1)
            if r is not None:
                return r
    return None


def _unwrap(d):
    """The XML encoding wraps a struct in a single `no_name` child (key "");
    the LV encoding does not. Return the struct either way."""
    while isinstance(d, dict) and len(d) == 1 and "" in d:
        d = d[""]
    return d


def _values(d) -> list:
    """Members of a Nikon container. The LV encoding lists members under
    repeated empty names (a list under key ""), the XML encoding under
    repeated `no_name` children (the same), and a one-member container is a
    single dict under "". Named members (a0, a1, ...) are a dict's values."""
    if isinstance(d, list):
        return d
    if not isinstance(d, dict):
        return []
    if len(d) == 1 and "" in d:
        inner = d[""]
        return inner if isinstance(inner, list) else [inner]
    return list(d.values())


def _int(v) -> Optional[int]:
    return int(v) if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) else None


# ── experiment loops ───────────────────────────────────────────────────────────

def _loop_count(etype: int, pars) -> Tuple[Optional[int], Dict[str, object]]:
    extra: Dict[str, object] = {}
    pars = _unwrap(pars)
    if not isinstance(pars, dict):
        return None, extra
    if etype == 8:  # NETimeLoop: sum of valid periods
        valid = _values(pars.get("pPeriodValid"))
        total = 0
        for i, period in enumerate(_values(pars.get("pPeriod"))):
            if i < len(valid) and not valid[i]:
                continue
            period = _unwrap(period)
            c = _int(period.get("uiCount")) if isinstance(period, dict) else None
            total += c or 0
        return (total or _int(pars.get("uiCount"))), extra
    if etype == 6:  # SpectLoop: planes
        planes = pars.get("pPlanes")
        c = _int(_find(planes, "uiCount")) if planes is not None else None
        return c, extra
    count = _int(pars.get("uiCount"))
    if etype == 4:
        step = _finite(pars.get("dZStep"))
        lo, hi = _finite(pars.get("dZLow")), _finite(pars.get("dZHigh"))
        if (not step) and count and count > 1 and lo is not None and hi is not None:
            step = abs(hi - lo) / (count - 1)
        extra["zStepUm"] = round(step, 6) if step else None
    if etype == 2 and count:
        valid = _values(pars.get("pItemValid"))
        if valid:
            extra["validPositions"] = sum(1 for v in valid if v)
    return count, extra


def _walk_loops(exp, depth: int, out: List[Dict[str, object]]) -> None:
    if depth > _MAX_DEPTH or not isinstance(exp, dict) or len(out) > 32:
        return
    if "eType" not in exp:
        # A wrapper (the XML root, a no_name container): look one level down.
        for v in exp.values():
            _walk_loops(v, depth + 1, out)
        return
    etype = _int(exp.get("eType"))
    if etype:
        count, extra = _loop_count(etype, exp.get("uLoopPars"))
        if count:
            rec: Dict[str, object] = {"type": _LOOP_NAMES.get(etype, "Loop%d" % etype), "count": count}
            rec.update(extra)
            out.append(rec)
    nxt = exp.get("ppNextLevelEx")
    for child in _values(nxt):
        _walk_loops(child, depth + 1, out)


# ── decode ─────────────────────────────────────────────────────────────────────

def _read_chunk(src: Source, pos: int, cap: int = _MAX_CHUNK_BYTES):
    """(name, data, data_start, data_len) of the chunk at `pos`, data capped
    to `cap` bytes, or None."""
    hdr = src.read_at(pos, _CHUNK_HDR.size)
    if len(hdr) < _CHUNK_HDR.size:
        return None
    magic, name_len, data_len = _CHUNK_HDR.unpack(hdr)
    if magic != _MAGIC or name_len > 65536:  # Ver3 writers pad the chunk-map name to 4128 bytes
        return None
    body = src.read_at(pos + _CHUNK_HDR.size, name_len + min(data_len, cap))
    if len(body) < name_len:
        return None
    return body[:name_len].split(b"\x00", 1)[0], body[name_len:], pos + _CHUNK_HDR.size + name_len, data_len


# The per-frame metadata chunk (ImageMetadataSeq|0) can run to hundreds of KB
# of camera settings. The channel names sit near its start and the objective
# fields near its end, so read a head and a tail rather than the whole chunk.
_SEQ_HEAD = 32 * 1024
_SEQ_TAIL = 8 * 1024


def _scan_scalars(buf: bytes) -> Dict[str, object]:
    """wsObjectiveName / dObjectiveMag / dObjectiveNA out of raw chunk bytes,
    either encoding, last occurrence wins (the top-level record is written
    after the per-sample copies)."""
    out: Dict[str, object] = {}
    if b"runtype=" in buf:
        for key in ("wsObjectiveName", "dObjectiveMag", "dObjectiveNA"):
            hits = re.findall(rb"<" + key.encode() + rb' runtype="[^"]*" value="([^"]*)"', buf)
            if hits:
                v = hits[-1].decode("utf-8", "replace")
                out[key] = v if key.startswith("ws") else _finite(v)
        return out
    for key, kind in (("wsObjectiveName", 8), ("dObjectiveMag", 6), ("dObjectiveNA", 6)):
        pat = bytes([kind, len(key) + 1]) + key.encode("utf-16-le") + b"\x00\x00"
        i = buf.rfind(pat)
        if i < 0:
            continue
        j = i + len(pat)
        if kind == 6:
            if j + 8 <= len(buf):
                out[key] = _finite(struct.unpack_from("<d", buf, j)[0])
        else:
            k = j
            while k + 2 <= len(buf) and buf[k:k + 2] != b"\x00\x00":
                k += 2
            out[key] = buf[j:k].decode("utf-16-le", "replace")
    return out


def _chunk_map(src: Source) -> Optional[Dict[bytes, Tuple[int, int]]]:
    tail = src.read_tail(40)
    if len(tail) < 40 or tail[:32] != _MAP_SIG:
        return None
    (pos,) = struct.unpack_from("<Q", tail, 32)
    chunk = _read_chunk(src, pos, _MAX_MAP_BYTES)
    if chunk is None or chunk[0] != _FILEMAP_SIG:
        return None
    data = chunk[1]
    out: Dict[bytes, Tuple[int, int]] = {}
    i = 0
    while len(out) < _MAX_MAP_ENTRIES:
        j = data.find(b"!", i)
        if j < 0:
            break
        name = data[i:j + 1]
        if name == _MAP_SIG:
            break
        if j + 17 > len(data):
            break
        off, size = struct.unpack_from("<QQ", data, j + 1)
        out[name] = (off, size)
        i = j + 17
    return out


def _text_items(info) -> Dict[str, Optional[str]]:
    def item(n):
        v = _find(info, "TextInfoItem_%d" % n)
        v = v.strip() if isinstance(v, str) else None
        return v or None
    return {"description": item(5), "capturing": item(6), "date": item(9), "optics": item(13)}


def decode_nd2(src: Source) -> Optional[DecodedHeader]:
    head = src.read_at(0, 64)
    if len(head) >= 12 and head[:12] == _JP2_SIG:
        return DecodedHeader(format="nd2-legacy", summary="Nikon ND2 (legacy JPEG 2000 container)",
                             fields={"legacy": True})
    if len(head) < 64:
        return None
    magic, name_len, data_len = _CHUNK_HDR.unpack_from(head, 0)
    if magic != _MAGIC or name_len != 32 or head[16:48] != _FILE_SIG:
        return None
    ver = head[48:64].split(b"\x00", 1)[0].decode("ascii", "replace")  # "Ver3.0"
    m = re.match(r"Ver(\d+)\.(\d+)", ver)
    version = "%s.%s" % (m.group(1), m.group(2)) if m else None

    cmap = _chunk_map(src)
    if cmap is None:
        return None
    frames_in_map = sum(1 for k in cmap if k.startswith(b"ImageDataSeq|"))

    trees: Dict[str, Dict[str, object]] = {}
    seq_tail = b""
    for base in _WANTED:
        # Ver3 names carry "LV" before any "|index": ImageCalibrationLV|0!
        stem, _, index = base.partition("|")
        lv = (stem + "LV" + ("|" + index if index else "")).encode() + b"!"
        for key in (lv, base.encode() + b"!"):
            if key in cmap:
                is_seq = base.startswith("ImageMetadataSeq")
                chunk = _read_chunk(src, cmap[key][0], _SEQ_HEAD if is_seq else _MAX_CHUNK_BYTES)
                if chunk is not None:
                    trees[base] = _decode_chunk(chunk[1])
                    if is_seq:
                        _name, data, start, length = chunk
                        if length > len(data):
                            tail_len = min(_SEQ_TAIL, length - len(data))
                            seq_tail = src.read_at(start + length - tail_len, tail_len)
                        else:
                            seq_tail = data
                break

    attrs = trees.get("ImageAttributes", {})
    width = _int(_find(attrs, "uiWidth"))
    height = _int(_find(attrs, "uiHeight"))
    comp = _int(_find(attrs, "uiComp"))
    bpc = _int(_find(attrs, "uiBpcInMemory"))
    bpc_sig = _int(_find(attrs, "uiBpcSignificant"))
    frames = _int(_find(attrs, "uiSequenceCount"))
    compression = None
    for k in ("eCompression", "uiCompression", "Compression"):
        c = _int(_find(attrs, k))
        if c is not None:
            compression = c
            break
    if width is None or height is None:
        return None

    cal = _finite(_find(trees.get("ImageCalibration|0", {}), "dCalibration"))
    seq = trees.get("ImageMetadataSeq|0", {})
    planes = _find(seq, "sPicturePlanes")
    channels: List[str] = []
    plane_count = None
    if isinstance(planes, dict):
        plane_count = _int(planes.get("uiCount"))
        for p in _values(planes.get("sPlaneNew") or planes.get("sPlane"))[:_MAX_CHANNELS]:
            if isinstance(p, dict):
                name = p.get("sDescription") or p.get("sOpticalConfigName")
                channels.append(name if isinstance(name, str) and name else "Ch%d" % (len(channels) + 1))
    scal = _scan_scalars(seq_tail)
    objective = scal.get("wsObjectiveName")
    objective = objective.strip() or None if isinstance(objective, str) else None
    mag = _finite(scal.get("dObjectiveMag"))
    na = _finite(scal.get("dObjectiveNA"))
    mag = mag if mag and mag > 0 else None
    na = na if na and na > 0 else None

    loops: List[Dict[str, object]] = []
    exp = _find(trees.get("ImageMetadata", {}), "SLxExperiment")
    _walk_loops(exp if exp is not None else trees.get("ImageMetadata", {}), 0, loops)
    size_t = size_z = positions = None
    z_step = None
    for lp in loops:
        if lp["type"] in ("TimeLoop", "NETimeLoop") and size_t is None:
            size_t = lp["count"]
        elif lp["type"] == "ZStackLoop" and size_z is None:
            size_z = lp["count"]
            z_step = lp.get("zStepUm")
        elif lp["type"] == "XYPosLoop" and positions is None:
            positions = lp.get("validPositions") or lp["count"]

    text = _text_items(trees.get("ImageTextInfo", {}))
    if bpc == 8:
        dtype = "uint8"
    elif bpc == 16:
        dtype = "uint16"
    elif bpc == 32:
        dtype = "float32"
    else:
        dtype = None
    size_c = plane_count or comp
    fields: Dict[str, object] = {
        "fileVersion": version,
        "sizeX": width, "sizeY": height,
        "sizeZ": size_z if size_z else (1 if loops or frames else None),
        "sizeT": size_t if size_t else (1 if loops or frames else None),
        "sizeC": size_c, "channelCount": size_c,
        "componentsPerFrame": comp,
        "positionCount": positions,
        "frameCount": frames,
        "framesInChunkMap": frames_in_map,
        "dtype": dtype,
        "bitsPerComponent": bpc,
        "significantBits": bpc_sig,
        "compression": compression,
        "voxelSizeUm": [round(cal, 6) if cal else None, round(cal, 6) if cal else None, z_step],
        "channels": channels,
        "loops": loops,
        "objective": objective or text["optics"],
        "magnification": mag, "numericalAperture": na,
        "acquisitionDate": text["date"],
        "description": text["description"][:500] if text["description"] else None,
        "camera": text["capturing"].splitlines()[0].strip() if text["capturing"] else None,
    }
    parts = ["Nikon ND2 %dx%d" % (width, height)]
    if (size_c or 0) > 1:
        parts.append("%d channels" % size_c)
    if (size_z or 0) > 1:
        parts.append("%d z-planes" % size_z)
    if (size_t or 0) > 1:
        parts.append("%d timepoints" % size_t)
    if (positions or 0) > 1:
        parts.append("%d positions" % positions)
    if dtype:
        parts.append(dtype)
    vs = [v for v in fields["voxelSizeUm"] if v is not None]
    if len(vs) >= 2:
        parts.append("x".join("%g" % v for v in vs) + " um")
    if frames is not None:
        parts.append("%d frames" % frames)
    return DecodedHeader(format="nd2", summary=", ".join(parts), fields=fields)


register_decoder("nd2", decode_nd2, read=Read(seek=ND2_READ_BUDGET))
