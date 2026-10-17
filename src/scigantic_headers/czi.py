"""Decode a Zeiss CZI (Carl Zeiss Image) file header.

A CZI is a sequence of segments, each a 32-byte header followed by a payload:

    0   Id              char[16]  b"ZISRAWFILE", b"ZISRAWMETADATA", b"ZISRAWDIRECTORY",
                                  b"ZISRAWATTDIR", b"ZISRAWSUBBLOCK", b"ZISRAWATTACH"
    16  AllocatedSize   int64     payload bytes reserved
    24  UsedSize        int64     payload bytes written
    32  payload

The file header segment is always at offset 0 and its 512-byte payload points
at everything else:

    0   Major           int32     (1)
    4   Minor           int32     (0)
    8   PrimaryFileGuid 16 bytes
    24  FileGuid        16 bytes
    48  FilePart        int32     0 for a single-file image, N for part N of a split one
    52  DirectoryPosition         int64   offset of the ZISRAWDIRECTORY segment
    60  MetadataPosition          int64   offset of the ZISRAWMETADATA segment
    68  UpdatePending   int32
    72  AttachmentDirectoryPosition int64 offset of the ZISRAWATTDIR segment, 0 if none

The metadata segment payload is `XmlSize int32 @0, AttachmentSize int32 @4`, 248
reserved bytes, then UTF-8 XML (`<ImageDocument><Metadata>...`) of XmlSize
bytes. That XML carries the image shape (Information/Image/Size*), pixel type,
scaling (Scaling/Items/Distance, in metres), channels, instrument and software.
Where it sits varies by writer: a few KB into the file for ZEN blue, or past the
pixel data (hundreds of MB in) for ZEN black / Lightsheet Z.1, so this is a
seeking decoder: file header, then one read at MetadataPosition, then the
directory and attachment segment headers for their entry counts. Four reads,
a few hundred KB, regardless of file size.

The directory segment payload is `EntryCount int32 @0` then 124 reserved bytes
then the entries; only the count is read. Channel ExposureTime is in
nanoseconds in the file and reported here in milliseconds.
"""

from __future__ import annotations

import math
import re
import struct
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

from .decoders import DecodedHeader, Read, Source, register_decoder

_SEG_HEADER = 32
_FILE_HEADER_PAYLOAD = 512
_META_HEADER_PAYLOAD = 256
_ID_FILE = b"ZISRAWFILE"
_ID_META = b"ZISRAWMETADATA"
_ID_DIR = b"ZISRAWDIRECTORY"
_ID_ATTDIR = b"ZISRAWATTDIR"

# Bytes of metadata XML this will read. Real files run from 16 KB (Z.1) to
# 860 KB (LSM980 with a full experiment block); 3 MiB leaves headroom and keeps
# the whole decode under the 4 MiB budget registered below.
_XML_CAP = 3 * 1024 * 1024
_BUDGET = 4 * 1024 * 1024

# PixelType -> (dtype, samples per pixel)
_PIXEL_TYPE = {
    "Gray8": ("uint8", 1), "Gray16": ("uint16", 1), "Gray32": ("int32", 1),
    "Gray32Float": ("float32", 1), "Gray64": ("float64", 1),
    "Gray64ComplexFloat": ("complex64", 1), "Bgr24": ("uint8", 3), "Bgr48": ("uint16", 3),
    "Bgr96Float": ("float32", 3), "Bgra32": ("uint8", 4), "Bgr192ComplexFloat": ("complex64", 3),
}

# Dimension letters beyond X/Y/Z/C/T and the field names they report under.
_EXTRA_DIMS = {
    "S": "scenes", "M": "tiles", "B": "blocks", "V": "views",
    "I": "illuminations", "R": "rotations", "H": "phases",
}


def _finite(x: Optional[float]) -> Optional[float]:
    return x if x is not None and math.isfinite(x) else None


def _num(text) -> Optional[float]:
    try:
        return _finite(float(text))
    except (TypeError, ValueError):
        return None


def _int(text) -> Optional[int]:
    try:
        return int(str(text).strip())
    except (TypeError, ValueError):
        return None


def _segment_id(buf: bytes) -> bytes:
    return buf[:16].rstrip(b"\x00")


def _text(el, path: str) -> Optional[str]:
    """Stripped text of the first element at `path` under `el`, or None."""
    if el is None:
        return None
    found = el.find(path)
    if found is None or found.text is None:
        return None
    t = found.text.strip()
    return t or None


def _xml_fields(xml: bytes) -> Dict[str, object]:
    """Fields from the metadata XML. ElementTree first; if the XML was cut by
    the budget or is otherwise malformed, a regex fallback still recovers the
    shape, pixel type and scaling."""
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        root = None
    if root is not None and root.tag == "ImageDocument":
        return _fields_from_tree(root)
    return _fields_from_regex(xml)


def _fields_from_tree(root) -> Dict[str, object]:
    meta = root.find("Metadata")
    info = meta.find("Information") if meta is not None else None
    image = info.find("Image") if info is not None else None
    f: Dict[str, object] = {}

    for k in ("X", "Y", "Z", "C", "T"):
        f["size" + k] = _int(_text(image, "Size" + k))
    for letter, name in _EXTRA_DIMS.items():
        f[name] = _int(_text(image, "Size" + letter))
    f["pixelType"] = _text(image, "PixelType")
    f["componentBitCount"] = _int(_text(image, "ComponentBitCount"))
    f["acquisitionDate"] = _text(image, "AcquisitionDateAndTime")

    # Scaling: metres in the file, micrometres here. A 0 (ZEN black writes 0 for
    # an unset Z) is unknown, not a zero-thickness voxel.
    scale = {}
    scaling = meta.find("Scaling") if meta is not None else None
    if scaling is not None:
        for d in scaling.findall("Items/Distance"):
            v = _num(_text(d, "Value"))
            if d.get("Id") in ("X", "Y", "Z") and v is not None and v > 0:
                scale[d.get("Id")] = round(v * 1e6, 6)
    f["voxelSizeUm"] = [scale.get("X"), scale.get("Y"), scale.get("Z")]

    channels: List[Optional[str]] = []
    exposures: List[Optional[float]] = []
    modes: List[str] = []
    if image is not None:
        for ch in image.findall("Dimensions/Channels/Channel")[:32]:
            channels.append(ch.get("Name") or ch.get("Id"))
            exp = _num(_text(ch, "ExposureTime"))
            exposures.append(round(exp / 1e6, 4) if exp is not None else None)  # ns -> ms
            mode = _text(ch, "AcquisitionMode")
            if mode and mode not in modes:
                modes.append(mode)
    f["channels"] = channels
    f["channelCount"] = f.get("sizeC") or (len(channels) or None)
    f["exposureTimesMs"] = exposures if any(e is not None for e in exposures) else None
    f["acquisitionMode"] = ", ".join(modes) if modes else None

    inst = info.find("Instrument") if info is not None else None
    mic = inst.find("Microscopes/Microscope") if inst is not None else None
    f["microscope"] = None
    if mic is not None:
        f["microscope"] = mic.get("Name") or _text(mic, "System")
        f["microscopeSystem"] = _text(mic, "System")
    obj = inst.find("Objectives/Objective") if inst is not None else None
    f["objective"] = f["magnification"] = f["numericalAperture"] = f["immersion"] = None
    if obj is not None:
        f["objective"] = obj.get("Name") or _text(obj, "Manufacturer/Model")
        mag = _num(_text(obj, "NominalMagnification")) or _num(_text(obj, "CalibratedMagnification"))
        na = _num(_text(obj, "LensNA"))
        f["magnification"] = round(mag, 6) if mag is not None else None
        f["numericalAperture"] = round(na, 6) if na is not None else None
        f["immersion"] = _text(obj, "Immersion")

    app = info.find("Application") if info is not None else None
    name, version = _text(app, "Name"), _text(app, "Version")
    f["software"] = " ".join(x for x in (name, version) if x) or None
    return f


def _fields_from_regex(xml: bytes) -> Dict[str, object]:
    def tag(name: str) -> Optional[str]:
        m = re.search(rb"<" + name.encode() + rb">([^<]*)</" + name.encode() + rb">", xml)
        return m.group(1).decode("utf-8", "replace").strip() if m else None

    f: Dict[str, object] = {}
    for k in ("X", "Y", "Z", "C", "T"):
        f["size" + k] = _int(tag("Size" + k))
    for letter, name in _EXTRA_DIMS.items():
        f[name] = _int(tag("Size" + letter))
    f["pixelType"] = tag("PixelType")
    f["componentBitCount"] = _int(tag("ComponentBitCount"))
    f["acquisitionDate"] = tag("AcquisitionDateAndTime")
    scale = {}
    for m in re.finditer(rb'<Distance Id="([XYZ])">\s*<Value>([^<]*)</Value>', xml):
        v = _num(m.group(2).decode("ascii", "replace"))
        axis = m.group(1).decode()
        if v is not None and v > 0 and axis not in scale:
            scale[axis] = round(v * 1e6, 6)
    f["voxelSizeUm"] = [scale.get("X"), scale.get("Y"), scale.get("Z")]
    channels = [m.group(1).decode("utf-8", "replace")
                for m in re.finditer(rb'<Channel Id="[^"]*" Name="([^"]*)"', xml)][:32]
    f["channels"] = channels
    f["channelCount"] = f.get("sizeC") or (len(channels) or None)
    f["exposureTimesMs"] = None
    f["acquisitionMode"] = None
    f["microscope"] = f["objective"] = f["magnification"] = f["numericalAperture"] = f["immersion"] = None
    f["software"] = None
    f["xmlTruncated"] = True
    return f


_IMAGE_KEYS = (
    "sizeX", "sizeY", "sizeZ", "sizeC", "sizeT", "scenes", "tiles", "blocks", "views",
    "illuminations", "rotations", "phases", "pixelType", "componentBitCount", "acquisitionDate",
    "voxelSizeUm", "channels", "channelCount", "exposureTimesMs", "acquisitionMode", "microscope",
    "objective", "magnification", "numericalAperture", "immersion", "software",
)


def _entry_count(src: Source, pos: int, want_id: bytes) -> Optional[int]:
    """EntryCount of a directory-like segment at `pos`, or None."""
    if pos <= 0 or (src.size is not None and pos + _SEG_HEADER + 4 > src.size):
        return None
    buf = src.read_at(pos, _SEG_HEADER + 4)
    if len(buf) < _SEG_HEADER + 4 or _segment_id(buf) != want_id:
        return None
    (n,) = struct.unpack_from("<i", buf, _SEG_HEADER)
    return n if 0 <= n <= 100_000_000 else None


def decode_czi(src: Source) -> Optional[DecodedHeader]:
    head = src.read_at(0, _SEG_HEADER + _FILE_HEADER_PAYLOAD)
    if len(head) < _SEG_HEADER + _FILE_HEADER_PAYLOAD or _segment_id(head) != _ID_FILE:
        return None
    p = _SEG_HEADER
    major, minor = struct.unpack_from("<ii", head, p)
    (file_part,) = struct.unpack_from("<i", head, p + 48)
    (dir_pos,) = struct.unpack_from("<q", head, p + 52)
    (meta_pos,) = struct.unpack_from("<q", head, p + 60)
    (update_pending,) = struct.unpack_from("<i", head, p + 68)
    (att_pos,) = struct.unpack_from("<q", head, p + 72)

    fields: Dict[str, object] = {
        "version": "%d.%d" % (major, minor),
        "filePart": file_part,
        "updatePending": bool(update_pending),
        "subBlockCount": _entry_count(src, dir_pos, _ID_DIR),
        "attachmentCount": _entry_count(src, att_pos, _ID_ATTDIR),
        "hasAttachments": att_pos > 0,
        "metadataPosition": meta_pos if meta_pos > 0 else None,
    }

    # The image keys are always present, None when the metadata segment is
    # missing or unreadable, so consumers see one shape of record.
    xml_fields: Dict[str, object] = {k: None for k in _IMAGE_KEYS}
    xml_fields["channels"] = []
    meta_ok = False
    if meta_pos > 0 and (src.size is None or meta_pos + _SEG_HEADER + _META_HEADER_PAYLOAD <= src.size):
        mh = src.read_at(meta_pos, _SEG_HEADER + _META_HEADER_PAYLOAD)
        if len(mh) == _SEG_HEADER + _META_HEADER_PAYLOAD and _segment_id(mh) == _ID_META:
            (xml_size,) = struct.unpack_from("<i", mh, _SEG_HEADER)
            if xml_size > 0:
                xml = src.read_at(meta_pos + _SEG_HEADER + _META_HEADER_PAYLOAD, min(xml_size, _XML_CAP))
                if xml:
                    meta_ok = True
                    xml_fields.update(_xml_fields(xml))
                    if len(xml) < xml_size:
                        xml_fields["xmlTruncated"] = True
                    fields["metadataBytes"] = xml_size
    fields["hasMetadata"] = meta_ok
    fields.update(xml_fields)

    pt = fields.get("pixelType")
    dtype, spp = _PIXEL_TYPE.get(pt or "", (None, None))
    fields["dtype"] = dtype
    fields["samplesPerPixel"] = spp

    return DecodedHeader(format="czi", summary=_summary(fields), fields=fields)


def _summary(f: Dict[str, object]) -> str:
    parts = ["CZI image"]
    x, y, z = f.get("sizeX"), f.get("sizeY"), f.get("sizeZ")
    if x and y:
        dims = [x, y] + ([z] if z and z > 1 else [])
        parts[0] += " " + "x".join(str(d) for d in dims)
    elif not f.get("hasMetadata"):
        parts[0] = "CZI image (no metadata segment)"
    for key, label in (("sizeC", "channels"), ("sizeT", "timepoints"), ("scenes", "scenes"),
                       ("tiles", "tiles"), ("views", "views"), ("illuminations", "illuminations")):
        v = f.get(key)
        if v and v > 1:
            parts.append("%d %s" % (v, label))
    if f.get("dtype"):
        parts.append(str(f["dtype"]))
    vs = f.get("voxelSizeUm") or []
    if vs and vs[0] is not None and vs[1] is not None:
        parts.append("x".join("%.4g" % v for v in vs if v is not None) + " um")
    if f.get("filePart"):
        parts.append("part %d" % f["filePart"])
    return ", ".join(parts)


register_decoder("czi", decode_czi, read=Read(seek=_BUDGET))
