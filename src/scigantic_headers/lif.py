"""Decode a Leica LIF (Leica Image File) container header.

A LIF file is a small binary preamble, one UTF-16 XML document describing every
image in the experiment, then the pixel memory blocks. Everything this library
reports is in that XML, which sits at a fixed offset, so one seek-free read of
the leading bytes is enough; the read is bounded because an experiment with
hundreds of images writes megabytes of XML.

    0   int32   0x70                marker
    4   int32   chunk length        bytes that follow this field
    8   uint8   0x2A                test byte
    9   int32   xmlCharCount        UTF-16 code units in the XML
    13  bytes   XML                 xmlCharCount * 2 bytes, UTF-16-LE
    ... memory blocks: int32 0x70, int32 length, 0x2A, memory size, 0x2A,
        int32 description length, UTF-16 description (MemBlock_N), pixel data

The XML root is `LMSDataContainerHeader` (attribute Version). Each `Element`
has a Name; an image element carries `Data/Image/ImageDescription` with
`Channels/ChannelDescription` (DataType 0 integer / 1 float, Resolution bits,
LUTName) and `Dimensions/DimensionDescription` (DimID 1 X, 2 Y, 3 Z, 4 T,
5 lambda, 10 mosaic tile; NumberOfElements; Length = physical extent in
metres across the dimension, so the voxel pitch is Length / (N - 1)). The
hardware attachment (`Attachment[@Name="HardwareSetting"]`) holds an
`ATLConfocalSettingDefinition` with the objective, magnification and NA.

Only the first image is detailed; the count and the first few names are
reported so a multi-image experiment is visible as such.
"""

from __future__ import annotations

import math
import re
import struct
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

from .decoders import DecodedHeader, Read, Source, register_decoder

# The XML of a large LIF experiment runs to megabytes; read what fits this.
LIF_READ_BUDGET = 8 * 1024 * 1024
_PREAMBLE = 13
_MAX_IMAGE_NAMES = 16
_MAX_CHANNELS = 32

_DIM_X, _DIM_Y, _DIM_Z, _DIM_T, _DIM_LAMBDA, _DIM_TILE = 1, 2, 3, 4, 5, 10


def _finite(x: Optional[float]) -> Optional[float]:
    return x if x is not None and math.isfinite(x) else None


def _num(text) -> Optional[float]:
    try:
        return _finite(float(text))
    except (TypeError, ValueError):
        return None


def _int(text) -> Optional[int]:
    try:
        return int(text)
    except (TypeError, ValueError):
        return None


def _parse_xml(text: str):
    """ET root, or None. A truncated document (the XML outgrew the budget) is
    cut back to its last complete `Element` and closed up so the images that
    did fit still parse."""
    try:
        return ET.fromstring(text)
    except ET.ParseError:
        pass
    cut = text.rfind("</Element>")
    if cut < 0:
        return None
    text = text[:cut + len("</Element>")]
    # Close whatever is still open: count opens minus closes per tag, innermost last.
    opens = re.findall(r"<([A-Za-z_][\w.-]*)(?=[\s>])(?![^>]*/>)", text)
    closes = re.findall(r"</([A-Za-z_][\w.-]*)>", text)
    stack: List[str] = []
    ci = 0
    # A cheap reconstruction: walk opens, pop when the next close matches.
    for tag in opens:
        stack.append(tag)
        while ci < len(closes) and stack and stack[-1] == closes[ci]:
            stack.pop()
            ci += 1
    try:
        return ET.fromstring(text + "".join("</%s>" % t for t in reversed(stack)))
    except ET.ParseError:
        return None


def _lif_dtype(ch) -> Optional[str]:
    bits = _int(ch.get("Resolution"))
    kind = _int(ch.get("DataType"))
    if bits is None:
        return None
    if kind == 1:
        return "float32" if bits <= 32 else "float64"
    if bits <= 8:
        return "uint8"
    if bits <= 16:
        return "uint16"
    return "uint32"


def _image_fields(el) -> Dict[str, object]:
    """Fields of one image `Element`."""
    img = el.find("Data/Image")
    desc = img.find("ImageDescription") if img is not None else None
    f: Dict[str, object] = {}
    channels = desc.findall("Channels/ChannelDescription") if desc is not None else []
    dims = desc.findall("Dimensions/DimensionDescription") if desc is not None else []

    sizes: Dict[int, int] = {}
    voxel: Dict[int, Optional[float]] = {}
    for d in dims:
        dim_id = _int(d.get("DimID"))
        n = _int(d.get("NumberOfElements"))
        if dim_id is None or n is None:
            continue
        sizes[dim_id] = n
        length = _num(d.get("Length"))
        if length is not None and n > 1 and (d.get("Unit") or "m") == "m":
            voxel[dim_id] = round(abs(length) / (n - 1) * 1e6, 6)
        else:
            voxel[dim_id] = None

    f["sizeX"] = sizes.get(_DIM_X)
    f["sizeY"] = sizes.get(_DIM_Y)
    f["sizeZ"] = sizes.get(_DIM_Z, 1 if sizes else None)
    f["sizeT"] = sizes.get(_DIM_T, 1 if sizes else None)
    f["sizeC"] = len(channels) or None
    f["channelCount"] = f["sizeC"]
    f["tileCount"] = sizes.get(_DIM_TILE)
    f["lambdaCount"] = sizes.get(_DIM_LAMBDA)
    f["dtype"] = _lif_dtype(channels[0]) if channels else None
    f["voxelSizeUm"] = [voxel.get(_DIM_X), voxel.get(_DIM_Y), voxel.get(_DIM_Z)]
    names = []
    for i, ch in enumerate(channels[:_MAX_CHANNELS]):
        names.append(ch.get("LUTName") or ch.get("NameOfMeasuredQuantity") or "Ch%d" % (i + 1))
    f["channels"] = names

    hw = None
    for att in (img.findall("Attachment") if img is not None else []):
        if att.get("Name") == "HardwareSetting":
            hw = att
            break
    f["software"] = hw.get("Software") if hw is not None else None
    f["microscope"] = None
    f["objective"] = f["magnification"] = f["numericalAperture"] = None
    f["immersion"] = f["scanMode"] = None
    if hw is not None:
        f["microscope"] = hw.get("SystemTypeName") or None
        for s in hw.iter("ATLConfocalSettingDefinition"):
            if not s.get("ObjectiveName"):
                continue
            f["objective"] = s.get("ObjectiveName")
            f["magnification"] = _num(s.get("Magnification"))
            f["numericalAperture"] = _num(s.get("NumericalAperture"))
            f["immersion"] = s.get("Immersion") or None
            f["scanMode"] = s.get("ScanMode") or None
            f["microscope"] = s.get("MicroscopeModel") or f["microscope"]
            break
    return f


def _summary(f: Dict[str, object]) -> str:
    parts = ["Leica LIF"]
    if f.get("sizeX") and f.get("sizeY"):
        dims = "%dx%d" % (f["sizeX"], f["sizeY"])
        if (f.get("sizeZ") or 1) > 1:
            dims += "x%d" % f["sizeZ"]
        parts[0] += " " + dims
    if (f.get("sizeC") or 0) > 1:
        parts.append("%d channels" % f["sizeC"])
    if (f.get("sizeT") or 0) > 1:
        parts.append("%d timepoints" % f["sizeT"])
    if f.get("dtype"):
        parts.append(str(f["dtype"]))
    vs = [v for v in (f.get("voxelSizeUm") or []) if v is not None]
    if len(vs) >= 2:
        parts.append("x".join("%g" % v for v in vs) + " um")
    if (f.get("tileCount") or 0) > 1:
        parts.append("%d tiles" % f["tileCount"])
    if (f.get("imageCount") or 0) > 1:
        parts.append("%d images" % f["imageCount"])
    return ", ".join(parts)


def decode_lif(src: Source) -> Optional[DecodedHeader]:
    head = src.read_at(0, _PREAMBLE)
    if len(head) < _PREAMBLE:
        return None
    marker, _chunk_len, test, n_chars = struct.unpack_from("<iiBi", head, 0)
    if marker != 0x70 or test != 0x2A or n_chars <= 0:
        return None
    want = min(2 * n_chars, LIF_READ_BUDGET - _PREAMBLE)
    raw = src.read_at(_PREAMBLE, want)
    if len(raw) < 32:
        return None
    text = raw[: len(raw) - (len(raw) % 2)].decode("utf-16-le", "replace")
    root = _parse_xml(text)
    if root is None or root.tag != "LMSDataContainerHeader":
        return None

    images = [el for el in root.iter("Element") if el.find("Data/Image") is not None]
    fields: Dict[str, object] = {
        "lifVersion": _int(root.get("Version")),
        "imageCount": len(images),
        "images": [el.get("Name") for el in images[:_MAX_IMAGE_NAMES]],
        "xmlTruncated": len(raw) < 2 * n_chars,
        "imageName": images[0].get("Name") if images else None,
    }
    if images:
        fields.update(_image_fields(images[0]))
    return DecodedHeader(format="lif", summary=_summary(fields), fields=fields)


register_decoder("lif", decode_lif, read=Read(seek=LIF_READ_BUDGET))
