"""Shared OME-XML reader: the metadata model behind OME-TIFF (in the TIFF
ImageDescription), OME-XML sidecars (`.ome.xml`) and companion files.

One parse serves every carrier, so an OME-TIFF and its sidecar describe an
image with the same keys. Namespace-agnostic (the schema URL changes by
release year) and total: anything that is not an OME document yields None.

Only the first Image is detailed; the count of images is reported so a
multi-image file (a plate, a multi-position acquisition) is visible as such.

Physical sizes are returned in micrometres whatever unit the file uses
(PhysicalSizeXUnit defaults to the micron sign in the schema).
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

# OME Pixels Type -> the dtype names used across this library.
OME_PIXEL_TYPE = {
    "int8": "int8", "int16": "int16", "int32": "int32",
    "uint8": "uint8", "uint16": "uint16", "uint32": "uint32",
    "float": "float32", "double": "float64", "bit": "bit",
    "complex": "complex64", "double-complex": "complex128",
}

# Length units the schema allows for PhysicalSize*Unit, as a factor to micrometres.
_TO_UM = {
    "µm": 1.0, "um": 1.0, "μm": 1.0, "micron": 1.0, "microns": 1.0, "micrometer": 1.0,
    "nm": 1e-3, "mm": 1e3, "cm": 1e4, "m": 1e6, "pm": 1e-6, "Å": 1e-4, "angstrom": 1e-4,
}


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def _child(el, name: str):
    for c in el:
        if _local(c.tag) == name:
            return c
    return None


def _children(el, name: str):
    return [c for c in el if _local(c.tag) == name]


def _num(text) -> Optional[float]:
    try:
        v = float(text)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _int(text) -> Optional[int]:
    try:
        return int(text)
    except (TypeError, ValueError):
        return None


def _um(value, unit) -> Optional[float]:
    v = _num(value)
    if v is None:
        return None
    factor = _TO_UM.get((unit or "µm").strip(), None)
    if factor is None:
        return None
    return round(v * factor, 6)


def is_ome_xml(data: bytes) -> bool:
    """Cheap sniff: an OME root element somewhere in the leading bytes."""
    head = data[:4096]
    return b"<OME" in head and b"openmicroscopy.org" in head


def parse_ome_xml(data) -> Optional[Dict[str, object]]:
    """OME-XML (bytes or str) -> flat fields for the first Image, or None if
    the document is not OME. Keys: imageCount, imageName, sizeX..sizeT, dtype,
    dimensionOrder, voxelSizeUm [x, y, z], channels, channelCount,
    acquisitionDate, microscope, objective, magnification, numericalAperture,
    significantBits, creator."""
    if isinstance(data, str):
        data = data.encode("utf-8", "replace")
    # Trim anything after the closing root tag (a TIFF ImageDescription may be
    # NUL-padded) and tolerate a truncated tail by cutting at the last '>'.
    end = data.rfind(b"</OME>")
    if end >= 0:
        data = data[:end + len(b"</OME>")]
    try:
        root = ET.fromstring(data)
    except ET.ParseError:
        return None
    if _local(root.tag) != "OME":
        return None

    images = _children(root, "Image")
    fields: Dict[str, object] = {
        "imageCount": len(images),
        "creator": root.get("Creator"),
    }
    if not images:
        return fields

    img = images[0]
    fields["imageName"] = img.get("Name")
    acq = _child(img, "AcquisitionDate")
    fields["acquisitionDate"] = (acq.text or "").strip() or None if acq is not None else None

    px = _child(img, "Pixels")
    channels: List[Optional[str]] = []
    if px is not None:
        for k in ("X", "Y", "Z", "C", "T"):
            fields["size" + k] = _int(px.get("Size" + k))
        fields["dtype"] = OME_PIXEL_TYPE.get(px.get("Type") or "", px.get("Type"))
        fields["dimensionOrder"] = px.get("DimensionOrder")
        fields["significantBits"] = _int(px.get("SignificantBits"))
        fields["voxelSizeUm"] = [
            _um(px.get("PhysicalSize" + k), px.get("PhysicalSize" + k + "Unit")) for k in ("X", "Y", "Z")
        ]
        for ch in _children(px, "Channel")[:32]:
            channels.append(ch.get("Name") or ch.get("ID"))
        emission = [_num(ch.get("EmissionWavelength")) for ch in _children(px, "Channel")[:32]]
        if any(e is not None for e in emission):
            fields["emissionWavelengthsNm"] = emission
    fields["channels"] = channels
    fields["channelCount"] = fields.get("sizeC") or (len(channels) or None)

    # Instrument: the one the image references, else the first.
    instruments = _children(root, "Instrument")
    ref = _child(img, "InstrumentRef")
    inst = None
    if ref is not None:
        inst = next((i for i in instruments if i.get("ID") == ref.get("ID")), None)
    if inst is None and instruments:
        inst = instruments[0]
    microscope = objective = None
    mag = na = None
    if inst is not None:
        m = _child(inst, "Microscope")
        if m is not None:
            microscope = " ".join(x for x in (m.get("Manufacturer"), m.get("Model")) if x) or m.get("Type")
        objs = _children(inst, "Objective")
        oref = _child(img, "ObjectiveSettings")
        obj = None
        if oref is not None:
            obj = next((o for o in objs if o.get("ID") == oref.get("ID")), None)
        if obj is None and objs:
            obj = objs[0]
        if obj is not None:
            objective = " ".join(x for x in (obj.get("Manufacturer"), obj.get("Model")) if x) or obj.get("ID")
            mag = _num(obj.get("NominalMagnification"))
            na = _num(obj.get("LensNA"))
    fields["microscope"] = microscope
    fields["objective"] = objective
    fields["magnification"] = mag
    fields["numericalAperture"] = na
    return fields


def ome_summary(kind: str, f: Dict[str, object]) -> str:
    """One line: '<kind> 512x512x30, 3 channels, 10 timepoints, uint16, 0.1x0.1x0.5 um'."""
    dims = [f.get("sizeX"), f.get("sizeY")]
    if (f.get("sizeZ") or 1) > 1:
        dims.append(f.get("sizeZ"))
    parts = [kind]
    if dims[0] and dims[1]:
        parts[0] += " " + "x".join(str(d) for d in dims)
    if (f.get("sizeC") or 0) > 1:
        parts.append("%d channels" % f["sizeC"])
    if (f.get("sizeT") or 0) > 1:
        parts.append("%d timepoints" % f["sizeT"])
    if f.get("dtype"):
        parts.append(str(f["dtype"]))
    vs = f.get("voxelSizeUm") or []
    if vs and all(v is not None for v in vs[:2]):
        parts.append("x".join("%g" % v for v in vs if v is not None) + " um")
    if (f.get("imageCount") or 0) > 1:
        parts.append("%d images" % f["imageCount"])
    return ", ".join(parts)
