"""Decode the XML sidecars that describe a light-microscopy dataset, by root
element.

`.xml` is too generic an extension for any one format to own, so this module
claims it once and routes on the root element found in the leading bytes:

    <SpimData ...>         BigDataViewer / BigStitcher dataset XML   -> bdv-xml
    <OME xmlns=...openmicroscopy.org...>   OME-XML sidecar / companion -> ome-xml
    <MicroscopeImage ...>  Thermo Fisher EPU image XML (epu.py)       -> epu-xml
    anything else          None (no context)

Exact-name registrations (Illumina RunInfo.xml / RunParameters.xml) win over
the extension in the registry, so they never reach this dispatcher. mzML keeps
its own `.mzml` extension.

BigDataViewer XML (SpimData version 0.2) layout:

    <SpimData version="0.2">
      <BasePath type="relative">.</BasePath>
      <SequenceDescription>
        <ImageLoader format="bdv.hdf5">            bdv.n5, bdv.multimg.zarr, ...
          <hdf5 type="relative">data.h5</hdf5>     or <n5>, <zarr>
        </ImageLoader>
        <ViewSetups>
          <ViewSetup>
            <id>0</id> <name>Angle 0</name>
            <size>1388 1040 81</size>                 x y z, pixels
            <voxelSize><unit>um</unit><size>0.73 0.73 2.0</size></voxelSize>
            <attributes><illumination>0</illumination><channel>0</channel>
                        <angle>0</angle><tile>0</tile></attributes>
          </ViewSetup> ...
          <Attributes name="channel"><Channel><id>0</id><name>GFP</name></Channel></Attributes>
          <Attributes name="angle"> ... <Attributes name="illumination"> ... <Attributes name="tile">
        </ViewSetups>
        <Timepoints type="range"><first>0</first><last>2</last></Timepoints>
                    type="list"><Timepoint>..</Timepoint>...
                    type="pattern"><integerpattern>0-10:2, 18</integerpattern>
        <MissingViews><MissingView timepoint=".." setup=".."/></MissingViews>
      </SequenceDescription>
      <ViewRegistrations> ... </ViewRegistrations>
    </SpimData>

Older exports (the OME drosophila sample) omit <size>, <voxelSize> and the
<Attributes> blocks entirely; those fields are then None.
"""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

from .decoders import DecodedHeader, Read, register_decoder
from .epu import decode_epu_xml
from .ome_xml import ome_summary, parse_ome_xml

# How far into the file the root element must appear. A declaration, a BOM,
# comments and a DOCTYPE precede it; 4 KiB is generous for all of those.
_SNIFF_BYTES = 4096

# Skips a BOM, whitespace, the XML declaration, comments and processing
# instructions, and a DOCTYPE, then captures the root element's local name.
_ROOT = re.compile(
    rb"\A(?:\xef\xbb\xbf)?\s*(?:<\?[^>]*\?>\s*|<!--.*?-->\s*|<!DOCTYPE[^>]*>\s*)*<(?:[\w.-]+:)?([\w.-]+)",
    re.S,
)

# Voxel-size units BDV writers emit, as a factor to micrometres.
_TO_UM = {
    "um": 1.0, "µm": 1.0, "μm": 1.0, "micron": 1.0, "microns": 1.0,
    "micrometer": 1.0, "micrometers": 1.0, "micrometre": 1.0, "micrometres": 1.0,
    "nm": 1e-3, "nanometer": 1e-3, "nanometers": 1e-3,
    "mm": 1e3, "millimeter": 1e3, "millimeters": 1e3,
}


def root_element(data: bytes) -> Optional[str]:
    """Local name of the root element in the leading bytes, or None."""
    m = _ROOT.match(data[:_SNIFF_BYTES])
    return m.group(1).decode("ascii", "replace") if m else None


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def _child(el, name: str):
    if el is None:
        return None
    for c in el:
        if _local(c.tag) == name:
            return c
    return None


def _children(el, name: str):
    return [] if el is None else [c for c in el if _local(c.tag) == name]


def _text(el) -> Optional[str]:
    if el is None or el.text is None:
        return None
    t = el.text.strip()
    return t or None


def _ints(text: Optional[str]) -> Optional[List[int]]:
    if not text:
        return None
    try:
        return [int(float(x)) for x in text.split()]
    except ValueError:
        return None


def _floats(text: Optional[str]) -> Optional[List[float]]:
    if not text:
        return None
    try:
        vals = [float(x) for x in text.split()]
    except ValueError:
        return None
    return vals if all(math.isfinite(v) for v in vals) else None


def _timepoint_count(tp) -> Optional[int]:
    """Number of timepoints from a <Timepoints> element of type range, list or
    pattern. A pattern is comma-separated items 'a', 'a-b' or 'a-b:step'."""
    if tp is None:
        return None
    kind = (tp.get("type") or "").lower()
    if kind == "range":
        first, last = _ints(_text(_child(tp, "first"))), _ints(_text(_child(tp, "last")))
        if first and last and last[0] >= first[0]:
            return last[0] - first[0] + 1
        return None
    if kind == "list":
        return len(_children(tp, "Timepoint")) or None
    if kind == "pattern":
        pat = _text(_child(tp, "integerpattern"))
        if not pat:
            return None
        seen = set()
        for item in pat.split(","):
            item = item.strip()
            if not item:
                continue
            m = re.fullmatch(r"(-?\d+)(?:-(-?\d+)(?::(\d+))?)?", item)
            if not m:
                return None
            a = int(m.group(1))
            b = int(m.group(2)) if m.group(2) is not None else a
            step = int(m.group(3)) if m.group(3) else 1
            if step <= 0 or b < a or (b - a) // step > 1_000_000:
                return None
            seen.update(range(a, b + 1, step))
        return len(seen) or None
    return None


def decode_bdv_xml(data: bytes) -> Optional[DecodedHeader]:
    """A BigDataViewer / BigStitcher dataset XML (root <SpimData>) -> fields.
    None for anything else."""
    try:
        root = ET.fromstring(data)
    except ET.ParseError:
        return None
    if _local(root.tag) != "SpimData":
        return None

    seq = _child(root, "SequenceDescription")
    loader = _child(seq, "ImageLoader")
    loader_format = loader.get("format") if loader is not None else None
    data_file = None
    if loader is not None:
        for name in ("hdf5", "n5", "zarr", "zgroup", "ome.zarr"):
            data_file = _text(_child(loader, name))
            if data_file:
                break
        if data_file is None:
            # Other loaders (e.g. spimreconstruction.stack.ij) name a directory
            # and a file pattern instead.
            data_file = _text(_child(loader, "imagedirectory")) or _text(_child(loader, "filePattern"))

    setups_el = _child(seq, "ViewSetups")
    setups = _children(setups_el, "ViewSetup")
    first = setups[0] if setups else None
    size = _ints(_text(_child(first, "size")))
    vs_el = _child(first, "voxelSize")
    unit = _text(_child(vs_el, "unit"))
    vs = _floats(_text(_child(vs_el, "size")))
    factor = _TO_UM.get(unit.lower()) if unit else None
    voxel_um: Optional[List[Optional[float]]] = None
    if vs and factor is not None:
        voxel_um = [round(v * factor, 6) for v in vs[:3]] + [None] * (3 - len(vs[:3]))
    elif vs:
        voxel_um = [None, None, None]  # a size without a length unit ("pixel") is not a measurement

    def attribute_entries(name: str):
        for a in _children(setups_el, "Attributes"):
            if (a.get("name") or "").lower() == name:
                return list(a)
        return None

    channel_entries = attribute_entries("channel")
    channels = None
    if channel_entries is not None:
        channels = [(_text(_child(c, "name")) or _text(_child(c, "id")) or "") for c in channel_entries[:32]]
    illum = attribute_entries("illumination")
    angles = attribute_entries("angle")
    tiles = attribute_entries("tile")

    timepoints = _timepoint_count(_child(seq, "Timepoints"))
    missing = len(_children(_child(seq, "MissingViews"), "MissingView"))
    regs = _child(root, "ViewRegistrations")
    has_regs = bool(regs is not None and len(_children(regs, "ViewRegistration")))

    fields: Dict[str, object] = {
        "bdvVersion": root.get("version"),
        "imageLoaderFormat": loader_format,
        "dataFile": data_file,
        "viewSetupCount": len(setups) or None,
        "timepointCount": timepoints,
        "sizeX": size[0] if size and len(size) > 0 else None,
        "sizeY": size[1] if size and len(size) > 1 else None,
        "sizeZ": size[2] if size and len(size) > 2 else None,
        "voxelSizeUm": voxel_um,
        "voxelSizeUnit": unit,
        "channels": channels,
        "channelCount": len(channel_entries) if channel_entries is not None else None,
        "illuminationCount": len(illum) if illum is not None else None,
        "angleCount": len(angles) if angles is not None else None,
        "tileCount": len(tiles) if tiles is not None else None,
        "hasViewRegistrations": has_regs,
        "missingViewCount": missing,
    }

    parts = ["BDV dataset"]
    if size and len(size) >= 3:
        parts[0] += " %dx%dx%d" % tuple(size[:3])
    elif fields["viewSetupCount"]:
        parts[0] += " (%d view setups)" % fields["viewSetupCount"]
    if (fields["channelCount"] or 0) > 1:
        parts.append("%d channels" % fields["channelCount"])
    if (fields["angleCount"] or 0) > 1:
        parts.append("%d angles" % fields["angleCount"])
    if (fields["tileCount"] or 0) > 1:
        parts.append("%d tiles" % fields["tileCount"])
    if (timepoints or 0) > 1:
        parts.append("%d timepoints" % timepoints)
    if voxel_um and all(v is not None for v in voxel_um):
        parts.append("x".join("%g" % v for v in voxel_um) + " um")
    if loader_format:
        parts.append(loader_format)
    return DecodedHeader(format="bdv-xml", summary=", ".join(parts), fields=fields)


def decode_ome_xml(data: bytes) -> Optional[DecodedHeader]:
    """An OME-XML document (sidecar or companion file) -> the shared OME fields."""
    fields = parse_ome_xml(data)
    if fields is None:
        return None
    return DecodedHeader(format="ome-xml", summary=ome_summary("OME-XML", fields), fields=fields)


def decode_xml_sidecar(data: bytes) -> Optional[DecodedHeader]:
    """Route an `.xml` file on its root element. None when the root is not one
    this library describes, so an arbitrary XML is 'no context', never an error."""
    root = root_element(data)
    if root == "SpimData":
        return decode_bdv_xml(data)
    if root == "OME":
        return decode_ome_xml(data)
    if root == "MicroscopeImage":
        return decode_epu_xml(data)
    return None


# A BDV XML with hundreds of view registrations runs to a few hundred KB; an
# OME-XML companion for a plate can be larger still. 1 MiB covers both; the
# metadata of interest sits at the front in either case.
register_decoder("xml", decode_xml_sidecar, read=Read(leading=1024 * 1024))
