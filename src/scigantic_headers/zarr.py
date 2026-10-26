"""Decode the JSON sidecars of chunked array stores: Zarr v2 (`.zattrs`,
`.zarray`, `.zgroup`), Zarr v3 (`zarr.json`), OME-Zarr / OME-NGFF on either,
and N5 (`attributes.json`).

A store is a directory, not a file, so there is no header to read: the
metadata is a small JSON document at the root of the store or of each array.
These decoders are dispatched by exact file name (see
decoders.register_decoder_for_name), and `sources.decode_file` on a directory
tries the registered names in the order they are registered here, so pointing
it at `image.zarr/` or `volume.n5/` decodes the store.

What each file carries:

  zarr.json       Zarr v3: zarr_format 3, node_type "group" | "array"; an array
                  node has shape, data_type, chunk_grid.configuration.chunk_shape,
                  codecs; OME-NGFF 0.5 keeps its model under attributes.ome, and
                  a store written with consolidated metadata keeps the whole
                  tree here as a flat path map.
  .zmetadata      Zarr v2 consolidated metadata: every node's `.zattrs`,
                  `.zgroup` and `.zarray` in one document. Claimed here only
                  for an NWB store, which it describes in full.
  .zattrs         Zarr v2 group attributes. OME-NGFF 0.1-0.4 keeps multiscales
                  (axes, datasets with coordinateTransformations), omero
                  (channels), plus plate / well / labels / bioformats2raw.layout.
  .zarray         Zarr v2 array: shape, chunks, dtype ("<u2"), compressor,
                  order, dimension_separator, fill_value.
  .zgroup         Zarr v2 group marker: zarr_format 2.
  attributes.json N5: the root holds {"n5": "2.0.0"}; a dataset holds
                  dimensions, blockSize, dataType, compression.type, and often
                  pixelResolution {dimensions, unit} or transform {scale, units}.

Formats returned: "nwb-zarr" when the root group is an NWB file (see
nwb_zarr.py), "ome-zarr" when an OME multiscales model is present,
"zarr-array" for an array node, "zarr-group" for a bare group, "n5" for N5.

Units: OME-NGFF axis units are UDUNITS names ("micrometer", "nanometer");
voxelSizeUm converts the level-0 scale of each space axis to micrometres when
the unit is a known length, else leaves None.
"""

from __future__ import annotations

import json
import math
from typing import Dict, List, Optional

from . import nwb_zarr
from .decoders import DecodedHeader, Read, register_decoder_for_name

_READ = Read(leading=1024 * 1024)
# A v2 store's `.zmetadata` holds every node's metadata in one document and
# must parse whole, so it gets a bigger read than the per-node sidecars: a
# real 3170-node NWB store from DANDI is 1.4 MB. A store past this declines
# here and falls back to the root `.zattrs`, which still names the file.
_ZMETADATA_READ = Read(leading=8 * 1024 * 1024)

# UDUNITS-2 length names (OME-NGFF axes) and the common abbreviations, to micrometres.
_TO_UM = {
    "micrometer": 1.0, "micrometre": 1.0, "micron": 1.0, "um": 1.0, "µm": 1.0, "μm": 1.0,
    "nanometer": 1e-3, "nanometre": 1e-3, "nm": 1e-3,
    "millimeter": 1e3, "millimetre": 1e3, "mm": 1e3,
    "centimeter": 1e4, "centimetre": 1e4, "cm": 1e4,
    "meter": 1e6, "metre": 1e6, "m": 1e6,
    "picometer": 1e-6, "picometre": 1e-6, "pm": 1e-6,
    "angstrom": 1e-4, "Å": 1e-4,
}

# Zarr v2 dtype strings ("<u2", "|u1", ">f4") -> names used across this library.
_V2_KIND = {"u": "uint", "i": "int", "f": "float", "b": "bool", "c": "complex"}

# Zarr v3 / N5 type names are already the names this library uses, bar these.
_TYPE_ALIAS = {"float": "float32", "double": "float64", "bool": "bool"}


def _finite(x) -> Optional[float]:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _int(x) -> Optional[int]:
    return x if isinstance(x, int) and not isinstance(x, bool) and 0 <= x < 1 << 62 else None


def _int_list(xs) -> Optional[List[int]]:
    if not isinstance(xs, list) or len(xs) > 32:
        return None
    out = [_int(x) for x in xs]
    return out if all(v is not None for v in out) else None


def _load(data: bytes):
    """The JSON object in `data`, or None if it is not a JSON object. The read
    is bounded, so a document longer than the read is cut at the end of the
    last complete top-level object where possible."""
    text = data.decode("utf-8", "replace").lstrip("﻿")
    try:
        obj = json.loads(text)
    except ValueError:
        return None
    return obj if isinstance(obj, dict) else None


def v2_dtype_name(descr) -> Optional[str]:
    """'<u2' -> 'uint16', '|b1' -> 'bool', '<f8' -> 'float64'. None if unknown."""
    if not isinstance(descr, str) or len(descr) < 2:
        return None
    s = descr[1:] if descr[0] in "<>|=" else descr
    kind, size = s[0], s[1:]
    if kind in _V2_KIND and size.isdigit():
        return "bool" if kind == "b" else "%s%d" % (_V2_KIND[kind], int(size) * 8)
    return descr


def _fill_value(v):
    """A JSON-safe fill value: numbers pass through (NaN/inf become the strings
    Zarr itself uses), strings pass through, anything else is stringified."""
    if v is None or isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return v if _finite(v) is not None else str(v)
    return v if isinstance(v, str) else str(v)


def _type_name(t) -> Optional[str]:
    if not isinstance(t, str) or not t:
        return None
    return _TYPE_ALIAS.get(t, t)


def _dims_summary(shape: Optional[List[int]], axes: Optional[List[str]]) -> str:
    if not shape:
        return ""
    dims = " x ".join(str(d) for d in shape)
    if axes and len(axes) == len(shape):
        dims += " (%s)" % "".join(a if len(a) == 1 else "[%s]" % a for a in axes)
    return dims


# ── OME-NGFF multiscales (shared by .zattrs and zarr.json) ─────────────────────

def _ome_fields(ome: dict) -> Optional[Dict[str, object]]:
    """Fields from an OME-NGFF attribute block ({multiscales, omero, plate, ...}),
    or None if it carries no multiscales. Only the first multiscale is detailed."""
    ms_list = ome.get("multiscales")
    if not isinstance(ms_list, list) or not ms_list or not isinstance(ms_list[0], dict):
        return None
    ms = ms_list[0]
    version = ome.get("version") or ms.get("version")

    axes_raw = ms.get("axes")
    axes: List[str] = []
    types: List[Optional[str]] = []
    units: List[Optional[str]] = []
    if isinstance(axes_raw, list):
        for a in axes_raw[:16]:
            if isinstance(a, dict):
                axes.append(str(a.get("name", "")))
                types.append(a.get("type") if isinstance(a.get("type"), str) else None)
                units.append(a.get("unit") if isinstance(a.get("unit"), str) else None)
            elif isinstance(a, str):  # NGFF 0.3 wrote axes as plain names
                axes.append(a)
                types.append("space" if a in ("x", "y", "z") else "channel" if a == "c" else "time" if a == "t" else None)
                units.append(None)

    datasets = ms.get("datasets") if isinstance(ms.get("datasets"), list) else []
    datasets = [d for d in datasets if isinstance(d, dict)]
    level0 = datasets[0] if datasets else {}
    scale = None
    for ct in level0.get("coordinateTransformations", []) if isinstance(level0.get("coordinateTransformations"), list) else []:
        if isinstance(ct, dict) and ct.get("type") == "scale" and isinstance(ct.get("scale"), list):
            scale = [_finite(v) for v in ct["scale"][:16]]
            break

    # Voxel size per space axis x, y, z in micrometres, when the unit is a length.
    voxel: List[Optional[float]] = [None, None, None]
    if scale is not None and axes:
        for i, name in enumerate(axes):
            if i < len(scale) and name in ("x", "y", "z") and scale[i] is not None:
                unit = units[i] if i < len(units) else None
                factor = _TO_UM.get(unit) if unit else None
                if factor is not None:
                    voxel["xyz".index(name)] = round(scale[i] * factor, 6)

    channels: List[Optional[str]] = []
    omero = ome.get("omero")
    if isinstance(omero, dict) and isinstance(omero.get("channels"), list):
        for ch in omero["channels"][:32]:
            if isinstance(ch, dict):
                label = ch.get("label")
                channels.append(str(label) if label is not None else None)

    fields: Dict[str, object] = {
        "ngffVersion": str(version) if version is not None else None,
        "axes": axes,
        "axisTypes": types,
        "axisUnits": units,
        "levels": len(datasets),
        "level0Path": str(level0.get("path")) if level0.get("path") is not None else None,
        "level0Scale": scale,
        "voxelSizeUm": voxel,
        "channels": channels,
        "channelCount": len(channels) or None,
        "multiscaleName": ms.get("name") if isinstance(ms.get("name"), str) else None,
        "isPlate": "plate" in ome,
        "isWell": "well" in ome,
        "isLabels": "labels" in ome or "image-label" in ome,
        "isBioformats2raw": isinstance(ome.get("bioformats2raw.layout"), int),
    }
    creator = ome.get("_creator")
    if isinstance(creator, dict):
        fields["creator"] = " ".join(str(creator[k]) for k in ("name", "version") if creator.get(k)) or None
    return fields


def _ome_summary(f: Dict[str, object], shape=None, dtype=None) -> str:
    parts = ["OME-Zarr"]
    nch = f.get("channelCount")
    body = []
    if shape:
        body.append(_dims_summary(shape, f.get("axes")))
    elif nch:
        body.append("%d channels" % nch)
    if body:
        parts[0] += " " + body[0]
    if shape and nch:
        parts.append("%d channels" % nch)
    if f.get("levels"):
        parts.append("%d levels" % f["levels"])
    if dtype:
        parts.append(dtype)
    vs = f.get("voxelSizeUm") or []
    known = [v for v in vs if v is not None]
    if len(known) >= 2:
        parts.append("x".join("%.4g" % v for v in known) + " um")
    return ", ".join(parts)


def with_level0_array(group: DecodedHeader, array: Optional[DecodedHeader]) -> DecodedHeader:
    """An OME-Zarr group header completed with its level-0 array's shape, dtype
    and chunks (the group's own metadata never states them). Used when a whole
    store directory is decoded, so the result reads like one image. Returns the
    group unchanged when the array is missing or not an array."""
    if array is None or array.format != "zarr-array" or group.format != "ome-zarr":
        return group
    fields = dict(group.fields)
    shape = array.fields.get("shape")
    fields.update({
        "shape": shape,
        "dtype": array.fields.get("dtype"),
        "chunks": array.fields.get("chunks"),
        "compressor": array.fields.get("compressor") or array.fields.get("codecs"),
    })
    axes = fields.get("axes")
    if isinstance(shape, list) and isinstance(axes, list) and len(axes) == len(shape):
        for ax, n in zip(axes, shape):
            key = {"x": "sizeX", "y": "sizeY", "z": "sizeZ", "c": "sizeC", "t": "sizeT"}.get(str(ax).lower())
            if key:
                fields[key] = n
    return DecodedHeader(format=group.format, summary=_ome_summary(fields, shape, fields.get("dtype")), fields=fields)


# ── Zarr v2 ───────────────────────────────────────────────────────────────────

def decode_zattrs(data: bytes) -> Optional[DecodedHeader]:
    """A Zarr v2 group's `.zattrs`: NWB when it carries `nwb_version`,
    OME-Zarr when it holds multiscales, else a plain group whose attribute
    keys are listed. A v2 store keeps its consolidated metadata in a separate
    `.zmetadata` beside this file, so an NWB record from here names the file
    and its version but not its structure."""
    obj = _load(data)
    if obj is None:
        return None
    if nwb_zarr.is_nwb(obj):
        return nwb_zarr.decode(obj, None, 2)
    ome = _ome_fields(obj)
    if ome is not None:
        ome["zarrFormat"] = 2
        return DecodedHeader(format="ome-zarr", summary=_ome_summary(ome), fields=ome)
    keys = sorted(str(k) for k in obj)[:64]
    fields = {"zarrFormat": 2, "attributeKeys": keys, "isPlate": "plate" in obj,
              "isLabels": "labels" in obj or "image-label" in obj,
              "isBioformats2raw": isinstance(obj.get("bioformats2raw.layout"), int)}
    return DecodedHeader(format="zarr-group", summary="Zarr group, %d attributes" % len(obj), fields=fields)


def decode_zarray(data: bytes) -> Optional[DecodedHeader]:
    """A Zarr v2 array's `.zarray`."""
    obj = _load(data)
    if obj is None or obj.get("zarr_format") != 2 or "shape" not in obj:
        return None
    shape = _int_list(obj.get("shape"))
    if shape is None:
        return None
    chunks = _int_list(obj.get("chunks"))
    dtype = v2_dtype_name(obj.get("dtype"))
    comp = obj.get("compressor")
    comp_name = None
    if isinstance(comp, dict):
        cid = comp.get("id")
        comp_name = str(cid) if cid is not None else None
        if cid == "blosc" and comp.get("cname"):
            comp_name = "blosc/%s" % comp["cname"]
    filters = obj.get("filters")
    fields = {
        "zarrFormat": 2,
        "shape": shape,
        "ndim": len(shape),
        "chunks": chunks,
        "dtype": dtype,
        "compressor": comp_name,
        "filters": [str(f.get("id")) for f in filters if isinstance(f, dict)][:8] if isinstance(filters, list) else None,
        "order": obj.get("order") if obj.get("order") in ("C", "F") else None,
        "dimensionSeparator": obj.get("dimension_separator") if isinstance(obj.get("dimension_separator"), str) else ".",
        "fillValue": _fill_value(obj.get("fill_value")),
    }
    summary = "Zarr array %s, %s%s" % (
        " x ".join(str(d) for d in shape) or "scalar", dtype or "unknown dtype",
        (", chunks %s" % " x ".join(str(c) for c in chunks)) if chunks else "")
    return DecodedHeader(format="zarr-array", summary=summary, fields=fields)


def decode_zgroup(data: bytes) -> Optional[DecodedHeader]:
    """A Zarr v2 group marker `.zgroup` ({"zarr_format": 2})."""
    obj = _load(data)
    if obj is None or obj.get("zarr_format") != 2:
        return None
    return DecodedHeader(format="zarr-group", summary="Zarr v2 group", fields={"zarrFormat": 2})


# ── Zarr v3 ───────────────────────────────────────────────────────────────────

def _consolidated(obj: dict) -> Optional[Dict[str, dict]]:
    """A v3 store written with consolidated metadata keeps the whole tree in
    the root's own file, as a flat map of slash-separated path to node. That
    is what lets one bounded read describe a session with hundreds of
    objects; without it, only the root's own attributes are here."""
    block = obj.get("consolidated_metadata")
    if not isinstance(block, dict):
        return None
    metadata = block.get("metadata")
    if not isinstance(metadata, dict):
        return None
    return {str(k): v for k, v in metadata.items() if isinstance(v, dict)}

def _v3_chunk_shape(obj: dict) -> Optional[List[int]]:
    grid = obj.get("chunk_grid")
    if isinstance(grid, dict) and isinstance(grid.get("configuration"), dict):
        return _int_list(grid["configuration"].get("chunk_shape"))
    return None


def _v3_codecs(obj: dict) -> List[str]:
    codecs = obj.get("codecs")
    if not isinstance(codecs, list):
        return []
    names = []
    for c in codecs[:8]:
        if isinstance(c, dict) and isinstance(c.get("name"), str):
            name = c["name"]
            cfg = c.get("configuration")
            if name == "blosc" and isinstance(cfg, dict) and cfg.get("cname"):
                name = "blosc/%s" % cfg["cname"]
            names.append(name)
    return names


def decode_zarr_json(data: bytes) -> Optional[DecodedHeader]:
    """A Zarr v3 node's `zarr.json`: array or group, OME-Zarr when attributes.ome
    (OME-NGFF 0.5) or legacy top-level attributes carry multiscales."""
    obj = _load(data)
    if obj is None or obj.get("zarr_format") != 3:
        return None
    node = obj.get("node_type")
    attrs = obj.get("attributes") if isinstance(obj.get("attributes"), dict) else {}
    ome_block = attrs.get("ome") if isinstance(attrs.get("ome"), dict) else attrs
    ome = _ome_fields(ome_block) if isinstance(ome_block, dict) else None

    if node == "array":
        shape = _int_list(obj.get("shape"))
        if shape is None:
            return None
        dtype = _type_name(obj.get("data_type"))
        chunks = _v3_chunk_shape(obj)
        codecs = _v3_codecs(obj)
        fields: Dict[str, object] = {
            "zarrFormat": 3, "shape": shape, "ndim": len(shape), "chunks": chunks,
            "dtype": dtype, "codecs": codecs,
            "compressor": next((c for c in codecs if c.split("/")[0] in ("blosc", "gzip", "zstd", "lz4", "bz2", "zlib")), None),
            "dimensionNames": [str(n) for n in obj["dimension_names"][:16]] if isinstance(obj.get("dimension_names"), list) else None,
            "fillValue": _fill_value(obj.get("fill_value")),
        }
        if ome is not None:
            fields.update(ome)
            return DecodedHeader(format="ome-zarr", summary=_ome_summary(ome, shape, dtype), fields=fields)
        summary = "Zarr v3 array %s, %s%s" % (
            " x ".join(str(d) for d in shape) or "scalar", dtype or "unknown dtype",
            (", chunks %s" % " x ".join(str(c) for c in chunks)) if chunks else "")
        return DecodedHeader(format="zarr-array", summary=summary, fields=fields)

    if node != "group":
        return None
    if nwb_zarr.is_nwb(attrs):
        return nwb_zarr.decode(attrs, _consolidated(obj), 3)
    if ome is not None:
        ome["zarrFormat"] = 3
        return DecodedHeader(format="ome-zarr", summary=_ome_summary(ome), fields=ome)
    keys = sorted(str(k) for k in attrs)[:64]
    return DecodedHeader(format="zarr-group", summary="Zarr v3 group, %d attributes" % len(attrs),
                         fields={"zarrFormat": 3, "attributeKeys": keys})


def decode_zmetadata(data: bytes) -> Optional[DecodedHeader]:
    """A Zarr v2 store's consolidated `.zmetadata`. Only an NWB store is
    claimed here: for anything else the per-node sidecars already say
    everything this document would, and the root `.zattrs` is what the
    directory dispatch falls through to."""
    obj = _load(data)
    if obj is None:
        return None
    folded = nwb_zarr.normalize_v2(obj)
    if folded is None:
        return None
    attrs, nodes = folded
    if not nwb_zarr.is_nwb(attrs):
        return None
    return nwb_zarr.decode(attrs, nodes, 2)


# ── N5 ────────────────────────────────────────────────────────────────────────

def decode_n5_attributes(data: bytes) -> Optional[DecodedHeader]:
    """An N5 `attributes.json`: the container root ({"n5": version}) or a
    dataset (dimensions, blockSize, dataType, compression), with the voxel
    size from pixelResolution or a transform when present."""
    obj = _load(data)
    if obj is None:
        return None
    is_root = isinstance(obj.get("n5"), str)
    dims = _int_list(obj.get("dimensions")) if "dimensions" in obj else None
    if not is_root and dims is None:
        return None

    fields: Dict[str, object] = {"n5Version": obj.get("n5") if is_root else None}
    if dims is not None:
        dtype = _type_name(obj.get("dataType"))
        comp = obj.get("compression")
        comp_name = None
        if isinstance(comp, dict) and comp.get("type") is not None:
            comp_name = str(comp["type"])
        elif isinstance(obj.get("compressionType"), str):  # N5 v1 wrote a string
            comp_name = obj["compressionType"]
        fields.update({
            "shape": dims, "ndim": len(dims),
            "blockSize": _int_list(obj.get("blockSize")),
            "dtype": dtype, "compression": comp_name,
        })
        # N5 dimensions are in Fortran (x fastest) order: dimensions[0] is x.
        voxel: List[Optional[float]] = [None, None, None]
        axes = None
        pr = obj.get("pixelResolution")
        tr = obj.get("transform")
        if isinstance(pr, dict) and isinstance(pr.get("dimensions"), list):
            unit = pr.get("unit") if isinstance(pr.get("unit"), str) else None
            factor = _TO_UM.get(unit) if unit else None
            for i, v in enumerate(pr["dimensions"][:3]):
                fv = _finite(v)
                if fv is not None and factor is not None:
                    voxel[i] = round(fv * factor, 6)
            fields["pixelResolution"] = [_finite(v) for v in pr["dimensions"][:3]]
            fields["pixelResolutionUnit"] = unit
        elif isinstance(tr, dict) and isinstance(tr.get("scale"), list) and isinstance(tr.get("axes"), list):
            axes = [str(a) for a in tr["axes"][:8]]
            units = tr.get("units") if isinstance(tr.get("units"), list) else []
            for i, name in enumerate(axes):
                if name in ("x", "y", "z") and i < len(tr["scale"]):
                    fv = _finite(tr["scale"][i])
                    unit = units[i] if i < len(units) and isinstance(units[i], str) else None
                    factor = _TO_UM.get(unit) if unit else None
                    if fv is not None and factor is not None:
                        voxel["xyz".index(name)] = round(fv * factor, 6)
        fields["voxelSizeUm"] = voxel
        if isinstance(tr, dict) and isinstance(tr.get("axes"), list):
            fields["axes"] = [str(a) for a in tr["axes"][:8]]
        if isinstance(obj.get("downsamplingFactors"), list):
            fields["downsamplingFactors"] = _int_list(obj["downsamplingFactors"])

    # BDV-N5 / multiscale markers at the root or on a setup group.
    for key in ("multiScale", "resolution", "setups", "scales"):
        if key in obj and key not in fields and isinstance(obj[key], (list, str, bool, int)):
            fields[key] = obj[key] if not isinstance(obj[key], list) else obj[key][:16]
    other = sorted(k for k in obj if k not in ("n5", "dimensions", "blockSize", "dataType", "compression",
                                                "pixelResolution", "transform", "downsamplingFactors"))[:32]
    fields["attributeKeys"] = other

    if dims is not None:
        summary = "N5 dataset %s, %s" % (" x ".join(str(d) for d in dims), fields.get("dtype") or "unknown dtype")
        bs = fields.get("blockSize")
        if bs:
            summary += ", blocks %s" % " x ".join(str(b) for b in bs)
        vs = fields.get("voxelSizeUm") or []
        if all(v is not None for v in vs):
            summary += ", %s um" % "x".join("%.4g" % v for v in vs)
    else:
        summary = "N5 container v%s" % fields["n5Version"]
    return DecodedHeader(format="n5", summary=summary, fields=fields)


# Registration order is the order sources.decode_file tries them in a directory:
# the v3 document first, then the v2 group attributes (OME model), the array,
# the bare group marker, and the N5 sidecar.
register_decoder_for_name("zarr.json", decode_zarr_json, read=_READ)
# Registered before `.zattrs`: `sources._decode_directory` tries the sidecar
# names in registration order, so a consolidated NWB store is described in
# full and anything else falls through to the per-node files below.
register_decoder_for_name(".zmetadata", decode_zmetadata, read=_ZMETADATA_READ)
register_decoder_for_name(".zattrs", decode_zattrs, read=_READ)
register_decoder_for_name(".zarray", decode_zarray, read=_READ)
register_decoder_for_name(".zgroup", decode_zgroup, read=_READ)
register_decoder_for_name("attributes.json", decode_n5_attributes, read=_READ)
