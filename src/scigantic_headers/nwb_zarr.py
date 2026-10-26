"""NWB (Neurodata Without Borders) stored in a Zarr store, not HDF5.

`hdmf-zarr` writes the same NWB model into a Zarr directory: the root group's
attributes carry `nwb_version`, `namespace`, `neurodata_type: NWBFile` and
`object_id`, and every neurodata object below is a Zarr group or array whose
own attributes carry its type, unit and rate. A store written with
consolidated metadata keeps that whole tree in the root's own metadata file,
so one bounded read describes the session.

**What this cannot report, and why.** In HDF5 a scalar string sits in the
object header, so `hdf5.py` reads the identifier, the session description and
the subject straight out of the metadata. In Zarr a string dataset's value is
chunk data, usually behind a compressor this library does not carry (the
sample this was built against uses zstd). No chunk is ever opened here, the
same rule the rest of `zarr.py` follows, so those values are absent rather
than guessed and `stringValuesUnread` names the paths a reader would have to
open to get them. Everything that lives in metadata, which is the structure,
the shapes and dtypes, the units and the sampling rates, is reported in full.

Dispatched from `zarr.py`, which owns the `zarr.json` and `.zattrs` names.
"""

import math
from typing import Dict, List, Optional, Tuple

from .decoders import DecodedHeader

# Caps mirroring hdf5.py's NWB record, so the two are comparable.
NWB_LIST_CAP = 32
NWB_SERIES_DETAIL = 16
_PATH_CAP = 32

_CORE_NAMESPACES = ("core", "hdmf-common", "hdmf-experimental")


def normalize_v2(consolidated: dict) -> Optional[Tuple[dict, Dict[str, dict]]]:
    """A Zarr v2 store keeps its consolidated metadata in `.zmetadata`, keyed
    by `<path>/.zattrs`, `<path>/.zgroup` and `<path>/.zarray` rather than by
    node. Fold those back into the one-node-per-path map the v3 reader uses,
    so both container versions walk the same way. Returns the root attributes
    and the map, or None when the document is not a consolidated v2 store."""
    if not isinstance(consolidated, dict):
        return None
    metadata = consolidated.get("metadata")
    if not isinstance(metadata, dict):
        return None
    nodes: Dict[str, dict] = {}
    for key, value in metadata.items():
        if not isinstance(key, str) or not isinstance(value, dict):
            continue
        head, slash, tail = key.rpartition("/")
        if tail not in (".zattrs", ".zgroup", ".zarray"):
            continue
        path = head if slash else ""
        node = nodes.setdefault(path, {})
        if tail == ".zattrs":
            node["attributes"] = value
        elif tail == ".zgroup":
            node.setdefault("node_type", "group")
        else:
            node["node_type"] = "array"
            node["shape"] = value.get("shape")
            node["data_type"] = _v2_dtype(value)
    root = nodes.pop("", {})
    return (root.get("attributes") or {}), nodes


def _v2_dtype(zarray: dict):
    """A v2 dtype descriptor ('<i8', '|O') as a name. An object dtype with a
    vlen-utf8 filter is how zarr v2 stores strings, which is what the v3
    `"string"` data type says directly."""
    from .zarr import v2_dtype_name  # imported here: zarr.py imports this module
    descr = zarray.get("dtype")
    if descr in ("|O", "O", "object"):
        for f in zarray.get("filters") or []:
            if isinstance(f, dict) and str(f.get("id", "")).startswith("vlen-utf8"):
                return "string"
        return "object"
    return v2_dtype_name(descr)


def is_nwb(attributes: dict) -> bool:
    """Whether a Zarr group's attributes are an NWB file's root."""
    if not isinstance(attributes, dict):
        return False
    return "nwb_version" in attributes or attributes.get("neurodata_type") == "NWBFile"


def _children(metadata: Dict[str, dict], prefix: str) -> List[str]:
    """Immediate child names of a path in the flat consolidated map, in the
    order the writer stored them."""
    out: List[str] = []
    head = prefix + "/" if prefix else ""
    for path in metadata:
        if not path.startswith(head):
            continue
        rest = path[len(head):]
        if not rest or "/" in rest:
            continue
        out.append(rest)
    return out


def _node(metadata: Dict[str, dict], path: str) -> Optional[dict]:
    node = metadata.get(path)
    return node if isinstance(node, dict) else None


def _safe(value):
    """An attribute value that is valid JSON. Zarr metadata is written by
    Python writers and carries bare `NaN` tokens, which `json.load` accepts
    as floats: an NWB TimeSeries writes `"resolution": NaN` whenever the
    resolution is unspecified, 278 times over in one real DANDI store. A
    non-finite number cannot be written back out, so it becomes None, which
    is also what the Go port sees after it rewrites the token to null."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, list):
        return [_safe(v) for v in value]
    return value


def _attrs(node: Optional[dict]) -> dict:
    if not node:
        return {}
    a = node.get("attributes")
    return a if isinstance(a, dict) else {}


def _shape(node: Optional[dict]) -> Optional[List[int]]:
    if not node:
        return None
    shape = node.get("shape")
    if not isinstance(shape, list):
        return None
    out = []
    for d in shape:
        if not isinstance(d, int) or isinstance(d, bool) or d < 0:
            return None
        out.append(d)
    return out


def _rows(metadata: Dict[str, dict], path: str) -> Optional[int]:
    """Row count of a DynamicTable: the length of its `id` column."""
    shape = _shape(_node(metadata, path + "/id"))
    return shape[0] if shape else None


def _columns(node: Optional[dict]) -> Optional[List[str]]:
    """A DynamicTable names its columns in a `colnames` attribute, which HDF5
    does not have to be asked for either."""
    names = _attrs(node).get("colnames")
    if not isinstance(names, list):
        return None
    return [str(n) for n in names][:NWB_LIST_CAP]


def _is_string_array(node: Optional[dict]) -> bool:
    return bool(node) and node.get("node_type") == "array" and node.get("data_type") == "string"


def _series(metadata: Dict[str, dict], base: str, names: List[str], cap: int) -> List[Dict[str, object]]:
    """Describe the neurodata objects in a container: name, type, and for a
    TimeSeries its data shape, dtype, unit and rate."""
    out = []
    for name in names[:cap]:
        path = base + "/" + name if base else name
        attrs = _attrs(_node(metadata, path))
        entry: Dict[str, object] = {"name": name, "neurodataType": _safe(attrs.get("neurodata_type"))}
        data = _node(metadata, path + "/data")
        if data is not None:
            entry["shape"] = _shape(data)
            dtype = data.get("data_type")
            entry["dtype"] = str(dtype) if isinstance(dtype, str) else None
            unit = _safe(_attrs(data).get("unit"))
            if unit is not None:
                entry["unit"] = unit
        start = _node(metadata, path + "/starting_time")
        if start is not None:
            rate = _safe(_attrs(start).get("rate"))
            if isinstance(rate, (int, float)) and not isinstance(rate, bool):
                entry["rate"] = rate
        elif _node(metadata, path + "/timestamps") is not None:
            entry["timestamps"] = True
        out.append(entry)
    return out


def _specifications(metadata: Dict[str, dict]) -> Tuple[Dict[str, object], List[str]]:
    """`specifications/<namespace>/<version>/...` gives both the namespace and
    the version it was cached at, which the HDF5 layout also carries but this
    one states in the path."""
    specs: Dict[str, object] = {}
    for namespace in _children(metadata, "specifications")[:NWB_LIST_CAP]:
        versions = _children(metadata, "specifications/" + namespace)
        specs[namespace] = versions[:8]
    extensions = sorted(n for n in specs if n not in _CORE_NAMESPACES)
    return specs, extensions


def _unread_strings(metadata: Dict[str, dict]) -> List[str]:
    """Paths of the string datasets whose values are chunk data this never
    opens: the session's own fields, and the subject's."""
    out = []
    for path, node in metadata.items():
        if not _is_string_array(node):
            continue
        depth = path.count("/")
        if depth == 0 or path.startswith("general/"):
            out.append(path)
    return sorted(out)[:_PATH_CAP]


def decode(attributes: dict, metadata: Optional[Dict[str, dict]], zarr_format: int) -> DecodedHeader:
    """An NWB-on-Zarr record from the root group's attributes and, when the
    store was written with consolidated metadata, the flat path map that comes
    with them."""
    fields: Dict[str, object] = {
        "nwbVersion": attributes.get("nwb_version"),
        "namespace": attributes.get("namespace"),
        "objectId": attributes.get("object_id"),
        "container": "zarr",
        "zarrFormat": zarr_format,
        "consolidatedMetadata": metadata is not None,
    }
    if metadata is None:
        # Without consolidated metadata the root file names no children, and
        # reading them would mean one request per node. The record says what
        # it is and stops.
        fields["nodeCount"] = None
        return DecodedHeader(
            format="nwb-zarr",
            summary=_summary(fields, 0, 0, None, None, []),
            fields=fields,
        )

    fields["nodeCount"] = len(metadata)
    fields["stringValuesUnread"] = _unread_strings(metadata)

    devices = []
    for name in _children(metadata, "general/devices")[:NWB_LIST_CAP]:
        a = _attrs(_node(metadata, "general/devices/" + name))
        devices.append({"name": name, "description": _safe(a.get("description")),
                        "manufacturer": _safe(a.get("manufacturer"))})
    fields["devices"] = devices

    ephys = _children(metadata, "general/extracellular_ephys")
    if ephys:
        fields["electrodeGroups"] = [n for n in ephys if n != "electrodes"][:NWB_LIST_CAP]
        fields["electrodeCount"] = _rows(metadata, "general/extracellular_ephys/electrodes")
    else:
        fields["electrodeGroups"] = None
        fields["electrodeCount"] = None
    planes = _children(metadata, "general/optophysiology")[:NWB_LIST_CAP]
    fields["imagingPlanes"] = planes or None
    subject = _node(metadata, "general/subject")
    fields["subjectPresent"] = subject is not None

    acquisition = _children(metadata, "acquisition")
    fields["acquisitionNames"] = acquisition[:NWB_LIST_CAP * 2]
    fields["acquisitionCount"] = len(acquisition)
    fields["acquisition"] = _series(metadata, "acquisition", acquisition, NWB_SERIES_DETAIL)

    modules = _children(metadata, "processing")[:NWB_LIST_CAP]
    fields["processingModules"] = [
        {
            "name": name,
            "description": _safe(_attrs(_node(metadata, "processing/" + name)).get("description")),
            "interfaceNames": _children(metadata, "processing/" + name)[:NWB_LIST_CAP],
            "interfaces": _series(metadata, "processing/" + name,
                                  _children(metadata, "processing/" + name), 8),
        }
        for name in modules
    ]

    units = _node(metadata, "units")
    unit_count = _rows(metadata, "units") if units is not None else None
    fields["unitCount"] = unit_count
    fields["unitColumns"] = _columns(units) if units is not None else None

    fields["intervals"] = [
        {"name": name, "rows": _rows(metadata, "intervals/" + name),
         "columns": _columns(_node(metadata, "intervals/" + name))}
        for name in _children(metadata, "intervals")[:NWB_LIST_CAP]
    ]

    stimulus = _children(metadata, "stimulus/presentation")
    fields["stimulusPresentationNames"] = stimulus[:NWB_LIST_CAP * 2]
    fields["stimulusPresentation"] = _series(metadata, "stimulus/presentation", stimulus, 8)
    fields["stimulusTemplates"] = _children(metadata, "stimulus/templates")[:NWB_LIST_CAP]
    fields["analysis"] = _children(metadata, "analysis")[:NWB_LIST_CAP]

    specs, extensions = _specifications(metadata)
    fields["specifications"] = specs
    fields["extensions"] = extensions

    return DecodedHeader(
        format="nwb-zarr",
        summary=_summary(fields, len(acquisition), len(modules), unit_count,
                         fields["electrodeCount"], extensions),
        fields=fields,
    )


def _summary(fields: Dict[str, object], acquisition: int, modules: int,
             units: Optional[int], electrodes: Optional[int], extensions: List[str]) -> str:
    version = fields.get("nwbVersion")
    head = "NWB %s on Zarr v%s" % (version, fields.get("zarrFormat")) if version \
        else "NWB on Zarr v%s" % fields.get("zarrFormat")
    parts = [head]
    if acquisition:
        parts.append("%d acquisition" % acquisition)
    if modules:
        parts.append("%d processing module%s" % (modules, "" if modules == 1 else "s"))
    if units:
        parts.append("%d units" % units)
    if electrodes:
        parts.append("%d electrodes" % electrodes)
    if extensions:
        parts.append("extensions " + ", ".join(extensions))
    if not fields.get("consolidatedMetadata"):
        parts.append("no consolidated metadata")
    return ", ".join(parts)
