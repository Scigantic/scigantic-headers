"""Chromatography, electrophoresis and surface plasmon resonance exports.

Liquid chromatography: Agilent OpenLab CDS result sets (.rslt, a zip holding
the ACAML sequence document), Cytiva Unicorn result zips, and the JSON that
Benchling's Waters Empower and Thermo Fisher Scientific Chromeleon adapters
write. Electrophoresis: the Agilent TapeStation Analysis XML export. Binding
affinity: Cytiva Biacore T200 Control result files (.blr) and T200 Evaluation
files (.bme), both OLE2 compound files, and Biacore Insight Evaluation .xlsx
exports.

The JSON and XML exports are read from their leading bytes with the same
truncation-tolerant parsing as scigantic-edge, so both libraries report the
same fields for the same bytes. The zip, OLE2 and .xlsx containers are read
through the bounded Source (zip directory, one or two members, a handful of
OLE2 sectors) and are Python only.
"""

from __future__ import annotations

import datetime
import io
import re
import struct
import xml.etree.ElementTree as ET
import zipfile
from typing import Dict, List, Optional, Tuple

from . import (
    head,
    head_text,
    header,
    instrument,
    open_zip,
    register_instrument_sniffer,
    xlsx_rows,
    xlsx_sheet_names,
    zip_member,
)

OPENLAB = "Agilent OpenLab CDS"
UNICORN = "Cytiva Unicorn"
EMPOWER = "Benchling Waters Empower Adapter"
CHROMELEON = "Benchling Thermo Fisher Scientific Chromeleon"
TAPESTATION = "Agilent TapeStation Analysis"
T200_CONTROL = "Cytiva Biacore T200 Control"
T200_EVALUATION = "Cytiva Biacore T200 Evaluation"
INSIGHT = "Cytiva Biacore Insight"

instrument(OPENLAB, "Liquid Chromatography", "Agilent", "agilent_openlab_cds",
           ["Absorbance", "Fluorescence"])
instrument(UNICORN, "Liquid Chromatography", "Cytiva", "cytiva_unicorn", ["Absorbance"])
instrument(EMPOWER, "Liquid Chromatography", "Waters", "benchling_empower", ["Absorbance"])
instrument(CHROMELEON, "Liquid Chromatography", "Thermo Fisher Scientific",
           "benchling_chromeleon", ["Absorbance", "Fluorescence", "Conductivity"])
instrument(TAPESTATION, "Electrophoresis", "Agilent", "agilent_tapestation_analysis",
           ["Fluorescence"])
instrument(T200_CONTROL, "Binding Affinity Analyzer", "Cytiva", "cytiva_biacore_t200_control",
           ["Surface Plasmon Resonance"])
instrument(T200_EVALUATION, "Binding Affinity Analyzer", "Cytiva",
           "cytiva_biacore_t200_evaluation", ["Surface Plasmon Resonance"])
instrument(INSIGHT, "Binding Affinity Analyzer", "Cytiva", "cytiva_biacore_insight",
           ["Surface Plasmon Resonance"])

# Leading bytes the JSON and XML sniffers parse. Equal to scigantic-edge's
# InstrumentHeadBytes so both libraries see the same bytes.
JSON_HEAD = 256 * 1024
# Decompressed bytes read from one zip member (the ACAML document, a Unicorn
# chromatogram) and the largest nested zip opened in memory.
MEMBER_LIMIT = 8 * 1024 * 1024
NESTED_ZIP_LIMIT = 16 * 1024 * 1024


# Whitespace trimmed from values. Spelled out (not str.strip()) so the Go port
# trims exactly the same characters.
_STRIP = " \t\r\n\x0b\x0c"


def _plural(n: int, word: str) -> str:
    return "%d %s%s" % (n, word, "" if n == 1 else "s")


def _num(v: float):
    """220.0 -> 220, 254.5 -> 254.5, so integral values print like scigantic-edge."""
    return int(v) if float(v).is_integer() and abs(v) < 1e15 else v


def _unique(items) -> List:
    out: List = []
    for x in items:
        if x not in out:
            out.append(x)
    return out


# ── XML text helpers (shared with scigantic-edge) ────────────────────────────

_ENTITY = re.compile(r"&(#[0-9]+|#x[0-9A-Fa-f]+|lt|gt|amp|quot|apos);")
_NAMED = {"lt": "<", "gt": ">", "amp": "&", "quot": '"', "apos": "'"}


def _unescape(s: str) -> str:
    def rep(m):
        e = m.group(1)
        if e in _NAMED:
            return _NAMED[e]
        try:
            cp = int(e[2:], 16) if e[1] in "xX" else int(e[1:])
            return chr(cp) if 0 < cp <= 0x10FFFF and not 0xD800 <= cp <= 0xDFFF else m.group(0)
        except ValueError:
            return m.group(0)

    return _ENTITY.sub(rep, s)


def _tag_text(s: str, tag: str) -> Optional[str]:
    """Unescaped, stripped text of the first <tag>...</tag> in `s`, or None
    when absent, empty, or holding markup."""
    m = re.search(r"<%s>([^<]*)</%s>" % (tag, tag), s)
    if not m:
        return None
    v = _unescape(m.group(1)).strip(_STRIP)
    return v or None


# ── Truncation-tolerant JSON (shared with scigantic-edge) ────────────────────
#
# A JSON export is parsed from its leading bytes. Containers cut by the end of
# the read come back marked incomplete; a string, number or literal cut by the
# end is dropped. Counts are only reported from complete containers.


class PDict(dict):
    complete = False


class PList(list):
    complete = False


_NUMBER = re.compile(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?")
_WS_RE = re.compile(r"[ \t\r\n]*")
_STR_RUN = re.compile(r'[^"\\\x00-\x1f]*')
_HEX4 = re.compile(r"[0-9A-Fa-f]{4}")
_MAX_DEPTH = 64
_ESC = {'"': '"', "\\": "\\", "/": "/", "b": "\b", "f": "\f", "n": "\n", "r": "\r", "t": "\t"}


class _Stop(Exception):
    """The input ended (or stopped being JSON). `partial` is the innermost
    container read so far, already holding its partial children."""

    def __init__(self):
        super().__init__()
        self.partial = None


class _PartialJSON:
    def __init__(self, s: str):
        self.s = s
        self.i = 0
        self.n = len(s)

    def ws(self):
        self.i = _WS_RE.match(self.s, self.i).end()

    def value(self, depth: int):
        self.ws()
        if self.i >= self.n or depth > _MAX_DEPTH:
            raise _Stop()
        c = self.s[self.i]
        if c == "{":
            return self.obj(depth)
        if c == "[":
            return self.arr(depth)
        if c == '"':
            return self.string()
        for lit, val in (("true", True), ("false", False), ("null", None)):
            if self.s.startswith(lit, self.i):
                self.i += len(lit)
                if self.i >= self.n:
                    raise _Stop()
                return val
        m = _NUMBER.match(self.s, self.i)
        # A number touching the end of the read may have been cut.
        if not m or m.end() >= self.n:
            raise _Stop()
        self.i = m.end()
        t = m.group(0)
        if "." in t or "e" in t or "E" in t or len(t) > 15:
            f = float(t)
            if f != f or f in (float("inf"), float("-inf")):
                raise _Stop()
            return f
        return int(t)

    def string(self) -> str:
        s, n = self.s, self.n
        i = self.i + 1
        out = []
        while True:
            j = _STR_RUN.match(s, i).end()
            out.append(s[i:j])
            if j >= n:
                raise _Stop()
            c = s[j]
            if c == '"':
                self.i = j + 1
                return "".join(out)
            if c != "\\" or j + 1 >= n:
                raise _Stop()  # a raw control character, or a cut escape
            e = s[j + 1]
            if e in _ESC:
                out.append(_ESC[e])
                i = j + 2
                continue
            if e != "u" or not _HEX4.fullmatch(s, j + 2, j + 6):
                raise _Stop()
            cp = int(s[j + 2:j + 6], 16)
            i = j + 6
            if 0xD800 <= cp < 0xDC00 and s.startswith("\\u", i) and _HEX4.fullmatch(s, i + 2, i + 6):
                lo = int(s[i + 2:i + 6], 16)
                if 0xDC00 <= lo < 0xE000:
                    cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00)
                    i += 6
            out.append("\ufffd" if 0xD800 <= cp < 0xE000 else chr(cp))

    def obj(self, depth: int) -> PDict:
        d = PDict()
        self.i += 1
        try:
            self.ws()
            if self.i < self.n and self.s[self.i] == "}":
                self.i += 1
                d.complete = True
                return d
            while True:
                self.ws()
                if self.i >= self.n or self.s[self.i] != '"':
                    raise _Stop()
                k = self.string()
                self.ws()
                if self.i >= self.n or self.s[self.i] != ":":
                    raise _Stop()
                self.i += 1
                try:
                    d[k] = self.value(depth + 1)
                except _Stop as e:
                    if e.partial is not None:
                        d[k] = e.partial
                    raise
                self.ws()
                if self.i >= self.n:
                    raise _Stop()
                c = self.s[self.i]
                self.i += 1
                if c == "}":
                    d.complete = True
                    return d
                if c != ",":
                    raise _Stop()
        except _Stop as e:
            e.partial = d
            raise

    def arr(self, depth: int) -> PList:
        a = PList()
        self.i += 1
        try:
            self.ws()
            if self.i < self.n and self.s[self.i] == "]":
                self.i += 1
                a.complete = True
                return a
            while True:
                try:
                    a.append(self.value(depth + 1))
                except _Stop as e:
                    if e.partial is not None:
                        a.append(e.partial)
                    raise
                self.ws()
                if self.i >= self.n:
                    raise _Stop()
                c = self.s[self.i]
                self.i += 1
                if c == "]":
                    a.complete = True
                    return a
                if c != ",":
                    raise _Stop()
        except _Stop as e:
            e.partial = a
            raise


def partial_json(text: Optional[str]):
    """The JSON object or array at the start of `text`, tolerating truncation:
    containers cut by the end of the text come back with `complete` False and
    the members read so far. None when the text does not start with { or [."""
    if not text:
        return None
    p = _PartialJSON(text)
    p.ws()
    if p.i >= p.n or p.s[p.i] not in "{[":
        return None
    try:
        return p.value(0)
    except _Stop as e:
        return e.partial
    except (RecursionError, ValueError, OverflowError):
        return None


def _json_head(src):
    key = "_chrom_json"
    if not hasattr(src, key):
        t = head_text(src, JSON_HEAD)
        setattr(src, key, partial_json(t.lstrip("\ufeff") if t else t))
    return getattr(src, key)


def _s(v) -> Optional[str]:
    """A non-empty string value, stripped, or None."""
    if isinstance(v, str):
        v = v.strip(_STRIP)
        return v or None
    return None


def _complete_len(v) -> Optional[int]:
    return len(v) if isinstance(v, PList) and v.complete else None


_NM = re.compile(r"([0-9]+(?:\.[0-9]+)?)[ \t]*nm(?:[^A-Za-z0-9_]|$)")


# ── Benchling Waters Empower adapter (JSON) ──────────────────────────────────

_EMPOWER_KEYS = ("SampleSetId", "SampleSetName", "SystemName", "SampleSetStartDate",
                 "SampleSetMethod", "SampleSetType", "Node")


def sniff_empower(src) -> Optional[object]:
    doc = _json_head(src)
    if not isinstance(doc, dict):
        return None
    values = doc.get("values")
    if not isinstance(values, dict):
        return None
    fields = values.get("fields")
    if not isinstance(fields, dict) or sum(1 for k in _EMPOWER_KEYS if k in fields) < 3:
        return None
    out: Dict[str, object] = {
        "systemName": _s(fields.get("SystemName")),
        "sampleSetName": _s(fields.get("SampleSetName")),
        "sampleSetMethod": _s(fields.get("SampleSetMethod")),
        "runDate": _s(fields.get("SampleSetStartDate")),
    }
    channels = values.get("channels")
    chan_fields = [c.get("fields") for c in channels if isinstance(c, dict)] if isinstance(channels, list) else []
    chan_fields = [f for f in chan_fields if isinstance(f, dict)]
    for f in chan_fields:
        v = _s(f.get("AcqSWVersion"))
        if v:
            out["softwareVersion"] = v
            break
    out["channelCount"] = _complete_len(channels)
    if out["channelCount"] is not None:
        desc = _unique(d for d in (_s(f.get("ChannelDescription")) for f in chan_fields) if d)
        if desc:
            out["channels"] = desc
            wl = _unique(_num(float(m)) for d in desc for m in _NM.findall(d))
            if wl:
                out["wavelengthsNm"] = wl
    injections = values.get("injections")
    out["injectionCount"] = _complete_len(injections)
    if out["injectionCount"] is not None:
        methods = _unique(m for m in (
            _s(i["fields"].get("InstrumentMethodName")) for i in injections
            if isinstance(i, dict) and isinstance(i.get("fields"), dict)) if m)
        if methods:
            out["instrumentMethods"] = methods
    out["exportFormat"] = "json"
    parts = [EMPOWER + " export"]
    if out.get("sampleSetName"):
        parts.append("sample set %s" % out["sampleSetName"])
    if out.get("injectionCount") is not None:
        parts.append(_plural(out["injectionCount"], "injection"))
    return header(EMPOWER, ", ".join(parts), out)


# ── Benchling Thermo Fisher Scientific Chromeleon adapter (JSON) ─────────────

_MODULES = (("pump", "pump model number"), ("uv", "uv model number"),
            ("sampler", "sampler model number"))


def sniff_chromeleon(src) -> Optional[object]:
    doc = _json_head(src)
    if not isinstance(doc, dict) or not isinstance(doc.get("injections"), list):
        return None
    seq = doc.get("sequence")
    dev = doc.get("device information")
    has_seq = isinstance(seq, dict) and any(
        k in seq for k in ("sequence name", "sequence creation time", "sequence directory"))
    has_dev = isinstance(dev, dict) and any(k.endswith(" model number") for k in dev)
    if not (has_seq or has_dev):
        return None
    injections = doc["injections"]
    first = injections[0] if injections and isinstance(injections[0], dict) else None
    if first is not None and not any(
            k in first for k in ("injection identifier", "injection number", "signals")):
        return None
    if first is None and not has_seq:
        return None
    out: Dict[str, object] = {}
    if isinstance(seq, dict):
        out["sequenceName"] = _s(seq.get("sequence name"))
        out["runDate"] = _s(seq.get("sequence creation time"))
    if isinstance(dev, dict):
        mods = [{"type": t, "model": _s(dev.get(k))} for t, k in _MODULES if _s(dev.get(k))]
        if mods:
            out["modules"] = mods
    out["injectionCount"] = _complete_len(injections)
    if first is not None and isinstance(first.get("signals"), PList) and first["signals"].complete:
        sigs = [x for x in first["signals"] if isinstance(x, dict)]
        names = _unique(n for n in (_s(x.get("signal name")) for x in sigs) if n)
        if names:
            out["signals"] = names
        wl = _unique(_num(float(w)) for w in (x.get("wavelength setting") for x in sigs)
                     if isinstance(w, (int, float)) and not isinstance(w, bool) and w > 0)
        if wl:
            out["wavelengthsNm"] = wl
    out["exportFormat"] = "json"
    parts = [CHROMELEON + " export"]
    if out.get("sequenceName"):
        parts.append("sequence %s" % out["sequenceName"])
    if out.get("injectionCount") is not None:
        parts.append(_plural(out["injectionCount"], "injection"))
    return header(CHROMELEON, ", ".join(parts), out)


# ── Agilent TapeStation Analysis (XML) ───────────────────────────────────────

_XML_ROOT = re.compile(
    r"\A[ \t\r\n]*(?:<\?[^>]*\?>[ \t\r\n]*|<!--.*?-->[ \t\r\n]*)*<([A-Za-z_][A-Za-z0-9_.-]*)[ \t\r\n>/]",
    re.S)


def _root_tag(text: str) -> Optional[str]:
    m = _XML_ROOT.match(text)
    return m.group(1) if m else None


def sniff_tapestation(src) -> Optional[object]:
    t = head_text(src, JSON_HEAD)
    if not t or _root_tag(t) != "File":
        return None
    a = t.find("<FileInformation>")
    b = t.find("</FileInformation>", a)
    if a < 0 or b < 0 or "<ScreenTapes>" not in t:
        return None
    info = t[a:b]
    assay = _tag_text(info, "Assay")
    if assay is None:
        return None
    out: Dict[str, object] = {"assay": assay, "runDate": _tag_text(info, "RunEndDate")}
    for tag, key in (("RINeVersion", "rineVersion"), ("DINVersion", "dinVersion")):
        out[key] = _tag_text(info, tag)
    e1 = t.find("<Environment>")
    e2 = t.find("</Environment>", e1)
    if e1 >= 0 and e2 >= 0:
        env = t[e1:e2]
        out["softwareVersion"] = _tag_text(env, "AnalysisVersion")
        out["instrumentType"] = _tag_text(env, "InstrumentType")
        out["serialNumber"] = _tag_text(env, "InstrumentSerialNumber")
    s1 = t.find("<ScreenTapes>")
    s2 = t.find("</ScreenTapes>", s1)
    if s2 >= 0:
        out["screenTapeCount"] = len(re.findall(r"<ScreenTape>", t[s1:s2]))
    m1 = t.find("<Samples>")
    m2 = t.find("</Samples>", m1) if m1 >= 0 else -1
    if m2 >= 0:
        out["sampleCount"] = len(re.findall(r"<Sample>", t[m1:m2]))
    out["exportFormat"] = "xml"
    parts = [TAPESTATION + " export", "%s assay" % assay]
    if out.get("sampleCount") is not None:
        parts.append(_plural(out["sampleCount"], "sample"))
    return header(TAPESTATION, ", ".join(parts), out)


# ── Zip helpers ──────────────────────────────────────────────────────────────


def _visible(name: str) -> bool:
    base = name.rsplit("/", 1)[-1]
    return bool(base) and not base.startswith(".") and not name.startswith("__MACOSX/")


def _zip_of(src) -> Optional[zipfile.ZipFile]:
    zf = getattr(src, "_instrument_zip", None)
    if zf is None:
        zf = open_zip(src)
        src._instrument_zip = zf  # type: ignore[attr-defined]
    return zf


def _names(zf: zipfile.ZipFile) -> List[str]:
    try:
        return [n for n in zf.namelist() if _visible(n)][:100000]
    except Exception:
        return []


def _inner_zip(zf: zipfile.ZipFile, name: str, limit: int) -> Optional[zipfile.ZipFile]:
    try:
        if zf.getinfo(name).file_size > limit:
            return None
    except KeyError:
        return None
    data = zip_member(zf, name, limit)
    if not data or not data.startswith(b"PK\x03\x04"):
        return None
    try:
        return zipfile.ZipFile(io.BytesIO(data))
    except Exception:
        return None


# ── Agilent OpenLab CDS result set (.rslt zip) ───────────────────────────────


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _child(e, name: str):
    if e is None:
        return None
    for c in e:
        if _local(c.tag) == name:
            return c
    return None


def _children(e, name: str) -> list:
    return [c for c in e if _local(c.tag) == name] if e is not None else []


def _path(e, *names):
    for n in names:
        e = _child(e, n)
    return e


def _text(e) -> Optional[str]:
    if e is None or e.text is None:
        return None
    return e.text.strip() or None


def _parse_xml_partial(data: bytes, ends: Tuple[str, ...]):
    """Parse XML that may be cut short. Returns (root, set of local names of
    `ends` elements whose end tag was read)."""
    if b"<!DOCTYPE" in data[:4096] or b"<!ENTITY" in data[:8192]:
        return None, set()
    parser = ET.XMLPullParser(events=("start", "end"))
    root = None
    closed = set()
    try:
        for i in range(0, len(data), 1 << 20):
            parser.feed(data[i:i + (1 << 20)])
            for ev, el in parser.read_events():
                if ev == "start" and root is None:
                    root = el
                elif ev == "end" and _local(el.tag) in ends:
                    closed.add(_local(el.tag))
        parser.close()
    except ET.ParseError:
        pass
    return root, closed


_SIG = re.compile(r"Sig=(\d+(?:\.\d+)?)")


def _acaml_member(zf: zipfile.ZipFile, depth: int = 0):
    names = _names(zf)
    for n in names:
        if n.lower().endswith(".acaml"):
            return zf, n, depth
    if depth >= 2:
        return None
    for n in names:
        if n.lower().endswith((".rslt", ".zip")):
            inner = _inner_zip(zf, n, NESTED_ZIP_LIMIT)
            if inner is not None:
                return _acaml_member(inner, depth + 1)
    return None


def sniff_openlab(src) -> Optional[object]:
    zf = _zip_of(src)
    if zf is None:
        return None
    found = _acaml_member(zf)
    if found is None:
        return None
    zf, name, depth = found
    data = zip_member(zf, name, MEMBER_LIMIT)
    if not data:
        return None
    if data.startswith(b"\xef\xbb\xbf"):
        data = data[3:]
    root, closed = _parse_xml_partial(data, ("Injections", "Samples"))
    if root is None or _local(root.tag) != "ACAML" or not root.tag.startswith("{urn:schemas-agilent-com:acaml"):
        return None
    content = _path(root, "Doc", "Content")
    instr = _path(content, "Resources", "Instrument")
    out: Dict[str, object] = {
        "instrumentName": _text(_child(instr, "Name")),
        "technique": _text(_child(instr, "Technique")),
    }
    mods = []
    for m in _children(instr, "Module"):
        d = {
            "type": _text(_child(m, "Type")),
            "name": _text(_child(m, "Name")),
            "model": _text(_child(m, "PartNo")),
            "serialNumber": _text(_child(m, "SerialNo")),
            "firmwareVersion": _text(_child(m, "FirmwareRevision")),
        }
        d = {k: v for k, v in d.items() if v is not None}
        if d:
            mods.append(d)
    if mods:
        out["modules"] = mods
    seq = _path(content, "SampleContexts", "Setup", "IdentParam", "Name")
    out["sequenceName"] = _text(seq)
    methods = _unique(n for n in (_text(_child(m, "Name")) for m in _children(content, "Method")) if n)
    if methods:
        out["methodNames"] = methods
    samples = _child(content, "Samples")
    if samples is not None and "Samples" in closed:
        out["sampleCount"] = len(_children(samples, "Setup"))
    inj_el = _child(content, "Injections")
    injections = _children(inj_el, "MeasData")
    if inj_el is not None and "Injections" in closed:
        out["injectionCount"] = len(injections)
    if injections:
        first = injections[0]
        app = _path(first, "AcquisitionApplication", "AgilentApp")
        out["softwareVersion"] = _text(_child(app, "Version"))
        out["runDate"] = _text(_path(first, "Info", "CreatedDate"))
        sigs = [s for s in _children(first, "Signal") if _text(_child(s, "Type")) == "Chromatogram"]
        descs = _unique(d for d in (_text(_child(s, "Description")) for s in sigs) if d)
        if descs:
            out["signals"] = descs
        dets = _unique(d for d in (_text(_child(s, "DetectorName")) for s in sigs) if d)
        if dets:
            out["detectors"] = dets
        wl = _unique(_num(float(m.group(1))) for d in descs for m in [_SIG.search(d)] if m)
        if wl:
            out["wavelengthsNm"] = wl
    out["exportFormat"] = "rslt result set" if depth == 0 else "zip holding an rslt result set"
    parts = [OPENLAB + " result set"]
    if out.get("instrumentName"):
        parts.append(str(out["instrumentName"]))
    if out.get("injectionCount") is not None:
        parts.append(_plural(out["injectionCount"], "injection"))
    return header(OPENLAB, ", ".join(parts), out)


# ── Cytiva Unicorn result zip ────────────────────────────────────────────────

_UNICORN_ROOT = re.compile(rb"<(Result|ResultData)\b[^>]*\bUNICORNVersion=\"([^\"]*)\"")
_CURVE = re.compile(r"<Curve\s+CurveDataType=\"([^\"]*)\"[^>]*>\s*<Name>([^<]*)</Name>")
_UV_NAME = re.compile(r"^UV(?:\s*\d+)?_(\d+(?:\.\d+)?)$")


def _xml_from(data: Optional[bytes]) -> str:
    """The XML text inside a Unicorn part: from the first '<' (the inner Xml
    parts carry a binary serialization prefix) to the last '>'."""
    if not data:
        return ""
    a = data.find(b"<")
    b = data.rfind(b">")
    if a < 0 or b < a:
        return ""
    return data[a:b + 1].decode("utf-8", "replace")


def _unicorn_part(zf, names, base: str) -> Optional[str]:
    for n in names:
        b = n.rsplit("/", 1)[-1]
        if b.endswith(base) or b.endswith(base + ".zip"):
            inner = _inner_zip(zf, n, 1 << 20)
            if inner is None:
                return None
            return _xml_from(zip_member(inner, "Xml", 1 << 20))
    return None


def sniff_unicorn(src) -> Optional[object]:
    zf = _zip_of(src)
    if zf is None:
        return None
    names = _names(zf)
    res = [n for n in names if n.endswith("Result.xml")]
    if len(res) != 1:
        return None
    rdata = zip_member(zf, res[0], 1 << 20)
    m = _UNICORN_ROOT.search(rdata[:4096] if rdata else b"")
    if not m:
        return None
    chroms = [n for n in names if n.endswith(".Xml")]
    if not chroms and not any(n.rsplit("/", 1)[-1].endswith(("SystemData", "SystemData.zip")) for n in names):
        return None
    r = _xml_from(rdata)
    # The run's own elements come before the nested documentation, search
    # criteria and chromatogram lists, which carry <Name> elements of their own.
    cut = re.search(r"<(ResultDocumentation|ResultRunInformation|ResultSearchCriterias|Chromatograms)\b", r)
    top = r[:cut.start()] if cut else r
    out: Dict[str, object] = {
        "softwareVersion": _unescape(m.group(2).decode("utf-8", "replace")).strip() or None,
        "systemName": _tag_text(top, "SystemName"),
        "systemType": _tag_text(top, "SystemTypeName"),
        "runName": _tag_text(top, "Name"),
        "runType": _tag_text(top, "RunType"),
        "resultState": _tag_text(top, "ResultState"),
    }
    sd = re.search(r"RunInformationType=\"RunStartDate\"><RunInformation>([^<]*)</RunInformation>", r)
    out["runDate"] = (sd.group(1).strip() if sd else None) or _tag_text(r, "Created")
    sysxml = _unicorn_part(zf, names, "SystemData")
    if sysxml:
        systems = re.findall(r"<System\b[^>]*>(.*?)</System>", sysxml, re.S)
        pick = None
        for s in systems:
            if out["systemName"] and _tag_text(s, "Name") == out["systemName"]:
                pick = s
                break
        if pick is None and len(systems) == 1:
            pick = systems[0]
        if pick is not None:
            d = re.search(r"<InstrumentConfiguration\b[^>]*\bDescription=\"([^\"]*)\"", pick)
            if d and d.group(1).strip():
                out["instrumentModel"] = _unescape(d.group(1)).strip()
    icd = _unicorn_part(zf, names, "InstrumentConfigurationData")
    if icd:
        out["firmwareVersion"] = _tag_text(icd, "FirmwareVersion")
    chrom_names = []
    for n in chroms[:64]:
        part = zip_member(zf, n, 4096)
        if part and b"<Chromatogram" in part:
            chrom_names.append(n)
    if len(chroms) <= 64:
        out["chromatogramCount"] = len(chrom_names)
    if chrom_names:
        first = next((n for n in chrom_names if n.endswith("Chrom.1.Xml")), chrom_names[0])
        text = _xml_from(zip_member(zf, first, MEMBER_LIMIT))
        curves = _CURVE.findall(text)
        names_ = _unique(_unescape(nm).strip() for _, nm in curves if nm.strip())
        if names_:
            out["curves"] = names_
        wl = _unique(_num(float(mm.group(1))) for nm in names_ for mm in [_UV_NAME.match(nm)]
                     if mm and float(mm.group(1)) > 0)
        if wl:
            out["wavelengthsNm"] = wl
    out["exportFormat"] = "result zip"
    parts = [UNICORN + " result"]
    if out.get("runName"):
        parts.append(str(out["runName"]))
    if out.get("instrumentModel"):
        parts.append(str(out["instrumentModel"]))
    return header(UNICORN, ", ".join(parts), out)


# ── OLE2 compound files (Biacore .blr / .bme) ────────────────────────────────

OLE_MAGIC = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
_END = 0xFFFFFFFA  # sector ids at or above this end a chain
_OLE_CHUNK = 4096


class _Short(Exception):
    """A read came back short: outside the file or past the Source budget."""


class _Ole:
    """A lazy OLE2 reader: directory entries, FAT and mini FAT sectors are read
    only when a lookup reaches them, through a small chunk cache, so finding a
    few named streams costs a few reads even in a large file."""

    def __init__(self, src):
        self.src = src
        self._cache: Dict[int, bytes] = {}
        h = self._read(0, 512)
        if h[:8] != OLE_MAGIC or h[28:30] != b"\xfe\xff":
            raise _Short
        shift = struct.unpack_from("<H", h, 30)[0]
        if shift not in (9, 12) or struct.unpack_from("<H", h, 32)[0] != 6:
            raise _Short
        self.ssz = 1 << shift
        self.nfat, self.dir0 = struct.unpack_from("<II", h, 44)
        self.cutoff, self.minifat0, _, self.difat_next, self.ndifat = struct.unpack_from("<IIIII", h, 56)
        if self.cutoff != 4096:
            raise _Short
        self.difat = list(struct.unpack_from("<109I", h, 76))
        self._difat_read = 0
        self._fat: Dict[int, tuple] = {}
        self._dir = [self.dir0]
        self._ents: Dict[int, dict] = {}
        self._root_chain: Optional[List[int]] = None
        self._minifat_chain: Optional[List[int]] = None

    def _read(self, off: int, n: int) -> bytes:
        out = b""
        c = off // _OLE_CHUNK
        while len(out) < n:
            if c not in self._cache:
                self._cache[c] = self.src.read_at(c * _OLE_CHUNK, _OLE_CHUNK)
            chunk = self._cache[c]
            lo = off - c * _OLE_CHUNK if not out else 0
            out += chunk[lo:lo + n - len(out)]
            if len(chunk) < _OLE_CHUNK:
                break
            c += 1
        if len(out) < n:
            raise _Short
        return out

    def sector(self, sid: int) -> bytes:
        if sid >= _END:
            raise _Short
        return self._read((sid + 1) * self.ssz, self.ssz)

    def fat(self, sid: int) -> int:
        per = self.ssz // 4
        i = sid // per
        if i >= self.nfat or i > 1 << 22:
            raise _Short
        while i >= len(self.difat):
            if self.difat_next >= _END or self._difat_read >= self.ndifat:
                raise _Short
            v = struct.unpack("<%dI" % per, self.sector(self.difat_next))
            self.difat.extend(v[:-1])
            self.difat_next = v[-1]
            self._difat_read += 1
        if i not in self._fat:
            self._fat[i] = struct.unpack("<%dI" % per, self.sector(self.difat[i]))
        return self._fat[i][sid % per]

    def _grow(self, chain: List[int], k: int):
        while len(chain) <= k:
            nxt = self.fat(chain[-1])
            if nxt >= _END or len(chain) > 1 << 22:
                raise _Short
            chain.append(nxt)

    def entry(self, i: int) -> dict:
        if i in self._ents:
            return self._ents[i]
        per = self.ssz // 128
        self._grow(self._dir, i // per)
        d = self.sector(self._dir[i // per])
        e = d[(i % per) * 128:(i % per) * 128 + 128]
        nl = struct.unpack_from("<H", e, 64)[0]
        name = e[:max(0, min(nl, 64) - 2)].decode("utf-16-le", "replace").lstrip("\x01\x02\x03\x04\x05")
        left, right, child = struct.unpack_from("<III", e, 68)
        size = struct.unpack_from("<Q", e, 120)[0]
        if self.ssz == 512:
            size &= 0xFFFFFFFF
        ent = {"name": name, "type": e[66], "left": left, "right": right, "child": child,
               "start": struct.unpack_from("<I", e, 116)[0], "size": size}
        self._ents[i] = ent
        return ent

    def children(self, i: int, limit: int = 4096) -> Optional[List[dict]]:
        """The entries directly under storage `i`, or None if the listing could
        not be completed within the budget."""
        try:
            out: List[dict] = []
            stack = [self.entry(i)["child"]]
            seen = set()
            while stack:
                j = stack.pop()
                if j >= _END:
                    continue
                if j in seen or len(out) >= limit:
                    return None
                seen.add(j)
                e = self.entry(j)
                out.append(e)
                stack.append(e["left"])
                stack.append(e["right"])
            return out
        except _Short:
            return None

    def find(self, i: int, name: str) -> Optional[dict]:
        """The child of storage `i` called `name`, walking the red-black tree
        by name order so only one path of entries is read."""
        try:
            j = self.entry(i)["child"]
            key = (len(name), name.upper())
            for _ in range(64):
                if j >= _END:
                    return None
                e = self.entry(j)
                k = (len(e["name"]), e["name"].upper())
                if k == key:
                    return e
                j = e["left"] if key < k else e["right"]
        except _Short:
            return None
        return None

    def read(self, e: dict, limit: int = 1 << 20) -> Optional[bytes]:
        """The stream's bytes (up to `limit`), or None if unreadable."""
        try:
            size = min(e["size"], limit)
            out = []
            got = 0
            if e["size"] < self.cutoff:
                root = self.entry(0)
                if self._root_chain is None:
                    self._root_chain = [root["start"]]
                    self._minifat_chain = [self.minifat0]
                per = self.ssz // 4
                m = e["start"]
                for _ in range((size + 63) // 64):
                    off = m * 64
                    self._grow(self._root_chain, off // self.ssz)
                    sec = self.sector(self._root_chain[off // self.ssz])
                    out.append(sec[off % self.ssz:off % self.ssz + 64])
                    got += 64
                    if got >= size:
                        break
                    self._grow(self._minifat_chain, m // per)
                    mf = self.sector(self._minifat_chain[m // per])
                    m = struct.unpack_from("<I", mf, (m % per) * 4)[0]
                    if m >= _END:
                        break
            else:
                chain = [e["start"]]
                for k in range((size + self.ssz - 1) // self.ssz):
                    self._grow(chain, k)
                    out.append(self.sector(chain[k]))
            data = b"".join(out)
            return data[:size] if len(data) >= size else None
        except _Short:
            return None


def _kv(data: Optional[bytes]) -> Dict[str, str]:
    """A Biacore key=value stream (Environment, Chip, DataManager)."""
    if not data:
        return {}
    text = data.rstrip(b"\x00").decode("utf-8", "replace")
    out: Dict[str, str] = {}
    for line in text.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out.setdefault(k.strip(), v.strip())
    return out


def _ole_date(v: Optional[str]) -> Optional[str]:
    """An OLE automation date (days since 1899-12-30) as ISO 8601, no zone."""
    try:
        days = float(v or "")
    except ValueError:
        return None
    if not 1.0 <= days <= 2958465.0:
        return None
    t = datetime.datetime(1899, 12, 30) + datetime.timedelta(seconds=int(days * 86400))
    return t.isoformat()


def _basename_win(p: Optional[str]) -> Optional[str]:
    if not p:
        return None
    b = re.split(r"[\\/]", p)[-1].strip()
    return b or None


def _ole_of(src) -> Optional[_Ole]:
    if head(src, 8) != OLE_MAGIC:
        return None
    o = getattr(src, "_chrom_ole", None)
    if o is None:
        try:
            o = _Ole(src)
        except (_Short, struct.error):
            o = False
        src._chrom_ole = o  # type: ignore[attr-defined]
    return o or None


def _top_stream(o: _Ole, name: str) -> Optional[dict]:
    """A top-level stream by name. Biacore writes its file-type streams as the
    first directory entries, so the first directory sector is checked before
    walking the root's tree, which in a large file spans the whole file."""
    try:
        for i in range(1, o.ssz // 128):
            e = o.entry(i)
            if e["name"] == name and e["type"] == 2:
                return e
    except _Short:
        pass
    e = o.find(0, name)
    return e if e is not None and e["type"] == 2 else None


def _biacore_file_type(o: _Ole) -> Optional[Tuple[str, Dict[str, str], Optional[str]]]:
    comp = _top_stream(o, "BIA compability info")
    if comp is None:
        return None
    kv = _kv(o.read(comp, 4096))
    app = _top_stream(o, "BIA application info")
    app_info = None
    if app is not None:
        raw = o.read(app, 4096)
        if raw:
            app_info = raw.split(b"\x00", 1)[0].decode("latin-1").strip() or None
    return kv.get("FileType", ""), kv, app_info


def _run_fields(env: Dict[str, str], chip: Dict[str, str]) -> Dict[str, object]:
    out: Dict[str, object] = {
        "instrumentModel": env.get("ProcessingUnit") or None,
        "instrumentId": env.get("InstrumentId") or None,
        "runType": env.get("RunType") or None,
        "methodFile": _basename_win(env.get("TemplateFile") or env.get("MethodFile")),
        "runDate": _ole_date(env.get("Timestamp")),
        "chipType": chip.get("Name") or None,
        "chipLotNumber": chip.get("LotNo") or None,
    }
    if (chip.get("NoFcs") or "").isdigit():
        out["flowCellCount"] = int(chip["NoFcs"])
    if out["methodFile"] == "(Untitled)":
        out["methodFile"] = None
    return out


def _cycle_count(o: _Ole, storage: int) -> Optional[int]:
    kids = o.children(storage)
    if kids is None:
        return None
    return sum(1 for e in kids if e["type"] == 1 and re.fullmatch(r"_Cycle \d+", e["name"]))


def _entry_index(o: _Ole, e: dict) -> Optional[int]:
    for i, x in o._ents.items():
        if x is e:
            return i
    return None


def sniff_t200_control(src) -> Optional[object]:
    o = _ole_of(src)
    if o is None:
        return None
    ft = _biacore_file_type(o)
    if ft is None or ft[0] != "Result File":
        return None
    env_e = o.find(0, "Environment")
    env = _kv(o.read(env_e)) if env_e is not None and env_e["type"] == 2 else {}
    if "T200" not in env.get("Application", "") or "Control" not in env.get("Application", ""):
        return None
    chip_e = o.find(0, "Chip")
    chip = _kv(o.read(chip_e)) if chip_e is not None and chip_e["type"] == 2 else {}
    out: Dict[str, object] = {"softwareVersion": env.get("Version") or None}
    out.update(_run_fields(env, chip))
    out["cycleCount"] = _cycle_count(o, 0)
    out["exportFormat"] = "blr result file"
    parts = [T200_CONTROL + " result"]
    if out.get("runType"):
        parts.append(str(out["runType"]))
    if out.get("cycleCount") is not None:
        parts.append(_plural(out["cycleCount"], "cycle"))
    return header(T200_CONTROL, ", ".join(parts), out)


def sniff_t200_evaluation(src) -> Optional[object]:
    o = _ole_of(src)
    if o is None:
        return None
    ft = _biacore_file_type(o)
    if ft is None or ft[0] != "T200 Evaluation File" or not (ft[2] or "").startswith("T200Eval"):
        return None
    out: Dict[str, object] = {"applicationInfo": ft[2]}
    env_e = o.find(0, "Environment")
    env = _kv(o.read(env_e)) if env_e is not None and env_e["type"] == 2 else {}
    if "Evaluation" in env.get("Application", ""):
        out["softwareVersion"] = env.get("Version") or None
    dm = o.find(0, "_DataManager 1")
    if dm is not None and dm["type"] == 1:
        dmi = _entry_index(o, dm)
        if dmi is not None:
            denv_e = o.find(dmi, "Environment")
            denv = _kv(o.read(denv_e)) if denv_e is not None and denv_e["type"] == 2 else {}
            chip_e = o.find(dmi, "Chip")
            chip = _kv(o.read(chip_e)) if chip_e is not None and chip_e["type"] == 2 else {}
            if denv.get("Application"):
                out["acquisitionSoftware"] = denv["Application"]
                out["acquisitionSoftwareVersion"] = denv.get("Version") or None
            out.update(_run_fields(denv, chip))
            out["cycleCount"] = _cycle_count(o, dmi)
    dmgr = o.find(0, "DataManager")
    if out.get("cycleCount") is None and dmgr is not None and dmgr["type"] == 2:
        total = _kv(o.read(dmgr, 1 << 16)).get("TotalNumberOfCycles", "")
        if total.isdigit():
            out["cycleCount"] = int(total)
    ev = o.find(0, "Evaluation")
    if ev is not None and ev["type"] == 1:
        evi = _entry_index(o, ev)
        kids = o.children(evi) if evi is not None else None
        if kids is not None:
            out["evaluationItemCount"] = sum(
                1 for e in kids if re.fullmatch(r"EvaluationItem\d+", e["name"]))
    out["exportFormat"] = "bme evaluation file"
    parts = [T200_EVALUATION + " file"]
    if out.get("cycleCount") is not None:
        parts.append(_plural(out["cycleCount"], "cycle"))
    if out.get("evaluationItemCount") is not None:
        parts.append(_plural(out["evaluationItemCount"], "evaluation item"))
    return header(T200_EVALUATION, ", ".join(parts), out)


# ── Cytiva Biacore Insight (.xlsx) ───────────────────────────────────────────


def _cell(row: List[str], i: int) -> str:
    return row[i].strip() if i < len(row) and row[i] is not None else ""


def _intlike(v: str) -> str:
    """Excel stores 1234 as 1234 or 1234.0; report 1234."""
    m = re.fullmatch(r"(-?\d+)\.0+", v)
    return m.group(1) if m else v


def sniff_insight(src) -> Optional[object]:
    names = xlsx_sheet_names(src)
    if not names or "Properties" not in names or "Report point table" not in names:
        return None
    rows = xlsx_rows(src, "Properties", 2000)
    if not rows:
        return None
    labels = [_cell(r, 0) for r in rows]
    if "Software" not in labels or "Evaluation" not in labels:
        return None
    out: Dict[str, object] = {}
    i = labels.index("Software")
    for r in rows[i + 1:i + 3]:
        if _cell(r, 0) == "Name":
            out["softwareName"] = _cell(r, 1) or None
        elif _cell(r, 0) == "Version":
            out["softwareVersion"] = _intlike(_cell(r, 1)) or None
    if not str(out.get("softwareName", "")).startswith("Biacore Insight"):
        return None
    i = labels.index("Evaluation")
    for r in rows[i + 1:i + 5]:
        if _cell(r, 0) in ("Software", "Properties") or re.fullmatch(r"Run \d+", _cell(r, 0)):
            break
        if _cell(r, 0) == "Name":
            out["evaluationName"] = _cell(r, 1) or None
            break
    runs = [k for k, lab in enumerate(labels) if re.fullmatch(r"Run \d+", lab)]
    out["runCount"] = len(runs) if len(rows) < 2000 else None
    if runs:
        start = runs[0]
        end = runs[1] if len(runs) > 1 else len(rows)
        info: Dict[str, str] = {}
        chip: Dict[str, str] = {}
        for r in rows[start:end]:
            for lc, vc in ((3, 4), (6, 7)):
                lab = _cell(r, lc)
                if lab and lab not in info:
                    info[lab] = _intlike(_cell(r, vc))
            if _cell(r, 0) in ("Chip type", "Chip id", "Lot number") and _cell(r, 0) not in chip:
                chip[_cell(r, 0)] = _intlike(_cell(r, 1))
        out["instrumentModel"] = info.get("Instrument type") or None
        out["serialNumber"] = info.get("Instrument id") or None
        out["methodName"] = info.get("Method") or None
        out["runningBuffer"] = info.get("Running buffer") or None
        cyc = info.get("Cycles", "")
        if cyc.isdigit():
            out["cycleCount"] = int(cyc)
        rate = info.get("Data collection rate", "")
        try:
            out["dataCollectionRateHz"] = _num(float(rate)) if rate else None
        except ValueError:
            pass
        out["chipType"] = chip.get("Chip type") or None
        out["chipLotNumber"] = chip.get("Lot number") or None
    others = [n for n in names if n not in ("Properties", "Report point table")]
    if others:
        out["evaluationSheets"] = others
    out["exportFormat"] = "xlsx"
    parts = [INSIGHT + " export"]
    if out.get("instrumentModel"):
        parts.append(str(out["instrumentModel"]))
    if out.get("cycleCount") is not None:
        parts.append(_plural(out["cycleCount"], "cycle"))
    return header(INSIGHT, ", ".join(parts), out)


register_instrument_sniffer("json", sniff_empower)
register_instrument_sniffer("json", sniff_chromeleon)
register_instrument_sniffer("xml", sniff_tapestation)
register_instrument_sniffer(["rslt", "zip"], sniff_openlab)
register_instrument_sniffer("zip", sniff_unicorn)
register_instrument_sniffer("blr", sniff_t200_control)
register_instrument_sniffer("bme", sniff_t200_evaluation)
register_instrument_sniffer(["xlsx", "xlsm"], sniff_insight)
