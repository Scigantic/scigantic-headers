"""Flow cytometry, multi-analyte profiling and liquid handler exports.

    FlowJo                       .wsp / .xml workspace (root <Workspace flowJoVersion=...>)
    BD Biosciences FACSDiva      .xml experiment export (root <bdfacs version=...>)
    Luminex xPONENT              .csv batch export ("Program","xPONENT",,<model>)
    Luminex INTELLIFLEX          .csv / .xlsx flat results table (INSTRUMENT TYPE,
                                 ..., WELL LOCATION, SAMPLE ID columns), and the
                                 xPONENT-style batch export an INTELLIFLEX writes
                                 ("Program","xPONENT",,"Intelliflex")
    Bio-Rad Bio-Plex Manager     .xml export (root <bpm:BioPlexManager ...>)
    Beckman Coulter Biomek       .log / .csv transfer log ("Method = ..." header)
    Beckman Echo Plate Reformat  .csv transfer report (Application Name ... Plate Reformat)
    Beckman Echo Cherry Pick     .csv transfer report (Application Name ... Cherry Pick)

FCS list-mode files are decoded by the FCS decoder, not here. The XML sniffers
match on the root element only, so any other XML (sidecars, EPU, OME) is left
to the ordinary .xml decoder.

Counts (samples, wells, transfers) are reported only when the part of the file
they describe was read to its end within the leading bytes, so a large export
never reports a count cut short by the read. Dates are converted to ISO 8601
only when written year first; month/day order is locale dependent, so other
forms are kept as written. Sample names are kept as the strings in the file.

The Go port (scigantic-edge header/instrument_flow_map_lh.go) reports the same
fields for every text and XML export; the INTELLIFLEX .xlsx export is read
only here.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

from . import (
    TEXT_HEAD_BYTES,
    head,
    head_text,
    header,
    instrument,
    register_instrument_sniffer,
    xlsx_rows,
    xlsx_sheet_names,
    zip_member,
)

FLOWJO = "FlowJo"
FACSDIVA = "BD Biosciences FACSDiva"
INTELLIFLEX = "Luminex INTELLIFLEX"
XPONENT = "Luminex xPONENT"
BIOPLEX = "Bio-Rad Bio-Plex Manager"
BIOMEK = "Beckman Coulter Biomek"
ECHO_REFORMAT = "Beckman Echo Plate Reformat"
ECHO_CHERRY = "Beckman Echo Cherry Pick"

instrument(FLOWJO, "Flow Cytometry", "BD Biosciences", "flowjo", ["Fluorescence"])
instrument(FACSDIVA, "Flow Cytometry", "BD Biosciences", "bd_biosciences_facsdiva", ["Fluorescence"])
instrument(INTELLIFLEX, "Multi Analyte Profiling", "Luminex", "luminex_intelliflex", ["Fluorescence"])
instrument(XPONENT, "Multi Analyte Profiling", "Luminex", "luminex_xponent", ["Fluorescence"])
instrument(BIOPLEX, "Multi Analyte Profiling", "Bio-Rad", "biorad_bioplex_manager", ["Fluorescence"])
instrument(BIOMEK, "Liquid Handler", "Beckman Coulter", "beckman_coulter_biomek", [])
instrument(ECHO_REFORMAT, "Liquid Handler", "Beckman Coulter", "beckman_echo_plate_reformat", [])
instrument(ECHO_CHERRY, "Liquid Handler", "Beckman Coulter", "beckman_echo_cherry_pick", [])

# Sample name lists longer than this are left out (the count is still given).
MAX_SAMPLE_NAMES = 384


# --- shared text helpers (mirrored one for one in the Go port) -------------

def _re(pattern):
    """Patterns are ASCII only (\\d, \\s), the same classes as Go's regexp."""
    return re.compile(pattern, re.ASCII)


def _text(src) -> Tuple[Optional[str], bool]:
    """The decoded leading text and whether it is the whole file."""
    t = head_text(src)
    if t is None:
        return None, False
    return t, len(head(src)) < TEXT_HEAD_BYTES


def _lines(text: str, complete: bool) -> List[str]:
    """Lines split on LF with CR stripped. When the read stopped before the end
    of the file, the last line may be cut and is dropped."""
    lines = text.split("\n")
    if not complete and lines:
        lines.pop()
    return [ln[:-1] if ln.endswith("\r") else ln for ln in lines]


def _cells(line: str) -> List[str]:
    """One CSV line split into cells: comma separated, a cell that starts with
    a double quote runs to the closing quote with "" as an escaped quote; text
    after the closing quote is kept as is."""
    out: List[str] = []
    cur: List[str] = []
    i, n = 0, len(line)
    quoted = False
    at_start = True
    while i < n:
        c = line[i]
        if quoted:
            if c == '"':
                if i + 1 < n and line[i + 1] == '"':
                    cur.append('"')
                    i += 2
                    continue
                quoted = False
            else:
                cur.append(c)
        elif c == ",":
            out.append("".join(cur))
            cur = []
            at_start = True
            i += 1
            continue
        elif c == '"' and at_start:
            quoted = True
        else:
            cur.append(c)
        at_start = False
        i += 1
    out.append("".join(cur))
    return out


_BLANK = _re(r'^[,"\s]*$')


def _blank(line: str) -> bool:
    return _BLANK.match(line) is not None


_ISO = _re(
    r"^(\d{4})-(\d{2})-(\d{2})[ T](\d{1,2}):(\d{2})(?::(\d{2}))?(?:\.\d+)?(?:\s*([AaPp][Mm]))?$"
)


def _date(s: Optional[str]) -> Optional[str]:
    """ISO 8601 for a year-first date-time, else the string as written."""
    if s is None:
        return None
    s = s.strip()
    if not s:
        return None
    m = _ISO.match(s)
    if not m:
        return s
    y, mo, d, h, mi, sec, ampm = m.groups()
    h_i, mo_i, d_i, mi_i = int(h), int(mo), int(d), int(mi)
    sec_i = int(sec) if sec else 0
    if ampm:
        if not 1 <= h_i <= 12:
            return s
        if ampm.lower() == "pm" and h_i != 12:
            h_i += 12
        elif ampm.lower() == "am" and h_i == 12:
            h_i = 0
    if not (1 <= mo_i <= 12 and 1 <= d_i <= 31 and h_i <= 23 and mi_i <= 59 and sec_i <= 59):
        return s
    return "%s-%s-%sT%02d:%02d:%02d" % (y, mo, d, h_i, mi_i, sec_i)


def _int(s: Optional[str]) -> Optional[int]:
    if s is None:
        return None
    s = s.strip()
    if not s.isascii() or not s.isdigit() or len(s) > 9:
        return None
    return int(s)


def _value(s: Optional[str]) -> Optional[str]:
    """A stripped non-empty value, or None for empty, "None" and "N/A"."""
    if s is None:
        return None
    s = s.strip()
    if not s or s in ("None", "N/A", "<None>"):
        return None
    return s


def _distinct(names: List[str]) -> List[str]:
    seen = set()
    out = []
    for n in names:
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out


def _summary(parts: List[Optional[str]]) -> str:
    return ", ".join(p for p in parts if p)


def _plural(n: Optional[int], word: str) -> Optional[str]:
    if n is None:
        return None
    return "%d %s%s" % (n, word, "" if n == 1 else "s")


# --- XML helpers ------------------------------------------------------------

_ATTR = _re(r"""([A-Za-z_:][A-Za-z0-9_:.\-]*)\s*=\s*(?:"([^"]*)"|'([^']*)')""")
_ENTITY = _re(r"&(amp|lt|gt|quot|apos|#[0-9]{1,7}|#x[0-9A-Fa-f]{1,6});")
_NAMED = {"amp": "&", "lt": "<", "gt": ">", "quot": '"', "apos": "'"}


def _unescape(s: str) -> str:
    def rep(m):
        e = m.group(1)
        if e in _NAMED:
            return _NAMED[e]
        cp = int(e[2:], 16) if e[1] == "x" else int(e[1:])
        if cp == 0 or cp > 0x10FFFF or 0xD800 <= cp <= 0xDFFF:
            return m.group(0)
        return chr(cp)

    return _ENTITY.sub(rep, s)


def _attrs(tag: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for m in _ATTR.finditer(tag):
        if m.group(1) not in out:
            v = m.group(2) if m.group(2) is not None else m.group(3)
            out[m.group(1)] = _unescape(v)
    return out


def _tag_end(text: str, start: int) -> int:
    """Index just past the '>' closing the tag that starts at `start`, quote
    aware, or -1."""
    q = ""
    i, n = start, min(len(text), start + 65536)
    while i < n:
        c = text[i]
        if q:
            if c == q:
                q = ""
        elif c in "\"'":
            q = c
        elif c == ">":
            return i + 1
        i += 1
    return -1


_NAME = _re(r"[A-Za-z_][A-Za-z0-9_:.\-]*")


def _root(text: str) -> Optional[Tuple[str, Dict[str, str], int]]:
    """(name, attributes, end offset) of the document's root start tag, after
    the XML declaration, comments, processing instructions and DOCTYPE."""
    i, n = 0, min(len(text), 65536)
    while i < n:
        c = text[i]
        if c in " \t\r\n\ufeff":
            i += 1
            continue
        if c != "<":
            return None
        if text.startswith("<?", i):
            j = text.find("?>", i)
            if j < 0:
                return None
            i = j + 2
        elif text.startswith("<!--", i):
            j = text.find("-->", i)
            if j < 0:
                return None
            i = j + 3
        elif text.startswith("<!", i):
            j = text.find(">", i)
            if j < 0:
                return None
            i = j + 1
        else:
            m = _NAME.match(text, i + 1)
            if not m:
                return None
            end = _tag_end(text, i)
            if end < 0:
                return None
            return m.group(0), _attrs(text[m.end():end]), end
    return None


def _first_tag(text: str, name: str, start: int = 0, stop: int = -1) -> Optional[Dict[str, str]]:
    """Attributes of the first <name ...> start tag at or after `start`."""
    m = _re("<" + re.escape(name) + r"[\s/>]").search(text, start, stop if stop >= 0 else len(text))
    if not m:
        return None
    end = _tag_end(text, m.start())
    if end < 0:
        return None
    return _attrs(text[m.start() + 1 + len(name):end])


def _first_text(text: str, name: str, start: int = 0, stop: int = -1) -> Optional[str]:
    """Text of the first <name>text</name> element at or after `start`."""
    m = _re("<" + re.escape(name) + r"\s*>([^<]*)</" + re.escape(name) + r"\s*>").search(
        text, start, stop if stop >= 0 else len(text))
    if not m:
        return None
    return _unescape(m.group(1)).strip()


def _count(text: str, name: str, start: int, stop: int) -> int:
    return len(_re("<" + re.escape(name) + r"[\s/>]").findall(text, start, stop))


# --- FlowJo -----------------------------------------------------------------

def sniff_flowjo(src):
    text, _ = _text(src)
    if text is None:
        return None
    root = _root(text)
    if root is None:
        return None
    name, attrs, end = root
    if name != "Workspace" or "flowJoVersion" not in attrs:
        return None
    cyt = _first_tag(text, "Cytometer", end) or {}
    samples = None
    sl = _re(r"<SampleList[\s/>]").search(text, end)
    if sl:
        close = text.find("</SampleList>", sl.end())
        if close >= 0:
            samples = _count(text, "Sample", sl.end(), close)
    version = _value(attrs.get("flowJoVersion"))
    model = _value(cyt.get("cyt"))
    fields = {
        "softwareVersion": version,
        "workspaceVersion": _value(attrs.get("version")),
        "instrumentModel": model,
        "serialNumber": _value(cyt.get("serialnumber")),
        "sampleCount": samples,
        "exportFormat": "workspace xml",
    }
    summary = _summary(["FlowJo workspace" + (" " + version if version else ""), model,
                        _plural(samples, "sample")])
    return header(FLOWJO, summary, fields)


# --- BD FACSDiva ------------------------------------------------------------

def sniff_facsdiva(src):
    text, _ = _text(src)
    if text is None:
        return None
    root = _root(text)
    if root is None:
        return None
    name, attrs, end = root
    if name != "bdfacs" or "version" not in attrs:
        return None
    version = _value(attrs.get("version"))
    if version and version.startswith("Version "):
        version = _value(version[len("Version "):])
    exp = _re(r"<experiment[\s>]").search(text, end)
    exp_name = exp_date = export_time = None
    if exp:
        exp_attrs = _first_tag(text, "experiment", exp.start()) or {}
        exp_name = _value(exp_attrs.get("name"))
        exp_date = _date(_first_text(text, "date", exp.end(), exp.end() + 4096))
        export_time = _date(_first_text(text, "export_time", exp.end(), exp.end() + 4096))
    tubes = None
    close = text.find("</bdfacs>", end)
    if close >= 0:
        tubes = _count(text, "tube", end, close)
    model = _value(_first_text(text, "data_instrument_name", end))
    fields = {
        "softwareVersion": version,
        "instrumentModel": model,
        "serialNumber": _value(_first_text(text, "data_instrument_serial_number", end)),
        "experimentName": exp_name,
        "experimentDate": exp_date,
        "exportDate": export_time,
        "sampleCount": tubes,
        "exportFormat": "experiment xml",
    }
    summary = _summary(["FACSDiva experiment" + (" " + version if version else ""), model,
                        _plural(tubes, "tube")])
    return header(FACSDIVA, summary, fields)


# --- Bio-Rad Bio-Plex Manager -----------------------------------------------

def sniff_bioplex(src):
    text, _ = _text(src)
    if text is None:
        return None
    root = _root(text)
    if root is None:
        return None
    name, attrs, end = root
    if name.split(":")[-1] != "BioPlexManager" or "BioPlexManagerVersion" not in attrs:
        return None
    dims = _first_tag(text, "PlateDimensions", end) or {}
    mi = text.find("<MachineInfo>", end)
    serial = firmware = None
    if mi >= 0:
        mi_end = text.find("</MachineInfo>", mi)
        if mi_end >= 0:
            serial = _value(_first_text(text, "SerialNumber", mi, mi_end))
            firmware = _value(_first_text(text, "MicroControllerVersion", mi, mi_end))
    analytes = None
    roi = text.find("<RegionsOfInterest>", end)
    if roi >= 0:
        roi_end = text.find("</RegionsOfInterest>", roi)
        if roi_end >= 0:
            analytes = _count(text, "Region", roi, roi_end)
    wells = None
    w = _re(r"<Wells[\s>]").search(text, end)
    if w:
        w_end = text.find("</Wells>", w.end())
        if w_end >= 0:
            wells = _count(text, "Well", w.end(), w_end)
    version = _value(attrs.get("BioPlexManagerVersion"))
    plate = _int(dims.get("TotalWells"))
    fields = {
        "softwareVersion": version,
        "serialNumber": serial,
        "firmwareVersion": firmware,
        "plateFormat": plate,
        "analyteCount": analytes,
        "measurementCount": wells,
        "exportDate": _value(attrs.get("XmlExportDateTime")),
        "exportFormat": "xml",
    }
    summary = _summary(["Bio-Plex Manager" + (" " + version if version else "") + " export",
                        "%d-well plate" % plate if plate else None,
                        _plural(analytes, "analyte"), _plural(wells, "well")])
    return header(BIOPLEX, summary, fields)


# --- Luminex ----------------------------------------------------------------

_REGION_COL = _re(r"^R(\d+):")
_BUNDLE = _re(r"INTELLIFLEX Bundle - ([^/]+)")
_FIRMWARE = _re(r"Firmware Executable Version - ([^/]+)")
_PLATE_GEOMETRY = _re(r"^Plate(\d+)Well$")
_PLATE_TYPE = _re(r"^(\d+)[ -]?[Ww]ell")


def _xponent(lines: List[str], complete: bool):
    first = None
    for idx, ln in enumerate(lines[:20]):
        if not _blank(ln):
            first = idx
            break
    if first is None:
        return None
    cells = _cells(lines[first])
    if len(cells) < 2 or cells[0].strip() != "Program" or cells[1].strip() != "xPONENT":
        return None
    model = _value(cells[3]) if len(cells) > 3 else None
    kv: Dict[str, List[str]] = {}
    data_at = None
    for idx in range(first, len(lines)):
        ln = lines[idx]
        if _blank(ln):
            continue
        c = _cells(ln)
        key = c[0].strip()
        if key.startswith("DataType:"):
            data_at = idx
            break
        if key and key not in kv:
            kv[key] = c[1:]
    if "Build" not in kv and "SN" not in kv:
        return None

    def get(key, i=0):
        v = kv.get(key)
        return _value(v[i]) if v and len(v) > i else None

    plate = plate_count = None
    pp = kv.get("ProtocolPlate")
    if pp:
        for i in range(0, len(pp) - 1, 2):
            if pp[i].strip() == "Type":
                plate = _int(pp[i + 1])
            elif pp[i].strip() == "Plates":
                plate_count = _int(pp[i + 1])
    analytes = rows = samples = None
    names: List[str] = []
    if data_at is not None and data_at + 1 < len(lines):
        hdr = [h.strip() for h in _cells(lines[data_at + 1])]
        if len(hdr) >= 2 and hdr[0] == "Location" and hdr[1] == "Sample":
            cols = [h for h in hdr[2:] if h]
            if cols and cols[-1] == "Total Events":
                cols = cols[:-1]
            analytes = len(cols)
            n = 0
            ended = False
            for ln in lines[data_at + 2:]:
                if _blank(ln):
                    ended = True
                    break
                c = _cells(ln)
                names.append(c[1] if len(c) > 1 else "")
                n += 1
            if ended or complete:
                rows = n
                samples = len(_distinct(names))
            else:
                names = []
    return {
        "model": model,
        "fields": {
            "softwareVersion": get("Build"),
            "instrumentModel": model,
            "serialNumber": get("SN"),
            "batchName": get("Batch"),
            "protocolName": get("ProtocolName"),
            "runDate": _date(get("BatchStartTime")),
            "plateFormat": plate,
            "plateCount": plate_count,
            "analyteCount": analytes,
            "measurementCount": rows,
            "sampleCount": samples,
            "sampleNames": _distinct(names) if rows is not None and samples <= MAX_SAMPLE_NAMES else None,
            "exportFormat": "batch csv",
        },
    }


def _flat(rows: List[List[str]], complete: bool, export: str, date_of=None):
    """The INTELLIFLEX flat results table: one header row, one row per well."""
    if not rows:
        return None
    if date_of is None:
        date_of = _date
    hdr = [h.strip() for h in rows[0]]
    if not {"INSTRUMENT TYPE", "WELL LOCATION", "SAMPLE ID"}.issubset(hdr):
        return None
    data: List[List[str]] = []
    signed = False
    ended = False
    for r in rows[1:]:
        if r and r[0].strip().startswith("--DOCUMENT SIGNATURE--"):
            signed = True
            ended = True
            break
        if all(not c.strip() for c in r):
            ended = True
            break
        data.append(r)
    if not data:
        return None
    col = {h: i for i, h in reversed(list(enumerate(hdr)))}

    def cell(row, key):
        i = col.get(key)
        return row[i] if i is not None and i < len(row) else None

    first = data[0]
    itype = _value(cell(first, "INSTRUMENT TYPE"))
    sw = _value(cell(first, "SOFTWARE VERSION"))
    if not ((itype and "INTELLIFLEX" in itype.upper()) or (sw and "INTELLIFLEX" in sw.upper())):
        return None
    version = firmware = None
    if sw:
        m = _BUNDLE.search(sw)
        version = m.group(1).strip() if m else sw
        m = _FIRMWARE.search(sw)
        firmware = m.group(1).strip() if m else None
    plate = None
    geo = _value(cell(first, "PLATE GEOMETRY"))
    if geo and _PLATE_GEOMETRY.match(geo):
        plate = _int(_PLATE_GEOMETRY.match(geo).group(1))
    ptype = _value(cell(first, "PLATE TYPE NAME"))
    if plate is None and ptype and _PLATE_TYPE.match(ptype):
        plate = _int(_PLATE_TYPE.match(ptype).group(1))
    regions = len({m.group(1) for m in (_REGION_COL.match(h) for h in hdr) if m})
    rows_n = samples = None
    names: List[str] = []
    if (complete or ended) and all(len(r) == len(rows[0]) for r in data):
        rows_n = len(data)
        names = _distinct([cell(r, "SAMPLE ID") or "" for r in data])
        samples = len(names)
    model = _value(cell(first, "MODEL NAME"))
    return {
        "model": model or itype,
        "fields": {
            "softwareVersion": version,
            "firmwareVersion": firmware,
            "instrumentModel": model,
            "instrumentType": itype,
            "serialNumber": _value(cell(first, "SERIAL NUMBER")),
            "plateName": _value(cell(first, "PLATE NAME")),
            "protocolName": _value(cell(first, "PROTOCOL NAME")),
            "runDate": date_of(cell(first, "PLATE START")),
            "plateFormat": plate,
            "regionCount": regions or None,
            "measurementCount": rows_n,
            "sampleCount": samples,
            "sampleNames": names if rows_n is not None and samples <= MAX_SAMPLE_NAMES else None,
            "documentSigned": True if signed else None,
            "exportFormat": export,
        },
    }


def _luminex_header(r):
    f = r["fields"]
    model = r["model"]
    name = XPONENT
    if f.get("exportFormat") != "batch csv" or (model and model.lower() == "intelliflex"):
        name = INTELLIFLEX
    version = f.get("softwareVersion")
    label = "xPONENT" if name == XPONENT else "INTELLIFLEX"
    summary = _summary([label + " " + f["exportFormat"] + (" " + version if version else ""), model,
                        _plural(f.get("analyteCount"), "analyte") or _plural(f.get("regionCount"), "region"),
                        _plural(f.get("measurementCount"), "well")])
    return header(name, summary, f)


def sniff_luminex_csv(src):
    text, complete = _text(src)
    if text is None:
        return None
    lines = _lines(text, complete)
    r = _xponent(lines, complete)
    if r is None:
        first = next((ln for ln in lines[:1]), "")
        if "INSTRUMENT TYPE" not in first:
            return None
        rows = [_cells(ln) for ln in lines]
        r = _flat(rows, complete, "results csv")
    if r is None:
        return None
    return _luminex_header(r)


XLSX_MAX_ROWS = 4000
_SERIAL = _re(r"^\d{1,6}(\.\d+)?$")


def _excel_date(v: Optional[str], date1904: bool) -> Optional[str]:
    """A worksheet cell holding a date: an Excel serial day number becomes ISO
    8601 (1900 or 1904 date system as the workbook declares), text is read the
    same way as in the CSV export."""
    if v is None or not _SERIAL.match(v.strip()):
        return _date(v)
    import datetime

    try:
        days = float(v)
        if days < 61 and not date1904:
            return None  # the 1900 leap-year bug makes early serials ambiguous
        base = datetime.datetime(1904, 1, 1) if date1904 else datetime.datetime(1899, 12, 30)
        t = base + datetime.timedelta(seconds=round(days * 86400))
        return t.strftime("%Y-%m-%dT%H:%M:%S")
    except (OverflowError, ValueError):
        return None


def sniff_luminex_xlsx(src):
    names = xlsx_sheet_names(src)
    if not names or "Results" not in names:
        return None
    rows = xlsx_rows(src, "Results", max_rows=XLSX_MAX_ROWS)
    if not rows:
        return None
    width = len(rows[0])
    padded = [r + [""] * (width - len(r)) if len(r) < width else r for r in rows]
    wb = zip_member(src._instrument_zip, "xl/workbook.xml", 256 * 1024) or b""
    date1904 = re.search(rb"date1904\s*=\s*[\"'](1|true)[\"']", wb) is not None
    r = _flat(padded, len(rows) < XLSX_MAX_ROWS, "results xlsx",
              lambda v: _excel_date(v, date1904))
    if r is None:
        return None
    return _luminex_header(r)


# --- Beckman Coulter Biomek ---------------------------------------------------

_BIOMEK_TS = _re(r"^\d{1,2}[/-]\d{1,2}[/-]\d{2,4}\s+\d{1,2}:\d{2}")
_POD_SERIAL = _re(r"^Pod(\d+) head serial number$")
_BIOMEK_LAYOUT = {9: "pipetting", 11: "unified pipetting", 14: "unified transfer"}


def sniff_biomek(src):
    text, complete = _text(src)
    if text is None:
        return None
    lines = _lines(text, complete)
    if not lines or not lines[0].startswith("Method ="):
        return None
    kv: Dict[str, str] = {}
    started = None
    pods: List[Tuple[int, Optional[str]]] = []
    data_from = None
    for idx, ln in enumerate(lines[:40]):
        first = _cells(ln)[0].strip()
        if _BIOMEK_TS.match(first):
            data_from = idx
            break
        if first.startswith("Time Stamp") or "Well Index" in ln:
            data_from = idx + 1
            break
        if first.startswith("Started "):
            started = first[len("Started "):].strip()
            continue
        if "=" in first:
            k, v = first.split("=", 1)
            k = k.strip()
            if k not in kv:
                kv[k] = v.strip()
                m = _POD_SERIAL.match(k)
                if m:
                    pods.append((int(m.group(1)), _value(v)))
    if "Unit serial number" not in kv or "Logged in user" not in kv:
        return None
    rows = layout = None
    if data_from is not None:
        data = [ln for ln in lines[data_from:] if ln.strip() and _BIOMEK_TS.match(ln)]
        if data:
            layout = _BIOMEK_LAYOUT.get(len(data[0].split(",")))
        if complete:
            rows = len(data)
    serials = [s for _, s in pods if s]
    method = _value(kv.get("Method"))
    fields = {
        "methodName": method,
        "serialNumber": _value(kv.get("Unit serial number")),
        "podCount": len(pods) or None,
        "podHeadSerialNumbers": serials or None,
        "runDate": _date(started),
        "measurementCount": rows,
        "logLayout": layout,
        "exportFormat": "transfer log",
    }
    summary = _summary(["Biomek transfer log", method, _plural(rows, "transfer step")])
    return header(BIOMEK, summary, fields)


# --- Beckman Echo -------------------------------------------------------------

_SECTION = _re(r"^\[(.+)\]\s*$")


def sniff_echo(src):
    text, complete = _text(src)
    if text is None:
        return None
    lines = _lines(text, complete)
    if not lines or not lines[0].startswith("Run ID,"):
        return None
    kv: Dict[str, str] = {}
    sections: List[str] = []
    counts: Dict[str, int] = {}
    i, n = 0, len(lines)
    while i < n:
        ln = lines[i]
        m = _SECTION.match(ln)
        if m:
            title = m.group(1).strip()
            sections.append(title)
            j = i + 2
            ended = False
            while j < n:
                if not lines[j].strip():
                    ended = True
                    break
                j += 1
            if (ended or complete) and i + 1 < n and title not in counts:
                counts[title] = max(0, j - i - 2)
            i = j + 1
            continue
        if ln.strip():
            c = _cells(ln)
            k = c[0].strip().replace("*", "")
            if k and k not in kv:
                kv[k] = c[1] if len(c) > 1 else ""
        i += 1
    app = kv.get("Application Name", "")
    if "Plate Reformat" in app:
        name = ECHO_REFORMAT
    elif "Cherry Pick" in app:
        name = ECHO_CHERRY
    else:
        return None
    if not sections:
        return None
    version = _value(kv.get("Application Version"))
    model = _value(kv.get("Instrument Model"))
    transfers = counts.get("DETAILS")
    exceptions = counts.get("EXCEPTIONS")
    fields = {
        "applicationName": _value(app),
        "softwareVersion": version,
        "instrumentModel": model,
        "serialNumber": _value(kv.get("Instrument Serial Number")),
        "instrumentSoftwareVersion": _value(kv.get("Instrument Software Version")),
        "runId": _value(kv.get("Run ID")),
        "runDate": _date(kv.get("Run Date/Time")),
        "protocolName": _value(kv.get("Protocol Name")),
        "sections": sections,
        "transferCount": transfers,
        "exceptionCount": exceptions,
        "exportFormat": "transfer report csv",
    }
    label = "Echo Plate Reformat" if name == ECHO_REFORMAT else "Echo Cherry Pick"
    summary = _summary([label + " report" + (" " + version if version else ""), model,
                        _plural(transfers, "transfer"), _plural(exceptions, "exception")])
    return header(name, summary, fields)


register_instrument_sniffer(["wsp", "xml"], sniff_flowjo)
register_instrument_sniffer("xml", sniff_facsdiva)
register_instrument_sniffer("xml", sniff_bioplex)
register_instrument_sniffer("csv", sniff_luminex_csv)
register_instrument_sniffer("xlsx", sniff_luminex_xlsx)
register_instrument_sniffer(["log", "csv"], sniff_biomek)
register_instrument_sniffer("csv", sniff_echo)
