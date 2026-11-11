"""Plate readers, family B: PerkinElmer Envision, Revvity Kaleido, MSD
Methodical Mind, MSD Discovery Workbench, CTL ImmunoSpot, Mabtech Apex and
Unchained Labs Lunatic & Stunner.

Text and CSV exports are read from their leading bytes; the Go port
(scigantic-edge header/instrument_plate_b.go) decodes the same files to the
same fields. The .xlsx exports (Workbench, Apex, Lunatic & Stunner) are read
here only.

Counts (plates, samples, measurements) are reported only when the whole file
fit in the leading read, so a truncated read never reports a partial count.
Dates are reported as ISO 8601 only when the file writes an unambiguous date
(YYYY-MM-DD, or a month name); a numeric day/month date is kept as written.

Telling the two MSD products apart: Methodical Mind (the MSD reader software,
"MMPR") writes a tab-separated .txt with a "Stack ID :" / "Plate # :" header
block and a "=====Data=====" rule. Discovery Workbench writes a CSV whose first
line is "Plate_<barcode>" and whose second line is the
"Sample,Assay,Well,Spot,..." table header, or an .xlsx with a
"Workbench data" sheet in that shape. allotropy's msd_workbench parser also
accepts the Methodical Mind .txt; those files are reported as Methodical Mind.
"""

from __future__ import annotations

import csv
import datetime as _dt
import io
import re
from typing import Dict, List, Optional, Tuple

from . import (
    TEXT_HEAD_BYTES,
    header,
    head,
    head_text,
    instrument,
    register_instrument_sniffer,
    xlsx_rows,
    xlsx_sheet_names,
    zip_member,
)

ENVISION = "PerkinElmer Envision"
KALEIDO = "Revvity Kaleido"
METHODICAL_MIND = "MSD Methodical Mind"
WORKBENCH = "MSD Discovery Workbench"
IMMUNOSPOT = "CTL ImmunoSpot"
APEX = "Mabtech Apex"
LUNATIC = "Unchained Labs Lunatic & Stunner"

instrument(ENVISION, "Plate Reader", "PerkinElmer", "perkin_elmer_envision",
           ["Absorbance", "Fluorescence", "Luminescence"])
instrument(KALEIDO, "Plate Reader", "Revvity", "revvity_kaleido",
           ["Absorbance", "Fluorescence", "Luminescence", "Optical Imaging"])
instrument(METHODICAL_MIND, "Plate Reader", "Meso Scale Discovery", "methodical_mind",
           ["Luminescence"])
instrument(WORKBENCH, "Plate Reader", "Meso Scale Discovery", "msd_workbench",
           ["Luminescence"])
instrument(IMMUNOSPOT, "Plate Reader", "CTL", "ctl_immunospot", ["Optical Imaging"])
instrument(APEX, "Plate Reader", "Mabtech", "mabtech_apex", ["Optical Imaging"])
instrument(LUNATIC, "Plate Reader", "Unchained Labs", "unchained_labs_lunatic_stunner",
           ["Absorbance", "Dynamic Light Scattering"])

# A list longer than this is summarized rather than listed.
MAX_LIST = 32


# ---- shared helpers (mirrored in Go) -------------------------------------

_ISO = re.compile(r"^(\d{4}-\d{2}-\d{2})[ T](\d{2}:\d{2}(?::\d{2}(?:\.\d+)?)?)(Z|[+-]\d{2}:\d{2})?$")
_MONTHS = {m: i + 1 for i, m in enumerate(
    ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])}
_CTIME = re.compile(r"^(?:[A-Z][a-z]{2} )?([A-Z][a-z]{2}) +(\d{1,2}) (\d{2}):(\d{2}):(\d{2}) (\d{4})$")


def _date(s: Optional[str]) -> Optional[str]:
    """ISO 8601 for an unambiguous date-time, else the string as written."""
    if s is None:
        return None
    s = s.strip()
    if not s:
        return None
    m = _ISO.match(s)
    if m:
        return "%sT%s%s" % (m.group(1), m.group(2), m.group(3) or "")
    m = _CTIME.match(s)
    if m and m.group(1) in _MONTHS:
        mon, day, hh, mm, ss, year = _MONTHS[m.group(1)], int(m.group(2)), m.group(3), m.group(4), m.group(5), m.group(6)
        if 1 <= day <= 31 and int(hh) < 24 and int(mm) < 60 and int(ss) < 61:
            return "%s-%02d-%02dT%s:%s:%s" % (year, mon, day, hh, mm, ss)
    return s


def _int(s: Optional[str]) -> Optional[int]:
    if s is None:
        return None
    s = s.strip()
    if re.match(r"^\d{1,9}$", s):
        return int(s)
    if re.match(r"^\d{1,9}\.0+$", s):
        return int(s.split(".")[0])
    return None


def _lead_int(s: Optional[str]) -> Optional[int]:
    m = re.match(r"^\s*(\d+)", s or "")
    return int(m.group(1)) if m and len(m.group(1)) <= 6 else None


def _plural(n: int, word: str) -> str:
    return "%d %s%s" % (n, word, "" if n == 1 else "s")


def _summary(name: str, f: Dict[str, object]) -> str:
    parts = [name + " export"]
    if f.get("softwareVersion"):
        parts.append("software " + str(f["softwareVersion"]))
    if f.get("instrumentModel"):
        parts.append("model " + str(f["instrumentModel"]))
    if f.get("plateFormat"):
        parts.append("%d-well plate" % f["plateFormat"])
    if f.get("plateCount"):
        parts.append(_plural(f["plateCount"], "plate"))  # type: ignore[arg-type]
    if f.get("sampleCount"):
        parts.append(_plural(f["sampleCount"], "sample"))  # type: ignore[arg-type]
    if f.get("wavelengthsNm"):
        parts.append("/".join(str(w) for w in f["wavelengthsNm"]) + " nm")  # type: ignore[union-attr]
    return ", ".join(parts)


def _result(name: str, f: Dict[str, object]) -> object:
    f = {k: v for k, v in f.items() if v is not None and v != "" and v != []}
    return header(name, _summary(name, f), f)


def _complete(src) -> bool:
    """True when the leading read held the whole file."""
    return src.size is not None and src.size <= TEXT_HEAD_BYTES


def _first(f: Dict[str, object], k: str, v) -> None:
    """Set f[k] to the first non-empty value seen."""
    if v is not None and v != "" and k not in f:
        f[k] = v


def _uniq(xs):
    out = []
    for x in xs:
        if x not in out:
            out.append(x)
    return out


def _lines(text: str) -> List[str]:
    """Lines split on CR LF, LF or CR only (as the Go port does)."""
    out = re.split(r"\r\n|\r|\n", text)
    if out and out[-1] == "":
        out.pop()
    return out


def _cells(line: str) -> List[str]:
    return [c.strip() for c in line.split(",")]


def _last_cell(line: str) -> str:
    return line.rsplit(",", 1)[-1].strip()


# ---- PerkinElmer Envision (.csv) ------------------------------------------

_ENV_FILTER_KEYS = ("Exc. filter", "2nd exc. filter", "Ems. filter", "2nd ems. filter")


def sniff_envision(src):
    text = head_text(src)
    if not text:
        return None
    lines = _lines(text)
    ok = False
    for i, line in enumerate(lines[:20]):
        if line.startswith("Plate information") and i + 1 < len(lines) \
                and lines[i + 1].startswith("Plate,Repeat,Barcode,"):
            ok = True
            break
    if not ok:
        return None
    f: Dict[str, object] = {}
    plates: List[str] = []
    filter_names: List[str] = []
    filters: Dict[str, int] = {}
    section = ""
    block_name = None
    prev_blank = True
    for i, line in enumerate(lines):
        s = line.strip()
        if line.startswith("Plate,Repeat,Barcode,") and i + 1 < len(lines):
            p = _cells(lines[i + 1])[0]
            if p and p not in plates:
                plates.append(p)
        if s.endswith(":") and "," not in s:
            section = s[:-1]
            block_name = None
            prev_blank = True
            continue
        if not s:
            prev_blank = True
            continue
        cells = _cells(line)
        if section == "Labels" and cells[0] in _ENV_FILTER_KEYS:
            v = _last_cell(line)
            if v and v != "N/A":
                filter_names.append(v)
        elif section == "Filters":
            if prev_blank:
                block_name = cells[0]
            elif cells[0] == "Description" and block_name:
                m = re.search(r"Longpass=(\d+)nm", line) or re.search(r"CWL=(\d+)nm", line)
                if m and block_name not in filters:
                    filters[block_name] = int(m.group(1))
        elif section == "Instrument" and cells[0] == "Serial number":
            _first(f, "serialNumber", _last_cell(line) or None)
        if cells[0] == "Number of the wells in the plate":
            _first(f, "plateFormat", _int(_last_cell(line)))
        elif cells[0] == "Assay Started:":
            _first(f, "runDate", _date(_last_cell(line)))
        elif cells[0] == "Protocol Name:":
            _first(f, "protocolName", _last_cell(line) or None)
        elif cells[0] == "Export format" and section == "Auto export parameters":
            _first(f, "exportFormat", _last_cell(line) or None)
        m = re.match(r"^Exported with (.+) version (\S+)", s)
        if m:
            f["softwareVersion"] = m.group(2)
        prev_blank = False
    waves = _uniq([filters[n] for n in filter_names if n in filters])
    if waves:
        f["wavelengthsNm"] = waves
    if _complete(src) and plates:
        f["plateCount"] = len(plates)
    return _result(ENVISION, f)


# ---- Revvity Kaleido (.csv) -----------------------------------------------

_KAL_SOFTWARE = re.compile(r"^Software [Vv]ersion:?,+Kaleido (\S+?),*$")
_KAL_EXPORTED = re.compile(r"^Exported with Kaleido (\d[^\s,]*)")


def _kv(line: str) -> Tuple[str, str]:
    """'Key:,,value,,,' or 'Key,,value' -> (Key, value)."""
    cells = _cells(line)
    key = cells[0].rstrip(":").strip()
    vals = [c for c in cells[1:] if c]
    return key, (vals[0] if vals else "")


def sniff_kaleido(src):
    text = head_text(src)
    if not text:
        return None
    # Kaleido 2 writes a no-break space after "Results for" and "Barcode:".
    lines = _lines(text.replace("\u00a0", " "))
    first = next((l.strip() for l in lines if l.strip().strip(",")), "")
    v3 = first.startswith("EnSight Results from")
    if not (v3 or re.match(r"^Results? for ", first)):
        return None
    version = None
    for line in lines:
        s = line.strip()
        m = _KAL_SOFTWARE.match(s)
        if m:
            version = m.group(1)
            break
    if version is None:
        for line in lines:
            m = _KAL_EXPORTED.match(line.strip())
            if m:
                version = m.group(1)
                break
    if version is None:
        return None
    f: Dict[str, object] = {"softwareVersion": version}
    if v3:
        f["instrumentModel"] = "EnSight"
    rows = cols = None
    waves: List[int] = []
    modes: List[str] = []
    for line in lines:
        s = line.strip()
        m = re.match(r"^Results? for (.+?) \d+,*$", s)
        if m and "(" not in m.group(1):
            if m.group(1) not in modes:
                modes.append(m.group(1))
            continue
        if "," not in s:
            continue
        key, val = _kv(s)
        if key == "Instrument Serial Number":
            _first(f, "serialNumber", val or None)
        elif key == "Measurement Started":
            _first(f, "runDate", _date(val))
        elif key == "Protocol Name":
            _first(f, "protocolName", val or None)
        elif key in ("Number of rows", "Number of Rows"):
            rows = rows if rows is not None else _int(val)
        elif key in ("Number of columns", "Number of Columns"):
            cols = cols if cols is not None else _int(val)
        elif key in ("Excitation wavelength [nm]", "Excitation Wavelength [nm]", "Emission Wavelength [nm]"):
            w = _lead_int(val)
            if w and w not in waves:
                waves.append(w)
        elif key == "Export Format":
            _first(f, "exportFormat", val or None)
    if rows and cols:
        f["plateFormat"] = rows * cols
    f["wavelengthsNm"] = waves
    f["readModes"] = modes
    return _result(KALEIDO, f)


# ---- MSD Methodical Mind (.txt) ------------------------------------------

_MM_DATA = re.compile(r"^=+Data=+\s*$")


def _mm_kv(line: str) -> Optional[Tuple[str, str]]:
    if ":" not in line:
        return None
    k, v = line.split(":", 1)
    v = v.replace("\\t", "\t")
    v = re.sub(r"\s+", " ", v).strip()
    return k.strip(), v


def sniff_methodical_mind(src):
    text = head_text(src)
    if not text:
        return None
    lines = _lines(text)
    stack = plate = data = False
    for line in lines[:60]:
        if re.match(r"^Stack ID\s*:", line):
            stack = True
        elif re.match(r"^Plate #\s*:", line):
            plate = True
        elif _MM_DATA.match(line):
            data = True
            break
    if not (stack and plate and data):
        return None
    f: Dict[str, object] = {}
    plates: List[Tuple[str, str]] = []
    cur_plate = ""
    in_header = True
    for line in lines:
        if _MM_DATA.match(line):
            in_header = False
            continue
        if re.match(r"^Stack ID\s*:", line):
            in_header = True
        if not in_header:
            continue
        kv = _mm_kv(line)
        if not kv:
            continue
        k, v = kv
        if k == "Plate #":
            cur_plate = v
        elif k == "Barcode1":
            pair = (cur_plate, v)
            if pair not in plates:
                plates.append(pair)
        if k in ("Version", "Serial #", "Model", "Type", "Spots Per Well", "Read Time"):
            if k not in f and v:
                f[k] = v
    out: Dict[str, object] = {
        "softwareVersion": f.get("Version"),
        "instrumentModel": f.get("Model"),
        "serialNumber": f.get("Serial #"),
        "plateFormat": _lead_int(f.get("Type")),  # type: ignore[arg-type]
        "runDate": _date(f.get("Read Time")),  # type: ignore[arg-type]
    }
    spots = _int(f.get("Spots Per Well"))  # type: ignore[arg-type]
    if spots is None:
        m = re.search(r"Multi-Spot (\d+)", str(f.get("Type") or ""))
        spots = int(m.group(1)) if m else None
    out["spotsPerWell"] = spots
    if _complete(src) and plates:
        out["plateCount"] = len(plates)
    return _result(METHODICAL_MIND, out)


# ---- MSD Discovery Workbench (.csv, .xlsx) --------------------------------

_WB_HEADER = ["Sample", "Assay", "Well", "Spot"]


def _workbench_rows(rows: List[List[str]], complete: bool, fmt: str):
    if len(rows) < 2:
        return None
    first = [c.strip() for c in rows[0]]
    if not first or not re.match(r"^Plate_\S", first[0]) \
            or any(c and not re.match(r"^Unnamed: \d+$", c) for c in first[1:]):
        return None
    hdr = [c.strip() for c in rows[1]]
    if hdr[:4] != _WB_HEADER:
        return None
    f: Dict[str, object] = {"exportFormat": fmt}
    if complete:
        samples: List[str] = []
        assays: List[str] = []
        n = 0
        for r in rows[2:]:
            r = [c.strip() for c in r]
            if not any(r):
                continue
            n += 1
            if r[0] and r[0] not in samples:
                samples.append(r[0])
            if len(r) > 1 and r[1] and r[1] not in assays:
                assays.append(r[1])
        f["measurementCount"] = n
        f["sampleCount"] = len(samples)
        if len(assays) <= MAX_LIST:
            f["assays"] = assays
    return _result(WORKBENCH, f)


def sniff_workbench_csv(src):
    text = head_text(src)
    if not text or not text.startswith("Plate_"):
        return None
    lines = _lines(text)
    complete = _complete(src)
    if not complete and lines:
        lines = lines[:-1]
    return _workbench_rows([l.split(",") for l in lines], complete, "csv")


_XLSX_MAX_ROWS = 50000


def sniff_workbench_xlsx(src):
    names = xlsx_sheet_names(src)
    if not names or "Workbench data" not in names:
        return None
    rows = xlsx_rows(src, "Workbench data", _XLSX_MAX_ROWS)
    if rows is None:
        return None
    return _workbench_rows(rows, len(rows) < _XLSX_MAX_ROWS, "xlsx")


# ---- CTL ImmunoSpot (.txt) -----------------------------------------------

def _ctl_grid(lines: List[str]) -> Optional[int]:
    for i, line in enumerate(lines):
        toks = line.split()
        if len(toks) < 2 or not all(t.isdigit() for t in toks):
            continue
        n = 0
        while n < len(toks) and toks[n] == str(n + 1):
            n += 1
        if n < 2:
            continue
        rows = 0
        for row in lines[i + 1:]:
            rt = row.split()
            if not rt or all(set(t) <= {"_"} for t in rt):
                continue
            if rows < 26 and rt[0] == chr(ord("A") + rows):
                rows += 1
                continue
            break
        if rows:
            return rows * n
        return None
    return None


def sniff_immunospot(src):
    text = head_text(src)
    if not text:
        return None
    lines = _lines(text)
    m = re.match(r"^Software version:\s*ImmunoSpot\s+(\d[\d.]*)\s*$", lines[0] if lines else "")
    if not m:
        return None
    f: Dict[str, object] = {"softwareVersion": m.group(1), "exportFormat": "counts txt"}
    for line in lines[:12]:
        s = line.strip()
        for part in s.split(";"):
            part = part.strip()
            mm = re.match(r"^Analyzer Serial number:\s*(\S+)", part)
            if mm:
                _first(f, "serialNumber", mm.group(1))
            mm = re.match(r"^Counted:\s*(.+)$", part)
            if mm:
                _first(f, "runDate", _date(mm.group(1)))
            mm = re.match(r"^QC review last updated on (.+?) by ", part)
            if mm:
                _first(f, "qcReviewDate", _date(mm.group(1)))
                f["exportFormat"] = "QC review txt"
    f["plateFormat"] = _ctl_grid(lines)
    return _result(IMMUNOSPOT, f)


# ---- .xlsx helpers --------------------------------------------------------

def _date1904(src) -> bool:
    zf = getattr(src, "_instrument_zip", None)
    if zf is None:
        return False
    wb = zip_member(zf, "xl/workbook.xml", 256 * 1024) or b""
    return bool(re.search(rb'date1904="(1|true)"', wb))


def _xlsx_date(src, v: Optional[str], t: Optional[str] = None) -> Optional[str]:
    """A date cell: an Excel serial (with an optional time-of-day serial) as
    ISO 8601, or a text date through `_date`."""
    if v is None or not v.strip():
        return None
    v = v.strip()
    try:
        days = float(v)
    except ValueError:
        return _date(v)
    frac = 0.0
    if t:
        try:
            frac = float(t)
        except ValueError:
            return None
        if not 0 <= frac < 1:
            return None
    if not 1 <= days < 2958466:
        return None
    base = _dt.datetime(1904, 1, 1) if _date1904(src) else _dt.datetime(1899, 12, 30)
    secs = int(round((days + frac) * 86400))
    return (base + _dt.timedelta(seconds=secs)).strftime("%Y-%m-%dT%H:%M:%S")


def _id_text(v: Optional[str]) -> Optional[str]:
    if v is None:
        return None
    v = v.strip()
    if re.match(r"^\d+\.0+$", v):
        v = v.split(".")[0]
    return v or None


# ---- Mabtech Apex (.xlsx) -------------------------------------------------

def sniff_apex(src):
    names = xlsx_sheet_names(src)
    if not names or "Plate Database" not in names:
        return None
    info = "Plate Information" if "Plate Information" in names else ("Plate Info" if "Plate Info" in names else None)
    if info is None:
        return None
    rows = xlsx_rows(src, info, 60)
    if not rows:
        return None
    kv: Dict[str, str] = {}
    for r in rows:
        cells = [c for c in r if c.strip()]
        if len(cells) >= 2 and cells[0].strip().endswith(":"):
            kv.setdefault(cells[0].strip()[:-1].strip(), cells[1].strip())
    if "Machine ID" not in kv and "Plate Read Date" not in kv:
        return None
    f: Dict[str, object] = {
        "softwareVersion": kv.get("Software Version"),
        "machineId": _id_text(kv.get("Machine ID")),
        "runDate": _xlsx_date(src, kv.get("Plate Read Date")),
        "exportFormat": "xlsx",
    }
    m = re.match(r"^([A-Za-z]+)\s+\d+$", str(f["machineId"] or ""))
    if m:
        f["instrumentModel"] = m.group(1)
    return _result(APEX, f)


# ---- Unchained Labs Lunatic & Stunner (.csv, .xlsx) ------------------------

_LUN_WAVE = re.compile(r"^a(\d{3})(?: \((\d+)?mm\))?$")
_LUN_DLS = ("z-avg dia. (nm)", "z ave. dia (nm)", "pdi", "diffusion coefficient (um^2/s)")


def _norm(c: str) -> str:
    return re.sub(r"\s+", " ", c).strip().lower()


def _lunatic_rows(src, rows: List[List[str]], complete: bool, fmt: str, xlsx: bool):
    hi = None
    for i, r in enumerate(rows[:80]):
        n = [_norm(c) for c in r]
        if "plate id" in n and "plate position" in n and "sample name" in n:
            hi = i
            break
    if hi is None:
        return None
    hdr = [_norm(c) for c in rows[hi]]
    waves: List[int] = []
    for c in hdr:
        m = _LUN_WAVE.match(c)
        if m and int(m.group(1)) not in waves:
            waves.append(int(m.group(1)))
    dls = any(c in _LUN_DLS for c in hdr)
    if not waves and not dls:
        return None
    meta: Dict[str, str] = {}
    report = False
    for r in rows[:hi]:
        cells = [c.strip() for c in r]
        nz = [c for c in cells if c]
        if nz and nz[0] == "Report" and not meta:
            report = True
        if len(cells) >= 2 and cells[0] and cells[1]:
            meta.setdefault(cells[0], cells[1])
    data = rows[hi + 1:]
    first: Dict[str, str] = {}
    if data:
        for k, v in zip(hdr, data[0]):
            first.setdefault(k, v.strip())

    def pick(meta_key: str, col: str) -> Optional[str]:
        return meta.get(meta_key) or first.get(col) or None

    f: Dict[str, object] = {
        "exportFormat": fmt + (" report" if report else " table"),
        "softwareVersion": meta.get("Software version"),
        "clientVersion": meta.get("Client version"),
        "instrumentId": _id_text(pick("Instrument", "instrument id")),
        "application": pick("Application", "application"),
        "plateType": pick("Plate type", "plate type"),
        "plateCount": _int(meta.get("Nr of Plates")),
    }
    if meta.get("Date"):
        f["runDate"] = _xlsx_date(src, meta["Date"]) if xlsx else _date(meta["Date"])
    else:
        d, t = first.get("date") or "", first.get("time") or ""
        if xlsx and re.match(r"^\d+(\.\d+)?$", d) and re.match(r"^0?\.\d+$|^0$", t):
            f["runDate"] = _xlsx_date(src, d, t)
        else:
            f["runDate"] = _date(" ".join(x for x in (d, t) if x))
    if complete:
        si = hdr.index("sample name")
        f["sampleCount"] = sum(1 for r in data if len(r) > si and r[si].strip())
    if len(waves) <= MAX_LIST:
        f["wavelengthsNm"] = waves
    else:
        f["wavelengthRangeNm"] = [min(waves), max(waves)]
        f["wavelengthCount"] = len(waves)
    modes = []
    if waves:
        modes.append("Absorbance")
    if dls:
        modes.append("Dynamic Light Scattering")
    f["readModes"] = modes
    return _result(LUNATIC, f)


def sniff_lunatic_csv(src):
    text = head_text(src)
    if not text:
        return None
    low = text.lower()
    if "plate id" not in low or "sample name" not in low:
        return None
    complete = _complete(src)
    if not complete:
        text = text[: text.rfind("\n") + 1]
    try:
        rows = [r for r in csv.reader(io.StringIO(text, newline="")) if r]
    except csv.Error:
        return None
    return _lunatic_rows(src, rows, complete, "csv", False)


def sniff_lunatic_xlsx(src):
    names = xlsx_sheet_names(src)
    if not names:
        return None
    rows = xlsx_rows(src, names[0], 80)
    if not rows or not any("sample name" in _norm(c) for r in rows for c in r):
        return None
    rows = xlsx_rows(src, names[0], _XLSX_MAX_ROWS) or rows
    return _lunatic_rows(src, rows, len(rows) < _XLSX_MAX_ROWS, "xlsx", True)


register_instrument_sniffer("csv", sniff_envision)
register_instrument_sniffer("csv", sniff_kaleido)
register_instrument_sniffer("csv", sniff_workbench_csv)
register_instrument_sniffer("csv", sniff_lunatic_csv)
register_instrument_sniffer("txt", sniff_methodical_mind)
register_instrument_sniffer("txt", sniff_immunospot)
register_instrument_sniffer("xlsx", sniff_workbench_xlsx)
register_instrument_sniffer("xlsx", sniff_apex)
register_instrument_sniffer("xlsx", sniff_lunatic_xlsx)
