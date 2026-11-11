"""Cell counters and solution analyzers: Beckman Coulter Vi-CELL XR and
Vi-CELL BLU, ChemoMetec NucleoView and NC-View, Revvity Matrix (Cellaca),
Roche Cedex HiRes and Cedex Bio HT, Beckman Coulter PharmSpec (HIAC) and Nova
Biomedical BioProfile FLEX2.

Signatures were read from allotropy's parsers for the same files:

  Vi-CELL XR     txt: first line "Vi-CELL XR <version>", then "Key : value"
                 lines. xlsx: "Vi-CELL XR <version>" in A1, "Unit S/N:" in A3,
                 a "Sample ID" table from row 4 or 5, or the 2.04 single-sample
                 report (A4 "Sample ID", A5 "File name"). Legacy .xls: the
                 "Vi-CELL XR <version>" string in the OLE2 workbook.
  Vi-CELL BLU    csv whose header starts "Sample ID" with the BLU result
                 columns, or the per-image detail export (settings header,
                 reagent table, "Image#" table per sample).
  NucleoView     NC-200 "Key;value;" csv, or the "Key:<TAB>value" csv.
  NC-View        csv table with NAME, SAMPLE ID, PROTOCOL, INSTRUMENT columns.
  Revvity Matrix csv/xlsx table "Well Name, Live Count, ...", optionally under
                 a "Plate Name / Instrument SN / Version" header block.
  Cedex HiRes    csv/xlsx table with "Cedex ID", "System name", "Sample
                 identifer" (sic) columns.
  Cedex Bio HT   tab-separated txt whose first row is
                 "0 <time> #ARC-FILE# <ver> ... CEDEX BIO HT <serial> <sw> <user>".
  PharmSpec      xlsx "HIAC PharmSpec vX Summary Report" with "Key : value"
                 cells, the HIAC "Run Counter Test" HTML saved as .xls, and the
                 binary .xls of the same report (identity only).
  NovaBio Flex2  csv whose header starts "Date & Time, Sample ID, Sample Type".

Counts (sampleCount, measurementCount) are reported only when the whole file
fit in the leading read, so a large export never reports a partial count.
Sample identifiers are strings exactly as written; from an .xlsx, a number
cell holding an integer is written without the ".0" the cell XML carries.
"""

from __future__ import annotations

import csv
import datetime as _dt
import html
import io
import re
from typing import Dict, List, Optional, Tuple

from . import (
    TEXT_HEAD_BYTES,
    head,
    head_text,
    header,
    instrument,
    open_zip,
    register_instrument_sniffer,
    xlsx_rows,
    zip_member,
)

XR = "Beckman Coulter Vi-Cell XR"
BLU = "Beckman Coulter Vi-Cell BLU"
NUCLEOVIEW = "ChemoMetec Nucleoview"
NCVIEW = "ChemoMetec NC View"
MATRIX = "Revvity Matrix"
HIRES = "Roche Cedex HiRes"
BIOHT = "Roche Cedex BioHT"
PHARMSPEC = "Beckman Coulter PharmSpec"
FLEX2 = "NovaBio Flex2"

instrument(XR, "Cell Counting", "Beckman Coulter", "beckman_vi_cell_xr", ["Brightfield"])
instrument(BLU, "Cell Counting", "Beckman Coulter", "beckman_vi_cell_blu", ["Brightfield"])
instrument(NUCLEOVIEW, "Cell Counting", "ChemoMetec", "chemometec_nucleoview", ["Dark Field"])
instrument(NCVIEW, "Cell Counting", "ChemoMetec", "chemometec_nc_view", ["Fluorescence"])
instrument(MATRIX, "Cell Counting", "Revvity", "revvity_matrix", ["Brightfield"])
instrument(HIRES, "Cell Counting", "Roche", "roche_cedex_hires", ["Brightfield"])
instrument(BIOHT, "Solution Analyzer", "Roche", "roche_cedex_bioht", ["Metabolite Detection"])
instrument(PHARMSPEC, "Solution Analyzer", "Beckman Coulter", "beckman_pharmspec", ["Light Obscuration"])
instrument(FLEX2, "Solution Analyzer", "Nova Biomedical", "novabio_flex2",
           ["Metabolite Detection", "Blood Gas", "pH", "Osmolality", "Cell Counting"])

# At most this many sample identifiers are listed; sampleCount has the total.
MAX_IDS = 20
XLSX_ROWS = 5000

OLE2_MAGIC = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"


# ---------------------------------------------------------------- helpers


def _complete(src) -> bool:
    """True when the leading read holds the whole file."""
    return len(head(src)) < TEXT_HEAD_BYTES


def _text(src) -> Optional[str]:
    t = head_text(src)
    if not t:
        return None
    return t[1:] if t.startswith("﻿") else t


def _csv(text: str, complete: bool, delimiter: str = ",") -> Optional[List[List[str]]]:
    try:
        rows = list(csv.reader(io.StringIO(text, newline=None), delimiter=delimiter))
    except csv.Error:
        return None
    rows = [r for r in rows if r]
    if not complete and rows:
        rows.pop()  # the last row may be cut by the leading read
    return rows


_NEWLINE = re.compile(r"\r\n|\r|\n")


def _lines(text: str) -> List[str]:
    """Lines split on CRLF, CR or LF (only those, as the Go port does)."""
    out = _NEWLINE.split(text)
    if out and out[-1] == "":
        out.pop()
    return out


def _at(row: List[str], i: int) -> str:
    return row[i] if 0 <= i < len(row) else ""


def _ids(values) -> List[str]:
    out: List[str] = []
    seen = set()
    for v in values:
        v = v.strip()
        if v and v not in seen:
            seen.add(v)
            out.append(v)
            if len(out) >= MAX_IDS:
                break
    return out


def _distinct(values) -> List[str]:
    out: List[str] = []
    seen = set()
    for v in values:
        v = v.strip()
        if v and v not in seen:
            seen.add(v)
            out.append(v)
    return out


def _only(values) -> Optional[str]:
    """The single distinct non-empty value, or None when there are several."""
    d = _distinct(values)
    return d[0] if len(d) == 1 else None


_INT_FLOAT = re.compile(r"^-?[0-9]+\.0+$")


def _xnum(s: str) -> str:
    """An .xlsx number cell holding an integer, as the integer ("3714019.0"
    -> "3714019"). Any other text is returned unchanged."""
    s = (s or "").strip()
    return s.split(".", 1)[0] if _INT_FLOAT.match(s) else s


def _iso(y, mo, d, h=0, mi=0, s=0, tz: str = "") -> Optional[str]:
    try:
        t = _dt.datetime(int(y), int(mo), int(d), int(h), int(mi), int(s))
    except (ValueError, OverflowError):
        return None
    return t.strftime("%Y-%m-%dT%H:%M:%S") + tz


_MONTHS = {m: i + 1 for i, m in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"])}
# Vi-CELL XR: "29 Sep 2023  5:21:32 AM"
_XR_DATE = re.compile(r"^([0-9]{1,2}) ([A-Za-z]{3}) ([0-9]{4})\s+([0-9]{1,2}):([0-9]{2}):([0-9]{2}) ?([AP]M)$", re.I)


def _xr_date(s: str) -> Optional[str]:
    m = _XR_DATE.match((s or "").strip())
    if not m:
        return None
    mon = _MONTHS.get(m.group(2).lower())
    hour = int(m.group(4))
    if mon is None or not 1 <= hour <= 12:
        return None
    hour = hour % 12 + (12 if m.group(7).upper() == "PM" else 0)
    return _iso(m.group(3), mon, m.group(1), hour, m.group(5), m.group(6))


def _plural(n: int, word: str) -> str:
    return "%d %s%s" % (n, word, "" if n == 1 else "s")


def _summary(name: str, fmt: str, samples: Optional[int]) -> str:
    s = "%s %s export" % (name, fmt)
    if samples is not None:
        s += ", " + _plural(samples, "sample")
    return s


def _strip_sn(v: Optional[str]) -> Optional[str]:
    """"S/N: 9002..." or "s/n 9000..." -> the serial itself."""
    if v is None:
        return None
    v = re.sub(r"^\s*s/n\s*:?\s*", "", v, flags=re.I).strip()
    return v or None


def _excel_serial(s: str, date1904: bool) -> Optional[str]:
    try:
        f = float(s)
    except ValueError:
        return None
    if not 1 <= f < 2958466:
        return None
    base = _dt.datetime(1904, 1, 1) if date1904 else _dt.datetime(1899, 12, 30)
    t = base + _dt.timedelta(seconds=round(f * 86400))
    return t.strftime("%Y-%m-%dT%H:%M:%S")


def _date1904(src) -> bool:
    zf = getattr(src, "_instrument_zip", None) or open_zip(src)
    data = zip_member(zf, "xl/workbook.xml", 1 << 20) if zf else None
    return bool(data and re.search(rb'date1904="(1|true)"', data))


def _ole2_has(src, *markers: str) -> bool:
    """Whether every marker string appears in the leading bytes of an OLE2
    (legacy .xls) file, as BIFF8 writes strings: 8-bit or UTF-16LE."""
    data = head(src)
    if not data.startswith(OLE2_MAGIC):
        return False
    for m in markers:
        if m.encode("latin-1") not in data and m.encode("utf-16-le") not in data:
            return False
    return True


# ---------------------------------------------------------- Vi-CELL XR

_XR_VERSION = re.compile(r"^Vi-CELL XR ([0-9]+\.[0-9]+(?:\.[0-9]+)*)\s*$")


def _xr_txt(src):
    text = _text(src)
    if not text or not text.lstrip().startswith("Vi-CELL XR"):
        return None
    lines = _lines(text)
    first = next((ln for ln in lines if ln.strip()), "")
    m = _XR_VERSION.match(first.strip())
    if not m:
        return None
    kv: Dict[str, str] = {}
    for ln in lines[1:200]:
        if " :" in ln:
            k, v = ln.split(" :", 1)
            kv.setdefault(k.strip(), v.strip())
    if "Sample ID" not in kv or not ("RunDate" in kv or "Results:" in text):
        return None
    sid = kv.get("Sample ID", "")
    fields = {
        "softwareVersion": m.group(1),
        "serialNumber": kv.get("Unit S/N") or None,
        "analysisVersion": kv.get("Analysis version") or None,
        "sampleCount": 1,
        "sampleIds": [sid] if sid else None,
        "runDate": _xr_date(kv.get("RunDate", "")),
        "exportFormat": "txt",
    }
    return header(XR, _summary(XR, "txt", 1), fields)


def _xr_xlsx(src):
    rows = xlsx_rows(src, None, XLSX_ROWS)
    if not rows or not rows[0]:
        return None
    m = _XR_VERSION.match(rows[0][0].strip())
    if not m:
        return None
    version = m.group(1)

    def cell(r, c):
        return rows[r][c].strip() if r < len(rows) and c < len(rows[r]) else ""

    # The worksheet reader skips rows with no cells, so rows are found by
    # content, not by their sheet row number.
    hdr = next((i for i in range(1, min(len(rows), 8)) if cell(i, 0) == "Sample ID"), None)
    if hdr is None:
        return None
    if cell(hdr + 1, 0) == "File name":
        sid = _xnum(cell(hdr, 2))
        date = next((cell(i, 2) for i in range(hdr, min(len(rows), hdr + 5))
                     if cell(i, 0) == "Sample date"), "")
        fields = {
            "softwareVersion": version,
            "sampleCount": 1,
            "sampleIds": [sid] if sid else None,
            "runDate": _xr_date(date),
            "exportFormat": "xlsx report",
        }
        return header(XR, _summary(XR, "xlsx report", 1), fields)

    date_col = next((i for i, h in enumerate(rows[hdr])
                     if re.match(r"^(Sample dat[ae]/time|Sample date|RunDate)\b", h.strip())), None)
    serial = None
    sn = next((cell(i, 0) for i in range(1, hdr) if cell(i, 0).startswith("Unit S/N")), "")
    if ":" in sn:
        serial = sn.split(":", 1)[1].strip()
    fields = {"softwareVersion": version, "serialNumber": serial or None, "exportFormat": "xlsx"}
    count = None
    if date_col is not None and len(rows) < XLSX_ROWS:
        ids = []
        for r in rows[hdr + 1:]:
            if date_col < len(r) and _xr_date(r[date_col]):
                ids.append(_xnum(r[0]))
        count = len(ids)
        fields["sampleCount"] = count
        fields["sampleIds"] = _ids(ids) or None
    return header(XR, _summary(XR, "xlsx", count), fields)


def _xr_xls(src):
    if not _ole2_has(src, "Vi-CELL XR "):
        return None
    data = head(src)
    version = None
    m = re.search(rb"Vi-CELL XR ([0-9]+\.[0-9]+(?:\.[0-9]+)*)", data)
    if m:
        version = m.group(1).decode()
    else:
        for off in (0, 1):
            u = re.search(r"Vi-CELL XR ([0-9]+\.[0-9]+(?:\.[0-9]+)*)", data[off:].decode("utf-16-le", "ignore"))
            if u:
                version = u.group(1)
                break
    fields = {"softwareVersion": version, "exportFormat": "xls"}
    return header(XR, _summary(XR, "xls", None), fields)


# --------------------------------------------------------- Vi-CELL BLU


def _blu(src):
    text = _text(src)
    if not text or not text.lstrip('"').startswith("Sample ID"):
        return None
    complete = _complete(src)
    rows = _csv(text, complete)
    if not rows or not rows[0] or rows[0][0] != "Sample ID":
        return None
    hdr = rows[0]
    summary_cols = ("Viability (%)" in hdr and "Viable (x10^6) cells/mL" in hdr
                    and any(h.startswith("Average viable diameter (") for h in hdr))
    if summary_cols:
        ids = [r[0] for r in rows[1:] if r and r[0].strip()]
        fmt = "summary csv"
    elif ("Images for analysis" in hdr and "Analysis date/time" in hdr
          and any(r and r[0] == "Image#" for r in rows[2:])):
        ids = [rows[i + 1][0] for i, r in enumerate(rows[:-1]) if r == hdr and rows[i + 1]]
        fmt = "detail csv"
    else:
        return None
    fields: Dict[str, object] = {"exportFormat": fmt}
    count = None
    if complete:
        count = len(ids)
        fields["sampleCount"] = count
    fields["sampleIds"] = _ids(ids) or None
    return header(BLU, _summary(BLU, fmt, count), fields)


# ------------------------------------------------------------ NucleoView

_NV_KEYS = ("Viability (%)", "Total (cells/ml)")
_NV_ANY = ("Live (cells/ml)", "Dead (cells/ml)")


def _nucleoview(src):
    text = _text(src)
    if not text:
        return None
    lines = _lines(text)[:200]
    first = lines[0] if lines else ""
    if ";" in first and "\t" not in first:
        kv: Dict[str, str] = {}
        for ln in lines:
            parts = ln.split(";")
            if len(parts) >= 2 and parts[0].strip():
                kv.setdefault(parts[0].strip(), parts[1].strip())
        if not all(k in kv for k in _NV_KEYS + ("Image",)) or not any(k in kv for k in _NV_ANY):
            return None
        run = None
        m = re.match(r"^([0-9]{4})([0-9]{2})([0-9]{2}) ([0-9]{2})([0-9]{2})([0-9]{2})$", kv.get("Date time", ""))
        if m:
            tz = ""
            off = kv.get("Time zone offset", "")
            if re.match(r"^[+-]?[0-9]{1,2}$", off) and abs(int(off)) <= 14:
                h = int(off)
                tz = "%s%02d:00" % ("-" if h < 0 else "+", abs(h))
            run = _iso(*m.groups(), tz=tz)
        fields = {
            "softwareVersion": kv.get("Application SW version") or None,
            "instrumentModel": kv.get("Instrument type") or None,
            "serialNumber": _strip_sn(kv.get("Instrument s/n")),
            "runDate": run,
            "protocol": kv.get("Protocol template title") or None,
            "exportVersion": kv.get("csv file version") or None,
            "sampleCount": 1,
            "exportFormat": "semicolon csv",
        }
        return header(NUCLEOVIEW, _summary(NUCLEOVIEW, "semicolon csv", 1), fields)
    if ":\t" in first:
        kv = {}
        for ln in lines:
            if ":\t" in ln:
                k, v = ln.split(":\t", 1)
                kv.setdefault(k.strip(), v.strip())
        if not all(k in kv for k in _NV_KEYS + ("Source ID", "Sample ID")) or not any(k in kv for k in _NV_ANY):
            return None
        run = None
        m = re.match(r"^([0-9]{4})-([0-9]{2})-([0-9]{2}) ([0-9]{2}):([0-9]{2}):([0-9]{2})$", kv.get("Time (UTC)", ""))
        if m:
            run = _iso(*m.groups(), tz="+00:00")
        sid = kv.get("Sample ID", "")
        fields = {
            "softwareVersion": kv.get("Software version") or None,
            "serialNumber": _strip_sn(kv.get("Instrument")),
            "runDate": run,
            "protocol": kv.get("Protocol") or None,
            "sampleCount": 1,
            "sampleIds": [sid] if sid else None,
            "exportFormat": "key-value csv",
        }
        return header(NUCLEOVIEW, _summary(NUCLEOVIEW, "key-value csv", 1), fields)
    return None


# --------------------------------------------------------------- NC-View


def _ncview(src):
    text = _text(src)
    if not text or not text.startswith("NAME,"):
        return None
    complete = _complete(src)
    rows = _csv(text, complete)
    if not rows:
        return None
    hdr = [h.strip() for h in rows[0]]
    need = ("NAME", "SAMPLE ID", "PROTOCOL", "INSTRUMENT", "VIABILITY (%)")
    if hdr[:2] != ["NAME", "SAMPLE ID"] or not all(k in hdr for k in need):
        return None
    data = [r for r in rows[1:] if len(r) == len(hdr) and r[0].strip()]
    col = hdr.index
    fields: Dict[str, object] = {
        "serialNumber": _strip_sn(_only(r[col("INSTRUMENT")] for r in data)),
        "protocol": _only(r[col("PROTOCOL")] for r in data),
        "sampleIds": _ids(r[col("SAMPLE ID")] for r in data) or None,
        "exportFormat": "csv",
    }
    count = len(data) if complete else None
    fields["sampleCount"] = count
    return header(NCVIEW, _summary(NCVIEW, "csv", count), fields)


# --------------------------------------------------------- Revvity Matrix

_MATRIX_COLS = ("Well Name", "Live Count", "Dead Count", "Total Count")


def _matrix_rows(rows: List[List[str]], complete: bool, xlsx: bool, src=None):
    rows = [[c.strip() for c in r] for r in rows]
    if not rows or not rows[0]:
        return None
    meta: Dict[str, str] = {}
    start = 0
    if rows[0][0] == "Plate Name":
        start = None
        for i, r in enumerate(rows[:60]):
            if "Row" in r and "Column" in r and "Well Name" in r:
                start = i
                break
            if len(r) >= 2 and r[0]:
                meta.setdefault(r[0], r[1])
        if start is None or "Instrument SN" not in meta:
            return None
    hdr = rows[start]
    if not all(c in hdr for c in _MATRIX_COLS) or not any(h.startswith("Viability") for h in hdr):
        return None
    w = hdr.index("Well Name")
    ids = [r[w] for r in rows[start + 1:] if w < len(r) and r[w]]
    if xlsx:
        ids = [_xnum(i) for i in ids]
    norm = _xnum if xlsx else (lambda s: s)
    run = meta.get("Results Timestamp") or None
    if run and xlsx:
        run = _excel_serial(run, _date1904(src)) or run
    fmt = "xlsx" if xlsx else "csv"
    fields: Dict[str, object] = {
        "serialNumber": norm(meta.get("Instrument SN", "")) or None,
        "exportVersion": norm(meta.get("Version", "")) or None,
        "plateName": norm(meta.get("Plate Name", "")) or None,
        "runDate": run,
        "sampleIds": _ids(ids) or None,
        "exportFormat": fmt,
    }
    count = len(ids) if complete else None
    fields["sampleCount"] = count
    return header(MATRIX, _summary(MATRIX, fmt, count), fields)


def _matrix_csv(src):
    text = _text(src)
    if not text or not (text.startswith("Well Name,") or text.startswith("Plate Name,")):
        return None
    complete = _complete(src)
    rows = _csv(text, complete)
    return _matrix_rows(rows, complete, False) if rows else None


def _matrix_xlsx(src):
    rows = xlsx_rows(src, None, XLSX_ROWS)
    if not rows or not rows[0] or rows[0][0].strip() not in ("Well Name", "Plate Name"):
        return None
    return _matrix_rows(rows, len(rows) < XLSX_ROWS, True, src)


# ----------------------------------------------------------- Cedex HiRes


def _hires_rows(rows: List[List[str]], complete: bool, xlsx: bool):
    if not rows:
        return None
    hdr = [h.strip().replace("identifer", "identifier") for h in rows[0]]
    if not all(c in hdr for c in ("Cedex ID", "System name", "Sample identifier", "Viable Cell Conc.")):
        return None
    norm = _xnum if xlsx else (lambda s: s.strip())
    data = [r for r in rows[1:] if len(r) > hdr.index("Sample identifier") and any(c.strip() for c in r)]
    col = hdr.index
    fmt = "xlsx" if xlsx else "csv"
    fields: Dict[str, object] = {
        "deviceIdentifier": _only(norm(_at(r, col("Cedex ID"))) for r in data),
        "systemName": _only(_at(r, col("System name")) for r in data),
        "sampleIds": _ids(norm(_at(r, col("Sample identifier"))) for r in data) or None,
        "exportFormat": fmt,
    }
    count = len(data) if complete else None
    fields["sampleCount"] = count
    return header(HIRES, _summary(HIRES, fmt, count), fields)


def _hires_csv(src):
    text = _text(src)
    if not text or not text.startswith("Cedex ID,"):
        return None
    complete = _complete(src)
    rows = _csv(text, complete)
    return _hires_rows(rows, complete, False) if rows else None


def _hires_xlsx(src):
    rows = xlsx_rows(src, None, XLSX_ROWS)
    if not rows or not rows[0] or rows[0][0].strip() != "Cedex ID":
        return None
    return _hires_rows(rows, len(rows) < XLSX_ROWS, True)


# ---------------------------------------------------------- Cedex Bio HT


def _bioht(src):
    text = _text(src)
    if not text or not text.startswith("0\t"):
        return None
    lines = _lines(text)
    first = lines[0].split("\t")
    if len(first) < 10 or first[2] != "#ARC-FILE#" or "CEDEX BIO HT" not in first[6].upper():
        return None
    complete = _complete(src)
    if not complete and lines:
        lines.pop()
    software = first[8].strip()
    sid_col, analyte_col = {"5": (5, 2), "6": (2, 6), "7": (2, 6)}.get(software[:1], (None, None))
    data = [ln.split("\t") for ln in lines[1:] if ln.startswith("40\t")]
    m = re.match(r"^([0-9]{4})-([0-9]{2})-([0-9]{2}) ([0-9]{2}):([0-9]{2}):([0-9]{2})$", first[1].strip())
    fields: Dict[str, object] = {
        "instrumentModel": first[6].strip() or None,
        "serialNumber": first[7].strip() or None,
        "softwareVersion": software or None,
        "exportVersion": first[3].strip() or None,
        "runDate": _iso(*m.groups()) if m else None,
        "exportFormat": "txt",
    }
    count = None
    if sid_col is not None:
        rows = [r for r in data if len(r) > max(sid_col, analyte_col)]
        fields["sampleIds"] = _ids(r[sid_col] for r in rows) or None
        fields["analytes"] = _distinct(r[analyte_col] for r in rows) or None
        if complete:
            count = len(_distinct(r[sid_col] for r in rows))
            fields["sampleCount"] = count
    if complete:
        fields["measurementCount"] = len(data)
    return header(BIOHT, _summary(BIOHT, "txt", count), fields)


# ------------------------------------------------------------- PharmSpec


def _kv_cells(rows) -> Dict[str, str]:
    """Header pairs laid out as key, ":", value (PharmSpec) or key, ": value"
    (HIAC Run Counter), several pairs per row."""
    kv: Dict[str, str] = {}
    for r in rows:
        r = [c.strip() for c in r]
        i = 0
        while i < len(r) - 1:
            k, nxt = r[i], r[i + 1]
            if k and not k.startswith(":") and nxt == ":" and i + 2 < len(r):
                kv.setdefault(k, r[i + 2])
                i += 3
            elif k and not k.startswith(":") and nxt.startswith(":") and nxt != ":":
                kv.setdefault(k, nxt[1:].strip())
                i += 2
            else:
                i += 1
    return kv


def _pharmspec_fields(kv: Dict[str, str], fmt: str, version: Optional[str], xlsx: bool):
    norm = _xnum if xlsx else (lambda s: s)
    sid = norm(kv.get("Probe") or kv.get("Sample Name") or "")
    runs = kv.get("No Of Runs", "")
    fields = {
        "softwareVersion": version,
        "sensorModel": kv.get("Sensor Model") or None,
        "sensorSerialNumber": norm(kv.get("Sensor Serial Number", "")) or None,
        "sampleIds": [sid] if sid else None,
        "runDate": kv.get("Sample Date") or None,
        "runCount": int(_xnum(runs)) if re.match(r"^[0-9]{1,4}(\.0+)?$", runs) else None,
        "exportFormat": fmt,
    }
    return header(PHARMSPEC, _summary(PHARMSPEC, fmt, None), fields)


def _pharmspec_xlsx(src):
    rows = xlsx_rows(src, None, 60)
    if not rows:
        return None
    title = next((c for r in rows[:3] for c in r if "PharmSpec" in c or "Run Counter Test" in c), None)
    if title is None:
        return None
    if not any("Particle Size" in c for r in rows for c in r):
        return None
    kv = _kv_cells(rows)
    if "Sensor Serial Number" not in kv:
        return None
    m = re.search(r"PharmSpec v([0-9]+(?:\.[0-9]+)*)", title)
    return _pharmspec_fields(kv, "xlsx", m.group(1) if m else None, True)


_TD = re.compile(r"<td\b[^>]*>(.*?)</td\s*>", re.I | re.S)
_TAG = re.compile(r"<[^>]+>")


def _pharmspec_html(src):
    text = _text(src)
    if not text or not text.lstrip()[:5].lower() == "<html":
        return None
    if "Sensor Serial Number" not in text or "Particle Size" not in text or "Run No." not in text:
        return None
    cells = []
    for m in _TD.finditer(text):
        c = html.unescape(_TAG.sub(" ", m.group(1)))
        c = " ".join(c.replace("\xa0", " ").split())
        if c:
            cells.append(c)
        if "Run No." in c:
            break
    kv = _kv_cells([cells])
    if "Sensor Serial Number" not in kv:
        return None
    return _pharmspec_fields(kv, "html xls", None, False)


def _pharmspec_xls(src):
    if not _ole2_has(src, "Sensor Serial Number", "Particle Size", "Run No."):
        return None
    return header(PHARMSPEC, _summary(PHARMSPEC, "xls", None), {"exportFormat": "xls"})


# ------------------------------------------------------------- NovaBio Flex2

_FLEX2_ANALYTES = ("Gln", "Glu", "Gluc", "Lac", "NH4+", "Na+", "K+", "Ca++", "HCO3",
                   "pH", "PO2", "PCO2", "Osm")
_FLEX2_MARKERS = ("Vessel ID", "Batch ID", "Cell Type", "Total Density", "Viable Density",
                  "Chemistry Dilution Ratio", "Vessel Temperature (°C)")


def _flex2(src):
    text = _text(src)
    if not text or not text.lstrip('"').startswith("Date & Time"):
        return None
    complete = _complete(src)
    rows = _csv(text, complete)
    if not rows:
        return None
    hdr = [h.strip() for h in rows[0]]
    if hdr[:3] != ["Date & Time", "Sample ID", "Sample Type"]:
        return None
    if not any(m in hdr for m in _FLEX2_MARKERS) or not any(a in hdr for a in _FLEX2_ANALYTES):
        return None
    data = [r for r in rows[1:] if len(r) >= 3 and r[0].strip() and r[1].strip()]
    fields: Dict[str, object] = {
        "analytes": [h for h in hdr if h in _FLEX2_ANALYTES],
        "sampleIds": _ids(r[1] for r in data) or None,
        "exportFormat": "csv",
    }
    count = len(data) if complete else None
    fields["sampleCount"] = count
    return header(FLEX2, _summary(FLEX2, "csv", count), fields)


# ------------------------------------------------------------ registration

register_instrument_sniffer("txt", _xr_txt)
register_instrument_sniffer("txt", _bioht)
register_instrument_sniffer("csv", _blu)
register_instrument_sniffer("csv", _nucleoview)
register_instrument_sniffer("csv", _ncview)
register_instrument_sniffer("csv", _matrix_csv)
register_instrument_sniffer("csv", _hires_csv)
register_instrument_sniffer("csv", _flex2)
register_instrument_sniffer("xlsx", _xr_xlsx)
register_instrument_sniffer("xlsx", _matrix_xlsx)
register_instrument_sniffer("xlsx", _hires_xlsx)
register_instrument_sniffer("xlsx", _pharmspec_xlsx)
register_instrument_sniffer("xls", _xr_xls)
register_instrument_sniffer("xls", _pharmspec_html)
register_instrument_sniffer("xls", _pharmspec_xls)

SNIFFERS: Tuple = (_xr_txt, _bioht, _blu, _nucleoview, _ncview, _matrix_csv, _hires_csv, _flex2,
                   _xr_xlsx, _matrix_xlsx, _hires_xlsx, _pharmspec_xlsx,
                   _xr_xls, _pharmspec_html, _pharmspec_xls)
