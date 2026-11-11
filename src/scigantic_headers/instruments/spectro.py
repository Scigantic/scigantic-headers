"""Spectrophotometer exports: Thermo Fisher VISIONlite, NanoDrop 8000, NanoDrop
Eight, NanoDrop One, Qubit 4, Qubit Flex, Genesys 30 and Genesys On-Board.

Markers, from the files allotropy's parsers read:

  VISIONlite        fixed-wavelength CSV with "<nm> nm   [A]" columns, quant CSV
                    with "Ordinate [A]" and "Concentration [unit]", or a scan /
                    rate CSV whose first line is file,sample,analyst,time,...
                    and whose second line is "nm,A" / "min,A"
  VISIONlite or     UTF-16 tab-separated fixed export with "Sample name" or
  Genesys On-Board  "Well Position" plus "Absorbance"/"A<nm>" and "Result"
                    columns. Both programs write this layout and allotropy's
                    two sniffers accept the same bytes, so the result names
                    both candidates instead of picking one.
  NanoDrop 8000     tab-separated report whose columns start Plate ID, Well,
                    Sample ID, User ID, Date, Time
  NanoDrop Eight    "Application:" / "Serial Number:" / "User Name:" header
                    lines over a tab-separated table
  NanoDrop One      comma CSV or .xlsx starting Date, Sample Name and a
                    "Nucleic Acid(...)" / "Protein A280" / "Protein & Label"
                    column, with no header block
  Qubit 4           Run ID, Assay Name, Test Name, Test Date, ... Green RFU /
                    Far Red RFU (CSV or .xlsx)
  Qubit Flex        Run ID, Test Date, Assay Name, Sample Name, ... Sample RFU,
                    Well (CSV, CSV with a "sep=," line, or .xlsx)
  Genesys 30        "Scan" / "Mode" / "Lower" / "Upper" / "Length" key-value
                    block over a "wavelength(nm)" table (CSV or TSV)

Dates are reported as the file writes them: these exports mix day-first and
month-first dates and do not say which. Qubit 4's Run ID is an ISO date and is
reported as runDate. Date cells in an .xlsx are Excel serials that were
sometimes parsed from text with the wrong day/month order (Qubit Flex and
Qubit 4 test data disagree with their own Run IDs), so .xlsx test dates are
not reported.
"""

from __future__ import annotations

import datetime as _dt
import re
from typing import Dict, List, Optional

from . import (
    INSTRUMENTS,
    DecodedHeader,
    Source,
    head,
    head_text,
    header,
    instrument,
    register_instrument_sniffer,
    xlsx_rows,
    xlsx_sheet_names,
    zip_member,
)

CATEGORY = "Spectrophotometry"
VENDOR = "Thermo Fisher Scientific"

VISIONLITE = "Thermo Fisher Scientific VISIONlite"
ND8000 = "Thermo Fisher Scientific NanoDrop 8000"
ND_EIGHT = "Thermo Fisher Scientific NanoDrop Eight"
ND_ONE = "Thermo Fisher Scientific NanoDrop One"
QUBIT4 = "Thermo Fisher Scientific Qubit 4"
QUBIT_FLEX = "Thermo Fisher Scientific Qubit Flex"
GENESYS30 = "Thermo Fisher Scientific Genesys30"
GENESYS_OB = "Thermo Fisher Scientific Genesys On-Board"

instrument(VISIONLITE, CATEGORY, VENDOR, "thermo_fisher_visionlite", ["Absorbance"])
instrument(ND8000, CATEGORY, VENDOR, "thermo_fisher_nanodrop_8000", ["Absorbance"])
instrument(ND_EIGHT, CATEGORY, VENDOR, "thermo_fisher_nanodrop_eight", ["Absorbance"])
instrument(ND_ONE, CATEGORY, VENDOR, "thermo_fisher_nanodrop_one", ["Absorbance"])
instrument(QUBIT4, CATEGORY, VENDOR, "thermo_fisher_qubit4", ["Fluorescence"])
instrument(QUBIT_FLEX, CATEGORY, VENDOR, "thermo_fisher_qubit_flex", ["Fluorescence"])
instrument(GENESYS30, CATEGORY, VENDOR, "thermo_fisher_genesys30", ["Absorbance"])
instrument(GENESYS_OB, CATEGORY, VENDOR, "thermo_fisher_genesys_on_board", ["Absorbance"])

XLSX_MAX_ROWS = 20000


# ---------------------------------------------------------------- helpers


def _total(fn):
    def wrapped(src):
        try:
            return fn(src)
        except Exception:
            return None

    wrapped.__name__ = fn.__name__
    wrapped.__doc__ = fn.__doc__
    return wrapped


def _text(src: Source):
    """(lines, complete) of a text export. When the leading read stops short
    of the end of the file, the last (possibly cut) line is dropped and
    complete is False, so counts are not reported."""
    t = head_text(src)
    if not t:
        return None, False
    data = head(src)
    complete = src.size is not None and len(data) >= src.size
    lines = t.splitlines()
    if not complete and lines:
        lines = lines[:-1]
    return lines, complete


def _num(s: str) -> Optional[float]:
    s = (s or "").strip().replace(",", ".")
    if len(s) > 40 or not re.fullmatch(r"-?[0-9]+(\.[0-9]+)?", s):
        return None
    v = float(s)
    return v if abs(v) < 1e12 else None


def _intish(v: Optional[float]):
    if v is None:
        return None
    return int(v) if abs(v) < 1e15 and v == int(v) else v


def _a_wavelengths(cols: List[str], exact: bool = True) -> List[int]:
    """Wavelengths named by absorbance columns: A260, A280 (and, when not
    exact, "A280 10mm"). Ratio columns (A260/A280) are not absorbances."""
    pat = r"A([0-9]{3})" if exact else r"A([0-9]{3})(?: .*)?"
    out = []
    for c in cols:
        m = re.fullmatch(pat, c.strip())
        if m:
            w = int(m.group(1))
            if w not in out:
                out.append(w)
    return sorted(out)


def _block(lines: List[str], start: int, complete: bool):
    """(rows, counted) of the table starting at lines[start]: the lines up to
    the first blank one (several exports append a second table after a blank
    line). counted is True when the block's end is inside the leading read."""
    rows = []
    for ln in lines[start:]:
        if not ln.strip():
            return rows, True
        rows.append(ln)
    return rows, complete


def _same(values: List[str]) -> Optional[str]:
    """The value every row agrees on, or None."""
    vs = [v.strip() for v in values if v is not None and v.strip()]
    if not vs or any(v != vs[0] for v in vs):
        return None
    return vs[0]


def _col(rows: List[List[str]], cols: List[str], name: str) -> List[str]:
    if name not in cols:
        return []
    i = cols.index(name)
    return [r[i] if i < len(r) else "" for r in rows]


def _excel_date(src: Source, serial: str) -> Optional[str]:
    """An Excel date serial as an ISO date, honoring the 1904 date system."""
    v = _num(serial)
    if v is None or not (1 <= v < 2958466):
        return None
    zf = getattr(src, "_instrument_zip", None)
    wb = zip_member(zf, "xl/workbook.xml", 256 * 1024) if zf else None
    base = _dt.date(1899, 12, 30)
    if wb and re.search(rb'date1904="(1|true)"', wb):
        base = _dt.date(1904, 1, 1)
    return (base + _dt.timedelta(days=int(v))).isoformat()


def _summary(name: str, parts: List[str]) -> str:
    short = name.replace("Thermo Fisher Scientific ", "")
    return " ".join([short] + [p for p in parts if p])


def _plural(n, word):
    if n is None:
        return ""
    return "%d %s%s" % (n, word, "" if n == 1 else "s")


# ---------------------------------------------------------------- VISIONlite / Genesys On-Board

_NM_COL = re.compile(r"([0-9]{1,5}) nm   \[A\]")
_VL_DATE = re.compile(r"[0-9]{1,2}/[0-9]{1,2}/[0-9]{4} [0-9]{1,2}:[0-9]{2}(:[0-9]{2})?( [AP]M)?")


def _fixed_wavelengths(cols: List[str]) -> Optional[List[int]]:
    """Wavelengths of a fixed export's absorbance columns, or None when a
    column is not one the fixed layout writes."""
    out = []
    n_abs = 0
    for c in cols[1:]:
        c = c.strip()
        m = _NM_COL.fullmatch(c) or re.fullmatch(r"A([0-9]{3,4})", c)
        if m:
            out.append(int(m.group(1)))
            n_abs += 1
        elif c == "Absorbance":
            n_abs += 1
        elif c != "Result":
            return None
    return out if n_abs else None


@_total
def sniff_visionlite(src: Source) -> Optional[DecodedHeader]:
    lines, complete = _text(src)
    if not lines or len(lines) < 2:
        return None
    first = lines[0]
    raw = head(src, 4)

    # Scan / rate: file,sample,analyst,time,... then nm,A or min,A
    second = lines[1].strip()
    if second in ("nm,A", "min,A"):
        hdr = first.split(",")
        if len(hdr) < 4 or not hdr[0].lower().endswith(".csv") or not _VL_DATE.fullmatch(hdr[3].strip()):
            return None
        data, complete = _block(lines, 2, complete)
        if not all(len(r.split(",")) == 2 and _num(r.split(",")[0]) is not None for r in data):
            return None
        scan = second == "nm,A"
        fields: Dict[str, object] = {
            "exportFormat": "scan csv" if scan else "rate csv",
            "sampleCount": 1,
            "runDate": hdr[3].strip(),
        }
        if complete:
            fields["measurementCount"] = len(data)
            if scan and data:
                fields["wavelengthRangeNm"] = [
                    _intish(_num(data[0].split(",")[0])),
                    _intish(_num(data[-1].split(",")[0])),
                ]
        if not scan:
            # allotropy rejects kinetic (rate) exports.
            fields["allotropyParser"] = None
            fields["allotropySupported"] = False
        kind = "scan" if scan else "rate (kinetic)"
        return header(VISIONLITE, _summary(VISIONLITE, [kind, "export"]), fields)

    # UTF-16 tab-separated fixed export: VISIONlite or Genesys On-Board.
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff") and "\t" in first:
        cols = [c.strip() for c in first.split("\t")]
        if cols[0].lower() not in ("sample name", "well position"):
            return None
        waves = _fixed_wavelengths(cols)
        if waves is None:
            return None
        data, complete = _block(lines, 1, complete)
        if not all(len(r.split("\t")) == len(cols) for r in data):
            return None
        return _ambiguous_fixed(cols, waves, data, complete)

    # Fixed (comma, "<nm> nm   [A]" columns) and quant.
    cols = first.split(",")
    if cols[0] != "Sample Name" or len(cols) < 2:
        return None
    data, complete = _block(lines, 1, complete)
    if len(cols) >= 4 and cols[1:3] == ["Dilution factor", "Ordinate [A]"]:
        m = re.fullmatch(r"Concentration \[(.+)\]", cols[3])
        if not m:
            return None
        fields = {
            "exportFormat": "quant csv",
            "concentrationUnit": m.group(1),
            "sampleCount": len(data) if complete else None,
        }
        n = fields["sampleCount"]
        return header(VISIONLITE, _summary(VISIONLITE, ["quant export", _plural(n, "sample")]), fields)
    if not any(_NM_COL.fullmatch(c.strip()) for c in cols[1:]):
        return None
    waves = _fixed_wavelengths(cols)
    if waves is None:
        return None
    n = len(data) if complete else None
    fields = {"exportFormat": "fixed csv", "wavelengthsNm": waves, "sampleCount": n}
    return header(
        VISIONLITE,
        _summary(VISIONLITE, ["fixed export", _plural(n, "sample"), "at " + ", ".join("%d nm" % w for w in waves)]),
        fields,
    )


def _ambiguous_fixed(cols, waves, data, complete) -> DecodedHeader:
    names = [VISIONLITE, GENESYS_OB]
    n = len(data) if complete else None
    fields: Dict[str, object] = {
        "instrumentCandidates": names,
        "instrumentCategory": CATEGORY,
        "vendor": VENDOR,
        "detectionModes": ["Absorbance"],
        "allotropyParserCandidates": [INSTRUMENTS[x][2] for x in names],
        "exportFormat": "fixed tab-separated UTF-16",
    }
    if waves:
        fields["wavelengthsNm"] = waves
    if n is not None:
        fields["sampleCount"] = n
    summary = "VISIONlite or Genesys On-Board fixed export"
    if n is not None:
        summary += " " + _plural(n, "sample")
    summary += " (the layout does not say which program wrote it)"
    return DecodedHeader(format="Thermo Fisher Scientific VISIONlite or Genesys On-Board", summary=summary, fields=fields)


# ---------------------------------------------------------------- NanoDrop 8000

_ND8000_LEAD = ["Plate ID", "Well", "Sample ID", "User ID", "Date", "Time"]


@_total
def sniff_nanodrop_8000(src: Source) -> Optional[DecodedHeader]:
    lines, complete = _text(src)
    if not lines:
        return None
    cols = [c.strip() for c in lines[0].split("\t")]
    if cols[:6] != _ND8000_LEAD or len(cols) < 8:
        return None
    block, complete = _block(lines, 1, complete)
    rows = [r.split("\t") for r in block]
    fields: Dict[str, object] = {"exportFormat": "report txt"}
    if "NA Type" in cols:
        fields["application"] = _same(_col(rows, cols, "NA Type"))
    elif "A280 10mm" in cols and "E 1%" in cols:
        fields["application"] = "Protein A280"
    if "nm 1" in cols and rows:
        waves = []
        for k in range(1, 10):
            for v in _col(rows[:1], cols, "nm %d" % k):
                w = _intish(_num(v))
                if w is not None and w not in waves:
                    waves.append(w)
        fields["wavelengthsNm"] = sorted(waves)
    else:
        fields["wavelengthsNm"] = _a_wavelengths(cols, exact=False)
    if not fields["wavelengthsNm"]:
        fields["wavelengthsNm"] = None
    if rows:
        r = rows[0]
        d, t = (r[4].strip() if len(r) > 4 else ""), (r[5].strip() if len(r) > 5 else "")
        if d:
            fields["runDate"] = (d + " " + t).strip()
    n = len(rows) if complete else None
    fields["sampleCount"] = n
    return header(ND8000, _summary(ND8000, [fields.get("application") or "", "report", _plural(n, "sample")]), fields)


# ---------------------------------------------------------------- NanoDrop Eight

_KV = re.compile(r"[^\t]+:[\t]*[^\t]+\s*")


@_total
def sniff_nanodrop_eight(src: Source) -> Optional[DecodedHeader]:
    lines, complete = _text(src)
    if not lines:
        return None
    kv: Dict[str, str] = {}
    i = 0
    while i < len(lines) and i < 10 and _KV.fullmatch(lines[i]):
        key, _, value = lines[i].partition(":")
        kv[key.strip().lower()] = value.strip()
        i += 1
    if "application" not in kv or not ({"serial number", "user name"} & set(kv)):
        return None
    if i >= len(lines):
        return None
    cols = [c.strip() for c in lines[i].split("\t")]
    if len(cols) < 4:
        return None
    block, complete = _block(lines, i + 1, complete)
    rows = [r.split("\t") for r in block]
    fields: Dict[str, object] = {
        "exportFormat": "tsv",
        "application": kv.get("application") or None,
        "serialNumber": kv.get("serial number") or None,
    }
    fields["wavelengthsNm"] = _a_wavelengths(cols) or None
    spectrum = [_num(c) for c in cols]
    spectrum = [v for v in spectrum if v is not None]
    if len(spectrum) >= 2:
        fields["spectrumRangeNm"] = [_intish(spectrum[0]), _intish(spectrum[-1])]
    for dcol in ("Date & Time", "Date"):
        if dcol in cols and rows:
            v = _col(rows[:1], cols, dcol)[0].strip()
            if v:
                fields["runDate"] = v
            break
    n = len(rows) if complete else None
    fields["sampleCount"] = n
    return header(ND_EIGHT, _summary(ND_EIGHT, [fields.get("application") or "", "export", _plural(n, "sample")]), fields)


# ---------------------------------------------------------------- NanoDrop One

_ND1_PREFIXES = ("Nucleic Acid", "Protein A280", "Protein & Label")
_ND1_SHEET_APPS = ("Nucleic Acid", "Protein A280", "Protein & Label", "dsDNA", "ssDNA", "RNA")


def _nd1_cols_ok(cols: List[str]) -> Optional[str]:
    if cols[:2] != ["Date", "Sample Name"]:
        return None
    for c in cols[2:]:
        for p in _ND1_PREFIXES:
            if c.startswith(p):
                return p
    return None


def _nd1_fields(cols, rows, n, fmt, application):
    fields: Dict[str, object] = {
        "exportFormat": fmt,
        "application": application,
        "wavelengthsNm": _a_wavelengths(cols) or None,
        "sampleCount": n,
    }
    if rows:
        d = rows[0][0].strip() if rows[0] else ""
        if d and _num(d) is None:
            fields["runDate"] = d
    return header(ND_ONE, _summary(ND_ONE, [application or "", fmt.split()[0], "export", _plural(n, "sample")]), fields)


@_total
def sniff_nanodrop_one_csv(src: Source) -> Optional[DecodedHeader]:
    lines, complete = _text(src)
    if not lines:
        return None
    cols = [c.strip() for c in lines[0].split(",")]
    family = _nd1_cols_ok(cols)
    if not family:
        return None
    block, complete = _block(lines, 1, complete)
    rows = [r.split(",") for r in block]
    return _nd1_fields(cols, rows, len(rows) if complete else None, "csv", family)


def _xlsx_table(src: Source, sheet=None):
    rows = xlsx_rows(src, sheet, XLSX_MAX_ROWS)
    if not rows:
        return None, None, None
    cols = [c.strip() for c in rows[0]]
    while cols and not cols[-1]:
        cols.pop()
    data = [r for r in rows[1:] if any(x.strip() for x in r)]
    n = len(data) if len(rows) < XLSX_MAX_ROWS else None
    return cols, data, n


@_total
def sniff_nanodrop_one_xlsx(src: Source) -> Optional[DecodedHeader]:
    names = xlsx_sheet_names(src)
    if not names:
        return None
    cols, data, n = _xlsx_table(src)
    if not cols:
        return None
    family = _nd1_cols_ok(cols)
    if not family:
        return None
    first = names[0].split(" ", 1)[0]
    app = first if first in _ND1_SHEET_APPS else None
    if app is None:
        for p in _ND1_SHEET_APPS:
            if names[0].startswith(p):
                app = p
                break
    return _nd1_fields(cols, data, n, "xlsx", app or family)


# ---------------------------------------------------------------- Qubit 4 / Qubit Flex

_Q4_LEAD = ["Run ID", "Assay Name", "Test Name", "Test Date"]
_QFLEX_LEAD = ["Run ID", "Test Date", "Assay Name", "Sample Name"]


def _qubit_kind(cols: List[str]) -> Optional[str]:
    if cols[:4] == _Q4_LEAD and ("Green RFU" in cols or "Far Red RFU" in cols) and any(
        c.startswith("Qubit") and "tube conc" in c for c in cols
    ):
        return QUBIT4
    if cols[:4] == _QFLEX_LEAD and "Sample RFU" in cols and "Well" in cols:
        return QUBIT_FLEX
    return None


def _qubit(src, name, cols, rows, n, fmt) -> DecodedHeader:
    fields: Dict[str, object] = {"exportFormat": fmt, "sampleCount": n}
    fields["application"] = _same(_col(rows, cols, "Assay Name"))
    fields["excitation"] = _same(_col(rows, cols, "Excitation"))
    if name == QUBIT4:
        fields["emission"] = _same(_col(rows, cols, "Emission"))
        run = _same(_col(rows, cols, "Run ID"))
        if run:
            if fmt == "csv" and re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}", run):
                fields["runDate"] = run
            elif fmt == "xlsx":
                fields["runDate"] = _excel_date(src, run)
    else:
        fields["softwareVersion"] = _same(_col(rows, cols, "Software Version"))
        if fmt == "csv" and rows:
            d = _col(rows[:1], cols, "Test Date")[0].strip()
            if d:
                fields["runDate"] = d
    return header(name, _summary(name, [fields.get("application") or "", fmt, "export", _plural(n, "sample")]), fields)


@_total
def sniff_qubit_csv(src: Source) -> Optional[DecodedHeader]:
    lines, complete = _text(src)
    if not lines:
        return None
    start = 1 if lines[0].strip().startswith("sep=") else 0
    if start >= len(lines):
        return None
    cols = [c.strip() for c in lines[start].split(",")]
    name = _qubit_kind(cols)
    if not name:
        return None
    block, complete = _block(lines, start + 1, complete)
    rows = [r.split(",") for r in block]
    rows = [r for r in rows if r and r[0].strip()]
    return _qubit(src, name, cols, rows, len(rows) if complete else None, "csv")


@_total
def sniff_qubit_xlsx(src: Source) -> Optional[DecodedHeader]:
    cols, data, n = _xlsx_table(src)
    if not cols:
        return None
    name = _qubit_kind(cols)
    if not name:
        return None
    data = [r for r in data if r and r[0].strip()]
    if n is not None:
        n = len(data)
    return _qubit(src, name, cols, data, n, "xlsx")


# ---------------------------------------------------------------- Genesys 30


@_total
def sniff_genesys30(src: Source) -> Optional[DecodedHeader]:
    lines, complete = _text(src)
    if not lines:
        return None
    sep = "\t" if "\t" in lines[0] else ","
    if lines[0].split(sep)[0].strip() != "Scan":
        return None
    kv: Dict[str, str] = {}
    i = 0
    while i < len(lines) and i < 20:
        parts = lines[i].split(sep)
        if not lines[i].strip() or not parts[0].strip():
            break
        kv[parts[0].strip()] = parts[1].strip() if len(parts) > 1 else ""
        i += 1
    if not {"Scan", "Mode", "Lower", "Upper", "Length"} <= set(kv):
        return None
    i += 1
    if i >= len(lines) or not lines[i].strip().lower().startswith("wavelength(nm)"):
        return None
    lo, hi, length = _num(kv["Lower"]), _num(kv["Upper"]), _num(kv["Length"])
    fields: Dict[str, object] = {
        "exportFormat": "scan " + ("tsv" if sep == "\t" else "csv"),
        "measurementMode": kv["Mode"] or None,
        "sampleCount": 1,
    }
    if lo is not None and hi is not None:
        fields["wavelengthRangeNm"] = [_intish(lo), _intish(hi)]
    if length is not None and 0 <= length < 1e9 and length == int(length):
        fields["measurementCount"] = int(length)
    parts = [kv["Mode"], "scan"]
    if "wavelengthRangeNm" in fields:
        parts.append("%s-%s nm" % tuple(fields["wavelengthRangeNm"]))
    return header(GENESYS30, _summary(GENESYS30, parts), fields)


register_instrument_sniffer("csv", sniff_visionlite)
register_instrument_sniffer("csv", sniff_nanodrop_one_csv)
register_instrument_sniffer("csv", sniff_qubit_csv)
register_instrument_sniffer(["csv", "tsv"], sniff_genesys30)
register_instrument_sniffer(["txt", "tsv"], sniff_nanodrop_8000)
register_instrument_sniffer(["txt", "tsv"], sniff_nanodrop_eight)
register_instrument_sniffer("xlsx", sniff_nanodrop_one_xlsx)
register_instrument_sniffer("xlsx", sniff_qubit_xlsx)

SNIFFERS = [
    sniff_visionlite,
    sniff_nanodrop_8000,
    sniff_nanodrop_eight,
    sniff_nanodrop_one_csv,
    sniff_nanodrop_one_xlsx,
    sniff_qubit_csv,
    sniff_qubit_xlsx,
    sniff_genesys30,
]
