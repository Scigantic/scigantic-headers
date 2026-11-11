"""qPCR and dPCR exports: Applied Biosystems QuantStudio (Real-Time PCR
Software and Design & Analysis), Applied Biosystems Absolute Q, QIAGEN
Qiacuity, Thermo Fisher Diomni and Bio-Rad CFX Maestro.

Formats, as allotropy's parsers read them:

* QuantStudio Real-Time PCR Software (QuantStudio 6/7/12K Flex, ViiA 7):
  a tab-separated .txt whose header lines read `* Key = Value`, followed by
  `[Sample Setup]`, `[Amplification Data]`, `[Results]`, ... sections; or an
  .xlsx with the same key/value block at the top of the first sheet
  (usually "Sample Setup").
* QuantStudio Design & Analysis (QuantStudio 1/3/5/6 Pro/7 Pro): an .xlsx
  whose "Results" sheet opens with a key/value block that names the software
  ("Software Name and Version = Design & Analysis Software v2.7.0") and the
  analysis plugins. Older desktop releases (v1.x) write the Real-Time PCR
  layout above; those are recognized by the "QuantStudio(TM) N System"
  instrument type written with a trademark sign, which the Real-Time PCR
  Software does not use.
* Thermo Fisher Diomni: an .xlsx with the Design & Analysis style key/value
  block and "Software name = Diomni", on a "Target_Call" (v4.2) or
  "Target Call" (v4.3) sheet, sometimes beside a "Primary_result" sheet.
* Absolute Q: CSV exports with Run, Date, Instrument, Plate columns (group
  results, per-partition fluorescence, multichannel counts), or the summary
  CSV whose first row names the dyes and second row starts "Run name"; also
  a .zip holding a `*_summary.csv`.
* Qiacuity: CSV that opens with an Excel `sep=,` line and has
  "Sample/NTC/Control" and "Partitions (valid)" columns.
* CFX Maestro: the Quantification Cq Results CSV, whose header starts
  ",Well,Fluor,Target,Content,Sample" and has a Cq column.

Run times are reported as written in the file (`runDate`); the files mix
12- and 24-hour clocks with AM/PM ("16:28:28 PM EDT") and use time zone
abbreviations, so they are not normalized.
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
    open_zip,
    register_instrument_sniffer,
    text_of,
    xlsx_rows,
    xlsx_sheet_names,
)

QS_RTPCR = "AppBio QuantStudio RT-PCR"
QS_DA = "AppBio QuantStudio Design & Analysis"
ABSOLUTE_Q = "AppBio AbsoluteQ"
QIACUITY = "Qiacuity dPCR"
DIOMNI = "Thermo Fisher Diomni"
CFX = "Bio-Rad CFX Maestro"

instrument(QS_DA, "qPCR", "Applied Biosystems", "appbio_quantstudio_designandanalysis", ["Fluorescence"])
instrument(QS_RTPCR, "qPCR", "Applied Biosystems", "appbio_quantstudio", ["Fluorescence"])
instrument(ABSOLUTE_Q, "dPCR", "Applied Biosystems", "appbio_absolute_q", ["Fluorescence"])
instrument(QIACUITY, "dPCR", "QIAGEN", "qiacuity_dpcr", ["Fluorescence"])
instrument(DIOMNI, "qPCR", "Thermo Fisher Scientific", "thermo_fisher_diomni", ["Fluorescence"])
instrument(CFX, "qPCR", "Bio-Rad", "cfxmaestro", ["Fluorescence"])

# Rows of an .xlsx results sheet read for the header block and well table.
XLSX_ROWS = 4000

_MISSING = {"", "na", "n/a", "none", "nan"}
_WELL = re.compile(r"^[A-Za-z]{1,2}[0-9]{1,3}$")
_INT_FLOAT = re.compile(r"^-?[0-9]+\.0+$")
_SCI_INT = re.compile(r"^[0-9]\.[0-9]+[eE]\+?[0-9]{1,2}$")


def _clean(v) -> Optional[str]:
    """A header value, or None when it is blank or a not-applicable marker. An
    integer stored as a spreadsheet number ("12345.0") loses the ".0"."""
    if v is None:
        return None
    s = str(v).strip()
    if s.lower() in _MISSING:
        return None
    if _INT_FLOAT.match(s):
        s = s.split(".", 1)[0]
    elif _SCI_INT.match(s):
        # A long integer stored as a spreadsheet number ("1.23456789E8").
        x = float(s)
        mantissa, exp = s.lower().split("e")
        if x.is_integer() and x < 1e15 and len(mantissa) - 2 <= int(exp):
            s = str(int(x))
    return s


_BLOCK_WELLS = re.compile(r"(?:^|[^0-9])(96|384)(?:[^0-9]|$)")


def _plate_format(block_type: Optional[str]) -> Optional[int]:
    """Wells per plate from a QuantStudio block type, the way allotropy reads
    it ("96-Well Block (0.2mL)", "Fast 96-Well Block (0.1mL)", "384-Well
    Block")."""
    if not block_type:
        return None
    m = _BLOCK_WELLS.search(block_type)
    return int(m.group(1)) if m else None


def _distinct(values) -> List[str]:
    out: List[str] = []
    seen = set()
    for v in values:
        c = _clean(v)
        if c is not None and c not in seen:
            seen.add(c)
            out.append(c)
    return out


def _col(header_row: List[str], *names: str) -> Optional[int]:
    norm = [h.strip() for h in header_row]
    for n in names:
        if n in norm:
            return norm.index(n)
    return None


def _cell(row: List[str], i: Optional[int]) -> str:
    if i is None or i >= len(row):
        return ""
    return row[i]


def _well_table(header_row: List[str], rows: List[List[str]]) -> Dict[str, object]:
    """Targets, reporter dyes and counts from a QuantStudio / Diomni well table
    (one row per well and target). Only rows that name a target count, as in
    allotropy, which makes one measurement per well and target."""
    well = _col(header_row, "Well Position", "Well")
    target = _col(header_row, "Target Name", "Target")
    snp = no_target = False
    if target is None and _col(header_row, "SNP Assay Name") is not None:
        # Genotyping: one SNP assay per well, two allele reporters.
        snp = True
        target = _col(header_row, "SNP Assay Name")
    elif target is None:
        # Design & Analysis "simple" primary analysis: no target column, one
        # row per well and reporter.
        no_target = True
        target = _col(header_row, "Reporter")
    sample = _col(header_row, "Sample Name", "Sample")
    reporters = [i for i in (_col(header_row, "Reporter"), _col(header_row, "Allele1 Reporter"),
                             _col(header_row, "Allele2 Reporter")) if i is not None]
    if well is None or target is None:
        return {}
    wells, targets, samples, dyes = [], [], [], []
    measured = 0
    for r in rows:
        t = _clean(_cell(r, target))
        if t is None:
            continue
        measured += 1
        wells.append(_cell(r, well))
        targets.append(t)
        samples.append(_cell(r, sample))
        for i in reporters:
            dyes.append(_cell(r, i))
    out: Dict[str, object] = {}
    targets_d = _distinct(targets)
    if targets_d and not no_target:
        out["snpAssays" if snp else "targets"] = targets_d
    dyes_d = _distinct(dyes)
    if dyes_d:
        out["reporterDyes"] = dyes_d
    if measured:
        out["sampleCount"] = len(_distinct(wells))
        if not snp:
            out["measurementCount"] = measured
        names = len(_distinct(samples))
        if names:
            out["sampleNameCount"] = names
    return out


# --------------------------------------------------------------------------
# QuantStudio key/value header (txt `* Key = Value`, or the top of a sheet)


def _qs_fields(kv: Dict[str, str]) -> Dict[str, object]:
    """Fields every QuantStudio-family header block shares."""
    block = _clean(kv.get("Block Type"))
    run_end = None
    for k in ("Experiment Run End Time", "Experiment Run Stop Time", "Run End Date/Time", "Run End Data/Time"):
        run_end = run_end or _clean(kv.get(k))
    return {
        "instrumentModel": _clean(kv.get("Instrument Type")),
        "instrumentName": _clean(kv.get("Instrument Name")),
        "serialNumber": _clean(kv.get("Instrument Serial Number")),
        "blockType": block,
        "plateFormat": _plate_format(block),
        "chemistry": _clean(kv.get("Chemistry")),
        "experimentType": _clean(kv.get("Experiment Type")),
        "experimentName": _clean(kv.get("Experiment Name")),
        "passiveReference": _clean(kv.get("Passive Reference")),
        "quantificationCycleMethod": _clean(kv.get("Quantification Cycle Method")),
        "runDate": run_end,
    }


# Keys the Real-Time PCR layout writes; a header block must carry several.
_QS_LEGACY_KEYS = ("Block Type", "Chemistry", "Experiment Name", "Experiment Type",
                   "Experiment Run End Time", "Experiment Run Stop Time",
                   "Quantification Cycle Method", "Passive Reference", "Experiment File Name")


def _is_qs_legacy(kv: Dict[str, str]) -> bool:
    if "Instrument Type" not in kv:
        return False
    model = kv.get("Instrument Type", "")
    if not re.search(r"QuantStudio|ViiA|StepOne", model, re.I):
        return False
    return sum(1 for k in _QS_LEGACY_KEYS if k in kv) >= 3


def _summary(name: str, f: Dict[str, object], variant: str) -> str:
    parts = [name, variant]
    if f.get("instrumentModel"):
        parts.append(str(f["instrumentModel"]))
    if f.get("plateFormat"):
        parts.append("%s-well" % f["plateFormat"])
    if f.get("targets"):
        t = f["targets"]
        parts.append("%d target%s" % (len(t), "" if len(t) == 1 else "s"))
    if f.get("sampleCount") is not None:
        parts.append("%s wells" % f["sampleCount"])
    return ", ".join(parts)


_TXT_KV = re.compile(r"^\*[ \t]*(.+?)[ \t]*=[ \t]?(.*)$")
_SECTION = re.compile(r"^\[(.+)\]$")


def sniff_quantstudio_txt(src) -> Optional[object]:
    text = head_text(src)
    if not text:
        return None
    lines = text.splitlines()
    if not lines or not lines[0].startswith("*"):
        return None
    kv: Dict[str, str] = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if _SECTION.match(line):
            break
        m = _TXT_KV.match(line)
        if not m:
            return None
        kv[m.group(1).strip()] = m.group(2).strip()
        i += 1
    if i >= len(lines) or not _is_qs_legacy(kv):
        return None
    # Sections and the well table of the first one that has one.
    truncated = len(head(src)) >= TEXT_HEAD_BYTES and (src.size is None or src.size > TEXT_HEAD_BYTES)
    sections: List[str] = []
    table: Dict[str, object] = {}
    while i < len(lines):
        m = _SECTION.match(lines[i].strip())
        if not m:
            i += 1
            continue
        title = m.group(1).strip()
        sections.append(title)
        i += 1
        body: List[List[str]] = []
        closed = False
        while i < len(lines):
            line = lines[i].rstrip("\r\n")
            if not line.strip():
                closed = True
                break
            if _SECTION.match(line.strip()):
                closed = True
                break
            body.append(line.split("\t"))
            i += 1
        if i >= len(lines) and not truncated:
            closed = True
        if not table and title in ("Sample Setup", "Results") and body and closed:
            table = _well_table(body[0], body[1:])
    f = _qs_fields(kv)
    f.update(table)
    f["sections"] = sections
    f["exportFormat"] = "txt"
    return header(QS_RTPCR, _summary(QS_RTPCR, f, "txt"), f)


# --------------------------------------------------------------------------
# .xlsx exports (Python only)


def _kv_block(rows: List[List[str]]) -> Tuple[Dict[str, str], int]:
    """The key/value block at the top of a sheet: column A keys, column B
    values, up to the first row whose column A is blank. Returns the block and
    the index of that row."""
    kv: Dict[str, str] = {}
    for i, r in enumerate(rows):
        k = r[0].strip() if r else ""
        if not k or k in ("Well", "Well Position"):
            # A blank row ends the block; the reader drops rows the sheet
            # does not store, so the well table's column row can follow
            # directly.
            return kv, i
        kv[k] = r[1].strip() if len(r) > 1 else ""
    return kv, len(rows)


def _table_after(rows: List[List[str]], start: int) -> Tuple[Optional[List[str]], List[List[str]], bool]:
    """The well table after a header block: its column row (the first row
    whose column A is "Well"), the data rows up to the next blank row, and
    whether the table ended inside the rows read."""
    j = start
    while j < len(rows) and (not rows[j] or rows[j][0].strip() not in ("Well", "Well Position")):
        j += 1
        if j - start > 8:
            return None, [], False
    if j >= len(rows):
        return None, [], False
    cols = rows[j]
    body: List[List[str]] = []
    for r in rows[j + 1:]:
        if not r or not r[0].strip():
            return cols, body, True
        body.append(r)
    return cols, body, len(rows) < XLSX_ROWS


def _xlsx_kv(src, sheet: Optional[str]):
    rows = xlsx_rows(src, sheet=sheet, max_rows=XLSX_ROWS)
    if not rows:
        return None, None, None
    kv, end = _kv_block(rows)
    return rows, kv, end


def sniff_quantstudio_xlsx(src) -> Optional[object]:
    names = xlsx_sheet_names(src)
    if not names:
        return None
    sheet = None
    if "Results" in names or "Primary_result" in names:
        sheet = "Results" if "Results" in names else "Primary_result"
    candidates = [sheet] if sheet else []
    if names[0] not in candidates:
        candidates.append(names[0])
    for sh in candidates:
        rows, kv, end = _xlsx_kv(src, sh)
        if not kv:
            continue
        hdr = _decode_qs_sheet(names, sh, rows, kv, end)
        if hdr is not None:
            return hdr
    # Diomni's single-sheet Target Call export.
    for sh in ("Target_Call", "Target Call"):
        if sh in names and sh not in candidates:
            rows, kv, end = _xlsx_kv(src, sh)
            if kv:
                hdr = _decode_qs_sheet(names, sh, rows, kv, end)
                if hdr is not None:
                    return hdr
    return None


def _decode_qs_sheet(names, sheet, rows, kv, end) -> Optional[object]:
    software = _clean(kv.get("Software Name and Version"))
    soft_name = _clean(kv.get("Software name"))
    version = None
    if soft_name and soft_name.lower() == "diomni":
        name = DIOMNI
        version = _clean(kv.get("Software version"))
    elif software and software.startswith("Design & Analysis"):
        name = QS_DA
        m = re.match(r"^(.*?)\s+v(\S+)$", software)
        version = m.group(2) if m else None
    elif _is_qs_legacy(kv):
        # Design & Analysis v1.x writes the trademark sign; the Real-Time PCR
        # Software writes "(TM)".
        name = QS_DA if "™" in kv.get("Instrument Type", "") else QS_RTPCR
    else:
        return None
    if "Instrument Type" not in kv or "Block Type" not in kv:
        return None
    f = _qs_fields(kv)
    f["softwareVersion"] = version
    cols, body, closed = _table_after(rows, end)
    if cols is not None and closed:
        f.update(_well_table(cols, body))
    f["sheets"] = list(names)
    if name == DIOMNI:
        f["assayName"] = _clean(kv.get("Assay name"))
        f["assayVersion"] = _clean(kv.get("Assay version"))
        has_call = "Target_Call" in names or "Target Call" in names
        f["exportFormat"] = "xlsx target call" if len(names) == 1 else "xlsx"
        if not has_call:
            # Without a Target Call sheet only the Design & Analysis parser
            # reads it (from Primary_result).
            f["allotropyParser"] = "appbio_quantstudio_designandanalysis"
    elif name == QS_DA:
        f["analysisPlugins"] = _clean(kv.get("Plugin Name and Version"))
        f["exportFormat"] = "xlsx" if software else "xlsx (Real-Time PCR layout)"
    else:
        f["exportFormat"] = "xlsx"
    if name in (QS_RTPCR, QS_DA) and _is_qs_legacy(kv) and not software:
        # The Real-Time PCR layout: both QuantStudio parsers in allotropy
        # convert it.
        other = "appbio_quantstudio_designandanalysis" if name == QS_RTPCR else "appbio_quantstudio"
        f["allotropyAlternateParser"] = other
    return header(name, _summary(name, f, str(f["exportFormat"])), f)


# --------------------------------------------------------------------------
# CSV exports


def _csv_line(line: str) -> List[str]:
    """One CSV line split the way Python's csv module (excel dialect) splits
    it, kept to a single line so a stray quote cannot swallow the rest of the
    file. The Go port uses the same state machine."""
    if not line:
        return []
    fields: List[str] = []
    cur: List[str] = []
    state = 0  # 0 field start, 1 unquoted, 2 quoted, 3 quote inside quoted
    for ch in line:
        if state == 0:
            if ch == '"':
                state = 2
            elif ch == ",":
                fields.append("")
            else:
                cur.append(ch)
                state = 1
        elif state == 1:
            if ch == ",":
                fields.append("".join(cur))
                cur = []
                state = 0
            else:
                cur.append(ch)
        elif state == 2:
            if ch == '"':
                state = 3
            else:
                cur.append(ch)
        else:
            if ch == '"':
                cur.append('"')
                state = 2
            elif ch == ",":
                fields.append("".join(cur))
                cur = []
                state = 0
            else:
                cur.append(ch)
                state = 1
    fields.append("".join(cur))
    return fields


def _csv_rows(text: str, truncated: bool) -> List[List[str]]:
    lines = text.splitlines()
    if truncated and lines:
        lines = lines[:-1]  # the last line may be cut
    return [_csv_line(ln) for ln in lines]


def _text_and_truncated(src) -> Tuple[Optional[str], bool]:
    text = head_text(src)
    raw = head(src)
    truncated = len(raw) >= TEXT_HEAD_BYTES and (src.size is None or src.size > TEXT_HEAD_BYTES)
    return text, truncated


_ABSQ_BASE = ("Run", "Date", "Instrument", "Plate")


def _absq_counts(f, rows, wcol, scol) -> None:
    """Distinct wells and sample names over the per-well rows (group and
    summary exports interleave "Average" / "Pooled" aggregate rows)."""
    per_well = [r for r in rows if _WELL.match(_cell(r, wcol).strip())]
    f["sampleCount"] = len(_distinct(_cell(r, wcol) for r in per_well))
    names = len(_distinct(_cell(r, scol) for r in per_well)) if scol is not None else 0
    if names:
        f["sampleNameCount"] = names


def _decode_absolute_q(text: str, truncated: bool, variant_prefix: str = "") -> Optional[object]:
    rows = _csv_rows(text, truncated)
    if len(rows) < 2:
        return None
    first = [c.strip() for c in rows[0]]
    f: Dict[str, object] = {}
    if first[:4] == list(_ABSQ_BASE):
        cols = first
        data = [r for r in rows[1:] if any(c.strip() for c in r)]
        run = _col(cols, "Run")
        date = _col(cols, "Date")
        inst = _col(cols, "Instrument")
        plate = _col(cols, "Plate")
        if "Dye" in cols and "Target" in cols and ("Name" in cols or "Sample" in cols) and "Well" in cols:
            variant = "group results csv"
            dyes = _distinct(_cell(r, _col(cols, "Dye")) for r in data)
            targets = _distinct(_cell(r, _col(cols, "Target")) for r in data)
            wcol, scol = _col(cols, "Well"), _col(cols, "Name", "Sample")
        elif "Well ID" in cols and "Index" in cols and any(c.endswith("_target") for c in cols):
            variant = "fluorescence csv"
            dye_names = [c[: -len("_target")] for c in cols if c.endswith("_target") and c[: -len("_target")] in cols]
            dyes = dye_names
            targets = _distinct(_cell(r, _col(cols, d + "_target")) for d in dye_names for r in data)
            wcol, scol = _col(cols, "Well ID"), _col(cols, "Sample")
        elif "Channels" in cols and "Sample" in cols and "Count" in cols:
            variant = "multichannel csv"
            chans = _distinct(_cell(r, _col(cols, "Channels")) for r in data)
            dyes = [c for c in chans if " " not in c]
            targets = []
            wcol, scol = None, _col(cols, "Sample")
        else:
            return None
        if not data:
            return None
        f["runName"] = _clean(_cell(data[0], run))
        f["runDate"] = _clean(_cell(data[0], date))
        f["instrumentName"] = _clean(_cell(data[0], inst))
        plates = _distinct(_cell(r, plate) for r in data)
        f["plateId"] = plates[0] if len(plates) == 1 else None
        if dyes:
            f["reporterDyes"] = dyes
        if targets:
            f["targets"] = targets
        if not truncated:
            if wcol is not None:
                _absq_counts(f, data, wcol, scol)
            elif scol is not None:
                names = len(_distinct(_cell(r, scol) for r in data))
                if names:
                    f["sampleNameCount"] = names
    elif len(rows) >= 2 and [c.strip() for c in rows[1][:4]] == ["Run name", "Group", "Sample", "Well"] \
            and "Lambda (cp/Rxn)" in [c.strip() for c in rows[1]]:
        variant = "summary csv"
        cols = [c.strip() for c in rows[1]]
        dyes = _distinct(rows[0])
        if not dyes:
            return None
        tcols = [i for i, c in enumerate(cols) if c == "Target"]
        data = [r for r in rows[2:] if any(c.strip() for c in r)]
        targets = _distinct(_cell(r, i) for i in tcols for r in data)
        f["reporterDyes"] = dyes
        if targets:
            f["targets"] = targets
        runs = _distinct(_cell(r, 0) for r in data)
        f["runName"] = runs[0] if len(runs) == 1 else None
        if not truncated:
            _absq_counts(f, data, _col(cols, "Well"), _col(cols, "Sample"))
    else:
        return None
    variant = variant_prefix or variant
    f["exportFormat"] = variant
    return header(ABSOLUTE_Q, _summary(ABSOLUTE_Q, f, variant), f)


def sniff_absolute_q_csv(src) -> Optional[object]:
    text, truncated = _text_and_truncated(src)
    if not text:
        return None
    if not (text.startswith("Run,Date,Instrument,Plate,")
            or (text.startswith(",") and "\nRun name,Group,Sample,Well," in text[:8192])):
        return None
    return _decode_absolute_q(text, truncated)


def sniff_absolute_q_zip(src) -> Optional[object]:
    zf = open_zip(src)
    if zf is None:
        return None
    try:
        names = [n for n in zf.namelist() if n.endswith("_summary.csv")]
        if len(names) != 1:
            return None
        with zf.open(names[0]) as fh:
            data = fh.read(TEXT_HEAD_BYTES + 1)
    except Exception:
        return None
    text = text_of(data[:TEXT_HEAD_BYTES])
    if not text:
        return None
    hdr = _decode_absolute_q(text, len(data) > TEXT_HEAD_BYTES, "summary zip")
    if hdr is not None and hdr.fields.get("exportFormat") != "summary zip":
        return None
    return hdr


_QIA_PLATE = re.compile(r"(?:^|[^0-9])(24|96)-well", re.I)


def sniff_qiacuity_csv(src) -> Optional[object]:
    text, truncated = _text_and_truncated(src)
    if not text:
        return None
    if text.startswith("\ufeff"):
        text = text[1:]
    if not text.startswith("sep="):
        return None
    rows = _csv_rows(text, truncated)
    if len(rows) < 2:
        return None
    cols = [c.strip() for c in rows[1]]
    if "Sample/NTC/Control" not in cols or "partitions (valid)" not in [c.lower() for c in cols]:
        return None
    data = [r for r in rows[2:] if any(c.strip() for c in r)]
    well = _col(cols, "Well")
    if well is None and cols[0] == "":
        well = 0
    target = _col(cols, "Target", "Target (Name)")
    f: Dict[str, object] = {}
    f["plateName"] = (_distinct(_cell(r, _col(cols, "Plate name")) for r in data) or [None])[0] \
        if _col(cols, "Plate name") is not None else None
    pids = _distinct(_cell(r, _col(cols, "Plate ID")) for r in data) if _col(cols, "Plate ID") is not None else []
    f["plateId"] = pids[0] if len(pids) == 1 else None
    ptypes = _distinct(_cell(r, _col(cols, "Plate type")) for r in data) if _col(cols, "Plate type") is not None else []
    if len(ptypes) == 1:
        f["plateType"] = ptypes[0]
        m = _QIA_PLATE.search(ptypes[0])
        f["plateFormat"] = int(m.group(1)) if m else None
    targets = _distinct(_cell(r, target) for r in data) if target is not None else []
    if targets:
        f["targets"] = targets
    mixes = _distinct(_cell(r, _col(cols, "Reaction Mix")) for r in data) if _col(cols, "Reaction Mix") is not None else []
    if mixes:
        f["reactionMixes"] = mixes
    if not truncated and well is not None:
        f["sampleCount"] = len(_distinct(_cell(r, well) for r in data))
        f["measurementCount"] = len(data)
        names = len(_distinct(_cell(r, _col(cols, "Sample/NTC/Control")) for r in data))
        if names:
            f["sampleNameCount"] = names
    f["exportFormat"] = "csv"
    return header(QIACUITY, _summary(QIACUITY, f, "csv"), f)


_CFX_HEAD = ["", "Well", "Fluor", "Target", "Content", "Sample"]


def sniff_cfx_csv(src) -> Optional[object]:
    text, truncated = _text_and_truncated(src)
    if not text:
        return None
    if text.startswith("\ufeff"):
        text = text[1:]
    if not text.startswith(",Well,Fluor,"):
        return None
    rows = _csv_rows(text, truncated)
    if not rows:
        return None
    cols = [c.strip() for c in rows[0]]
    if cols[:6] != _CFX_HEAD or "Cq" not in cols:
        return None
    data = [r for r in rows[1:] if len(r) > 1 and r[1].strip()]
    f: Dict[str, object] = {}
    dyes = _distinct(_cell(r, 2) for r in data)
    targets = _distinct(_cell(r, 3) for r in data)
    if dyes:
        f["reporterDyes"] = dyes
    if targets:
        f["targets"] = targets
    if not truncated:
        f["sampleCount"] = len(_distinct(_cell(r, 1) for r in data))
        f["measurementCount"] = len(data)
        names = len(_distinct(_cell(r, 5) for r in data))
        if names:
            f["sampleNameCount"] = names
    f["exportFormat"] = "quantification cq results csv"
    return header(CFX, _summary(CFX, f, "Cq results csv"), f)


register_instrument_sniffer("txt", sniff_quantstudio_txt)
register_instrument_sniffer("xlsx", sniff_quantstudio_xlsx)
register_instrument_sniffer("csv", sniff_qiacuity_csv)
register_instrument_sniffer("csv", sniff_cfx_csv)
register_instrument_sniffer("csv", sniff_absolute_q_csv)
register_instrument_sniffer("zip", sniff_absolute_q_zip)
