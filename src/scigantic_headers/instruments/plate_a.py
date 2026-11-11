"""Plate readers, family A: Agilent Gen5 and Gen5 Image, Molecular Devices
SoftMax Pro, BMG Labtech MARS and SMART Control, Tecan Magellan and Thermo
Fisher Scientific SkanIt.

Gen5, Gen5 Image, SoftMax Pro and MARS write text exports and are read from
the leading bytes (the Go port in scigantic-edge reads them too). SMART
Control, Magellan and SkanIt write .xlsx workbooks and are read here only.

Gen5 and Gen5 Image write the same text layout (a "Software Version" header
block, then "Procedure Details"). They are told apart by the read steps in
Procedure Details: imaging steps ("Image Montage", "Image Single Image",
"Image Z-Stack") with no photometric step mean Gen5 Image; any absorbance,
fluorescence, luminescence or Alpha step means Gen5.

Field conventions here: `wavelengthsNm` holds detection wavelengths
(absorbance, or emission for fluorescence and luminescence) and
`excitationWavelengthsNm` the excitation side. Dates are reported as the file
writes them unless the file's own format is unambiguous ISO (Magellan).
Placeholder values such as "n/a" and "Unknown" are dropped.
"""

from __future__ import annotations

import re as _re_unicode
from typing import Dict, List, Optional


class _AsciiRe:
    """`re` with ASCII classes (\\d is [0-9], \\s is ASCII space), so the
    patterns mean what Go's regexp means and the two ports agree."""

    A = _re_unicode.A

    @staticmethod
    def compile(p):
        return _re_unicode.compile(p, _re_unicode.A)

    @staticmethod
    def match(p, s):
        return _re_unicode.match(p, s, _re_unicode.A)

    @staticmethod
    def search(p, s):
        return _re_unicode.search(p, s, _re_unicode.A)

    @staticmethod
    def fullmatch(p, s):
        return _re_unicode.fullmatch(p, s, _re_unicode.A)

    @staticmethod
    def findall(p, s):
        return _re_unicode.findall(p, s, _re_unicode.A)

    @staticmethod
    def split(p, s):
        return _re_unicode.split(p, s, flags=_re_unicode.A)


re = _AsciiRe

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

GEN5 = "Agilent Gen5"
GEN5_IMAGE = "Agilent Gen5 Image"
SOFTMAX = "Molecular Devices SoftMax Pro"
MARS = "BMG Labtech MARS"
SMART = "BMG Labtech SMART Control"
MAGELLAN = "Tecan Magellan"
SKANIT = "Thermo Fisher Scientific SkanIt"

instrument(GEN5, "Plate Reader", "Agilent", "agilent_gen5", ["Absorbance", "Fluorescence", "Luminescence"])
instrument(GEN5_IMAGE, "Plate Reader", "Agilent", "agilent_gen5_image", ["Optical Imaging"])
instrument(SOFTMAX, "Plate Reader", "Molecular Devices", "moldev_softmax_pro", ["Absorbance", "Fluorescence", "Luminescence"])
instrument(MARS, "Plate Reader", "BMG Labtech", "bmg_mars", ["Absorbance", "Fluorescence", "Luminescence"])
instrument(SMART, "Plate Reader", "BMG Labtech", "bmg_labtech_smart_control", ["Fluorescence"])
instrument(MAGELLAN, "Plate Reader", "Tecan", "tecan_magellan", ["Absorbance"])
instrument(SKANIT, "Plate Reader", "Thermo Fisher Scientific", "thermo_skanit", ["Absorbance", "Fluorescence", "Luminescence"])

_MODES = ("Absorbance", "Fluorescence", "Luminescence")
_STANDARD_WELLS = (6, 12, 24, 48, 96, 384, 1536, 3456)
# (columns, rows) of a full plate grid -> wells.
_GRID_WELLS = {(3, 2): 6, (4, 3): 12, (6, 4): 24, (8, 6): 48, (12, 8): 96,
               (24, 16): 384, (48, 32): 1536, (72, 48): 3456}
_PLACEHOLDERS = {"", "n/a", "na", "none", "unknown", "-", "???"}


def _value(s: Optional[str]) -> Optional[str]:
    if s is None:
        return None
    s = s.strip().strip('"').strip()
    return None if s.lower() in _PLACEHOLDERS else s


def _uniq(xs) -> list:
    out = []
    for x in xs:
        if x not in out:
            out.append(x)
    return out


def _num(s: str):
    """An int (or float) from a wavelength token like '485', '485/20',
    '480-14', '450.5'; None if the token does not start with a number."""
    m = re.match(r"\s*(\d+(?:\.\d+)?)", s or "")
    if not m:
        return None
    v = float(m.group(1))
    return int(v) if v.is_integer() else v


def _plate_wells_from_name(name: Optional[str]) -> Optional[int]:
    """Wells per plate from a plate type name ('96 WELL PLATE', 'Nunc 384 flat
    bottom', 'Costar 6 well'): the one standard well count written as a
    standalone number, else None."""
    if not name:
        return None
    found = {int(n) for n in re.findall(r"(?<![\d.])(\d{1,4})(?![\d.])", name) if int(n) in _STANDARD_WELLS}
    return found.pop() if len(found) == 1 else None


def _complete(src) -> bool:
    """True when the leading read holds the whole file, so counts over it are
    counts over the file."""
    data = head(src)
    return len(data) < TEXT_HEAD_BYTES and (src.size is None or src.size <= len(data))


def _mode_word(s: str) -> Optional[str]:
    for m in _MODES:
        if s.strip().lower().startswith(m.lower()):
            return m
    return None


# ── Agilent Gen5 and Gen5 Image ───────────────────────────────────────────────

_GEN5_KEYS = ("Experiment File Path:", "Protocol File Path:", "Plate Number", "Date", "Time",
              "Reader Type:", "Reader Serial Number:", "Reading Type")
_GEN5_PHOTOMETRIC = re.compile(r"(?:^|\t)\s*(Absorbance|Fluorescence|Luminescence|Alpha) (Endpoint|Spectrum|Area Scan)\b")
_GEN5_IMAGE = re.compile(r"(?:^|\t)\s*Image (Single Image|Montage|Z-Stack)\b")
_GEN5_CHANNEL = re.compile(r"^\s*Channel \d+:\s+(.+?)\s*$")


def _basename(path: Optional[str]) -> Optional[str]:
    if not path:
        return None
    return re.split(r"[\\/]", path.strip())[-1] or None


def _gen5_parse(src) -> Optional[dict]:
    t = head_text(src)
    if not t or "Software Version" not in t[:4096]:
        return None
    lines = t.splitlines()
    i = 0
    while i < len(lines) and i < 20 and not lines[i].strip():
        i += 1
    if i >= len(lines) or not lines[i].startswith("Software Version\t"):
        return None
    hdr: Dict[str, str] = {}
    j = i
    proc = None
    while j < len(lines) and j < i + 80:
        line = lines[j]
        if line.strip() == "Procedure Details":
            proc = j
            break
        if "\t" in line:
            k, v = line.split("\t", 1)
            if k in ("Software Version",) + _GEN5_KEYS and k not in hdr:
                hdr[k] = v.strip("\t").strip()
        j += 1
    if proc is None or sum(1 for k in _GEN5_KEYS if k in hdr) < 3:
        return None
    version = hdr.get("Software Version", "").strip()
    if not re.match(r"^\d+(\.\d+)+$", version):
        return None
    # Procedure Details: a blank line after the title, then lines up to the next
    # blank line.
    k = proc + 1
    while k < len(lines) and not lines[k].strip():
        k += 1
    block: List[str] = []
    while k < len(lines) and lines[k].strip():
        block.append(lines[k])
        k += 1
    return {"lines": lines, "hdr": hdr, "version": version, "proc": block}


def _gen5_common(p: dict, src) -> dict:
    hdr = p["hdr"]
    plate_type = None
    for line in p["proc"]:
        if line.startswith("Plate Type\t"):
            plate_type = _value(line.split("\t", 1)[1])
            break
    date, time_ = _value(hdr.get("Date")), _value(hdr.get("Time"))
    fields = {
        "softwareVersion": p["version"],
        "instrumentModel": _value(hdr.get("Reader Type:")),
        "serialNumber": _value(hdr.get("Reader Serial Number:")),
        "plateType": plate_type,
        "plateFormat": _plate_wells_from_name(plate_type),
        "runDate": "%s %s" % (date, time_) if date and time_ else date,
        "readingType": _value(hdr.get("Reading Type")),
        "experimentFile": _basename(_value(hdr.get("Experiment File Path:"))),
        "protocolFile": _basename(_value(hdr.get("Protocol File Path:"))),
        "exportFormat": "text export",
    }
    if _complete(src):
        fields["plateCount"] = sum(1 for line in p["lines"] if line.startswith("Software Version\t"))
    return fields


def _gen5_read_type(block: List[str], kinds: List[str]) -> Optional[str]:
    if any(line.strip().startswith(("Start Kinetic", "Discontinuous Kinetics")) for line in block):
        return "Kinetic"
    kinds = _uniq(kinds)
    return kinds[0] if len(kinds) == 1 else None


def _gen5_filters(block: List[str]):
    ex, em = [], []
    for line in block:
        for part in line.strip().split(","):
            part = part.strip()
            if part.startswith("Excitation:"):
                v = _num(part[len("Excitation:"):])
                if v is not None:
                    ex.append(v)
            elif part.startswith("Emission:"):
                v = _num(part[len("Emission:"):])
                if v is not None:
                    em.append(v)
    return _uniq(ex), _uniq(em)


def sniff_gen5(src):
    p = _gen5_parse(src)
    if p is None:
        return None
    block = p["proc"]
    modes, kinds = [], []
    for line in block:
        for m in _GEN5_PHOTOMETRIC.finditer(line):
            modes.append(m.group(1))
            kinds.append(m.group(2))
    has_image = any(_GEN5_IMAGE.search(line) for line in block)
    if not modes and has_image:
        return None  # Gen5 Image
    fields = _gen5_common(p, src)
    abs_wl: List[object] = []
    for line in block:
        s = line.strip()
        if s.startswith("Wavelengths:"):
            toks = [x.strip() for x in s[len("Wavelengths:"):].split(",")]
            vals = [_num(x) for x in toks]
            if toks and all(v is not None and re.fullmatch(r"\d+(\.\d+)?", x) for v, x in zip(vals, toks)):
                abs_wl.extend(vals)
    ex, em = _gen5_filters(block)
    read_modes = _uniq(modes) + (["Imaging"] if has_image else [])
    fields.update({
        "readModes": read_modes or None,
        "readType": _gen5_read_type(block, kinds),
        "readCount": len(modes) or None,
        "wavelengthsNm": _uniq(abs_wl + em) or None,
        "excitationWavelengthsNm": ex or None,
    })
    bits = ["Gen5 %s" % p["version"]]
    if read_modes:
        bits.append(", ".join(read_modes))
    if fields.get("plateFormat"):
        bits.append("%d-well" % fields["plateFormat"])
    if fields.get("instrumentModel"):
        bits.append(fields["instrumentModel"])
    return header(GEN5, "Agilent Gen5 plate reader export: " + ", ".join(bits), fields)


def sniff_gen5_image(src):
    p = _gen5_parse(src)
    if p is None:
        return None
    block = p["proc"]
    if any(_GEN5_PHOTOMETRIC.search(line) for line in block):
        return None
    kinds = [m.group(1) for line in block for m in _GEN5_IMAGE.finditer(line)]
    if not kinds:
        return None
    fields = _gen5_common(p, src)
    channels, ex, em, objectives = [], [], [], []
    for line in block:
        m = _GEN5_CHANNEL.match(line)
        if m:
            name = m.group(1)
            channels.append(name)
            w = re.search(r"\s(\d{3}),(\d{3})$", name)
            if w:
                ex.append(int(w.group(1)))
                em.append(int(w.group(2)))
            continue
        s = line.strip()
        if s.startswith("Objective:"):
            objectives.append(s[len("Objective:"):].strip())
    kinetic = any(line.strip().startswith(("Start Kinetic", "Discontinuous Kinetics")) for line in block)
    fields.update({
        "readModes": ["Imaging"],
        "imagingModes": _uniq(kinds),
        "readType": "Kinetic" if kinetic else None,
        "readCount": len(kinds),
        "objectives": _uniq([o for o in objectives if o]) or None,
        "imagingChannels": _uniq(channels) or None,
        "excitationWavelengthsNm": _uniq(ex) or None,
        "wavelengthsNm": _uniq(em) or None,
    })
    bits = ["Gen5 %s" % p["version"], ", ".join(_uniq(kinds))]
    if fields.get("plateFormat"):
        bits.append("%d-well" % fields["plateFormat"])
    if fields.get("instrumentModel"):
        bits.append(fields["instrumentModel"])
    return header(GEN5_IMAGE, "Agilent Gen5 imaging export: " + ", ".join(bits), fields)


# ── Molecular Devices SoftMax Pro ─────────────────────────────────────────────

_SMP_BLOCKS = re.compile(r"^##BLOCKS=\s*(\d{1,9})\s*$")


def _ints(s: Optional[str]) -> Optional[list]:
    if not s or not s.strip():
        return None
    out = []
    for tok in s.split():
        if not re.fullmatch(r"\d+(\.\d+)?", tok):
            return None
        out.append(_num(tok))
    return out


def sniff_softmax_pro(src):
    t = head_text(src)
    if not t or "##BLOCKS=" not in t[:256]:
        return None
    lines = t.splitlines()
    first = next((ln for ln in lines[:5] if ln.strip()), "")
    m = _SMP_BLOCKS.match(first.strip())
    if not m:
        return None
    blocks = int(m.group(1))
    plates = [ln.split("\t") for ln in lines if ln.startswith("Plate:\t")]
    groups = [ln for ln in lines if ln.startswith("Group:")]
    if not plates and not groups and not any(ln.startswith("Note:") for ln in lines[:50]):
        return None
    versions, formats, modes, types, wells, wl, ex = [], [], [], [], [], [], []
    spectrum = None
    for h in plates:
        h = [c.strip() for c in h]
        if len(h) < 19:
            continue
        versions.append(h[2])
        formats.append(h[3])
        types.append(h[4])
        modes.append(h[5])
        if h[5] not in _MODES:
            continue  # other read modes lay the header out differently
        fl = h[5] == "Fluorescence"
        wl_i, wells_i = (16, 19) if fl else (15, 18)
        if len(h) > wells_i and re.fullmatch(r"\d{1,9}", h[wells_i]):
            wells.append(int(h[wells_i]))
        if h[4] == "Spectrum":
            s = [_num(x) for x in h[11:14]]
            if all(v is not None for v in s) and spectrum is None:
                spectrum = s
            continue
        vals = _ints(h[wl_i] if len(h) > wl_i else None)
        if vals:
            wl.extend(v for v in vals if v != 0)
        if fl and len(h) > 20:
            exv = _ints(h[20])
            if exv:
                ex.extend(exv)
    complete = _complete(src)
    fields = {
        "exportVersion": versions[0] if len(set(versions)) == 1 else None,
        "exportFormat": formats[0] if len(set(formats)) == 1 else None,
        "readModes": _uniq(modes) or None,
        "readType": types[0] if len(set(types)) == 1 else None,
        "plateFormat": wells[0] if wells and len(set(wells)) == 1 else None,
        "blockCount": blocks,
        "plateCount": len(plates) if complete else None,
        "groupCount": len(groups) if complete else None,
        "wavelengthsNm": _uniq(wl) or None,
        "excitationWavelengthsNm": _uniq(ex) or None,
    }
    if spectrum:
        fields["wavelengthRangeNm"] = spectrum[:2]
        fields["wavelengthStepNm"] = spectrum[2]
    bits = ["%d blocks" % blocks]
    if fields["readModes"]:
        bits.append(", ".join(fields["readModes"]))
    if fields["readType"]:
        bits.append(fields["readType"])
    if fields["plateFormat"]:
        bits.append("%d-well" % fields["plateFormat"])
    return header(SOFTMAX, "SoftMax Pro plate export: " + ", ".join(bits), fields)


# ── BMG Labtech MARS ──────────────────────────────────────────────────────────

def _kv_cells(line: str) -> Dict[str, str]:
    out = {}
    for cell in line.split(","):
        cell = cell.strip().strip('"')
        if ": " in cell:
            k, v = cell.split(": ", 1)
            out.setdefault(k.strip(), v.strip())
        elif cell.endswith(":"):
            out.setdefault(cell[:-1].strip(), "")
    return out


def _grid(lines: List[str], start: int, sep: str, complete: bool):
    """(columns, rows, non-empty cells) of the plate grid whose header row
    (',1,2,...') is at or after `start`, or None. `rows` is None when the grid
    runs past the leading read."""
    i = start
    while i < len(lines) and i < start + 4 and not lines[i].strip(sep + ' "'):
        i += 1
    if i >= len(lines):
        return None
    cells = [c.strip().strip('"') for c in lines[i].split(sep)]
    if not cells or cells[0] != "":
        return None
    labels = [c for c in cells[1:] if c != ""]
    nums = [_num(c) for c in labels]
    if not labels or any(not re.fullmatch(r"\d+(\.0)?", c) for c in labels):
        return None
    ncols = len(labels)
    if [int(n) for n in nums] != list(range(1, ncols + 1)):
        return None
    rows = 0
    filled = 0
    j = i + 1
    ended = False
    while j < len(lines):
        row = [c.strip().strip('"') for c in lines[j].split(sep)]
        if not row or not re.fullmatch(r"[A-Za-z]{1,2}", row[0]):
            ended = True
            break
        rows += 1
        filled += sum(1 for c in row[1:ncols + 1] if c != "")
        j += 1
    return ncols, (rows if ended or complete else None), filled


def sniff_mars(src):
    t = head_text(src)
    if not t or not t.lstrip('"').startswith("User: "):
        return None
    lines = t.splitlines()
    if len(lines) < 4:
        return None
    l1, l2 = _kv_cells(lines[0]), _kv_cells(lines[1])
    if "Path" not in l1 or "Test name" not in l2 or "Date" not in l2 or "Time" not in l2:
        return None
    raw_i = next((i for i, ln in enumerate(lines[:40]) if re.match(r'^"?\s*Raw Data\b', ln.lstrip(","))), None)
    if raw_i is None:
        return None
    hdr: Dict[str, str] = {}
    mode_line = None
    for ln in lines[:raw_i]:
        hdr.update({k: v for k, v in _kv_cells(ln).items() if k not in hdr})
        s = ln.strip().strip(",").strip().strip('"')
        if s and ": " not in s and _mode_word(s) and mode_line is None:
            mode_line = s
    mode = _mode_word(mode_line) if mode_line else None
    raw = lines[raw_i].strip().strip(",").strip().strip('"')
    wl = ex = None
    fm = re.search(r"Raw Data\s*\((.*?)\)", raw)
    if fm:
        inner = fm.group(1).strip()
        if re.fullmatch(r"\d+(\.\d+)?/\d+(\.\d+)?", inner):
            a, b = inner.split("/")
            ex, wl = [_num(a)], [_num(b)]
        elif re.fullmatch(r"\d+(\.\d+)?", inner):
            wl = [_num(inner)]
    g = _grid(lines, raw_i + 1, ",", _complete(src))
    fields = {
        "testName": _value(hdr.get("Test name")),
        "testRunNumber": _value(hdr.get("Test run no.")),
        "runDate": "%s %s" % (hdr["Date"], hdr["Time"]) if _value(hdr.get("Date")) and _value(hdr.get("Time")) else None,
        "plateIdentifier": _value(hdr.get("ID1")),
        "measurementMode": mode_line,
        "readModes": [mode] if mode else None,
        "wavelengthsNm": wl,
        "excitationWavelengthsNm": ex,
        "exportFormat": "csv",
    }
    if g:
        ncols, nrows, filled = g
        if nrows is not None:
            fields["plateFormat"] = _GRID_WELLS.get((ncols, nrows))
            fields["sampleCount"] = filled
    bits = [fields["measurementMode"] or "unknown mode"]
    if fields.get("plateFormat"):
        bits.append("%d-well" % fields["plateFormat"])
    if fields["testName"]:
        bits.append("test %s" % fields["testName"])
    return header(MARS, "BMG MARS plate reader export: " + ", ".join(bits), fields)


# ── xlsx helpers ──────────────────────────────────────────────────────────────

def _sheet(names: List[str], wanted: str) -> Optional[str]:
    for n in names:
        if n.strip().lower() == wanted:
            return n
    return None


def _cells(r) -> List[str]:
    return [c.strip() for c in (r or []) if c and c.strip()]


def _first(rows, key: str) -> Optional[str]:
    """The value of the first 'key: value' cell, or of a row whose first
    non-empty cell is the label 'key' (value: the next non-empty cell)."""
    for r in rows or []:
        cs = _cells(r)
        if not cs:
            continue
        c0 = cs[0]
        if c0.startswith(key + ": "):
            return c0[len(key) + 2:].strip()
        if c0.rstrip(":").strip() == key:
            return cs[1] if len(cs) > 1 else ""
    return None


def _excel_date(src, s: Optional[str]) -> Optional[str]:
    """A date-time cell: its text when the workbook stored text, or the ISO
    form of a stored serial date (honoring the 1904 date system)."""
    if not s or not re.fullmatch(r"\d+(\.\d+)?", s):
        return s
    try:
        import datetime

        wb = zip_member(src._instrument_zip, "xl/workbook.xml", 64 * 1024) or b""
        base = datetime.datetime(1904, 1, 1) if re.search(rb'date1904="(1|true)"', wb) else datetime.datetime(1899, 12, 30)
        v = float(s)
        if not 1 <= v < 2958466:
            return None
        dt = base + datetime.timedelta(days=v)
        dt = (dt + datetime.timedelta(microseconds=500000)).replace(microsecond=0)
        return dt.isoformat()
    except Exception:
        return None


def _xlsx_number_text(s: Optional[str]) -> Optional[str]:
    """'20.0' -> '20'; leaves other text alone."""
    if s and re.fullmatch(r"-?\d+\.0+", s):
        return s.split(".")[0]
    return s


# ── BMG Labtech SMART Control ─────────────────────────────────────────────────

def sniff_smart_control(src):
    names = xlsx_sheet_names(src)
    if not names:
        return None
    proto_n = _sheet(names, "protocol information")
    table_n = _sheet(names, "table end point")
    if not proto_n or not (table_n or _sheet(names, "microplate end point")):
        return None
    proto = xlsx_rows(src, proto_n, 120)
    if not proto or _first(proto, "Test Name") is None or _first(proto, "Measurement type") is None:
        return None
    path = _first(proto, "Path") or ""
    if "SMART Control" not in path and not table_n:
        return None
    mtype = _value(_first(proto, "Measurement type"))
    plate = _value(_first(proto, "Microplate name"))
    wells = None
    if plate and re.fullmatch(r"\d{1,4}", plate.split()[-1]) and int(plate.split()[-1]) in _STANDARD_WELLS:
        wells = int(plate.split()[-1])
    exs, ems = _value(_first(proto, "Excitation")), _value(_first(proto, "Emission"))
    samples = None
    if table_n:
        table = xlsx_rows(src, table_n, 60)
        for r in table or []:
            if r and r[0].strip() == "Well":
                samples = sum(1 for c in r[1:] if re.fullmatch(r"[A-Z]{1,2}\d{1,3}", c.strip()))
                break
    date, time_ = _value(_first(proto, "Date")), _value(_first(proto, "Time"))
    mode = _mode_word(mtype) if mtype else None
    fields = {
        "testName": _value(_first(proto, "Test Name")),
        "runDate": "%s %s" % (date, time_) if date and time_ else date,
        "measurementMode": mtype,
        "readModes": [mode] if mode else None,
        "plateType": plate,
        "plateFormat": wells,
        "excitationWavelengthsNm": [_num(exs)] if exs and _num(exs) is not None else None,
        "wavelengthsNm": [_num(ems)] if ems and _num(ems) is not None else None,
        "opticPreset": _value(_first(proto, "Presetname")),
        "sampleCount": samples or None,
        "exportFormat": "xlsx",
    }
    bits = [mtype or "unknown mode"]
    if wells:
        bits.append("%d-well" % wells)
    if fields["testName"]:
        bits.append("test %s" % fields["testName"])
    return header(SMART, "BMG SMART Control export: " + ", ".join(bits), fields)


# ── Tecan Magellan ────────────────────────────────────────────────────────────

_MAG_TIME = re.compile(r"^Date of measurement: (\d{4}-\d{2}-\d{2})/Time of measurement: (\d{1,2}:\d{2}:\d{2})")


def sniff_magellan(src):
    names = xlsx_sheet_names(src)
    if not names:
        return None
    rows = xlsx_rows(src, names[0], 6000)
    if not rows or "Well positions" not in [c.strip() for c in rows[0]]:
        return None
    col = [c.strip() for c in rows[0]].index("Well positions")
    t_i = next((i for i, r in enumerate(rows) if r and _MAG_TIME.match(r[0].strip())), None)
    if t_i is None:
        return None
    meta = [r[0].strip() if r else "" for r in rows[t_i:]]
    if not any(x.lower().endswith((".mth", ".wsp")) for x in meta[1:6]):
        return None
    m = _MAG_TIME.match(meta[0])
    serial_i = next((i for i, x in enumerate(meta) if x.startswith("Instrument serial number:")), None)
    serial = model = None
    if serial_i is not None:
        serial = _value(meta[serial_i].split(":", 1)[1])
        if serial_i >= 2 and not meta[serial_i - 1].lower().endswith((".mth", ".wsp", "nm")):
            model = _value(meta[serial_i - 1])
    mth = next((x for x in meta[1:6] if x.lower().endswith(".mth") and "\\" not in x and "/" not in x), None)
    wsp = next((x for x in meta[1:8] if x.lower().endswith(".wsp") and "\\" not in x and "/" not in x), None)
    wl, modes = [], []
    for x in meta:
        if x.startswith("Measurement wavelength:"):
            v = _num(x.split(":", 1)[1])
            if v is not None:
                wl.append(v)
        elif x.startswith("Measurement mode:"):
            mo = _mode_word(x.split(":", 1)[1])
            if mo:
                modes.append(mo)
        elif x in _MODES:
            modes.append(x)
    wells = sum(1 for r in rows[1:t_i] if len(r) > col and re.fullmatch(r"[A-Z]{1,2}\d{1,3}", r[col].strip()))
    fields = {
        "runDate": "%sT%s" % (m.group(1), m.group(2).zfill(8)),
        # The device name Magellan writes above the serial number (allotropy's
        # device identifier); often the reader type, but user-configurable.
        "instrumentName": model,
        "serialNumber": serial,
        "methodFile": mth,
        "workspaceFile": wsp,
        "readModes": _uniq(modes) or None,
        "wavelengthsNm": _uniq(wl) or None,
        "sampleCount": wells or None,
        "exportFormat": "xlsx",
    }
    bits = []
    if fields["readModes"]:
        bits.append(", ".join(fields["readModes"]))
    if fields["wavelengthsNm"]:
        bits.append("at " + ", ".join("%s nm" % w for w in fields["wavelengthsNm"]))
    if wells:
        bits.append("%d wells" % wells)
    if model:
        bits.append(model)
    return header(MAGELLAN, "Tecan Magellan export: " + (", ".join(bits) or "plate data"), fields)


# ── Thermo Fisher Scientific SkanIt ───────────────────────────────────────────

def sniff_skanit(src):
    names = xlsx_sheet_names(src)
    if not names:
        return None
    meas = [n for n in names if _mode_word(n) and re.search(r"_\d+$", n)]
    gen_n = _sheet(names, "general information")
    first = xlsx_rows(src, meas[0], 3) if meas else None
    is_results = bool(first and len(first) >= 2 and first[0] and first[0][0].strip() == "Measurement results"
                      and first[1] and first[1][0].strip().lower().endswith(".skax"))
    gen = xlsx_rows(src, gen_n, 20) if gen_n else None
    sw = _first(gen, "Report generated with SW version") if gen else None
    if not is_results and not (sw and "SkanIt" in sw):
        return None
    sess_n = _sheet(names, "session information")
    sess = xlsx_rows(src, sess_n, 30) if sess_n else None
    executed = _first(sess, "Executed with") if sess else None
    version = None
    for s in (sw, executed):
        vm = re.search(r"\bver\.?\s*(\d+(?:\.\d+)+)", s or "")
        if vm:
            version = vm.group(1)
            break
    inst_n = _sheet(names, "instrument information")
    inst = xlsx_rows(src, inst_n, 12) if inst_n else None
    model = serial = esw = None
    for r in inst or []:
        cs = _cells(r)
        if not cs:
            continue
        if cs[0] == "Instrument modules":
            break
        if len(cs) > 1:
            k, v = cs[0], cs[1]
            if k == "Name" and model is None:
                model = _value(v)
            elif k == "Serial number" and serial is None:
                serial = _value(_xlsx_number_text(v))
            elif k == "ESW version" and esw is None:
                esw = _value(v)
    wl: List[object] = []
    grid = None
    for n in meas[:6]:
        rows = xlsx_rows(src, n, 400) or []
        for i, r in enumerate(rows):
            c0 = r[0].strip() if r else ""
            wm = re.fullmatch(r"Wavelength:\s*(\d+(?:\.\d+)?)\s*nm", c0)
            if wm and _num(wm.group(1)):
                wl.append(_num(wm.group(1)))
            if grid is None and c0 in ("Abs", "RLU", "RFU") and len(r) > 1:
                labels = [x.strip() for x in r[1:] if x.strip()]
                if labels and all(re.fullmatch(r"\d+(\.0)?", x) for x in labels) and \
                        [int(float(x)) for x in labels] == list(range(1, len(labels) + 1)):
                    nrows = 0
                    for rr in rows[i + 1:]:
                        if rr and re.fullmatch(r"[A-Z]{1,2}", rr[0].strip()):
                            nrows += 1
                        else:
                            break
                    grid = _GRID_WELLS.get((len(labels), nrows))
    run = None
    if first and len(first) >= 3 and first[2]:
        run = _excel_date(src, _value(first[2][0]))
    if run is None and sess:
        run = _excel_date(src, _value(_first(sess, "Execution time")))
    session = _value(first[1][0]) if is_results else (_value(_first(sess, "Session name")) if sess else None)
    modes = _uniq(_mode_word(n) for n in meas)
    fields = {
        "softwareVersion": version,
        "instrumentModel": model,
        "serialNumber": serial,
        "firmwareVersion": esw,
        "sessionName": session,
        "runDate": run,
        "readModes": modes or None,
        "wavelengthsNm": _uniq(w for w in wl if w) or None,
        "plateFormat": grid,
        "exportFormat": "xlsx",
    }
    bits = []
    if version:
        bits.append("SkanIt %s" % version)
    if modes:
        bits.append(", ".join(modes))
    if grid:
        bits.append("%d-well" % grid)
    if model:
        bits.append(model)
    return header(SKANIT, "Thermo SkanIt export: " + (", ".join(bits) or "plate data"), fields)


register_instrument_sniffer("txt", sniff_gen5)
register_instrument_sniffer("txt", sniff_gen5_image)
register_instrument_sniffer("txt", sniff_softmax_pro)
register_instrument_sniffer("csv", sniff_mars)
register_instrument_sniffer("xlsx", sniff_smart_control)
register_instrument_sniffer("xlsx", sniff_magellan)
register_instrument_sniffer("xlsx", sniff_skanit)
