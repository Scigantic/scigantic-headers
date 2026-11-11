"""Minimal .xlsx worksheet reader on the standard library: workbook sheet
names, shared strings, and the first rows of one sheet. Enough to recognize an
instrument export and read its header block, not a spreadsheet library."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_PKG_REL = "{http://schemas.openxmlformats.org/package/2006/relationships}"
_MAX_PART = 8 * 1024 * 1024
_MAX_SHARED = 200_000


def _read(zf, name: str, limit: int = _MAX_PART) -> Optional[bytes]:
    try:
        with zf.open(name) as fh:
            return fh.read(limit)
    except Exception:
        return None


def _parse(data: Optional[bytes]):
    if not data or b"<!DOCTYPE" in data[:2048] or b"<!ENTITY" in data[:4096]:
        return None
    try:
        return ET.fromstring(data)
    except ET.ParseError:
        return None


def _sheets(zf) -> List[tuple]:
    """[(name, part path)] in workbook order."""
    wb = _parse(_read(zf, "xl/workbook.xml"))
    rels = _parse(_read(zf, "xl/_rels/workbook.xml.rels"))
    if wb is None:
        return []
    targets: Dict[str, str] = {}
    if rels is not None:
        for r in rels.iter(_PKG_REL + "Relationship"):
            t = r.get("Target", "")
            t = t.lstrip("/")
            if not t.startswith("xl/"):
                t = "xl/" + t
            targets[r.get("Id", "")] = t
    out = []
    for i, s in enumerate(wb.iter(_NS + "sheet")):
        rid = s.get(_REL_NS + "id", "")
        out.append((s.get("name", ""), targets.get(rid, "xl/worksheets/sheet%d.xml" % (i + 1))))
    return out


def sheet_names(zf) -> Optional[List[str]]:
    s = _sheets(zf)
    return [n for n, _ in s] if s else None


def _shared_strings(zf) -> List[str]:
    root = _parse(_read(zf, "xl/sharedStrings.xml"))
    if root is None:
        return []
    out: List[str] = []
    for si in root.iter(_NS + "si"):
        out.append("".join(t.text or "" for t in si.iter(_NS + "t")))
        if len(out) >= _MAX_SHARED:
            break
    return out


_COL = re.compile(r"([A-Z]+)")


def _col_index(ref: str) -> int:
    m = _COL.match(ref or "")
    if not m:
        return -1
    n = 0
    for ch in m.group(1):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_rows(zf, sheet: Optional[str], max_rows: int) -> Optional[List[List[str]]]:
    sheets = _sheets(zf)
    if not sheets:
        return None
    part = None
    if sheet is None:
        part = sheets[0][1]
    else:
        for n, p in sheets:
            if n == sheet:
                part = p
                break
    if part is None:
        return None
    root = _parse(_read(zf, part))
    if root is None:
        return None
    shared = None
    rows: List[List[str]] = []
    for row in root.iter(_NS + "row"):
        vals: List[str] = []
        for c in row.iter(_NS + "c"):
            idx = _col_index(c.get("r", ""))
            if idx < 0:
                idx = len(vals)
            if idx > 16384:
                continue
            t = c.get("t")
            v = c.find(_NS + "v")
            if t == "s" and v is not None:
                if shared is None:
                    shared = _shared_strings(zf)
                try:
                    text = shared[int(v.text or "")]
                except (ValueError, IndexError):
                    text = ""
            elif t == "inlineStr":
                text = "".join(x.text or "" for x in c.iter(_NS + "t"))
            else:
                text = (v.text or "") if v is not None else ""
            while len(vals) < idx:
                vals.append("")
            vals.append(text)
        rows.append(vals)
        if len(rows) >= max_rows:
            break
    return rows
