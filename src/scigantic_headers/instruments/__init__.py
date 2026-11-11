"""Lab instrument exports: plate readers, qPCR and dPCR, cell counters,
solution analyzers, chromatography, spectrophotometers, binding analyzers,
flow cytometry, multi-analyte profiling and liquid handlers.

These files are the text, XML, JSON or spreadsheet exports instrument software
writes (Gen5, SoftMax Pro, QuantStudio, Vi-CELL, OpenLab CDS, ...). Each family
module registers content sniffers with `register_instrument_sniffer`; the
registry tries them before the extension's ordinary decoder (see decoders.py).

Every sniffer returns a DecodedHeader whose `fields` use the shared keys below
where they apply, plus any instrument-specific keys:

    instrumentSoftware   the software that wrote the file, as listed in
                         INSTRUMENTS (e.g. "Agilent Gen5")
    instrumentCategory   e.g. "Plate Reader", "qPCR", "Cell Counting"
    vendor               e.g. "Agilent"
    softwareVersion      version string written in the file
    instrumentModel      instrument model / reader type written in the file
    serialNumber         instrument serial number
    detectionModes       list, e.g. ["Absorbance"]
    plateFormat          wells per plate (int), when stated or derivable
    plateCount           number of plates in the file
    sampleCount          number of samples / wells with data
    measurementCount     number of measurement rows or reads
    wavelengthsNm        list of wavelengths in nm
    runDate              ISO 8601 date-time when stated, else the file's string
    exportFormat         the export variant, e.g. "results txt", "xlsx"
    allotropyParser      the matching allotropy parser module name, so a user
                         knows `allotropy` can convert the file to ASM

A sniffer must be total (never raise, never NaN), cheap (leading bytes for text
exports; the zip directory plus one or two sheet parts for .xlsx), and strict:
it returns None unless the file is unmistakably that instrument's export.
"""

from __future__ import annotations

import io
import zipfile
from typing import Dict, List, Optional

from ..decoders import DecodedHeader, Source, _finite, register_instrument_sniffer  # noqa: F401

# Leading bytes a text export sniffer looks at. Enough for the header block and
# the first data section of every supported export.
TEXT_HEAD_BYTES = 256 * 1024

# Software name -> (category, vendor, allotropy parser, detection modes).
INSTRUMENTS: Dict[str, tuple] = {}


def instrument(name: str, category: str, vendor: str, parser: str, modes) -> None:
    INSTRUMENTS[name] = (category, vendor, parser, list(modes))


def base_fields(name: str) -> Dict[str, object]:
    """The shared identity fields for a supported instrument software."""
    category, vendor, parser, modes = INSTRUMENTS[name]
    return {
        "instrumentSoftware": name,
        "instrumentCategory": category,
        "vendor": vendor,
        "detectionModes": list(modes),
        "allotropyParser": parser,
    }


def header(name: str, summary: str, fields: Dict[str, object]) -> DecodedHeader:
    """Build the DecodedHeader for an instrument export. Drops None values and
    non-finite floats so the record stays JSON-safe."""
    out: Dict[str, object] = {}
    for k, v in {**base_fields(name), **fields}.items():
        if v is None:
            continue
        if isinstance(v, float):
            v = _finite(v)
            if v is None:
                continue
        out[k] = v
    return DecodedHeader(format=name, summary=summary, fields=out)


def head(src: Source, n: int = TEXT_HEAD_BYTES) -> bytes:
    """The leading bytes of the file, read once per Source and shared by every
    sniffer tried on it."""
    cache = getattr(src, "_instrument_head", None)
    if cache is None or (len(cache) < n and len(cache) == getattr(src, "_instrument_head_n", 0)):
        cache = src.read_at(0, n)
        src._instrument_head = cache  # type: ignore[attr-defined]
        src._instrument_head_n = n  # type: ignore[attr-defined]
    return cache[:n]


def text_of(data: bytes) -> Optional[str]:
    """Decode exported text. Handles UTF-8 (with or without BOM), UTF-16 LE/BE
    (with BOM, or without one when every other byte is NUL, as QuantStudio and
    Vi-CELL XR write), and falls back to cp1252/latin-1. A cut multibyte
    character at the end of a leading read is dropped, not an error."""
    if not data:
        return None
    try:
        if data.startswith(b"\xef\xbb\xbf"):
            return data[3:].decode("utf-8", "ignore")
        if data.startswith((b"\xff\xfe", b"\xfe\xff")):
            return data.decode("utf-16", "ignore")
        sample = data[:4096]
        if len(sample) >= 4:
            odd_nul = sample[1::2].count(0)
            even_nul = sample[0::2].count(0)
            half = len(sample) // 2
            if odd_nul > 0.4 * half and even_nul < 0.05 * half:
                return data[: len(data) // 2 * 2].decode("utf-16-le", "ignore")
            if even_nul > 0.4 * half and odd_nul < 0.05 * half:
                return data[: len(data) // 2 * 2].decode("utf-16-be", "ignore")
        if b"\x00" in sample:
            return None  # binary
        try:
            return data.decode("utf-8")
        except UnicodeDecodeError as e:
            if e.reason == "unexpected end of data":
                return data[: e.start].decode("utf-8")
            return data.decode("cp1252", "replace")
    except Exception:
        return None


def head_text(src: Source, n: int = TEXT_HEAD_BYTES) -> Optional[str]:
    """`text_of(head(src, n))`, cached per Source."""
    key = "_instrument_text_%d" % n
    if not hasattr(src, key):
        setattr(src, key, text_of(head(src, n)))
    return getattr(src, key)


class _SourceFile(io.RawIOBase):
    """A read-only seekable file over a Source, so zipfile can read an .xlsx
    through the same bounded reads as every other decoder."""

    def __init__(self, src: Source):
        self._src = src
        self._pos = 0

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self._pos

    def seek(self, offset: int, whence: int = 0) -> int:
        size = self._src.size or 0
        if whence == 0:
            self._pos = offset
        elif whence == 1:
            self._pos += offset
        else:
            self._pos = size + offset
        return self._pos

    def readinto(self, b) -> int:
        data = self._src.read_at(self._pos, len(b))
        b[: len(data)] = data
        self._pos += len(data)
        return len(data)


def open_zip(src: Source) -> Optional[zipfile.ZipFile]:
    """The .xlsx (or other zip) behind `src`, or None if it is not a zip or its
    size is unknown. Only the parts a sniffer opens are read."""
    try:
        if src.size is None or head(src, 4) != b"PK\x03\x04":
            return None
        return zipfile.ZipFile(io.BufferedReader(_SourceFile(src), 64 * 1024))
    except Exception:
        return None


def zip_member(zf: zipfile.ZipFile, name: str, limit: int = 4 * 1024 * 1024) -> Optional[bytes]:
    """Up to `limit` decompressed bytes of one zip member, or None."""
    try:
        with zf.open(name) as fh:
            return fh.read(limit)
    except Exception:
        return None


def xlsx_rows(src: Source, sheet: Optional[str] = None, max_rows: int = 200) -> Optional[List[List[str]]]:
    """The first `max_rows` rows of one worksheet of an .xlsx as lists of
    strings (shared strings resolved, empty cells as ""), or None. `sheet` is a
    sheet name; None means the first sheet. Cached per Source and sheet."""
    cache = getattr(src, "_instrument_xlsx", None)
    if cache is None:
        cache = {}
        src._instrument_xlsx = cache  # type: ignore[attr-defined]
    k = (sheet, max_rows)
    if k not in cache:
        from ._xlsx import read_rows

        zf = getattr(src, "_instrument_zip", None)
        if zf is None:
            zf = open_zip(src)
            src._instrument_zip = zf  # type: ignore[attr-defined]
        cache[k] = read_rows(zf, sheet, max_rows) if zf else None
    return cache[k]


def xlsx_sheet_names(src: Source) -> Optional[List[str]]:
    """Worksheet names of an .xlsx in workbook order, or None."""
    from ._xlsx import sheet_names

    zf = getattr(src, "_instrument_zip", None)
    if zf is None:
        zf = open_zip(src)
        src._instrument_zip = zf  # type: ignore[attr-defined]
    return sheet_names(zf) if zf else None


# Family modules register their sniffers on import. Each one lives in its own
# file; add new families here.
from . import _families  # noqa: E402,F401
