"""Decode an FCS (Flow Cytometry Standard) file's header and TEXT segment.

An FCS file opens with a HEADER: a 6-byte version ('FCS3.1'), four spaces, then
six 8-byte ASCII integers giving the byte offsets of the TEXT, DATA, and ANALYSIS
segments:

    0   'FCS3.1'    version (6 bytes)
    6   '    '      reserved (4 spaces)
    10  TEXT begin  8-byte ASCII int, right-justified
    18  TEXT end    8-byte ASCII int (inclusive offset of the last TEXT byte)
    26  DATA begin  ...
    34  DATA end
    42  ANALYSIS begin
    50  ANALYSIS end

The TEXT segment is a run of delimiter-separated keyword/value pairs. Its first
byte is the delimiter; the rest is `$KEY <delim> value <delim> $KEY <delim> ...`.
A doubled delimiter inside a value is an escaped literal delimiter (FCS 3.1
section 3.2.7). The keywords this reads: $PAR (parameters), $TOT (events),
$DATATYPE, $MODE, $CYT (cytometer), $NEXTDATA, $BEGINSTEXT/$ENDSTEXT, and
$PnN / $PnS (per-channel short and long names).

Legacy writers (BD CellQuest and other FCS 2.0 / 3.0 vendors) write an empty
value as a doubled delimiter: `KEY\\\\NEXTKEY\\value\\`. Read strictly, that
doubled delimiter is a literal, KEY's value swallows NEXTKEY, and every later
pair swaps keyword and value. The two readings can be told apart by position:
a keyword or value may not begin with the delimiter (section 3.2.7), so a
doubled delimiter that directly follows a keyword cannot start its value and
must be the empty value. Inside a value it is a literal. See `_split`.

The TEXT segment sits near the start but can run well past 1 KiB, so this
decoder registers a larger leading read (see the register_decoder call below).
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional

from .decoders import DecodedHeader, Read, register_decoder

_FCS_VERSIONS = (b"FCS2.0", b"FCS3.0", b"FCS3.1", b"FCS3.2")
_HEADER_MIN = 58  # version + reserved + six 8-byte offset fields
_DATATYPE = {"I": "int", "F": "float32", "D": "float64", "A": "ascii"}
_MAX_PAR = 100000
_PNN = re.compile(r"^\$P(\d+)N$")


def _offset(data: bytes, start: int) -> int:
    """One 8-byte ASCII offset field. -1 if it holds no non-negative integer.
    Blank means 0 (the spec's value for a segment that is absent or too far
    into the file for the HEADER). A stray trailing byte after the digits, seen
    in some LMD exports, is ignored."""
    s = data[start:start + 8].strip()
    if not s:
        return 0
    try:
        return int(s)
    except ValueError:
        m = re.match(rb"\d+", s)
        return int(m.group()) if m else -1


def _split(body: bytes, delim: int) -> List[bytes]:
    """Split a TEXT body (the bytes after the leading delimiter) into alternating
    keyword and value fields.

    Delimiters come in runs. A run that follows a value is `run // 2` literal
    delimiters, plus one separator if the run is odd (spec escaping: a doubled
    delimiter is a literal). A run of two or more that follows a keyword is the
    keyword's separator and an empty value: a value may not begin with the
    delimiter, so it cannot be a literal there. That reads spec-escaped values
    such as `$P3F/488//10/` (value `488/10`) and legacy empty values such as
    `KEY\\\\NEXT\\v\\` (KEY empty, NEXT = v) from the same text. A keyword that
    itself contains an escaped delimiter (a `CD GAG/100` label written with `/`
    as the delimiter) reads as two keywords around an empty value; that misreads
    only that pair, where the strict reading of a legacy empty value shifts
    every pair after it."""
    d = bytes((delim,))
    pieces = body.lstrip(d).split(d)  # a leading delimiter run is junk: no field may start with one
    out: List[bytes] = []
    cur = bytearray(pieces[0])
    npieces = len(pieces)
    k = 1
    while k < npieces:
        j = k
        while j < npieces and not pieces[j]:
            j += 1
        if j < npieces:
            run = j - k + 1  # delimiters between cur and pieces[j]
            nxt = pieces[j]
        else:
            run = j - k  # delimiters after cur at the end of the body
            nxt = b""
        while run:
            if len(out) % 2 == 0:  # cur is a keyword
                out.append(bytes(cur))
                cur = bytearray()
                if run >= 2:  # keyword, then an empty value (legacy writers)
                    out.append(b"")
                    run -= 2
                else:
                    run -= 1
            else:  # cur is a value: pairs are literals, an odd one ends it
                cur.extend(d * (run // 2))
                if run % 2:
                    out.append(bytes(cur))
                    cur = bytearray()
                run = 0
        cur.extend(nxt)
        k = j + 1
    if cur:
        out.append(bytes(cur))
    return out


def _text(v: bytes) -> str:
    """A keyword value as text. FCS 3.1 values are UTF-8; older files from
    Windows instruments carry Latin-1 bytes (CellQuest's '\\xaa'), which are
    kept rather than replaced."""
    try:
        return v.decode("utf-8")
    except UnicodeDecodeError:
        return v.decode("latin-1")


def _pairs(body: bytes, delim: int, *, drop_last: bool = False) -> Dict[str, str]:
    """Parse a TEXT body into a dict. Keys are upper-cased and stripped, empty
    keys are skipped, and a repeated key keeps its last value (as flowio and
    fcsparser do). drop_last discards the final field, which may have been cut
    by a short read."""
    parts = _split(body, delim)
    if drop_last and parts:
        parts = parts[:-1]
    kv: Dict[str, str] = {}
    for i in range(0, len(parts) - 1, 2):
        key = parts[i].decode("ascii", "replace").strip().upper()
        if key:
            kv[key] = _text(parts[i + 1])
    return kv


def _as_int(kv: Dict[str, str], key: str) -> Optional[int]:
    v = kv.get(key)
    if v is None:
        return None
    try:
        return int(v.strip())
    except ValueError:
        return None


def decode_fcs(data: bytes) -> Optional[DecodedHeader]:
    """Decode an FCS file from its leading bytes (version, HEADER offsets, and the
    TEXT segment). Returns None if the magic is wrong. A TEXT segment longer than
    `data` (many parameters, a large spillover matrix) is decoded as far as it
    goes: the last part may be cut mid-value, so it is dropped, and
    textTruncated says keywords may be missing. A supplemental TEXT segment
    that lies within `data` is read too; primary TEXT keywords win."""
    if len(data) < _HEADER_MIN or data[:6] not in _FCS_VERSIONS:
        return None
    version = data[:6].decode("ascii")

    text_begin = _offset(data, 10)
    text_end = _offset(data, 18)
    if text_begin <= 0 or text_end <= text_begin:
        return None
    truncated = False
    if text_end + 1 > len(data):  # TEXT runs past what we were handed
        if text_begin + 1 >= len(data):
            return None
        truncated = True
        text_end = len(data) - 1

    delim = data[text_begin]
    kv = _pairs(data[text_begin + 1:text_end + 1], delim, drop_last=truncated)

    # Supplemental TEXT: optional keywords only ($CYT, $PnS, spillover...), same
    # delimiter (3.2.5). Some writers start it with the delimiter, some do not.
    stext_begin = _as_int(kv, "$BEGINSTEXT") or 0
    stext_end = _as_int(kv, "$ENDSTEXT") or 0
    if (
        0 < stext_begin < stext_end < len(data)
        and (stext_begin, stext_end) != (text_begin, text_end)
        and not truncated
    ):
        seg = data[stext_begin:stext_end + 1]
        if seg[:1] == bytes((delim,)):
            seg = seg[1:]
        for key, val in _pairs(seg, delim).items():
            kv.setdefault(key, val)  # primary TEXT wins

    npar = _as_int(kv, "$PAR")
    ntot = _as_int(kv, "$TOT")

    # Channel count: $PAR, or when it is missing or unusable, the highest $PnN
    # index present. A channel without $PnN (seen in FCS 2.0) takes its $PnS.
    count = npar if npar is not None and 0 < npar <= _MAX_PAR else 0
    if count == 0:
        for key in kv:
            m = _PNN.match(key)
            if m:
                count = max(count, min(int(m.group(1)), _MAX_PAR))
    channels = []
    for i in range(1, count + 1):
        name = kv.get("$P%dN" % i)
        long = kv.get("$P%dS" % i)
        if name is None:
            if not long:
                continue
            name, long = long, None
        ch = {"name": name}
        if long:
            ch["label"] = long
        channels.append(ch)

    datatype = (kv.get("$DATATYPE") or "").strip().upper()
    fields = {
        "version": version,
        "numParameters": npar,
        "numEvents": ntot,
        "dataType": _DATATYPE.get(datatype, kv.get("$DATATYPE")),
        "mode": kv.get("$MODE"),
        "channels": channels,
    }
    cyt = (kv.get("$CYT") or "").strip()
    if cyt:
        fields["cytometer"] = cyt
    next_data = _as_int(kv, "$NEXTDATA")
    if next_data:
        fields["multipleDatasets"] = True  # a Coulter LMD holds an FCS 2.0 and an FCS 3.0 copy
    if truncated:
        fields["textTruncated"] = True

    evt = ntot if ntot is not None else "?"
    par = npar if npar is not None else "?"
    return DecodedHeader(
        format="fcs",
        summary="FCS %s, %s events, %s parameters" % (version[3:], evt, par),
        fields=fields,
    )


register_decoder("fcs", decode_fcs, read=Read(leading=1024 * 1024))  # TEXT can be long; a longer one is decoded as far as it goes
