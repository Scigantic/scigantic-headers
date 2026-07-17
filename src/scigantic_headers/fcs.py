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
A doubled delimiter is an escaped literal delimiter inside a value. The keywords
this reads: $PAR (parameters), $TOT (events), $DATATYPE, $MODE, $CYT (cytometer),
and $PnN / $PnS (per-channel short and long names).

The TEXT segment sits near the start but can run well past 1 KiB, so this
decoder registers a larger leading read (see the register_decoder call below).
"""

from __future__ import annotations

from typing import List, Optional

from .decoders import DecodedHeader, Read, register_decoder

_FCS_VERSIONS = (b"FCS2.0", b"FCS3.0", b"FCS3.1")
_HEADER_MIN = 58  # version + reserved + six 8-byte offset fields
_DATATYPE = {"I": "int", "F": "float32", "D": "float64", "A": "ascii"}


def _offset(data: bytes, start: int) -> int:
    """One 8-byte ASCII offset field. -1 if it is not a non-negative integer."""
    s = data[start:start + 8].strip()
    if not s:
        return 0
    try:
        return int(s)
    except ValueError:
        return -1


def _split(body: bytes, delim: int) -> List[bytes]:
    """Split the TEXT body on the delimiter, treating a doubled delimiter as an
    escaped literal delimiter within a value (the FCS convention)."""
    out: List[bytes] = []
    cur = bytearray()
    i, n = 0, len(body)
    while i < n:
        if body[i] == delim:
            if i + 1 < n and body[i + 1] == delim:  # escaped literal delimiter
                cur.append(delim)
                i += 2
                continue
            out.append(bytes(cur))
            cur = bytearray()
            i += 1
        else:
            cur.append(body[i])
            i += 1
    if cur:
        out.append(bytes(cur))
    return out


def decode_fcs(data: bytes) -> Optional[DecodedHeader]:
    """Decode an FCS file from its leading bytes (version, HEADER offsets, and the
    TEXT segment). Returns None if the magic is wrong or the TEXT segment is not
    fully present in `data` (read more leading bytes and retry)."""
    if len(data) < _HEADER_MIN or data[:6] not in _FCS_VERSIONS:
        return None
    version = data[:6].decode("ascii")

    text_begin = _offset(data, 10)
    text_end = _offset(data, 18)
    if text_begin <= 0 or text_end <= text_begin:
        return None
    if text_end + 1 > len(data):  # TEXT runs past what we were handed
        return None

    body = data[text_begin + 1:text_end + 1]  # skip the leading delimiter byte
    delim = data[text_begin]
    parts = _split(body, delim)

    kv = {}
    it = iter(parts)
    for key in it:
        try:
            val = next(it)
        except StopIteration:
            break
        kv[key.decode("ascii", "replace").strip().upper()] = val.decode("utf-8", "replace")

    def as_int(key: str) -> Optional[int]:
        v = kv.get(key)
        if v is None:
            return None
        try:
            return int(v.strip())
        except ValueError:
            return None

    npar = as_int("$PAR")
    ntot = as_int("$TOT")

    channels = []
    if npar is not None and 0 < npar <= 100000:
        for i in range(1, npar + 1):
            name = kv.get("$P%dN" % i)
            if name is None:
                continue
            ch = {"name": name}
            long = kv.get("$P%dS" % i)
            if long:
                ch["label"] = long
            channels.append(ch)

    fields = {
        "version": version,
        "numParameters": npar,
        "numEvents": ntot,
        "dataType": _DATATYPE.get(kv.get("$DATATYPE", ""), kv.get("$DATATYPE")),
        "mode": kv.get("$MODE"),
        "channels": channels,
    }
    cyt = kv.get("$CYT")
    if cyt:
        fields["cytometer"] = cyt

    evt = ntot if ntot is not None else "?"
    par = npar if npar is not None else "?"
    return DecodedHeader(
        format="fcs",
        summary="FCS %s, %s events, %s parameters" % (version[3:], evt, par),
        fields=fields,
    )


register_decoder("fcs", decode_fcs, read=Read(leading=256 * 1024))  # TEXT segment runs past 1 KiB
