"""Decode a GenBank flat-file header.

A GenBank record begins with a LOCUS line (name, length, molecule type, topology,
division, date) followed by DEFINITION, ACCESSION, VERSION, and SOURCE / ORGANISM,
before the FEATURES table and the ORIGIN sequence. This reads that header block.
"""

from __future__ import annotations

import re
from typing import Optional

from .decoders import DecodedHeader, register_decoder

_UNITS = (b"bp", b"aa", b"rc")


def _block(data: bytes, keyword: bytes) -> Optional[str]:
    """A keyword's value at column 0, joined across its indented continuation
    lines (which stop at the next column-0 keyword)."""
    m = re.search(rb"(?ms)^" + keyword + rb"\s+(.*?)(?=^\S)", data)
    if not m:
        m = re.search(rb"(?ms)^" + keyword + rb"\s+(.*)$", data)
    if not m:
        return None
    text = b" ".join(line.strip() for line in m.group(1).splitlines())
    return text.decode("utf-8", "replace").strip() or None


def decode_genbank(data: bytes) -> Optional[DecodedHeader]:
    if not data.startswith(b"LOCUS"):
        return None

    toks = data.split(b"\n", 1)[0].split()
    name = toks[1].decode("utf-8", "replace") if len(toks) > 1 else None
    length = units = moltype = topology = None
    for k, t in enumerate(toks):
        if t in _UNITS:
            units = t.decode()
            if k >= 1 and toks[k - 1].isdigit():
                length = int(toks[k - 1])
            # The token after bp/aa is the molecule type, unless it is the
            # topology (proteins have no molecule type on the LOCUS line).
            if k + 1 < len(toks) and toks[k + 1] not in (b"linear", b"circular"):
                moltype = toks[k + 1].decode("utf-8", "replace")
            break
    for t in toks:
        if t in (b"linear", b"circular"):
            topology = t.decode()

    organism = None
    om = re.search(rb"(?m)^\s+ORGANISM\s+([^\n]+)", data)
    if om:
        organism = om.group(1).decode("utf-8", "replace").strip() or None

    fields = {
        "locus": name,
        "length": length,
        "lengthUnits": units,
        "moleculeType": moltype,
        "topology": topology,
        "definition": _block(data, b"DEFINITION"),
        "accession": _block(data, b"ACCESSION"),
        "version": _block(data, b"VERSION"),
        "organism": organism,
    }
    size = "%d %s" % (length, units) if length and units else "?"
    summary = "GenBank %s, %s %s" % (name or "?", size, moltype or "")
    return DecodedHeader(format="genbank", summary=summary.strip(), fields=fields)


register_decoder(["gb", "gbk", "genbank", "gbff"], decode_genbank)
