"""Decode a PDB structure file header.

PDB is a fixed-column text format. The leading records name the structure: HEADER
(classification, deposition date, PDB id), TITLE, EXPDTA (experiment method), and
REMARK 2 (resolution). This reads those, not the atom records. mmCIF is a
different format and is not handled here.
"""

from __future__ import annotations

import re
from typing import Optional

from .decoders import DecodedHeader, Read, _finite, register_decoder

# Six-character, space-padded record names a PDB file can begin with.
_LEADING = {
    b"HEADER", b"OBSLTE", b"TITLE ", b"SPLIT ", b"CAVEAT", b"COMPND", b"SOURCE",
    b"KEYWDS", b"EXPDTA", b"NUMMDL", b"MDLTYP", b"AUTHOR", b"REVDAT", b"SPRSDE",
    b"REMARK",
}
_RESOLUTION = re.compile(rb"REMARK\s+2\s+RESOLUTION\.\s+([0-9]+\.[0-9]+)\s+ANGSTROM")


def _record(data: bytes, name: bytes, start: int, end: int) -> Optional[str]:
    m = re.search(rb"^" + name + rb"[^\n]*", data, re.MULTILINE)
    if not m:
        return None
    return m.group(0)[start:end].decode("utf-8", "replace").strip() or None


def decode_pdb(data: bytes) -> Optional[DecodedHeader]:
    if data[:6] not in _LEADING:
        return None

    resolution = None
    rm = _RESOLUTION.search(data)
    if rm:
        try:
            resolution = _finite(float(rm.group(1)))
        except ValueError:
            resolution = None

    titles = re.findall(rb"^TITLE [^\n]*", data, re.MULTILINE)
    title = None
    if titles:
        title = " ".join(t[10:80].decode("utf-8", "replace").strip() for t in titles).strip() or None

    pdb_id = _record(data, b"HEADER", 62, 66)
    classification = _record(data, b"HEADER", 10, 50)

    fields = {
        "pdbId": pdb_id,
        "classification": classification,
        "depositionDate": _record(data, b"HEADER", 50, 59),
        "title": title,
        "experimentMethod": _record(data, b"EXPDTA", 10, 79),
        "resolutionAngstrom": resolution,
    }
    res = ", %.2f A" % resolution if resolution else ""
    return DecodedHeader(
        format="pdb",
        summary="PDB %s, %s%s" % (pdb_id or "?", classification or "structure", res),
        fields=fields,
    )


register_decoder(["pdb", "ent"], decode_pdb, read=Read(leading=256 * 1024))  # title records precede the atoms
