"""Decode an mmCIF / CIF structure file header.

mmCIF is text organized as a data block of _category.item values. This reads the
block code and the single-value header items (entry id, title, keywords,
experiment method, resolution, deposition date) without reading the atom loop.
"""

from __future__ import annotations

import re
from typing import Optional

from .decoders import DecodedHeader, _finite, register_decoder


def _unquote(v: bytes) -> Optional[str]:
    v = v.strip()
    if len(v) >= 2 and ((v[:1] == b"'" and v[-1:] == b"'") or (v[:1] == b'"' and v[-1:] == b'"')):
        v = v[1:-1]
    return v.decode("utf-8", "replace").strip() or None


def _item(data: bytes, name: bytes) -> Optional[str]:
    # value on the same line: `_category.item   value`
    m = re.search(rb"(?m)^" + re.escape(name) + rb"[ \t]+(\S.*?)[ \t]*$", data)
    if m:
        return _unquote(m.group(1))
    # value in a following ';'-delimited text block
    m = re.search(rb"(?m)^" + re.escape(name) + rb"[ \t]*$", data)
    if m:
        tail = data[m.end():].lstrip(b"\r\n")
        if tail.startswith(b";"):
            stop = tail[1:].find(b"\n;")
            if stop != -1:
                return tail[1: 1 + stop].strip().decode("utf-8", "replace") or None
    return None


def decode_mmcif(data: bytes) -> Optional[DecodedHeader]:
    m = re.match(rb"\s*data_(\S+)", data)
    if not m:
        return None
    block = m.group(1).decode("utf-8", "replace")

    resolution = None
    for item in (b"_refine.ls_d_res_high", b"_reflns.d_resolution_high",
                 b"_em_3d_reconstruction.resolution"):
        rv = _item(data, item)
        if rv:
            try:
                resolution = _finite(float(rv))
                break
            except ValueError:
                pass

    classification = _item(data, b"_struct_keywords.pdbx_keywords")
    method = _item(data, b"_exptl.method")
    fields = {
        "entryId": _item(data, b"_entry.id") or block,
        "title": _item(data, b"_struct.title"),
        "classification": classification,
        "experimentMethod": method,
        "depositionDate": _item(data, b"_pdbx_database_status.recvd_initial_deposition_date"),
        "resolutionAngstrom": resolution,
    }
    res = ", %.2f A" % resolution if resolution else ""
    return DecodedHeader(
        format="mmcif",
        summary="mmCIF %s, %s%s" % (fields["entryId"] or "?", classification or method or "structure", res),
        fields=fields,
    )


register_decoder(["cif", "mmcif"], decode_mmcif)
