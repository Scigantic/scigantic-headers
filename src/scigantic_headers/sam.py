"""Decode a SAM alignment header.

SAM is text; the header is '@'-prefixed lines before the alignments: @HD (format
version, sort order), @SQ (one per reference sequence), @RG (read groups), @PG
(programs). This reads those, not the alignments. Binary BAM is not handled (its
BGZF framing is out of scope for a header read).
"""

from __future__ import annotations

import re
from typing import Optional

from .decoders import DecodedHeader, register_decoder

_HD_LINE = re.compile(rb"^@HD\t[^\n]*", re.MULTILINE)
_VN = re.compile(rb"\bVN:([^\t\n]+)")
_SO = re.compile(rb"\bSO:([^\t\n]+)")
_SN = re.compile(rb"\bSN:([^\t\n]+)")


def _count(data: bytes, record: bytes) -> int:
    return data.count(b"\n" + record) + (1 if data.startswith(record) else 0)


def decode_sam(data: bytes) -> Optional[DecodedHeader]:
    # A SAM header begins with an @-record (@HD is recommended, @SQ is common).
    if not (data.startswith(b"@HD\t") or data.startswith(b"@SQ\t")):
        return None

    version = sort_order = None
    hd = _HD_LINE.search(data)
    if hd:
        vn = _VN.search(hd.group(0))
        so = _SO.search(hd.group(0))
        version = vn.group(1).decode("ascii", "replace") if vn else None
        sort_order = so.group(1).decode("ascii", "replace") if so else None

    refs = []
    for line in data.split(b"\n"):
        if line.startswith(b"@SQ\t"):
            sn = _SN.search(line)
            if sn:
                refs.append(sn.group(1).decode("utf-8", "replace"))
            if len(refs) >= 64:
                break

    fields = {
        "version": version,
        "sortOrder": sort_order,
        "numReferences": _count(data, b"@SQ\t"),
        "referenceNames": refs,
        "numReadGroups": _count(data, b"@RG\t"),
        "numPrograms": _count(data, b"@PG\t"),
    }
    return DecodedHeader(
        format="sam",
        summary="SAM %s, %d references, sort %s"
        % (version or "?", fields["numReferences"], sort_order or "unsorted"),
        fields=fields,
    )


register_decoder("sam", decode_sam)
