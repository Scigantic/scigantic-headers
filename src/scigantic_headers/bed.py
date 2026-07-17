"""Decode a BED genome-interval file's structure.

BED has no header (optional 'track' / 'browser' lines). This reports the BED
flavor (column count: BED3, BED6, BED12) and whether a track line is present,
from the first data line, validating that start/end are integers.
"""

from __future__ import annotations

from typing import Optional

from .decoders import DecodedHeader, register_decoder

_BED_COLUMNS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts",
]


def decode_bed(data: bytes) -> Optional[DecodedHeader]:
    has_track = False
    for line in data.split(b"\n"):
        s = line.strip()
        if not s:
            continue
        if s.startswith((b"track", b"browser", b"#")):
            has_track = has_track or s.startswith(b"track")
            continue
        cols = line.rstrip(b"\r").split(b"\t")
        if len(cols) < 3 or not (cols[1].isdigit() and cols[2].isdigit()):
            return None  # not BED: needs chrom + integer start/end
        fields = {
            "numColumns": len(cols),
            "bedType": "BED%d" % len(cols),
            "columns": _BED_COLUMNS[: min(len(cols), len(_BED_COLUMNS))],
            "hasTrackLine": has_track,
            "firstChrom": cols[0].decode("utf-8", "replace"),
        }
        return DecodedHeader(format="bed", summary="BED%d interval file" % len(cols), fields=fields)
    return None


register_decoder("bed", decode_bed)
