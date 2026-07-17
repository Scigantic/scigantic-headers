"""Decode a GFF3 / GTF genome-annotation header.

GFF3 starts with '##gff-version 3'; GTF has no version pragma. Both are 9-column
tab-separated feature lines (seqid, source, type, start, end, score, strand,
phase, attributes). This reads the version and a preview of the sources and
feature types from the leading lines.
"""

from __future__ import annotations

from typing import Optional

from .decoders import DecodedHeader, Read, register_decoder


def decode_gff(data: bytes) -> Optional[DecodedHeader]:
    is_gff3 = data.startswith(b"##gff-version")
    version = None
    if is_gff3:
        parts = data.split(b"\n", 1)[0].split()
        version = parts[1].decode("utf-8", "replace") if len(parts) > 1 else "3"

    sources, types = set(), set()
    n_features = n_regions = 0
    for line in data.split(b"\n"):
        if line.startswith(b"##sequence-region"):
            n_regions += 1
            continue
        if line.startswith(b"#") or not line.strip():
            continue
        cols = line.split(b"\t")
        if len(cols) >= 8:
            n_features += 1
            sources.add(cols[1].decode("utf-8", "replace"))
            types.add(cols[2].decode("utf-8", "replace"))
            if len(sources) + len(types) > 400:
                break

    if not is_gff3 and n_features == 0:
        return None  # neither a GFF3 pragma nor recognizable feature lines

    fields = {
        "version": version,
        "sources": sorted(s for s in sources if s)[:32],
        "featureTypes": sorted(t for t in types if t)[:32],
        "sequenceRegionsInHeader": n_regions,
    }
    kind = "GFF3" if is_gff3 else "GFF/GTF"
    return DecodedHeader(
        format="gff",
        summary="%s, %d feature types, %d sources (preview)" % (kind, len(types), len(sources)),
        fields=fields,
    )


register_decoder(["gff", "gff3", "gtf"], decode_gff, read=Read(leading=256 * 1024))  # source + feature-type preview
