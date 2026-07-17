"""Decode a VCF (Variant Call Format) header.

VCF is text: '##' meta lines then a '#CHROM ... FORMAT sample1 sample2' column
line. This reads the fileformat version, the sample names, the reference, and
counts of the meta lines from the header, without reading variant records. `.gz`
(bgzip is gzip-compatible for the leading block) is handled by the reader.
"""

from __future__ import annotations

import re
from typing import Optional

from .decoders import DecodedHeader, Read, register_decoder

_FILEFORMAT = re.compile(rb"^##fileformat=(\S+)", re.MULTILINE)
_REFERENCE = re.compile(rb"^##reference=(\S+)", re.MULTILINE)


def decode_vcf(data: bytes) -> Optional[DecodedHeader]:
    if not data.startswith(b"##fileformat=VCF"):
        return None

    m = _FILEFORMAT.search(data)
    version = m.group(1).decode("ascii", "replace") if m else None
    ref = _REFERENCE.search(data)
    reference = ref.group(1).decode("utf-8", "replace") if ref else None

    samples = []
    for line in data.split(b"\n"):
        if line.startswith(b"#CHROM"):
            cols = line.rstrip(b"\r").split(b"\t")
            if len(cols) > 9:  # columns 10+ (after FORMAT) are sample names
                samples = [c.decode("utf-8", "replace") for c in cols[9:]]
            break

    fields = {
        "version": version,
        "numSamples": len(samples),
        "samples": samples[:64],
        "numContigs": data.count(b"\n##contig") + data.startswith(b"##contig"),
        "numInfoFields": data.count(b"\n##INFO"),
        "numFilterFields": data.count(b"\n##FILTER"),
        "numFormatFields": data.count(b"\n##FORMAT"),
        "reference": reference,
    }
    return DecodedHeader(
        format="vcf",
        summary="VCF %s, %d samples, %d contigs"
        % (version or "?", len(samples), fields["numContigs"]),
        fields=fields,
    )


register_decoder("vcf", decode_vcf, read=Read(leading=256 * 1024))  # header can carry many ##contig lines
