"""Decode a FASTQ file's first record.

FASTQ has no file header; the first record's name line carries the metadata. For
Illumina reads the name is

    @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filter>:<control>:<index>

so instrument, run, flowcell, and lane come from the first line, and the read
length from the first sequence line. `.gz` is handled by the reader.
"""

from __future__ import annotations

from typing import Optional

from .decoders import DecodedHeader, register_decoder


def decode_fastq(data: bytes) -> Optional[DecodedHeader]:
    if data[:1] != b"@":
        return None
    lines = data.split(b"\n")
    if len(lines) < 3 or not lines[2].startswith(b"+"):
        return None  # not a FASTQ record (line 1 '@name', line 3 '+')

    ident = lines[0][1:].split(b" ", 1)[0]  # drop '@', take up to first space
    read_length = len(lines[1])  # sequence is complete: a '+' line follows it

    platform = instrument = run = flowcell = None
    lane: Optional[int] = None
    parts = ident.split(b":")
    if len(parts) >= 7:  # the Illumina read-name convention
        platform = "Illumina"
        instrument = parts[0].decode("ascii", "replace")
        run = parts[1].decode("ascii", "replace")
        flowcell = parts[2].decode("ascii", "replace")
        try:
            lane = int(parts[3])
        except ValueError:
            lane = None

    fields = {
        "platform": platform,
        "instrument": instrument,
        "runNumber": run,
        "flowcellId": flowcell,
        "lane": lane,
        "firstReadLength": read_length,
        "readName": ident.decode("ascii", "replace"),
    }
    return DecodedHeader(
        format="fastq",
        summary="%s reads, first read %d bp" % (platform or "FASTQ", read_length),
        fields=fields,
    )


register_decoder(["fastq", "fq"], decode_fastq)
