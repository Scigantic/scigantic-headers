"""SAM header decoding."""

import json

from scigantic_headers import decode_sam
from scigantic_headers.decoders import decode_bytes

SAM = (
    b"@HD\tVN:1.6\tSO:coordinate\n"
    b"@SQ\tSN:chr1\tLN:248956422\n"
    b"@SQ\tSN:chr2\tLN:242193529\n"
    b"@RG\tID:sample1\tSM:s1\n"
    b"@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n"
    b"read1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tACGT\tIIII\n"
)


def test_core_fields():
    h = decode_sam(SAM)
    assert h.format == "sam"
    assert h.fields["version"] == "1.6"
    assert h.fields["sortOrder"] == "coordinate"
    assert h.fields["numReferences"] == 2
    assert h.fields["referenceNames"] == ["chr1", "chr2"]
    assert h.fields["numReadGroups"] == 1
    assert h.fields["numPrograms"] == 1


def test_header_starting_with_sq():
    s = b"@SQ\tSN:c1\tLN:100\nread\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
    h = decode_sam(s)
    assert h.fields["numReferences"] == 1
    assert h.fields["version"] is None
    assert h.fields["sortOrder"] is None


def test_rejects_non_sam():
    assert decode_sam(b"not a sam file") is None
    assert decode_sam(b"read1\t0\tchr1\t100\n") is None  # alignment, no @ header


def test_dispatch():
    assert decode_bytes("a.sam", SAM).format == "sam"


def test_json_safe():
    json.dumps(decode_sam(SAM).to_dict(), allow_nan=False)
