"""FASTQ first-record decoding."""

import gzip
import json

from scigantic_headers import decode_fastq, decode_file
from scigantic_headers.decoders import decode_bytes

ILLUMINA = (
    b"@M00123:7:000000000-ABCDE:1:1101:15589:1332 1:N:0:1\n"
    b"ACGTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIII\n"
)


def test_illumina_read_name():
    h = decode_fastq(ILLUMINA)
    assert h.format == "fastq"
    assert h.fields["platform"] == "Illumina"
    assert h.fields["instrument"] == "M00123"
    assert h.fields["runNumber"] == "7"
    assert h.fields["flowcellId"] == "000000000-ABCDE"
    assert h.fields["lane"] == 1
    assert h.fields["firstReadLength"] == 18


def test_non_illumina_name_still_decodes():
    fq = b"@read1 a description\nACGT\n+\n!!!!\n"
    h = decode_fastq(fq)
    assert h.fields["platform"] is None
    assert h.fields["firstReadLength"] == 4
    assert h.fields["readName"] == "read1"


def test_rejects_non_fastq():
    assert decode_fastq(b">seqid\nACGT\n") is None  # FASTA, not FASTQ
    assert decode_fastq(b"just some random bytes here") is None
    assert decode_fastq(b"@only one line") is None
    assert decode_fastq(b"@name\nACGT\nNOTPLUS\n") is None  # line 3 must start '+'


def test_dispatch_and_gz(tmp_path):
    assert decode_bytes("r.fastq", ILLUMINA).format == "fastq"
    p = tmp_path / "r.fq.gz"
    with gzip.open(p, "wb") as f:
        f.write(ILLUMINA)
    h = decode_file(str(p))
    assert h is not None
    assert h.fields["instrument"] == "M00123"


def test_json_safe():
    json.dumps(decode_fastq(ILLUMINA).to_dict(), allow_nan=False)
