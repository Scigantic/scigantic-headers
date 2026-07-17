"""VCF header decoding."""

import gzip
import json

from scigantic_headers import decode_file, decode_vcf
from scigantic_headers.decoders import decode_bytes

VCF = (
    b"##fileformat=VCFv4.2\n"
    b"##reference=GRCh38\n"
    b"##contig=<ID=chr1,length=248956422>\n"
    b"##contig=<ID=chr2,length=242193529>\n"
    b'##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">\n'
    b'##FILTER=<ID=q10,Description="q10">\n'
    b'##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n'
    b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA001\tNA002\n"
    b"chr1\t100\t.\tA\tT\t50\tPASS\tDP=30\tGT\t0/1\t1/1\n"
)


def test_core_fields():
    h = decode_vcf(VCF)
    assert h.format == "vcf"
    assert h.fields["version"] == "VCFv4.2"
    assert h.fields["reference"] == "GRCh38"
    assert h.fields["numSamples"] == 2
    assert h.fields["samples"] == ["NA001", "NA002"]
    assert h.fields["numContigs"] == 2
    assert h.fields["numInfoFields"] == 1
    assert h.fields["numFilterFields"] == 1
    assert h.fields["numFormatFields"] == 1


def test_sites_only_no_samples():
    v = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    h = decode_vcf(v)
    assert h.fields["version"] == "VCFv4.3"
    assert h.fields["numSamples"] == 0
    assert h.fields["samples"] == []


def test_rejects_non_vcf():
    assert decode_vcf(b"not a vcf file") is None
    assert decode_vcf(b"##something=1\n") is None
    assert decode_vcf(b"") is None


def test_dispatch_and_gz(tmp_path):
    assert decode_bytes("v.vcf", VCF).format == "vcf"
    p = tmp_path / "v.vcf.gz"
    with gzip.open(p, "wb") as f:
        f.write(VCF)
    h = decode_file(str(p))
    assert h is not None
    assert h.fields["numSamples"] == 2


def test_json_safe():
    json.dumps(decode_vcf(VCF).to_dict(), allow_nan=False)
