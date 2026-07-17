"""GFF3 / GTF annotation header decoding."""

import json

from scigantic_headers import decode_gff
from scigantic_headers.decoders import decode_bytes

GFF3 = (
    b"##gff-version 3\n"
    b"##sequence-region chr1 1 248956422\n"
    b"chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=gene1\n"
    b"chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tID=exon1\n"
    b"chr1\tENSEMBL\tCDS\t12010\t12057\t.\t+\t0\tID=cds1\n"
)


def test_gff3_core():
    h = decode_gff(GFF3)
    assert h.format == "gff"
    assert h.fields["version"] == "3"
    assert h.fields["sequenceRegionsInHeader"] == 1
    assert set(h.fields["featureTypes"]) == {"gene", "exon", "CDS"}
    assert set(h.fields["sources"]) == {"HAVANA", "ENSEMBL"}


def test_gtf_without_pragma():
    gtf = (
        b'chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id "g1";\n'
        b'chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "g1";\n'
    )
    h = decode_gff(gtf)
    assert h.fields["version"] is None
    assert set(h.fields["featureTypes"]) == {"gene", "transcript"}


def test_rejects_non_gff():
    assert decode_gff(b"just some text\nwith no tab columns\n") is None
    assert decode_gff(b"") is None


def test_dispatch():
    assert decode_bytes("a.gff3", GFF3).format == "gff"
    assert decode_bytes("a.gtf", GFF3).format == "gff"


def test_json_safe():
    json.dumps(decode_gff(GFF3).to_dict(), allow_nan=False)
