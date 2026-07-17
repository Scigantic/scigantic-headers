"""GenBank header decoding."""

import json

from scigantic_headers import decode_genbank
from scigantic_headers.decoders import decode_bytes

GB = b"""LOCUS       SCU49845     5028 bp    DNA     linear   PLN 21-JUN-1999
DEFINITION  Saccharomyces cerevisiae TCP1-beta gene, partial cds; and Axl2p
            (AXL2) and Rev7p (REV7) genes, complete cds.
ACCESSION   U49845
VERSION     U49845.1  GI:1293613
KEYWORDS    .
SOURCE      Saccharomyces cerevisiae (baker's yeast)
  ORGANISM  Saccharomyces cerevisiae
            Eukaryota; Fungi; Dikarya; Ascomycota.
FEATURES             Location/Qualifiers
     source          1..5028
ORIGIN
        1 gatcctccat atacaacggt
//
"""


def test_locus_fields():
    h = decode_genbank(GB)
    assert h.format == "genbank"
    assert h.fields["locus"] == "SCU49845"
    assert h.fields["length"] == 5028
    assert h.fields["lengthUnits"] == "bp"
    assert h.fields["moleculeType"] == "DNA"
    assert h.fields["topology"] == "linear"


def test_definition_accession_version_organism():
    h = decode_genbank(GB)
    assert "TCP1-beta" in h.fields["definition"]
    assert h.fields["accession"] == "U49845"
    assert h.fields["version"].startswith("U49845.1")
    assert h.fields["organism"] == "Saccharomyces cerevisiae"


def test_protein_locus_has_no_molecule_type():
    p = b"LOCUS       PROT001  360 aa  linear   PRI 01-JAN-2020\nDEFINITION  a protein.\n//\n"
    h = decode_genbank(p)
    assert h.fields["lengthUnits"] == "aa"
    assert h.fields["length"] == 360
    assert h.fields["moleculeType"] is None
    assert h.fields["topology"] == "linear"


def test_rejects_non_genbank():
    assert decode_genbank(b"not a genbank file") is None
    assert decode_genbank(b"") is None


def test_dispatch():
    assert decode_bytes("x.gb", GB).format == "genbank"
    assert decode_bytes("x.gbk", GB).format == "genbank"


def test_json_safe():
    json.dumps(decode_genbank(GB).to_dict(), allow_nan=False)
