"""PDB structure-header decoding."""

import json

from scigantic_headers import decode_pdb
from scigantic_headers.decoders import decode_bytes

# Build the HEADER line at exact PDB columns: classification 11-50, deposition
# date 51-59, PDB id 63-66.
_HEADER = ("HEADER" + " " * 4 + "HYDROLASE".ljust(40) + "01-JAN-98" + " " * 3 + "1ABC")
PDB = (
    _HEADER + "\n"
    "TITLE     CRYSTAL STRUCTURE OF SOMETHING INTERESTING\n"
    "EXPDTA    X-RAY DIFFRACTION\n"
    "REMARK   2 RESOLUTION.    2.00 ANGSTROMS.\n"
    "ATOM      1  N   MET A   1      11.104  13.207  10.567  1.00  0.00           N\n"
).encode("ascii")


def test_core_fields():
    h = decode_pdb(PDB)
    assert h.format == "pdb"
    assert h.fields["pdbId"] == "1ABC"
    assert h.fields["classification"] == "HYDROLASE"
    assert h.fields["depositionDate"] == "01-JAN-98"
    assert h.fields["experimentMethod"] == "X-RAY DIFFRACTION"
    assert h.fields["resolutionAngstrom"] == 2.00
    assert "CRYSTAL STRUCTURE" in h.fields["title"]


def test_nmr_has_no_resolution():
    p = (
        "HEADER" + " " * 4 + "DE NOVO PROTEIN".ljust(40) + "01-JAN-20" + " " * 3 + "2XYZ" + "\n"
        "EXPDTA    SOLUTION NMR\n"
    ).encode("ascii")
    h = decode_pdb(p)
    assert h.fields["resolutionAngstrom"] is None
    assert h.fields["experimentMethod"] == "SOLUTION NMR"
    assert h.fields["pdbId"] == "2XYZ"


def test_rejects_non_pdb():
    assert decode_pdb(b"not a pdb file at all") is None
    assert decode_pdb(b"ATOM      1  N   MET A   1  ") is None  # ATOM is not a leading record


def test_dispatch():
    assert decode_bytes("x.pdb", PDB).format == "pdb"
    assert decode_bytes("x.ent", PDB).format == "pdb"


def test_json_safe():
    json.dumps(decode_pdb(PDB).to_dict(), allow_nan=False)
