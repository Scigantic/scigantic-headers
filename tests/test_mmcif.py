"""mmCIF / CIF structure-header decoding."""

import json

from scigantic_headers import decode_mmcif
from scigantic_headers.decoders import decode_bytes

CIF = b"""data_1ABC
#
_entry.id   1ABC
#
_struct.title   'Crystal structure of a hydrolase'
#
_struct_keywords.pdbx_keywords   HYDROLASE
#
_exptl.method   'X-RAY DIFFRACTION'
#
_reflns.d_resolution_high   1.80
#
_pdbx_database_status.recvd_initial_deposition_date   1998-01-01
#
loop_
_atom_site.group_PDB
_atom_site.id
ATOM 1
"""


def test_core_fields():
    h = decode_mmcif(CIF)
    assert h.format == "mmcif"
    assert h.fields["entryId"] == "1ABC"
    assert h.fields["title"] == "Crystal structure of a hydrolase"
    assert h.fields["classification"] == "HYDROLASE"
    assert h.fields["experimentMethod"] == "X-RAY DIFFRACTION"
    assert h.fields["resolutionAngstrom"] == 1.80
    assert h.fields["depositionDate"] == "1998-01-01"


def test_cryoem_resolution_source():
    c = b"data_XYZ\n_exptl.method 'ELECTRON MICROSCOPY'\n_em_3d_reconstruction.resolution 3.2\n"
    h = decode_mmcif(c)
    assert h.fields["experimentMethod"] == "ELECTRON MICROSCOPY"
    assert h.fields["resolutionAngstrom"] == 3.2


def test_semicolon_delimited_title_block():
    c = b"data_ABC\n_struct.title\n;\nA very long title spanning\na text block\n;\n_entry.id ABC\n"
    h = decode_mmcif(c)
    assert "long title" in h.fields["title"]


def test_block_code_when_no_entry_id():
    c = b"data_WXYZ\n_exptl.method 'SOLUTION NMR'\n"
    h = decode_mmcif(c)
    assert h.fields["entryId"] == "WXYZ"


def test_rejects_non_cif():
    assert decode_mmcif(b"not a cif file") is None
    assert decode_mmcif(b"") is None


def test_dispatch():
    assert decode_bytes("x.cif", CIF).format == "mmcif"
    assert decode_bytes("x.mmcif", CIF).format == "mmcif"


def test_json_safe():
    json.dumps(decode_mmcif(CIF).to_dict(), allow_nan=False)
