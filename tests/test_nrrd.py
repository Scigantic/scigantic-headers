"""NRRD header decoding.

Golden: tests/fixtures/nrrd-cases.json. The first case is the real header of the
Allen Institute CCFv3 25 um annotation volume (download.alleninstitute.org,
mouse_ccf/annotation/ccf_2017/annotation_25.nrrd, fetched by Range request; the
file is 528 x 320 x 456 uint32, gzip, 25 um isotropic in LPS space). The others
were written by pynrrd 1.1.3: a gzip volume with space directions, units, labels
and a key/value pair; a detached `.nhdr` with `spacings`; an ascii RGB volume
whose colour axis has no space direction. Every expected value was checked
against pynrrd's read_header on the same bytes when the fixture was made.
"""

import json
import os

from scigantic_headers import decode_bytes, decode_nrrd, has_decoder_for
from scigantic_headers.sources import is_decodable
from scigantic_headers.decoders import read_for
from scigantic_headers.nrrd import NRRD_READ_BYTES

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "nrrd-cases.json")

with open(FIX) as _fh:
    CASES = json.load(_fh)["cases"]


def test_extensions_registered():
    assert has_decoder_for("atlas/annotation_25.nrrd") and has_decoder_for("vol.nhdr")
    assert is_decodable("vol.nrrd.gz")
    assert read_for("vol.nrrd").leading == NRRD_READ_BYTES


def test_golden_fixture_matches():
    assert len(CASES) >= 4
    for case in CASES:
        d = decode_bytes("fixture.nrrd", bytes.fromhex(case["hex"]))
        assert d is not None and d.to_dict() == case["expected"], case["name"]


def test_allen_ccf_annotation_header():
    d = decode_nrrd(bytes.fromhex(CASES[0]["hex"]))
    assert d.format == "nrrd"
    assert d.summary == "NRRD 3-D, 528x320x456, uint32, gzip, spacing 25x25x25"
    f = d.fields
    assert f["nrrdVersion"] == 4 and f["dimension"] == 3
    assert f["shape"] == [528, 320, 456] and f["numElements"] == 528 * 320 * 456
    assert f["dtype"] == "uint32" and f["encoding"] == "gzip" and f["endian"] == "little"
    assert f["space"] == "left-posterior-superior"
    assert f["spaceDirections"] == [[25.0, 0.0, 0.0], [0.0, 25.0, 0.0], [0.0, 0.0, 25.0]]
    assert f["spacing"] == [25.0, 25.0, 25.0] and f["spaceOrigin"] == [0.0, 0.0, 0.0]
    assert f["kinds"] == ["domain", "domain", "domain"]
    assert f["detached"] is False and f["headerComplete"] is True and f["headerBytes"] == 310


def test_units_labels_and_key_values():
    f = decode_nrrd(bytes.fromhex(CASES[1]["hex"])).fields
    assert f["dtype"] == "int16" and f["spacing"] == [0.7, 0.7, 1.2]
    assert f["spaceUnits"] == ["mm", "mm", "mm"] and f["labels"] == ["x", "y", "z"]
    assert f["spaceOrigin"] == [1.0, 2.0, 3.0]
    assert f["keyValues"] == {"custom key": "custom value"} and f["keyValueCount"] == 1
    assert decode_nrrd(bytes.fromhex(CASES[1]["hex"])).summary.endswith("spacing 0.7x0.7x1.2 mm")


def test_detached_header_names_its_data_file():
    d = decode_nrrd(bytes.fromhex(CASES[2]["hex"]))
    f = d.fields
    assert f["detached"] is True and f["dataFile"] == "t2.raw"
    assert f["dtype"] == "float64" and f["encoding"] == "raw" and f["spacing"] == [2.0, 2.0, 2.0]
    assert f["spaceDirections"] is None
    assert d.summary.endswith("detached header")


def test_non_spatial_axis_has_no_spacing():
    f = decode_nrrd(bytes.fromhex(CASES[3]["hex"])).fields
    assert f["kinds"] == ["RGB-color", "domain", "domain"]
    assert f["spaceDirections"] == [None, [1.0, 0.0], [0.0, 1.0]]
    assert f["spacing"] == [None, 1.0, 1.0] and f["spaceDimension"] == 2
    assert f["encoding"] == "ascii"


def test_type_aliases_and_encodings():
    def hdr(t, enc="raw"):
        return ("NRRD0001\ntype: %s\ndimension: 1\nsizes: 4\nencoding: %s\n\n" % (t, enc)).encode()
    for t, want in [("unsigned char", "uint8"), ("uchar", "uint8"), ("short int", "int16"),
                    ("unsigned  short", "uint16"), ("int", "int32"), ("uint64_t", "uint64"),
                    ("float", "float32"), ("double", "float64"), ("block", "block"), ("weird", "weird")]:
        assert decode_nrrd(hdr(t)).fields["dtype"] == want, t
    for e, want in [("gz", "gzip"), ("bz2", "bzip2"), ("txt", "ascii"), ("hex", "hex"), ("raw", "raw")]:
        assert decode_nrrd(hdr("int", e)).fields["encoding"] == want


def test_header_without_blank_line_is_marked_truncated():
    raw = b"NRRD0005\ntype: float\ndimension: 2\nsizes: 10 20\nencoding: raw\n"
    d = decode_nrrd(raw)
    assert d.fields["headerComplete"] is False and d.fields["headerBytes"] is None
    assert d.fields["shape"] == [10, 20] and "header truncated" in d.summary


def test_crlf_and_comments():
    raw = b"NRRD0003\r\n# a comment\r\ntype: short\r\ndimension: 2\r\nsizes: 3 4\r\nencoding: raw\r\nspacings: 1.5 nan\r\n\r\nDATA"
    f = decode_nrrd(raw).fields
    assert f["shape"] == [3, 4] and f["dtype"] == "int16"
    assert f["spacing"] == [1.5, None]
    assert f["headerBytes"] == len(raw) - 4


def test_rejects_non_nrrd_and_bad_dimension():
    assert decode_nrrd(b"NRRD0009\n") is None
    assert decode_nrrd(b"NRRD000") is None
    assert decode_nrrd(b"not a header\n\n") is None
    assert decode_nrrd(b"NRRD0004\ndimension: 99\n\n") is None
    assert decode_nrrd(b"") is None
