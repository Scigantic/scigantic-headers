"""ICS (Image Cytometry Standard) header decoding.

Golden: tests/fixtures/ics-huygens-benchmark.ics is the OME sample
ICS/jan/benchmark_v1_2018_x64y64z5c2s1t11_..._hrm.ics (downloads.openmicroscopy.org),
written by SVI Huygens through libics 1.6.3. Expected values are read from the
text itself (layout sizes 32 64 64 5 1 over order bits x y z p, real signed,
scale 3.6 3.6 10 micrometers); its 81,920-byte .ids = 64*64*5 float32 agrees.
"""

import json
import os

import pytest

from scigantic_headers import decode_bytes, has_decoder_for
import scigantic_headers.ics  # noqa: F401
from scigantic_headers.ics import decode_ics

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "ics-huygens-benchmark.ics")


def _ics(lines):
    return ("\t\n" + "\n".join("\t".join(l) for l in lines) + "\n").encode()


def test_extension_registered():
    assert has_decoder_for("a/b.ics") is True
    assert has_decoder_for("a/b.ids") is False   # the pixels file has no header


def test_huygens_sample():
    with open(FIX, "rb") as fh:
        d = decode_bytes("bench.ics", fh.read())
    assert d.format == "ics"
    assert d.summary == "ICS 1.0 image 64x64x5, float32, 3.6x3.6x10 um"
    f = d.fields
    assert f["icsVersion"] == "1.0"
    assert f["filename"].startswith("benchmark_v1_2018_x64y64z5c2s1t11")
    assert (f["sizeX"], f["sizeY"], f["sizeZ"], f["sizeC"], f["sizeT"]) == (64, 64, 5, None, None)
    assert f["axes"] == ["x", "y", "z", "p"] and f["sizes"] == [64, 64, 5, 1]
    assert f["dtype"] == "float32" and f["bits"] == 32 and f["significantBits"] == 32
    assert f["byteOrder"] == "little" and f["compression"] == "uncompressed"
    assert f["voxelSizeUm"] == [3.6, 3.6, 10.0]
    assert f["units"] == ["bits", "micrometers", "micrometers", "micrometers", ""]
    assert f["labels"] == ["intensity", "x", "y", "z", "p"]
    assert f["sensorType"] == "NipkowConfMicr"
    assert f["sensorParams"]["LambdaEx"] == "405.000000"
    assert f["numericalAperture"] == 0.963558
    assert f["historyCount"] == 4
    assert f["history"] == {"software": "SVI-Huygens", "libics_version": "1.6.3",
                            "version": "18.10.0p1", "date": "Fri 15 Mar 2019 15:13:13 UTC"}
    assert f["software"] == "SVI-Huygens"
    assert f["hasInlineData"] is False
    json.dumps(d.to_dict(), allow_nan=False)


def test_v2_inline_multichannel_timelapse():
    blob = _ics([
        ("ics_version", "2.0"),
        ("layout", "parameters", "6"),
        ("layout", "order", "bits", "x", "y", "z", "ch", "t"),
        ("layout", "sizes", "16", "512", "256", "30", "2", "11"),
        ("representation", "format", "integer"),
        ("representation", "sign", "unsigned"),
        ("representation", "byte_order", "1", "2"),
        ("representation", "compression", "gzip"),
        ("parameter", "scale", "1", "200", "200", "500", "1", "0.5"),
        ("parameter", "units", "bits", "nm", "nm", "nanometers", "", "s"),
        ("history", "software", "test"),
        ("history", "free text without a key"),
        ("end",),
    ]) + b"\x00\x01binary pixels follow"
    d = decode_ics(blob)
    assert d.summary == "ICS 2.0 image 512x256x30, 2 channels, 11 timepoints, uint16, 0.2x0.2x0.5 um"
    assert d.fields["sizeC"] == 2 and d.fields["sizeT"] == 11
    assert d.fields["voxelSizeUm"] == [0.2, 0.2, 0.5]
    assert d.fields["compression"] == "gzip"
    assert d.fields["hasInlineData"] is True
    assert d.fields["historyCount"] == 2
    assert d.fields["history"]["note2"] == "free text without a key"


def test_signed_integer_and_big_endian():
    d = decode_ics(_ics([("ics_version", "1.0"), ("layout", "order", "bits", "x", "y"),
                         ("layout", "sizes", "8", "4", "4"), ("representation", "format", "integer"),
                         ("representation", "sign", "signed"), ("representation", "byte_order", "2", "1")]))
    assert d.fields["dtype"] == "int8" and d.fields["byteOrder"] == "big"
    assert d.fields["voxelSizeUm"] == [None, None, None]   # no scale line: unknown


def test_history_is_bounded():
    lines = [("ics_version", "1.0"), ("layout", "order", "bits", "x"), ("layout", "sizes", "8", "4")]
    lines += [("history", "k%d" % i, "v") for i in range(100)]
    d = decode_ics(_ics(lines))
    assert d.fields["historyCount"] == 100 and len(d.fields["history"]) == 40


@pytest.mark.parametrize("blob", [b"", b"\x00" * 100, b"hello\tworld\n", b"\x89HDF\r\n\x1a\n" + b"\x00" * 50])
def test_rejects_non_ics(blob):
    assert decode_ics(blob) is None


def test_garbage_sizes_do_not_raise():
    d = decode_ics(_ics([("ics_version", "1.0"), ("layout", "order", "bits", "x", "y"),
                         ("layout", "sizes", "abc", "-5", "99999999999999"),
                         ("parameter", "scale", "nan", "inf", "x"), ("parameter", "units", "bits", "um", "um")]))
    assert d is not None and d.fields["sizeX"] is None and d.fields["voxelSizeUm"] == [None, None, None]
    json.dumps(d.to_dict(), allow_nan=False)
