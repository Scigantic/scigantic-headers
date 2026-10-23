"""KLB (Keller Lab Block) header decoding.

Golden: tests/fixtures/klb-cases.json holds the leading 1024 bytes of the OME
sample KLB/samples/img.klb (downloads.openmicroscopy.org). Its expected values
were read off the bytes with struct against the published 319-byte layout:
header version 2, 101 x 151 x 29 x 1 x 1, uint16, bzip2, pixel size 1.0 on
every axis, metadata "Testing metadata", blocks 96 x 96 x 8 x 1 x 1. No
second library was used (pyklb is not on PyPI for current Pythons).
"""

import json
import os
import struct

import pytest

from scigantic_headers import decode_bytes, has_decoder_for
import scigantic_headers.klb  # noqa: F401
from scigantic_headers.klb import KLB_HEADER_BYTES, decode_klb

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "klb-cases.json")


def make_klb(xyzct=(64, 64, 10, 2, 3), pixel=(0.4, 0.4, 2.0, 1.0, 1.0), version=2,
             dtype=1, compression=1, metadata=b"hello", block=(64, 64, 8, 1, 1), size=1024):
    buf = bytearray(size)
    struct.pack_into("<B5I5fBB256s5I", buf, 0, version, *xyzct, *pixel, dtype, compression,
                     metadata, *block)
    return bytes(buf)


def test_extension_registered():
    assert has_decoder_for("stack/img.klb") is True


def test_golden_fixture_matches():
    with open(FIX) as fh:
        cases = json.load(fh)["cases"]
    assert cases
    for case in cases:
        d = decode_bytes("fixture.klb", bytes.fromhex(case["hex"]))
        assert d is not None and d.to_dict() == case["expected"], case["name"]


def test_real_sample_values():
    with open(FIX) as fh:
        case = json.load(fh)["cases"][0]
    d = decode_klb(bytes.fromhex(case["hex"]))
    assert d.format == "klb"
    assert d.summary == "KLB volume 101x151x29, uint16, bzip2, 1x1x1 um"
    f = d.fields
    assert (f["sizeX"], f["sizeY"], f["sizeZ"], f["sizeC"], f["sizeT"]) == (101, 151, 29, 1, 1)
    assert f["headerVersion"] == 2 and f["dtype"] == "uint16" and f["compression"] == "bzip2"
    assert f["voxelSizeUm"] == [1.0, 1.0, 1.0]
    assert f["blockSize"] == [96, 96, 8, 1, 1]
    assert f["metadata"] == "Testing metadata"
    # 2 x 2 x 4 blocks, so a 16-entry uint64 offset table sits between the
    # header and the first block. In the full file a bzip2 stream ("BZh")
    # starts at byte 447, and 447 plus the table's last end offset is the
    # file size (382,141 bytes).
    assert f["dataOffset"] == KLB_HEADER_BYTES + 8 * 16 == 447


def test_data_offset_counts_partial_blocks():
    # 64x64x10x2x3 in 64x64x8x1x1 blocks: 1 * 1 * 2 * 2 * 3 = 12 blocks.
    assert decode_klb(make_klb()).fields["dataOffset"] == KLB_HEADER_BYTES + 8 * 12


def test_synthetic_multichannel_timelapse():
    d = decode_klb(make_klb())
    assert d.summary == "KLB volume 64x64x10, 2 channels, 3 timepoints, uint16, bzip2, 0.4x0.4x2 um"
    assert d.fields["voxelSizeUm"] == [0.4, 0.4, 2.0]
    assert d.fields["metadata"] == "hello"


@pytest.mark.parametrize("code,name", [(0, "uint8"), (3, "uint64"), (4, "int8"), (8, "float32"), (9, "float64")])
def test_dtype_table(code, name):
    assert decode_klb(make_klb(dtype=code)).fields["dtype"] == name


@pytest.mark.parametrize("kw", [
    {"version": 3}, {"version": 0}, {"dtype": 10}, {"compression": 3},
    {"xyzct": (0, 64, 10, 1, 1)}, {"xyzct": (1 << 31 | 1, 64, 10, 1, 1)},
])
def test_rejections(kw):
    assert decode_klb(make_klb(**kw)) is None


def test_too_short_and_garbage():
    assert decode_klb(b"\x02" + b"\x00" * 100) is None
    assert decode_klb(b"") is None


def test_nonfinite_pixel_size_is_dropped():
    d = decode_klb(make_klb(pixel=(float("nan"), float("inf"), 0.0, 1.0, 1.0)))
    assert d.fields["voxelSizeUm"] == [None, None, None]
    assert not d.summary.endswith(" um")
    json.dumps(d.to_dict(), allow_nan=False)
