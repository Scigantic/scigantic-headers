"""Parquet footer decode. Uses a real pyarrow-written fixture
(tests/fixtures/sample.parquet) so the suite needs no pyarrow at test time."""

import json
import os

import pytest

from scigantic_headers import decode_file, decode_parquet, has_decoder_for

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "sample.parquet")

EXPECTED_COLUMNS = [
    {"name": "name", "type": "BYTE_ARRAY"},
    {"name": "n64", "type": "INT64"},
    {"name": "n32", "type": "INT32"},
    {"name": "f64", "type": "DOUBLE"},
    {"name": "f32", "type": "FLOAT"},
    {"name": "flag", "type": "BOOLEAN"},
]


def test_parquet_extension_registered():
    assert has_decoder_for("data/table.parquet") is True


def test_decodes_real_fixture_schema_and_rows():
    hdr = decode_file(FIXTURE)
    assert hdr is not None
    assert hdr.format == "parquet"
    assert hdr.fields["numRows"] == 1000
    assert hdr.fields["numColumns"] == 6
    assert hdr.fields["columns"] == EXPECTED_COLUMNS
    assert hdr.summary == "Parquet, 1000 rows, 6 columns"


def test_result_is_json_safe():
    hdr = decode_file(FIXTURE)
    json.loads(json.dumps(hdr.to_dict(), allow_nan=False))


def test_decode_parquet_takes_trailing_bytes():
    # decode_parquet works on the file's trailing bytes directly.
    with open(FIXTURE, "rb") as fh:
        tail = fh.read()  # whole small file; the footer is at its end
    hdr = decode_parquet(tail)
    assert hdr.fields["numRows"] == 1000


@pytest.mark.parametrize("buf", [
    b"",
    b"not a parquet file",
    b"\x00" * 100 + b"PAR1",           # trailing magic, nonsense footer length
    b"PAR1" + b"\x00" * 50 + b"PAR1",  # both magics, no valid footer
])
def test_garbage_returns_none_not_raise(buf):
    assert decode_parquet(buf) is None


def test_footer_larger_than_tail_returns_none():
    # A valid trailing magic but a footer length pointing past the bytes we have.
    import struct
    buf = b"\x01\x02\x03\x04" + struct.pack("<I", 10_000) + b"PAR1"
    assert decode_parquet(buf) is None
