"""CryoSPARC .cs — structured-NPY header decode + first-record optics extraction.

Hand-built .npy bytes (verified live against numpy structured arrays), so the
suite stays dependency-free.
"""

import struct

import pytest

from scigantic_headers import (
    decode_bytes,
    decode_cryosparc_header,
    has_decoder_for,
    parse_cryosparc_optics,
)


def make_cs(descr_literal: str, n: int, record: bytes = b"") -> bytes:
    body = "{'descr': %s, 'fortran_order': False, 'shape': (%d,), }" % (descr_literal, n)
    preamble = 10 + len(body) + 1
    body = body + " " * ((64 - preamble % 64) % 64) + "\n"
    header = b"\x93NUMPY\x01\x00" + struct.pack("<H", len(body)) + body.encode("latin1")
    return header + record


# A CryoSPARC-like schema: uid, a fixed-width path string, a 2-elem subarray,
# then the optics float32 fields. Offsets: uid@0(8), path@8(66), shape@74(8),
# psize@82(4), kv@86(4), cs@90(4).
DESCR = ("[('uid', '<u8'), ('blob/path', 'S66'), ('blob/shape', '<u4', (2,)), "
         "('blob/psize_A', '<f4'), ('ctf/accel_kv', '<f4'), ('ctf/cs_mm', '<f4')]")


def _record(psize, kv, cs):
    return (struct.pack("<Q", 12345)
            + b"J1/movie_0001.mrc".ljust(66, b"\x00")
            + struct.pack("<2I", 4096, 4096)
            + struct.pack("<fff", psize, kv, cs))


def test_cs_extension_registered():
    assert has_decoder_for("run/particles.cs") is True


def test_cs_header_decode_reports_schema():
    d = decode_cryosparc_header(make_cs(DESCR, 50000))
    assert d.format == "cryosparc"
    assert d.fields["numRecords"] == 50000
    assert d.fields["numFields"] == 6
    assert "blob/psize_A" in d.fields["fields"]
    assert d.summary == "CryoSPARC dataset, 50000 records, 6 fields"


def test_cs_optics_extracted_from_first_record():
    data = make_cs(DESCR, 3, _record(1.06, 300.0, 2.7))
    optics = parse_cryosparc_optics(data)
    assert optics == {"pixelSizeA": 1.06, "voltageKv": 300.0, "sphericalAberrationMm": 2.7}


def test_cs_optics_handles_offset_after_fixed_string_and_subarray():
    # Different values to be sure the offsets (not defaults) are read.
    data = make_cs(DESCR, 1, _record(0.82, 200.0, 1.4))
    assert parse_cryosparc_optics(data) == {
        "pixelSizeA": 0.82, "voltageKv": 200.0, "sphericalAberrationMm": 1.4}


def test_cs_object_dtype_bails():
    # An object field ('|O') means .npy stores a pickle, not a flat buffer — the
    # optics reader must decline rather than read garbage. The header decode
    # still lists the schema (it never touches the data).
    descr = "[('uid', '<u8'), ('blob/path', '|O'), ('blob/psize_A', '<f4')]"
    data = make_cs(descr, 10)
    assert parse_cryosparc_optics(data) == {}
    assert decode_cryosparc_header(data).fields["numFields"] == 3  # schema still read


def test_cs_field_past_read_is_skipped():
    # Header present but the record bytes weren't read -> no fields extracted,
    # no error.
    assert parse_cryosparc_optics(make_cs(DESCR, 5)) == {}


def test_cs_non_structured_is_not_cryosparc():
    plain = make_cs("'<f8'", 100)  # a plain (non-structured) .npy renamed .cs
    assert decode_cryosparc_header(plain) is None
    assert parse_cryosparc_optics(plain) == {}


def test_cs_dispatch_via_decode_bytes():
    d = decode_bytes("x.cs", make_cs(DESCR, 7))
    assert d is not None and d.format == "cryosparc"
