"""read_session_optics: find a data file's session optics (STAR or .cs)."""

import os
import struct

from scigantic_headers import read_session_optics

STAR = ("data_optics\nloop_\n_rlnOpticsGroup #1\n_rlnImagePixelSize #2\n"
        "_rlnVoltage #3\n1  1.06  300.0\n")


def _cs_bytes(psize, kv, cs):
    descr = ("[('uid', '<u8'), ('blob/psize_A', '<f4'), "
             "('ctf/accel_kv', '<f4'), ('ctf/cs_mm', '<f4')]")
    body = "{'descr': %s, 'fortran_order': False, 'shape': (1,), }" % descr
    body = body + " " * ((64 - (10 + len(body) + 1) % 64) % 64) + "\n"
    return (b"\x93NUMPY\x01\x00" + struct.pack("<H", len(body)) + body.encode("latin1")
            + struct.pack("<Q", 1) + struct.pack("<fff", psize, kv, cs))


def _write(p, data, mode="w"):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, mode) as fh:
        fh.write(data)


def test_finds_star_next_to_movie(tmp_path):
    d = tmp_path / "Movies"
    _write(str(d / "frame.mrc"), "x")  # any file; we search by directory
    _write(str(d / "particles.star"), STAR)
    o = read_session_optics(str(d / "frame.mrc"))
    assert o["pixelSizeA"] == 1.06
    assert o["voltageKv"] == 300.0
    assert o["source"] == "relion-star:particles.star"


def test_finds_cs_next_to_movie(tmp_path):
    d = tmp_path / "extract"
    _write(str(d / "frame.mrc"), "x")
    _write(str(d / "particles.cs"), _cs_bytes(0.82, 200.0, 1.4), mode="wb")
    o = read_session_optics(str(d / "frame.mrc"))
    assert o["pixelSizeA"] == 0.82
    assert o["source"].startswith("cryosparc-cs:")


def test_searches_parent_directories(tmp_path):
    _write(str(tmp_path / "session" / "run.star"), STAR)
    deep = tmp_path / "session" / "Movies" / "gpu0"
    _write(str(deep / "frame.mrc"), "x")
    o = read_session_optics(str(deep / "frame.mrc"), depth=3)
    assert o["pixelSizeA"] == 1.06


def test_depth_limit_respected(tmp_path):
    _write(str(tmp_path / "run.star"), STAR)         # 3 levels up from the movie
    deep = tmp_path / "a" / "b" / "c"
    _write(str(deep / "frame.mrc"), "x")
    assert read_session_optics(str(deep / "frame.mrc"), depth=2) == {}


def test_none_found(tmp_path):
    _write(str(tmp_path / "frame.mrc"), "x")
    assert read_session_optics(str(tmp_path / "frame.mrc")) == {}


def test_accepts_a_directory(tmp_path):
    _write(str(tmp_path / "p.star"), STAR)
    assert read_session_optics(str(tmp_path))["pixelSizeA"] == 1.06
