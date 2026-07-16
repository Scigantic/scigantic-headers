"""I/O + bounded-parallel batch tests. No network — local files only."""

import os

import pytest

from scigantic_headers import (
    decode_file,
    decode_paths,
    iter_decodable_files,
    read_leading_bytes,
)
from scigantic_headers.decoders import HEADER_BYTES
from scigantic_headers.sources import _batch


def _write(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(data)


def test_read_is_bounded_regardless_of_file_size(tmp_path, mrc):
    # A large file — the read must return exactly HEADER_BYTES, not the file.
    big = tmp_path / "big.mrc"
    _write(str(big), mrc(4096, 4096, 16, 2, mx=4096, cella_x=4341.76) + b"\0" * (2 * 1024 * 1024))
    data = read_leading_bytes(str(big))
    assert len(data) == HEADER_BYTES


def test_decode_file(tmp_path, mrc):
    p = tmp_path / "frame.mrc"
    _write(str(p), mrc(2048, 2048, 8, 2, mx=2048, cella_x=2170.88))
    d = decode_file(str(p))
    assert d.fields["nz"] == 8
    assert d.summary == "MRC stack 2048x2048x8, float32, 1.06 A/px"


def test_decode_file_skips_unregistered(tmp_path):
    p = tmp_path / "notes.txt"
    _write(str(p), b"hello")
    assert decode_file(str(p)) is None


def test_decode_file_missing_returns_none():
    assert decode_file("/no/such/file.mrc") is None


def test_iter_decodable_files_filters(tmp_path, mrc):
    _write(str(tmp_path / "s1" / "a.mrc"), mrc(64, 64, 1, 2))
    _write(str(tmp_path / "s1" / "notes.txt"), b"x")
    _write(str(tmp_path / "s2" / "b.mrcs"), mrc(64, 64, 2, 2))
    found = sorted(os.path.basename(p) for p in iter_decodable_files(str(tmp_path)))
    assert found == ["a.mrc", "b.mrcs"]  # notes.txt excluded


def test_decode_paths_parallel_matches_serial(tmp_path, mrc):
    paths = []
    for i in range(12):
        p = tmp_path / f"f{i}.mrc"
        _write(str(p), mrc(128 + i, 128, 1 + (i % 3), 2, mx=128 + i, cella_x=(128 + i)))
        paths.append(str(p))
    # add a couple of undecodable files that must be dropped
    _write(str(tmp_path / "skip.txt"), b"x")

    serial = decode_paths(paths, workers=1)
    parallel = decode_paths(paths, workers=8)

    assert len(serial) == 12
    assert set(serial) == set(parallel)
    # identical results regardless of concurrency
    for p in paths:
        assert serial[p].to_dict() == parallel[p].to_dict()


def test_decode_paths_empty_and_all_skipped(tmp_path):
    assert decode_paths([]) == {}
    _write(str(tmp_path / "a.txt"), b"x")
    assert decode_paths([str(tmp_path / "a.txt")]) == {}


def test_batch_prefilters_before_pool(tmp_path, mrc):
    # Only registered extensions reach the pool; the .txt never causes a read.
    _write(str(tmp_path / "a.mrc"), mrc(64, 64, 1, 2, mx=64, cella_x=64))
    pairs = _batch([str(tmp_path / "a.mrc"), str(tmp_path / "b.txt")], decode_file, workers=4)
    assert [os.path.basename(k) for k, _ in pairs] == ["a.mrc"]


# ── gzip transparency (.nii.gz, .mrc.gz) ────────────────────────────────────

import gzip

from scigantic_headers import is_decodable
from scigantic_headers.sources import _gunzip_leading, _inner_key


@pytest.mark.parametrize("key,expected", [
    ("s/brain.nii.gz", True), ("s/frame.mrc.gz", True),
    ("s/notes.txt.gz", False), ("s/archive.gz", False), ("s/x.nii", True),
])
def test_is_decodable_sees_through_gz(key, expected):
    assert is_decodable(key) is expected


def test_inner_key_strips_one_gz():
    assert _inner_key("a/b.nii.gz") == "a/b.nii"
    assert _inner_key("a/B.NII.GZ") == "a/B.NII"
    assert _inner_key("a/b.nii") == "a/b.nii"


def test_decode_file_gzipped_mrc(tmp_path, mrc):
    p = tmp_path / "frame.mrc.gz"
    with gzip.open(str(p), "wb") as fh:
        fh.write(mrc(1024, 1024, 4, 2, mx=1024, cella_x=1085.44))  # 1.06 A/px
    d = decode_file(str(p))
    assert d is not None
    assert d.fields["nz"] == 4
    assert d.fields["pixelSizeA"] == 1.06


def test_decode_file_gzipped_undecodable_is_none(tmp_path):
    p = tmp_path / "notes.txt.gz"
    with gzip.open(str(p), "wb") as fh:
        fh.write(b"just text")
    assert decode_file(str(p)) is None


def test_iter_includes_gz(tmp_path, mrc):
    with gzip.open(str(tmp_path / "a.mrc.gz"), "wb") as fh:
        fh.write(mrc(64, 64, 1, 2))
    _write(str(tmp_path / "b.txt"), b"x")
    found = sorted(os.path.basename(p) for p in iter_decodable_files(str(tmp_path)))
    assert found == ["a.mrc.gz"]


def test_gunzip_leading_tolerates_truncation(mrc):
    full = gzip.compress(mrc(4096, 4096, 16, 2) + b"\x00" * 200000)
    truncated = full[: len(full) // 2]     # cut the stream mid-way
    head = _gunzip_leading(truncated, 1024)  # must not raise
    assert len(head) >= 212                  # enough to decode the MRC header
    assert head[208:212] == b"MAP "
