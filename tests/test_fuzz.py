"""Robustness fuzzing. A decoder is handed the leading bytes of arbitrary files,
so for ANY input it must either return None or a valid, JSON-safe DecodedHeader —
never raise, never emit NaN/inf. Seeded random for reproducibility.

Two input populations:
  * pure random bytes of every length (exercises the reject paths)
  * magic-seeded bytes (correct magic + random remainder), which get PAST the
    cheap magic check and into the field parsing where NaN/inf could leak.
"""

import json
import random
import struct

import pytest

from scigantic_headers import DecodedHeader, decode_bytes
from scigantic_headers.decoders import (
    decode_mrc_header,
    decode_nifti_header,
    decode_npy_header,
)
from scigantic_headers.decoders import _DECODERS_BY_EXT, _DECODERS_BY_NAME

SEED = 1729
ITERATIONS = 4000
# Every dispatch path: extension keys (registered + unregistered) and the
# file-name keys (Illumina RunInfo.xml etc.), so filename-dispatched decoders
# are held to the same never-raise / JSON-safe contract.
KEYS = (
    ["file.%s" % e for e in sorted(_DECODERS_BY_EXT)]
    + ["file.txt", "file.parquet", "file.unknown", "file"]
    + sorted(_DECODERS_BY_NAME)
)


def _assert_valid_or_none(result):
    if result is None:
        return
    assert isinstance(result, DecodedHeader)
    assert isinstance(result.format, str) and result.format
    assert isinstance(result.summary, str)
    assert isinstance(result.fields, dict)
    # allow_nan=False raises if any NaN/inf slipped through — the record must be
    # standards-valid JSON.
    text = json.dumps(result.to_dict(), allow_nan=False)
    json.loads(text)  # round-trips


def _seed_magic(rng, kind):
    n = rng.randint(0, 1200)
    buf = bytearray(rng.getrandbits(8) for _ in range(n))
    if kind == "mrc" and len(buf) >= 212:
        buf[208:212] = b"MAP "
    elif kind == "npy" and len(buf) >= 12:
        buf[0:6] = b"\x93NUMPY"
        buf[6] = rng.choice([1, 2])
        buf[7] = 0
    elif kind == "nifti" and len(buf) >= 348:
        struct.pack_into("<i", buf, 0, 348)   # sizeof_hdr passes
        buf[344:348] = b"n+1\x00"             # magic passes -> reaches pixdim
    return bytes(buf)


def test_decode_bytes_total_over_random():
    rng = random.Random(SEED)
    for _ in range(ITERATIONS):
        key = rng.choice(KEYS)
        buf = bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 2048)))
        try:
            _assert_valid_or_none(decode_bytes(key, buf))
        except Exception as e:  # noqa: BLE001 — the whole point is nothing escapes
            pytest.fail(f"decode_bytes raised on key={key!r} len={len(buf)}: {e!r}")


@pytest.mark.parametrize("kind,fn", [
    ("mrc", lambda b: decode_mrc_header(b)),
    ("mrc", lambda b: decode_mrc_header(b, strict=False)),
    ("npy", decode_npy_header),
    ("nifti", decode_nifti_header),
])
def test_decoders_total_on_magic_seeded(kind, fn):
    # Buffers that pass the magic check and drive the field-parsing paths where
    # NaN/inf could leak (esp. NIfTI pixdim, MRC CELLA).
    rng = random.Random(SEED + hash(kind) % 1000)
    for _ in range(ITERATIONS):
        buf = _seed_magic(rng, kind)
        try:
            _assert_valid_or_none(fn(buf))
        except Exception as e:  # noqa: BLE001
            pytest.fail(f"{kind} decoder raised on magic-seeded len={len(buf)}: {e!r}")


def test_nifti_nan_pixdim_is_sanitized():
    # A header that validates but carries NaN voxel spacing must not emit NaN.
    buf = bytearray(352)
    struct.pack_into("<i", buf, 0, 348)
    struct.pack_into("<8h", buf, 40, 3, 64, 64, 40, 1, 1, 1, 1)  # ndim=3
    struct.pack_into("<h", buf, 70, 16)  # float32
    struct.pack_into("<8f", buf, 76, 0.0, float("nan"), float("inf"), 2.0, 0, 0, 0, 0)
    buf[344:348] = b"n+1\x00"
    d = decode_nifti_header(bytes(buf))
    assert d is not None
    assert d.fields["voxelSizesMm"] == [None, None, 2.0]
    assert "mm" not in d.summary  # incomplete spacing -> omitted, not "NaN mm"
    json.dumps(d.to_dict(), allow_nan=False)  # standards-valid
