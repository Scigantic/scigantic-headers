import struct

import pytest


def make_mrc(nx, ny, nz, mode, mx=None, cella_x=0.0, nsymbt=0, stamp=b"MAP ", size=1024):
    """Synthetic MRC2014 header. Byte layout matches the decoder and the
    fixtures in backend/tests/headerDecoders.unit.test.ts, so the Python and
    TypeScript twins are held to the same expected values."""
    buf = bytearray(size)
    struct.pack_into("<4i", buf, 0, nx, ny, nz, mode)
    struct.pack_into("<i", buf, 28, mx if mx is not None else nx)
    struct.pack_into("<f", buf, 40, cella_x)
    struct.pack_into("<i", buf, 92, nsymbt)
    buf[208:212] = stamp
    return bytes(buf)


@pytest.fixture
def mrc():
    return make_mrc


# ── Sparse fixtures for seeking decoders ─────────────────────────────────────
#
# A seeking decoder reads a few KB out of a file that may be gigabytes. To hold
# a real file to a golden fixture, the test records exactly the byte ranges a
# decode touched and stores those (hex) with the file size; the fixture then
# replays as a BytesSource. tests/fixtures/*.sparse.json are such recordings.

import json
import os

from scigantic_headers.decoders import BytesSource, Source

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
SAMPLES_DIR = os.environ.get("SCIGANTIC_HEADERS_SAMPLES", "")


class RecordingSource(Source):
    """Wraps a Source, remembering every range actually read."""

    def __init__(self, inner: Source):
        super().__init__(1 << 40)
        self.inner = inner
        self.size = inner.size
        self.ranges = {}

    def _fetch(self, offset, n):
        data = self.inner.read_at(offset, n)
        if data:
            self.ranges[offset] = data
        return data

    def dump(self, path: str, note: str = ""):
        merged = _merge_ranges(self.ranges)
        with open(path, "w") as fh:
            json.dump({"note": note, "size": self.size,
                       "ranges": {str(o): b.hex() for o, b in merged.items()}}, fh, indent=1)


def _merge_ranges(ranges):
    out = {}
    cur_start, cur = None, b""
    for start in sorted(ranges):
        chunk = ranges[start]
        if cur_start is not None and start <= cur_start + len(cur):
            overlap = cur_start + len(cur) - start
            cur += chunk[overlap:]
        else:
            if cur_start is not None:
                out[cur_start] = cur
            cur_start, cur = start, chunk
    if cur_start is not None:
        out[cur_start] = cur
    return out


def load_sparse(name: str) -> BytesSource:
    """A BytesSource replaying tests/fixtures/<name>.sparse.json."""
    with open(os.path.join(FIXTURES_DIR, name + ".sparse.json")) as fh:
        rec = json.load(fh)
    return BytesSource(size=rec["size"],
                       ranges={int(o): bytes.fromhex(h) for o, h in rec["ranges"].items()})


def sample_path(name: str):
    """Path of a real sample file under $SCIGANTIC_HEADERS_SAMPLES, or None (the
    sparse fixture is the committed evidence; the sample is for re-recording)."""
    p = os.path.join(SAMPLES_DIR, name) if SAMPLES_DIR else ""
    return p if p and os.path.exists(p) else None
