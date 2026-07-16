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
