"""Pure-core decoder tests. Held to the same fixtures as the TS twin."""

import pytest

from scigantic_headers import (
    DecodedHeader,
    decode_bytes,
    decode_mrc_header,
    extension_of,
    has_decoder_for,
    register_decoder,
)


@pytest.mark.parametrize("key,ext", [
    ("a/b/frame.MRC", "mrc"),
    ("EMPIAR-10028/data/stack.mrcs", "mrcs"),
    ("tomo.REC", "rec"),
    ("README", ""),
    (".hidden", ""),
    ("no_ext_dir/", ""),
])
def test_extension_of(key, ext):
    assert extension_of(key) == ext


@pytest.mark.parametrize("key,expected", [
    ("x/foo.mrc", True), ("x/foo.rec", True), ("x/foo.map", True),
    ("x/foo.parquet", False), ("x/foo.csv", False), ("x/plain", False),
])
def test_has_decoder_for(key, expected):
    assert has_decoder_for(key) is expected


def test_single_image_with_pixel_size(mrc):
    d = decode_mrc_header(mrc(4096, 4096, 1, 2, mx=4096, cella_x=4341.76))
    assert isinstance(d, DecodedHeader)
    assert d.format == "mrc"
    assert d.fields["dtype"] == "float32"
    assert d.fields["pixelSizeA"] == 1.06
    assert d.fields["isStack"] is False
    assert d.fields["frameBytes"] == 4096 * 4096 * 4
    assert d.summary == "MRC image 4096x4096, float32, 1.06 A/px"


def test_stack_reports_frame_count(mrc):
    d = decode_mrc_header(mrc(5760, 4092, 40, 1, mx=5760, cella_x=6048))
    assert d.fields["isStack"] is True
    assert d.fields["nz"] == 40
    assert d.fields["dtype"] == "int16"
    assert d.fields["frameBytes"] == 5760 * 4092 * 2
    assert d.summary == "MRC stack 5760x4092x40, int16, 1.05 A/px"


def test_raw_stack_without_pixel_size(mrc):
    # EMPIAR-10002 case: dims/frames/dtype present, CELLA = 0 -> unknown apix
    d = decode_mrc_header(mrc(4096, 4096, 16, 2, mx=4096, cella_x=0))
    assert d.fields["pixelSizeA"] is None
    assert d.summary == "MRC stack 4096x4096x16, float32"


def test_pixel_size_null_when_grid_zero(mrc):
    d = decode_mrc_header(mrc(512, 512, 1, 2, mx=0, cella_x=500))
    assert d.fields["pixelSizeA"] is None


def test_extended_header_offset(mrc):
    d = decode_mrc_header(mrc(512, 512, 1, 2, nsymbt=640))
    assert d.fields["dataOffset"] == 1024 + 640


@pytest.mark.parametrize("mode,dtype", [(0, "int8"), (1, "int16"), (2, "float32"), (6, "uint16"), (12, "float16")])
def test_mode_dtype_table(mrc, mode, dtype):
    d = decode_mrc_header(mrc(64, 64, 1, mode, mx=64, cella_x=64))
    assert d.fields["dtype"] == dtype


def test_rejects_bad_stamp(mrc):
    assert decode_mrc_header(mrc(512, 512, 1, 2, stamp=b"JUNK")) is None


def test_strict_false_skips_stamp_check(mrc):
    # A caller that vouches for the format (e.g. scigantic_empiar) decodes even
    # a stamp-less pre-2014 MRC. Strict dispatch still rejects it.
    stampless = mrc(4096, 4096, 16, 2, mx=4096, cella_x=4341.76, stamp=b"\x00\x00\x00\x00")
    assert decode_mrc_header(stampless, strict=True) is None
    d = decode_mrc_header(stampless, strict=False)
    assert d is not None
    assert d.fields["pixelSizeA"] == 1.06


def test_strict_false_accepts_short_header(mrc):
    # Permissive callers only need the numeric fields (through byte 96), not the
    # 212 bytes the stamp requires.
    short = mrc(256, 256, 1, 2, mx=256, cella_x=256, size=128)
    assert decode_mrc_header(short, strict=True) is None
    assert decode_mrc_header(short, strict=False).fields["nx"] == 256


@pytest.mark.parametrize("dims", [(-1, 512, 1), (512, 512, 99_999_999), (0, 512, 1)])
def test_rejects_insane_dims(mrc, dims):
    assert decode_mrc_header(mrc(dims[0], dims[1], dims[2], 2)) is None


def test_rejects_truncated():
    assert decode_mrc_header(b"\x00" * 100) is None


def test_dispatch_by_extension(mrc):
    d = decode_bytes("EMPIAR-10028/data/img.mrc", mrc(1024, 1024, 1, 2, mx=1024, cella_x=1024))
    assert d.format == "mrc"
    assert d.fields["pixelSizeA"] == 1.0


def test_dispatch_unregistered_returns_none(mrc):
    assert decode_bytes("data/table.parquet", mrc(1024, 1024, 1, 2)) is None


def test_dispatch_registered_but_invalid_returns_none():
    assert decode_bytes("data/notreally.mrc", b"\x00" * 1024) is None


def test_to_dict_roundtrips(mrc):
    d = decode_bytes("x.mrc", mrc(256, 256, 1, 2, mx=256, cella_x=256))
    dd = d.to_dict()
    assert dd["format"] == "mrc"
    assert dd["fields"]["nx"] == 256
    assert set(dd) == {"format", "summary", "fields"}


def test_register_new_format_is_isolated():
    # Extending the registry is a pure function + one call; nothing else changes.
    sentinel = DecodedHeader(format="fake", summary="fake header", fields={"ok": True})
    register_decoder("fakefmt", lambda b: sentinel if b[:4] == b"FAKE" else None)
    assert has_decoder_for("x.fakefmt") is True
    assert decode_bytes("x.fakefmt", b"FAKE....") is sentinel
    assert decode_bytes("x.fakefmt", b"nope") is None
    # unrelated formats untouched
    assert has_decoder_for("x.mrc") is True


def test_decoded_header_is_frozen(mrc):
    d = decode_bytes("x.mrc", mrc(64, 64, 1, 2))
    with pytest.raises(Exception):
        d.format = "changed"  # frozen dataclass
