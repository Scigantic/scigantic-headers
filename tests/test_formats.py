"""NPY and NIfTI decoders — the non-cryo-EM formats that prove the framework is
format-general. Synthetic headers here; both are separately cross-checked
against numpy and nibabel (the reference libraries) in the benchmark/dev notes.
"""

import struct

import pytest

from scigantic_headers import decode_bytes, decode_nifti_header, decode_npy_header, has_decoder_for


# ── NPY ────────────────────────────────────────────────────────────────────

def make_npy(shape, descr, fortran=False, version=(1, 0)):
    body = "{'descr': %r, 'fortran_order': %s, 'shape': %r, }" % (descr, fortran, tuple(shape))
    if version[0] == 1:
        preamble = 10 + len(body) + 1
        pad = (64 - preamble % 64) % 64
        body = body + " " * pad + "\n"
        return b"\x93NUMPY" + bytes(version) + struct.pack("<H", len(body)) + body.encode("latin1")
    preamble = 12 + len(body) + 1
    pad = (64 - preamble % 64) % 64
    body = body + " " * pad + "\n"
    return b"\x93NUMPY" + bytes(version) + struct.pack("<I", len(body)) + body.encode("latin1")


def test_npy_extension_registered():
    assert has_decoder_for("x/arr.npy") is True


@pytest.mark.parametrize("descr,name", [
    ("<f8", "float64"), ("<f4", "float32"), ("<i4", "int32"),
    ("|u1", "uint8"), ("<i8", "int64"), ("|b1", "bool"), ("<c16", "complex128"),
])
def test_npy_dtype_names(descr, name):
    d = decode_npy_header(make_npy((10, 20), descr))
    assert d.fields["dtype"] == name


def test_npy_shape_and_summary():
    d = decode_npy_header(make_npy((100, 50, 3), "<f8"))
    assert d.format == "npy"
    assert d.fields["shape"] == [100, 50, 3]
    assert d.fields["ndim"] == 3
    assert d.fields["numElements"] == 15000
    assert d.summary == "NPY array 100x50x3, float64"


def test_npy_fortran_order_and_scalar():
    d = decode_npy_header(make_npy((), "<f8"))
    assert d.fields["shape"] == []
    assert d.summary == "NPY array scalar, float64"
    f = decode_npy_header(make_npy((4, 4), "<f4", fortran=True))
    assert f.fields["fortranOrder"] is True
    assert "Fortran-order" in f.summary


def test_npy_version_2():
    d = decode_npy_header(make_npy((8,), "<i2", version=(2, 0)))
    assert d.fields["dtype"] == "int16"
    assert d.fields["shape"] == [8]


def test_npy_rejects_bad_magic():
    assert decode_npy_header(b"NOTNPY" + b"\x00" * 100) is None
    assert decode_bytes("x.npy", b"\x00" * 200) is None


def _npy_raw(header_dict_str):
    body = header_dict_str
    preamble = 10 + len(body) + 1
    body = body + " " * ((64 - preamble % 64) % 64) + "\n"
    return b"\x93NUMPY\x01\x00" + struct.pack("<H", len(body)) + body.encode("latin1")


@pytest.mark.parametrize("bad", [
    "{'descr': '<f8', 'fortran_order': False, 'shape': 5, }",          # shape not a tuple
    "{'descr': '<f8', 'fortran_order': False, 'shape': ('a', 2), }",   # non-int in shape
    "{'descr': '<f8', 'fortran_order': False, 'shape': [1, 2], }",     # list, not tuple
    "{'descr': '<f8', 'fortran_order': False}",                        # no shape key
    "not even a dict",                                                  # not literal-eval-able to dict
])
def test_npy_malformed_header_returns_none_not_raise(bad):
    # Each of these reaches past the magic/length check into header parsing;
    # a total decoder returns None rather than raising.
    assert decode_npy_header(_npy_raw(bad)) is None


# ── NIfTI-1 ──────────────────────────────────────────────────────────────────

def make_nifti(shape, datatype=16, voxel=(1.0, 1.0, 1.0), order="<", magic=b"n+1\x00", sizeof=348):
    buf = bytearray(352)
    struct.pack_into(order + "i", buf, 0, sizeof)
    dim = [len(shape)] + list(shape)
    dim = (dim + [1] * 8)[:8]
    struct.pack_into(order + "8h", buf, 40, *dim)
    struct.pack_into(order + "h", buf, 70, datatype)
    pixdim = ([0.0] + list(voxel) + [0.0] * 8)[:8]
    struct.pack_into(order + "8f", buf, 76, *pixdim)
    buf[344:348] = magic
    return bytes(buf)


def test_nifti_extension_registered():
    assert has_decoder_for("scan.nii") is True


def test_nifti_3d_volume():
    d = decode_nifti_header(make_nifti((182, 218, 182), datatype=16, voxel=(1.0, 1.0, 1.0)))
    assert d.format == "nifti"
    assert d.fields["shape"] == [182, 218, 182]
    assert d.fields["dtype"] == "float32"
    assert d.fields["voxelSizesMm"] == [1.0, 1.0, 1.0]
    assert d.summary == "NIfTI-1 volume 182x218x182, float32, 1.0x1.0x1.0 mm"


def test_nifti_4d_fmri_anisotropic():
    d = decode_nifti_header(make_nifti((64, 64, 36, 200), datatype=4, voxel=(3.0, 3.0, 3.5)))
    assert d.fields["ndim"] == 4
    assert d.fields["shape"] == [64, 64, 36, 200]
    assert d.fields["dtype"] == "int16"
    assert d.fields["voxelSizesMm"] == [3.0, 3.0, 3.5]  # only x,y,z reported


def test_nifti_big_endian():
    d = decode_nifti_header(make_nifti((10, 10, 10), order=">"))
    assert d is not None
    assert d.fields["shape"] == [10, 10, 10]


def test_nifti_rejects_bad_sizeof_and_magic():
    assert decode_nifti_header(make_nifti((10, 10, 10), sizeof=999)) is None
    assert decode_nifti_header(make_nifti((10, 10, 10), magic=b"XXXX")) is None
    assert decode_nifti_header(b"\x00" * 100) is None  # too short
