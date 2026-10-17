"""HDF5 structure walker and the imaging flavors built on it (Imaris, BigDataViewer,
Luxendo, generic).

Real-file goldens are sparse recordings (tests/fixtures/hdf5-*.sparse.json) of OME
sample files: Imaris-IMS/davemason/Convallaria_3C_1T_confocal.ims (+ its IMS3 TIFF
twin), Imaris-IMS/eliana/retina_large.ims, BDV/samples/drosophila.h5. The expected
values in hdf5-cases.json were cross-checked against h5py 3.x on the same files
(dataset shapes, dtypes, chunks, compression, DataSetInfo attribute strings, BDV
resolutions/subdivisions tables). The synthetic .h5 fixtures were written by h5py
with libver 'earliest' (superblock v0, v1 object headers, symbol tables) and
'latest' (superblock v3, OHDR, dense links and attributes); when h5py is installed
the test re-derives their expectations live, otherwise the committed values stand.
"""

import json
import os
import struct

import pytest

from scigantic_headers import BytesSource, decode_bytes, decode_file, decode_hdf5, decode_source
from conftest import FIXTURES_DIR, load_sparse

with open(os.path.join(FIXTURES_DIR, "hdf5-cases.json")) as _fh:
    _CASES = json.load(_fh)


@pytest.mark.parametrize("name", sorted(_CASES))
def test_real_files_match_golden(name):
    case = _CASES[name]
    hdr = decode_source(case["key"], load_sparse(name))
    assert hdr is not None
    assert hdr.to_dict() == case["expected"]


def test_imaris_fields_are_the_physical_ones():
    f = _CASES["hdf5-imaris-retina"]["expected"]["fields"]
    assert (f["sizeX"], f["sizeY"], f["sizeZ"]) == (2048, 1567, 64)
    assert f["channels"] == ["CollagenIV (TxRed)", "GFAP (FITC)"]
    assert f["resolutionLevels"] == 4 and f["dtype"] == "uint8" and f["compression"] == ["deflate"]
    assert f["voxelSizeUm"] == [0.022906, 0.022899, 0.2]   # (ExtMax - ExtMin) / size, in um
    assert f["unit"] == "um" and f["objective"] == "63x, 1.3NA"
    g = _CASES["hdf5-imaris-convallaria-1t"]["expected"]["fields"]
    assert g["channels"] == ["Red", "Green", "Blue"] and g["microscopeMode"] == "SpinningDiskConfocal"
    assert g["numericalAperture"] == 0.45 and g["recordingDate"] == "2018-10-22 13:09:30.696"


def test_bdv_fields():
    f = _CASES["hdf5-bdv-drosophila"]["expected"]["fields"]
    assert (f["setupCount"], f["timepointCount"]) == (2, 3)
    assert (f["sizeX"], f["sizeY"], f["sizeZ"]) == (1920, 1060, 130)   # cells are stored z, y, x
    assert f["resolutions"] == [[1.0, 1.0, 1.0], [2.0, 2.0, 1.0], [4.0, 4.0, 2.0], [8.0, 8.0, 4.0]]
    assert f["subdivisions"][0] == [32, 32, 4]
    assert f["dtype"] == "int16" and f["compression"] == ["scaleoffset", "deflate"]


def test_ims3_is_a_tiff_container():
    hdr = _CASES["hdf5-imaris-ims3-tiff"]["expected"]
    assert hdr["format"] == "imaris-ims3"
    assert hdr["fields"]["container"] == "tiff"
    assert hdr["fields"]["width"] == 512


# ── synthetic h5py files (committed) ────────────────────────────────────────

def _fixture(name):
    return os.path.join(FIXTURES_DIR, name)


def test_generic_earliest_layout():
    h = decode_file(_fixture("hdf5-generic-earliest.h5"))
    assert h.format == "hdf5"
    assert h.fields["rootAttributes"] == {"title": "generic test file", "version": 3, "element_size_um": [2.0, 0.5, 0.5]}
    assert h.fields["groups"] == ["acq", "acq/raw", "empty"]
    assert h.fields["walkComplete"] is True
    ds = {d["path"]: d for d in h.fields["datasets"]}
    assert ds["acq/raw/stack"]["shape"] == [12, 64, 80]
    assert ds["acq/raw/stack"]["dtype"] == "uint16"
    assert ds["acq/raw/stack"]["chunks"] == [4, 32, 32]
    assert ds["acq/raw/stack"]["filters"] == ["shuffle", "deflate"]
    assert ds["acq/raw/timestamps"] == {"shape": [12], "dtype": "float64", "layout": "contiguous",
                                        "chunks": None, "filters": [], "path": "acq/raw/timestamps"}
    assert ds["notes"]["dtype"] == "string(5)"
    assert h.summary == "HDF5, acq/raw/stack 12x64x80, uint16, 3 groups, 3 datasets"


def test_generic_latest_layout_dense_links_and_attributes():
    h = decode_file(_fixture("hdf5-generic-latest.h5"))
    assert h.fields["rootAttributes"]["name"] == "latest layout"           # vlen string via the global heap
    assert [h.fields["rootAttributes"]["attr%02d" % i] for i in range(12)] == list(range(12))
    assert h.fields["groups"] == ["many"] + ["many/child%02d" % i for i in range(10)]
    assert h.fields["walkComplete"] is True
    assert h.fields["shape"] == [5, 6, 7] and h.fields["dtype"] == "float32"


def test_luxendo_by_extension():
    h = decode_file(_fixture("hdf5-luxendo-shaped.lux.h5"))
    assert h.format == "luxendo-h5"
    assert (h.fields["sizeX"], h.fields["sizeY"], h.fields["sizeZ"]) == (2048, 2048, 300)
    assert h.fields["firstDataset"] == "Data"
    # The same bytes under a plain .h5 name are generic HDF5.
    data = open(_fixture("hdf5-luxendo-shaped.lux.h5"), "rb").read()
    assert decode_bytes("x.h5", data).format == "hdf5"
    assert decode_bytes("x.lux.h5", data).format == "luxendo-h5"


h5py = pytest.importorskip("h5py", reason="live cross-check needs h5py") if os.environ.get("SCIGANTIC_HEADERS_H5PY") else None


@pytest.mark.skipif(h5py is None, reason="set SCIGANTIC_HEADERS_H5PY=1 with h5py installed to cross-check live")
def test_cross_check_against_h5py():
    for name in ("hdf5-generic-earliest.h5", "hdf5-generic-latest.h5", "hdf5-luxendo-shaped.lux.h5"):
        ours = decode_file(_fixture(name))
        with h5py.File(_fixture(name), "r") as f:
            theirs = {}
            f.visititems(lambda p, o: theirs.__setitem__(p, o) if isinstance(o, h5py.Dataset) else None)
            for d in ours.fields["datasets"]:
                ref = theirs[d["path"]]
                assert d["shape"] == list(ref.shape)
                assert d["chunks"] == (list(ref.chunks) if ref.chunks else None)


# ── rejection and totality ──────────────────────────────────────────────────

def test_rejects_non_hdf5():
    assert decode_hdf5(BytesSource(b"\x00" * 4096)) is None
    assert decode_hdf5(BytesSource(b"\x89HDF\r\n\x1a\n")) is None            # signature only
    assert decode_hdf5(BytesSource(b"\x89HDF\r\n\x1a\n" + b"\xff" * 2000)) is None  # garbage superblock
    assert decode_bytes("x.h5", b"") is None


def test_superblock_after_user_block():
    data = open(_fixture("hdf5-generic-earliest.h5"), "rb").read()
    # Relocating a file behind a 512-byte user block is how HDF5 reserves a
    # leading region; the superblock is then found at 512 and addresses are
    # relative to it.
    shifted = b"U" * 512 + data
    h = decode_hdf5(BytesSource(shifted))
    assert h is not None and h.fields["groups"] == ["acq", "acq/raw", "empty"]


def test_truncated_file_degrades_to_none_or_partial():
    data = open(_fixture("hdf5-generic-earliest.h5"), "rb").read()
    for cut in (9, 64, 200, 1024, 3000):
        h = decode_hdf5(BytesSource(data[:cut]))
        assert h is None or h.format == "hdf5"
        if h is not None:
            json.dumps(h.to_dict(), allow_nan=False)


def test_budget_exhaustion_is_reported_not_raised():
    data = open(_fixture("hdf5-generic-latest.h5"), "rb").read()
    h = decode_hdf5(BytesSource(data, budget=600))
    assert h is None or h.fields["walkComplete"] is False
