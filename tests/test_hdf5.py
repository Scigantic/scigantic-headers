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

NWB goldens (hdf5-nwb-*.sparse.json) are three DANDI archive files from three labs
and writers (000127 Miller lab NWB 2.4.0; 000004 Rutishauser lab NWB 2.1.0;
000020 Allen Institute patch-seq written by MIES, NWB 2.2.5 with the ndx-mies
extension, 38 acquisition and 38 stimulus series) plus a pynwb 4.2.0 file. Every
structural value in hdf5-cases.json (version, identifier, session start, subject,
unit and electrode counts, acquisition, processing and stimulus listings, interval
row counts, embedded namespaces) was checked against h5py on the full files. Two
of the four read wrong before the global-heap scan learned to widen past 8 KiB,
which is why they are here.
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


# ── NWB ──────────────────────────────────────────────────────────────────────

def test_nwb_extension_registered():
    from scigantic_headers import has_decoder_for
    assert has_decoder_for("sub-01/sub-01_ses-1_ecephys.nwb")


def test_nwb_session_fields_from_real_files():
    f = _CASES["hdf5-nwb-dandi-000127"]["expected"]["fields"]
    assert f["nwbVersion"] == "2.4.0" and f["namespace"] == "core"
    assert f["sessionStartTime"] == "2017-12-07T17:06:09.727000-05:00" and f["sessionId"] == "20171207"
    assert f["institution"] == "Northwestern University" and f["lab"] == "Miller"
    assert f["experimenter"] == ["Raeed H. Chowdhury"]
    assert f["subject"]["subjectId"] == "Han" and f["subject"]["species"] == "Macaca mulatta" and f["subject"]["sex"] == "M"
    assert f["unitCount"] == 49 and f["electrodeCount"] == 96
    assert f["intervals"] == [{"name": "trials", "rows": 98}]
    assert f["specifications"] == {"core": ["2.4.0"], "hdmf-common": ["1.5.0"], "hdmf-experimental": ["0.1.0"]}
    assert f["extensions"] == []
    assert _CASES["hdf5-nwb-dandi-000127"]["expected"]["summary"] == \
        "NWB 2.4.0, session 2017-12-07, Macaca mulatta, 49 units, 96 electrodes"


def test_nwb_strings_on_a_large_global_heap():
    # 000004 keeps its strings in a global heap collection far larger than the
    # 8 KiB the scan used to read; the subject and the acquisition types were
    # unreadable until the scan widened.
    f = _CASES["hdf5-nwb-dandi-000004"]["expected"]["fields"]
    assert f["subject"]["species"] == "Human" and f["subject"]["subjectId"] == "P10HMH" and f["subject"]["age"] == "37"
    assert f["identifier"] == "H10_7" and f["unitCount"] == 38 and f["electrodeCount"] == 18
    assert [a["neurodataType"] for a in f["acquisition"]] == ["AnnotationSeries", "TimeSeries"]
    assert f["acquisition"][0]["shape"] == [1004]
    assert f["specifications"] == {}  # pynwb of that era did not embed the schema


def test_nwb_detail_stops_before_the_read_cap_and_names_everything():
    # The patch-seq file has 38 acquisition and 38 stimulus series. Describing
    # each costs reads, so the detail pass stops short of the cap, but every
    # series is still named and nothing after it (modules, units, specs) is
    # starved: the ndx-mies extension is seen.
    f = _CASES["hdf5-nwb-dandi-000020"]["expected"]["fields"]
    assert f["nwbVersion"] == "2.2.5" and f["extensions"] == ["ndx-mies"]
    assert f["acquisitionCount"] == 38 and len(f["acquisitionNames"]) == 38
    assert 0 < len(f["acquisition"]) < 38
    # MIES writes scalar attributes as one-element arrays; rate and unit unwrap.
    assert f["acquisition"][0] == {"name": "data_00000_AD0", "neurodataType": "VoltageClampSeries",
                                   "shape": [66000], "dtype": "float32", "unit": "amperes", "rate": 200000.0}
    assert len(f["stimulusPresentationNames"]) == 38 and f["stimulusPresentation"] == []
    assert [m["name"] for m in f["processingModules"]] == ["spikes"]
    assert len(f["processingModules"][0]["interfaceNames"]) == 15
    assert f["subject"]["genotype"].startswith("Vipr2-IRES2-Cre") and f["subject"]["dateOfBirth"] == "2019-12-03T16:00:00-08:00"
    # A one-element string array (how MIES writes identifier) reads as the string.
    assert isinstance(f["identifier"], str) and len(f["identifier"]) == 64


def test_nwb_pynwb_file_is_read_in_full():
    f = _CASES["hdf5-nwb-pynwb-rich"]["expected"]["fields"]
    assert f["sessionId"] == "sess-42" and f["experimenter"] == ["Doe, Jane", "Roe, Rick"] and f["keywords"] == ["V1", "neuropixels"]
    assert f["subject"] == {"subjectId": "M123", "species": "Mus musculus", "sex": "M", "age": "P90D", "genotype": "wt",
                            "strain": "C57BL/6J", "dateOfBirth": None, "weight": None, "description": "a mouse"}
    assert f["devices"] == [{"name": "Camera", "description": "behavior camera", "manufacturer": None},
                            {"name": "Neuropixels 1.0", "description": "probe", "manufacturer": "imec"}]
    assert f["electrodeGroups"] == ["shank0"] and f["electrodeCount"] == 64
    assert f["acquisition"] == [
        {"name": "ElectricalSeries", "neurodataType": "ElectricalSeries", "shape": [1000, 64], "dtype": "int16", "unit": "volts", "rate": 30000.0},
        {"name": "running_speed", "neurodataType": "TimeSeries", "shape": [500], "dtype": "float32", "unit": "cm/s", "rate": 60.0},
    ]
    assert [(m["name"], m["interfaceNames"]) for m in f["processingModules"]] == [("behavior", []), ("ecephys", ["LFP"])]
    assert f["unitCount"] == 120 and f["unitColumns"] == ["quality", "spike_times"]
    assert f["intervals"] == [{"name": "trials", "rows": 30}]
    assert f["specifications"]["core"] == ["2.11.0"] and f["rootComplete"] is True


def test_nwb_is_not_claimed_by_a_plain_hdf5_file():
    # The synthetic generic fixtures have no nwb_version: they stay "hdf5".
    for name in ("hdf5-generic-earliest.h5", "hdf5-generic-latest.h5"):
        d = decode_file(os.path.join(FIXTURES_DIR, name))
        assert d is not None and d.format == "hdf5"


def test_nwb_repeated_lookups_are_memoized():
    # Resolving the same group twice must not spend a second read.
    from scigantic_headers.hdf5 import _H5, _lookup
    src = load_sparse("hdf5-nwb-pynwb-rich")
    h = _H5(src)
    assert h.open()
    root = h.object(h.root)
    _lookup(h, root, "general/subject/species")
    n = src._reads
    _lookup(h, root, "general/subject/species")
    assert src._reads == n


def test_soft_link_in_v1_group_does_not_end_the_listing():
    # hdf5-softlink-earliest.h5 was written by h5py 3.16 with libver 'earliest':
    # the group `series` holds a_first, data, a soft link `electrode`, then
    # starting_time. The soft link has no object address; before, the symbol
    # node listing stopped there and starting_time was never seen.
    d = decode_file(os.path.join(FIXTURES_DIR, "hdf5-softlink-earliest.h5"))
    assert d is not None and d.format == "hdf5"
    paths = [ds["path"] for ds in d.fields["datasets"]]
    assert paths == ["series/a_first", "series/data", "series/starting_time"]
    assert d.fields["groups"] == ["series", "table"]
    assert d.fields["walkComplete"] is True
