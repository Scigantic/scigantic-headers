"""NWB stored in a Zarr store rather than HDF5.

Golden: tests/fixtures/zarr-nwb-cases.json over two real metadata files, both
committed whole. `zarr-nwb-hdmf-v3.zarr.json` is the root of a store written
by hdmf-zarr 0.14.0 with pynwb 4.2.0; `zarr-nwb-dandi-000719.zmetadata` is the
real consolidated metadata of a DANDI 000719 asset, a Zarr v2 store of 3170
keys and 139 acquisition series, fetched from the DANDI S3 bucket.

What a Zarr NWB record cannot carry, and the HDF5 one can, is every scalar
string: in HDF5 those sit in the object header, in Zarr they are chunk data
behind a compressor this library does not carry. The tests below assert the
absence as deliberately as they assert the values.
"""

import json
import os

from scigantic_headers import decode_bytes, decode_file, has_decoder_for
from scigantic_headers.decoders import read_for
from scigantic_headers.nwb_zarr import is_nwb, normalize_v2

FIX = os.path.join(os.path.dirname(__file__), "fixtures")

with open(os.path.join(FIX, "zarr-nwb-cases.json")) as _fh:
    CASES = json.load(_fh)["cases"]


def _bytes(name):
    with open(os.path.join(FIX, name), "rb") as fh:
        return fh.read()


def test_sidecar_names_registered():
    assert has_decoder_for("session.nwb.zarr/.zmetadata")
    assert has_decoder_for("session.nwb.zarr/zarr.json")
    # The consolidated document must parse whole, so it reads more than the
    # per-node sidecars do.
    assert read_for(".zmetadata").leading > read_for(".zattrs").leading


def test_golden_fixtures_match():
    assert len(CASES) == 2
    for name, case in CASES.items():
        decoded = decode_bytes(case["key"], _bytes(name))
        assert decoded is not None and decoded.to_dict() == case["expected"], name


def test_zarr_v3_store_written_by_hdmf_zarr():
    d = decode_bytes("sample.nwb.zarr/zarr.json", _bytes("zarr-nwb-hdmf-v3.zarr.json"))
    assert d.format == "nwb-zarr"
    assert d.summary == ("NWB 2.11.0 on Zarr v3, 2 acquisition, 1 processing module, "
                         "2 units, 8 electrodes")
    f = d.fields
    assert f["nwbVersion"] == "2.11.0" and f["namespace"] == "core"
    assert f["container"] == "zarr" and f["zarrFormat"] == 3
    assert f["consolidatedMetadata"] is True and f["nodeCount"] == 84
    assert f["devices"] == [{"name": "Neuropixels", "description": "probe", "manufacturer": "imec"}]
    assert f["electrodeGroups"] == ["shank0"] and f["electrodeCount"] == 8
    assert f["acquisitionNames"] == ["ElectricalSeries", "running_speed"]
    first = f["acquisition"][0]
    assert first["neurodataType"] == "ElectricalSeries" and first["shape"] == [100, 8]
    assert first["dtype"] == "int16" and first["unit"] == "volts" and first["rate"] == 30000.0
    assert f["unitCount"] == 2 and f["unitColumns"] == ["spike_times"]
    assert f["intervals"][0]["name"] == "trials" and f["intervals"][0]["rows"] == 2
    assert f["processingModules"][0]["interfaces"][0]["rate"] == 50.0
    # The path carries the schema version, which the HDF5 layout does not.
    assert f["specifications"]["core"] == ["2.11.0"]
    assert f["extensions"] == []


def test_real_dandi_v2_store():
    d = decode_bytes("x.nwb.zarr/.zmetadata", _bytes("zarr-nwb-dandi-000719.zmetadata"))
    assert d.format == "nwb-zarr"
    assert d.summary == "NWB 2.6.0 on Zarr v2, 139 acquisition, 1 processing module"
    f = d.fields
    assert f["zarrFormat"] == 2 and f["consolidatedMetadata"] is True
    assert f["nodeCount"] == 1589 and f["acquisitionCount"] == 139
    series = f["acquisition"][0]
    assert series["neurodataType"] == "VoltageClampSeries"
    assert series["shape"] == [161000] and series["dtype"] == "float32"
    assert series["unit"] == "amperes" and series["rate"] == 200000.0
    # An icephys session has no units or electrode table, and this writer
    # cached no schema; all three are reported as absent, not invented.
    assert f["unitCount"] is None and f["electrodeCount"] is None
    assert f["specifications"] == {} and f["extensions"] == []
    assert f["intervals"][0] == {"name": "epochs", "rows": 1316,
                                 "columns": ["start_time", "stop_time", "tags",
                                             "timeseries", "treelevel"]}


def test_string_values_are_named_not_guessed():
    f = decode_bytes("x/zarr.json", _bytes("zarr-nwb-hdmf-v3.zarr.json")).fields
    # The fields the HDF5 reader gets from the object header are chunk data
    # here, so they are absent and the record says where they live.
    for absent in ("identifier", "sessionDescription", "sessionStartTime", "institution"):
        assert absent not in f
    assert "identifier" in f["stringValuesUnread"]
    assert "general/subject/species" in f["stringValuesUnread"]
    assert f["subjectPresent"] is True


def test_store_directory_prefers_consolidated_metadata(tmp_path):
    # A v2 store root holds both files; the directory dispatch must take the
    # consolidated one, which is the only one that names the structure.
    store = tmp_path / "session.nwb.zarr"
    store.mkdir()
    (store / ".zmetadata").write_bytes(_bytes("zarr-nwb-dandi-000719.zmetadata"))
    (store / ".zattrs").write_bytes(json.dumps(
        {"namespace": "core", "neurodata_type": "NWBFile", "nwb_version": "2.6.0"}).encode())
    d = decode_file(str(store))
    assert d is not None and d.format == "nwb-zarr"
    assert d.fields["consolidatedMetadata"] is True and d.fields["acquisitionCount"] == 139


def test_root_attributes_alone_still_name_the_file(tmp_path):
    store = tmp_path / "session.nwb.zarr"
    store.mkdir()
    (store / ".zattrs").write_bytes(json.dumps(
        {"namespace": "core", "neurodata_type": "NWBFile", "nwb_version": "2.6.0",
         "object_id": "abc"}).encode())
    d = decode_file(str(store))
    assert d is not None and d.format == "nwb-zarr"
    assert d.fields["consolidatedMetadata"] is False and d.fields["nodeCount"] is None
    assert "no consolidated metadata" in d.summary


def test_non_nwb_zarr_is_untouched():
    # A plain Zarr group and an OME-Zarr store keep their own formats.
    assert decode_bytes("x/.zattrs", b'{"foo": 1}').format == "zarr-group"
    assert decode_bytes("x/zarr.json", b'{"zarr_format": 3, "node_type": "group"}').format == "zarr-group"
    # A consolidated document that is not NWB declines, so the directory
    # dispatch falls through to the per-node sidecars as it always did.
    assert decode_bytes("x/.zmetadata", json.dumps(
        {"zarr_consolidated_format": 1,
         "metadata": {".zgroup": {"zarr_format": 2}, ".zattrs": {"foo": 1}}}).encode()) is None


def test_helpers_reject_garbage():
    assert is_nwb({}) is False and is_nwb(None) is False
    assert is_nwb({"nwb_version": "2.0.0"}) is True
    assert normalize_v2({}) is None
    assert normalize_v2({"metadata": "not a dict"}) is None
