"""Zarr v2 / v3, OME-Zarr (OME-NGFF 0.4 and 0.5) and N5 sidecar decoding.

Real goldens under tests/fixtures/:
  zarr-idr6001240.zattrs, -0.zarray, .zgroup   IDR idr0062A image 6001240.zarr
      (uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4, OME-NGFF 0.4, CC-BY). Cross-checked
      with zarr 3.4.1 opening the store over HTTP: multiscales version 0.4,
      datasets 0/1/2, level-0 scale [1, 0.5002, 0.3604, 0.3604], channels
      LaminB1/Dapi, array shape (2, 236, 275, 271), chunks (1, 1, 275, 271),
      uint16, Blosc lz4.
  zarr-idr6001240-v3.zarr.json   the same image as Zarr v3 / OME-NGFF 0.5 from
      the ome2024-ngff-challenge share (idr/share/ome2024-ngff-challenge/0.0.5).
  n5-jrc_hela-2-*.attributes.json   OpenOrganelle jrc_hela-2.n5 root and
      em/fibsem-uint16/s0 (janelia-cosem-datasets S3, CC-BY). The N5 goldens
      are read from the JSON itself (12000 x 1600 x 6368 uint16, 64^3 blocks,
      4 x 4 x 5.24 nm), no second library.
"""

import json
import os

import pytest

from scigantic_headers import decode_bytes, decode_file, has_decoder_for, is_decodable
import scigantic_headers.zarr  # noqa: F401  registers the sidecar names
from scigantic_headers.zarr import (
    decode_n5_attributes, decode_zarr_json, decode_zarray, decode_zattrs, decode_zgroup,
    v2_dtype_name,
)

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


def _fx(name):
    with open(os.path.join(FIX, name), "rb") as fh:
        return fh.read()


@pytest.mark.parametrize("key", ["img.zarr/.zattrs", "img.zarr/0/.zarray", "a/.zgroup",
                                 "img.zarr/zarr.json", "vol.n5/attributes.json",
                                 "x/.ZATTRS", "ZARR.JSON"])
def test_sidecar_names_dispatch(key):
    assert has_decoder_for(key) is True
    assert is_decodable(key) is True


def test_unrelated_json_content_is_not_claimed():
    # The .json extension is registered (BIDS sidecars are content-sniffed
    # under it), so dispatch reaches a decoder; no zarr name matches and the
    # content is not BIDS, so nothing is claimed.
    assert decode_bytes("x/config.json", b'{"debug": true, "port": 8080}') is None


# ── IDR 6001240 (OME-NGFF 0.4, Zarr v2) ───────────────────────────────────────

def test_idr_zattrs_is_ome_zarr():
    d = decode_bytes("6001240.zarr/.zattrs", _fx("zarr-idr6001240.zattrs"))
    assert d.format == "ome-zarr"
    assert d.summary == "OME-Zarr 2 channels, 3 levels, 0.3604x0.3604x0.5002 um"
    f = d.fields
    assert f["ngffVersion"] == "0.4"
    assert f["axes"] == ["c", "z", "y", "x"]
    assert f["axisTypes"] == ["channel", "space", "space", "space"]
    assert f["axisUnits"] == [None, "micrometer", "micrometer", "micrometer"]
    assert f["levels"] == 3 and f["level0Path"] == "0"
    assert f["level0Scale"] == [1.0, 0.5002025531914894, 0.3603981534640209, 0.3603981534640209]
    assert f["voxelSizeUm"] == [0.360398, 0.360398, 0.500203]   # x, y, z
    assert f["channels"] == ["LaminB1", "Dapi"] and f["channelCount"] == 2
    assert f["creator"] == "omero-zarr 0.4.0"
    assert f["zarrFormat"] == 2
    assert f["isPlate"] is False and f["isLabels"] is False and f["isBioformats2raw"] is False


def test_idr_zarray():
    d = decode_bytes("6001240.zarr/0/.zarray", _fx("zarr-idr6001240-0.zarray"))
    assert d.format == "zarr-array"
    assert d.summary == "Zarr array 2 x 236 x 275 x 271, uint16, chunks 1 x 1 x 275 x 271"
    assert d.fields == {
        "zarrFormat": 2, "shape": [2, 236, 275, 271], "ndim": 4, "chunks": [1, 1, 275, 271],
        "dtype": "uint16", "compressor": "blosc/lz4", "filters": None, "order": "C",
        "dimensionSeparator": "/", "fillValue": 0,
    }


def test_idr_zgroup():
    d = decode_bytes("6001240.zarr/.zgroup", _fx("zarr-idr6001240.zgroup"))
    assert d.format == "zarr-group" and d.fields == {"zarrFormat": 2}


# ── The same image as Zarr v3 / OME-NGFF 0.5 ─────────────────────────────────

def test_idr_v3_zarr_json_is_ome_zarr_05():
    d = decode_bytes("6001240.zarr/zarr.json", _fx("zarr-idr6001240-v3.zarr.json"))
    assert d.format == "ome-zarr"
    f = d.fields
    assert f["zarrFormat"] == 3 and f["ngffVersion"] == "0.5"
    assert f["axes"] == ["c", "z", "y", "x"]
    assert f["voxelSizeUm"] == [0.360398, 0.360398, 0.500203]
    assert f["channels"] == ["LaminB1", "Dapi"]
    assert f["levels"] == 3
    assert f["creator"] == "omero-zarr"
    # the v2 and v3 carriers describe the image identically where both speak
    v2 = decode_bytes("x/.zattrs", _fx("zarr-idr6001240.zattrs")).fields
    for k in ("axes", "axisTypes", "axisUnits", "levels", "level0Path", "level0Scale",
              "voxelSizeUm", "channels", "channelCount"):
        assert f[k] == v2[k], k


def test_v3_array_node():
    doc = {"zarr_format": 3, "node_type": "array", "shape": [10, 512, 512], "data_type": "uint16",
           "chunk_grid": {"name": "regular", "configuration": {"chunk_shape": [1, 256, 256]}},
           "chunk_key_encoding": {"name": "default"}, "fill_value": 0,
           "codecs": [{"name": "bytes", "configuration": {"endian": "little"}},
                      {"name": "blosc", "configuration": {"cname": "zstd", "clevel": 5}}],
           "dimension_names": ["z", "y", "x"], "attributes": {}}
    d = decode_zarr_json(json.dumps(doc).encode())
    assert d.format == "zarr-array"
    assert d.summary == "Zarr v3 array 10 x 512 x 512, uint16, chunks 1 x 256 x 256"
    assert d.fields["codecs"] == ["bytes", "blosc/zstd"]
    assert d.fields["compressor"] == "blosc/zstd"
    assert d.fields["dimensionNames"] == ["z", "y", "x"]


def test_v3_array_with_ome_attributes_folds_both():
    doc = {"zarr_format": 3, "node_type": "array", "shape": [2, 50, 100, 100], "data_type": "float32",
           "chunk_grid": {"name": "regular", "configuration": {"chunk_shape": [1, 1, 100, 100]}},
           "codecs": [], "fill_value": "NaN",
           "attributes": {"ome": {"version": "0.5", "multiscales": [{
               "axes": [{"name": "c", "type": "channel"}, {"name": "z", "type": "space", "unit": "nanometer"},
                        {"name": "y", "type": "space", "unit": "nanometer"}, {"name": "x", "type": "space", "unit": "nanometer"}],
               "datasets": [{"path": "0", "coordinateTransformations": [{"type": "scale", "scale": [1, 500, 100, 100]}]}]}]}}}
    d = decode_zarr_json(json.dumps(doc).encode())
    assert d.format == "ome-zarr"
    assert d.fields["shape"] == [2, 50, 100, 100] and d.fields["dtype"] == "float32"
    assert d.fields["voxelSizeUm"] == [0.1, 0.1, 0.5]       # nm -> um
    assert d.fields["fillValue"] == "NaN"
    assert d.summary == "OME-Zarr 2 x 50 x 100 x 100 (czyx), 1 levels, float32, 0.1x0.1x0.5 um"
    json.dumps(d.to_dict(), allow_nan=False)


def test_v3_bare_group():
    d = decode_zarr_json(b'{"zarr_format": 3, "node_type": "group", "attributes": {"a": 1}}')
    assert d.format == "zarr-group" and d.fields == {"zarrFormat": 3, "attributeKeys": ["a"]}


# ── Plain v2 attrs, dtype names, rejections ──────────────────────────────────

def test_plain_zattrs_without_multiscales_is_a_group():
    d = decode_zattrs(b'{"description": "x", "plate": {"rows": []}}')
    assert d.format == "zarr-group"
    assert d.fields["attributeKeys"] == ["description", "plate"] and d.fields["isPlate"] is True


def test_ngff_03_string_axes():
    doc = {"multiscales": [{"version": "0.3", "axes": ["t", "c", "z", "y", "x"],
                            "datasets": [{"path": "0"}, {"path": "1"}]}]}
    d = decode_zattrs(json.dumps(doc).encode())
    assert d.format == "ome-zarr"
    assert d.fields["axes"] == ["t", "c", "z", "y", "x"]
    assert d.fields["axisTypes"] == ["time", "channel", "space", "space", "space"]
    assert d.fields["voxelSizeUm"] == [None, None, None]   # no scale, no unit: unknown, not 1


@pytest.mark.parametrize("descr,name", [
    ("<u2", "uint16"), ("|u1", "uint8"), ("<f4", "float32"), (">f8", "float64"),
    ("<i8", "int64"), ("|b1", "bool"), ("<c16", "complex128"), ("|S10", "|S10"), (None, None),
])
def test_v2_dtype_names(descr, name):
    assert v2_dtype_name(descr) == name


@pytest.mark.parametrize("fn,blob", [
    (decode_zattrs, b"not json"), (decode_zattrs, b"[1, 2]"), (decode_zattrs, b""),
    (decode_zarray, b'{"zarr_format": 2}'),                 # no shape
    (decode_zarray, b'{"zarr_format": 3, "shape": [1]}'),   # wrong version for .zarray
    (decode_zarray, b'{"zarr_format": 2, "shape": ["a"]}'), # non-int dims
    (decode_zgroup, b'{"zarr_format": 3}'),
    (decode_zarr_json, b'{"zarr_format": 2}'),
    (decode_zarr_json, b'{"zarr_format": 3, "node_type": "array"}'),  # no shape
    (decode_zarr_json, b'{"zarr_format": 3, "node_type": "thing"}'),
    (decode_n5_attributes, b'{"hello": 1}'),                 # neither root nor dataset
    (decode_n5_attributes, b'{"dimensions": [1, -2]}'),
])
def test_rejections_return_none(fn, blob):
    assert fn(blob) is None


# ── N5 ────────────────────────────────────────────────────────────────────────

def test_n5_root():
    d = decode_bytes("jrc_hela-2.n5/attributes.json", _fx("n5-jrc_hela-2-root.attributes.json"))
    assert d.format == "n5" and d.summary == "N5 container v2.0.0"
    assert d.fields == {"n5Version": "2.0.0", "attributeKeys": []}


def test_n5_dataset_openorganelle():
    d = decode_bytes("s0/attributes.json", _fx("n5-jrc_hela-2-em-s0.attributes.json"))
    assert d.format == "n5"
    assert d.summary == "N5 dataset 12000 x 1600 x 6368, uint16, blocks 64 x 64 x 64, 0.004x0.004x0.00524 um"
    assert d.fields == {
        "n5Version": None, "shape": [12000, 1600, 6368], "ndim": 3, "blockSize": [64, 64, 64],
        "dtype": "uint16", "compression": "gzip", "pixelResolution": [4.0, 4.0, 5.24],
        "pixelResolutionUnit": "nm", "voxelSizeUm": [0.004, 0.004, 0.00524],
        "axes": ["z", "y", "x"], "attributeKeys": [],
    }


def test_n5_dataset_transform_only_and_v1_compression_string():
    doc = {"dimensions": [100, 200, 50], "blockSize": [32, 32, 32], "dataType": "float32",
           "compressionType": "raw",
           "transform": {"axes": ["x", "y", "z"], "scale": [0.5, 0.5, 2.0], "units": ["um", "um", "um"]},
           "downsamplingFactors": [1, 1, 1], "multiScale": True}
    d = decode_n5_attributes(json.dumps(doc).encode())
    assert d.fields["compression"] == "raw"
    assert d.fields["voxelSizeUm"] == [0.5, 0.5, 2.0]
    assert d.fields["downsamplingFactors"] == [1, 1, 1]
    assert d.fields["multiScale"] is True
    assert d.fields["attributeKeys"] == ["compressionType", "multiScale"]


# ── A store is a directory: decode_file resolves the sidecar ─────────────────

def test_decode_file_on_a_zarr_store_directory(tmp_path):
    store = tmp_path / "img.zarr"
    (store / "0").mkdir(parents=True)
    (store / ".zgroup").write_bytes(_fx("zarr-idr6001240.zgroup"))
    (store / ".zattrs").write_bytes(_fx("zarr-idr6001240.zattrs"))
    (store / "0" / ".zarray").write_bytes(_fx("zarr-idr6001240-0.zarray"))
    top = decode_file(str(store))
    assert top.format == "ome-zarr" and top.fields["channels"] == ["LaminB1", "Dapi"]
    arr = decode_file(str(store / "0"))
    assert arr.format == "zarr-array" and arr.fields["shape"] == [2, 236, 275, 271]


def test_decode_file_on_a_v3_store_prefers_zarr_json(tmp_path):
    store = tmp_path / "img.zarr"
    store.mkdir()
    (store / "zarr.json").write_bytes(_fx("zarr-idr6001240-v3.zarr.json"))
    (store / ".zattrs").write_bytes(b'{"stale": true}')
    d = decode_file(str(store))
    assert d.format == "ome-zarr" and d.fields["zarrFormat"] == 3


def test_decode_file_on_an_n5_dataset_directory(tmp_path):
    ds = tmp_path / "vol.n5" / "em" / "s0"
    ds.mkdir(parents=True)
    (ds / "attributes.json").write_bytes(_fx("n5-jrc_hela-2-em-s0.attributes.json"))
    (tmp_path / "vol.n5" / "attributes.json").write_bytes(_fx("n5-jrc_hela-2-root.attributes.json"))
    assert decode_file(str(ds)).fields["shape"] == [12000, 1600, 6368]
    assert decode_file(str(tmp_path / "vol.n5")).summary == "N5 container v2.0.0"


def test_decode_file_on_a_directory_without_sidecars(tmp_path):
    (tmp_path / "plain").mkdir()
    assert decode_file(str(tmp_path / "plain")) is None


def test_all_fixtures_are_json_safe():
    for name, key in [("zarr-idr6001240.zattrs", ".zattrs"), ("zarr-idr6001240-0.zarray", ".zarray"),
                      ("zarr-idr6001240-v3.zarr.json", "zarr.json"),
                      ("n5-jrc_hela-2-em-s0.attributes.json", "attributes.json")]:
        json.dumps(decode_bytes("x/" + key, _fx(name)).to_dict(), allow_nan=False)


# ── a whole store directory: group plus level-0 array ────────────────────────

def test_store_directory_merges_level0_array(tmp_path):
    import shutil
    store = tmp_path / "img.zarr"
    (store / "0").mkdir(parents=True)
    shutil.copy(os.path.join(FIX, "zarr-idr6001240.zattrs"), store / ".zattrs")
    shutil.copy(os.path.join(FIX, "zarr-idr6001240-0.zarray"), store / "0" / ".zarray")
    h = decode_file(str(store))
    assert h.format == "ome-zarr"
    assert h.fields["shape"] == [2, 236, 275, 271]
    assert (h.fields["sizeC"], h.fields["sizeZ"], h.fields["sizeY"], h.fields["sizeX"]) == (2, 236, 275, 271)
    assert h.fields["dtype"] == "uint16" and h.fields["chunks"] == [1, 1, 275, 271]
    assert h.summary.startswith("OME-Zarr ") and "uint16" in h.summary
    # Without the array directory the group still decodes, just without a shape.
    shutil.rmtree(store / "0")
    g = decode_file(str(store))
    assert g.format == "ome-zarr" and "shape" not in g.fields
