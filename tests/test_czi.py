"""Zeiss CZI decoding: synthetic segment layouts for the reject and
degraded paths, and sparse recordings of real files as golden fixtures.

Real-file expectations were cross-checked with czifile 2026.10.14 (header
pointers, sub-block and attachment directory counts, and every XML field used)
on the two whole files available locally, xt-scan-lsm980.czi and
Image_1_2023_08_18__14_32_31_964.czi (`czi-image1`). The two head-only
light-sheet recordings (`czi-mousebrain-z1`, `czi-lls7max`) were checked
against their raw metadata XML, since czifile needs the whole file.
"""

import json
import struct

import pytest

import scigantic_headers.czi  # noqa: F401  (registers the decoder)
from scigantic_headers.czi import decode_czi
from scigantic_headers.decoders import BytesSource, decode_bytes, has_decoder_for, read_for

from conftest import load_sparse


# ── synthetic files ───────────────────────────────────────────────────────────

def _segment(sid: bytes, payload: bytes, allocated=None) -> bytes:
    allocated = len(payload) if allocated is None else allocated
    return sid.ljust(16, b"\x00") + struct.pack("<qq", allocated, len(payload)) + payload.ljust(allocated, b"\x00")


def make_czi(xml: bytes, *, file_part=0, meta_pos=None, dir_entries=2, att_entries=1, magic=b"ZISRAWFILE"):
    """A minimal CZI: file header, then metadata, directory, attachment-dir segments."""
    meta_payload = struct.pack("<ii", len(xml), 0) + b"\x00" * 248 + xml
    meta = _segment(b"ZISRAWMETADATA", meta_payload)
    directory = _segment(b"ZISRAWDIRECTORY", struct.pack("<i", dir_entries) + b"\x00" * 124)
    attdir = _segment(b"ZISRAWATTDIR", struct.pack("<i", att_entries) + b"\x00" * 252)
    head_len = 32 + 512
    m_pos = head_len if meta_pos is None else meta_pos
    d_pos = head_len + len(meta)
    a_pos = d_pos + len(directory)
    hdr = bytearray(512)
    struct.pack_into("<ii", hdr, 0, 1, 0)
    struct.pack_into("<i", hdr, 48, file_part)
    struct.pack_into("<q", hdr, 52, d_pos)
    struct.pack_into("<q", hdr, 60, m_pos)
    struct.pack_into("<q", hdr, 72, a_pos)
    return _segment(magic, bytes(hdr)) + meta + directory + attdir


XML = b"""<ImageDocument>
  <Metadata>
    <Information>
      <Image>
        <SizeX>2048</SizeX><SizeY>2048</SizeY><SizeZ>200</SizeZ><SizeC>2</SizeC><SizeT>1</SizeT>
        <SizeI>2</SizeI><SizeV>4</SizeV>
        <PixelType>Gray16</PixelType><ComponentBitCount>16</ComponentBitCount>
        <AcquisitionDateAndTime>2026-01-02T03:04:05Z</AcquisitionDateAndTime>
        <Dimensions><Channels>
          <Channel Id="Channel:0" Name="GFP"><ExposureTime>20000000</ExposureTime><AcquisitionMode>SPIM</AcquisitionMode></Channel>
          <Channel Id="Channel:1" Name="RFP"><ExposureTime>50000000</ExposureTime><AcquisitionMode>SPIM</AcquisitionMode></Channel>
        </Channels></Dimensions>
      </Image>
      <Instrument>
        <Microscopes><Microscope Id="Microscope:0"><System>Lightsheet 7</System></Microscope></Microscopes>
        <Objectives><Objective Id="Objective:0" Name="W Plan-Apochromat 20x/1.0"><LensNA>1.0</LensNA><NominalMagnification>20</NominalMagnification><Immersion>Water</Immersion></Objective></Objectives>
      </Instrument>
      <Application><Name>ZEN</Name><Version>3.5</Version></Application>
    </Information>
    <Scaling><Items>
      <Distance Id="X"><Value>2.2E-07</Value></Distance>
      <Distance Id="Y"><Value>2.2E-07</Value></Distance>
      <Distance Id="Z"><Value>5E-07</Value></Distance>
    </Items></Scaling>
  </Metadata>
</ImageDocument>"""


def test_registered_as_seeking():
    assert has_decoder_for("x/img.czi") and has_decoder_for("IMG.CZI")
    assert read_for("img.czi").seek == 4 * 1024 * 1024


def test_synthetic_lightsheet_file():
    d = decode_bytes("a.czi", make_czi(XML))
    assert d.format == "czi"
    f = d.fields
    assert (f["sizeX"], f["sizeY"], f["sizeZ"], f["sizeC"], f["sizeT"]) == (2048, 2048, 200, 2, 1)
    assert f["illuminations"] == 2 and f["views"] == 4 and f["tiles"] is None
    assert f["dtype"] == "uint16" and f["samplesPerPixel"] == 1 and f["pixelType"] == "Gray16"
    assert f["voxelSizeUm"] == [0.22, 0.22, 0.5]
    assert f["channels"] == ["GFP", "RFP"] and f["channelCount"] == 2
    assert f["exposureTimesMs"] == [20.0, 50.0]            # nanoseconds in the file
    assert f["acquisitionMode"] == "SPIM"
    assert f["microscope"] == "Lightsheet 7" and f["objective"] == "W Plan-Apochromat 20x/1.0"
    assert f["magnification"] == 20.0 and f["numericalAperture"] == 1.0 and f["immersion"] == "Water"
    assert f["software"] == "ZEN 3.5" and f["acquisitionDate"] == "2026-01-02T03:04:05Z"
    assert f["subBlockCount"] == 2 and f["attachmentCount"] == 1 and f["hasAttachments"] is True
    assert f["filePart"] == 0 and "xmlTruncated" not in f
    assert d.summary == "CZI image 2048x2048x200, 2 channels, 4 views, 2 illuminations, uint16, 0.22x0.22x0.5 um"


def test_rejects_non_czi_and_short():
    assert decode_czi(BytesSource(b"\x00" * 600)) is None
    assert decode_czi(BytesSource(make_czi(XML, magic=b"ZISRAWSUBBLOCK"))) is None
    assert decode_czi(BytesSource(make_czi(XML)[:300])) is None
    assert decode_bytes("a.czi", b"") is None


def test_metadata_beyond_file_end_still_decodes_header():
    d = decode_czi(BytesSource(make_czi(XML, meta_pos=10 ** 12)))
    assert d is not None and d.fields["hasMetadata"] is False
    assert d.fields["sizeX"] is None and d.fields["subBlockCount"] == 2
    assert d.summary == "CZI image (no metadata segment)"


def test_file_part_is_reported():
    d = decode_czi(BytesSource(make_czi(XML, file_part=3)))
    assert d.fields["filePart"] == 3 and d.summary.endswith(", part 3")


def test_truncated_xml_falls_back_to_regex():
    # The budget cuts the XML mid-document: ElementTree fails, the regex path
    # still recovers shape, pixel type, scaling and channel names.
    whole = make_czi(XML)
    cut = whole[: 32 + 512 + 32 + 256 + len(XML) - 40]  # lose the closing tags
    d = decode_czi(BytesSource(cut))
    f = d.fields
    assert f["xmlTruncated"] is True
    assert (f["sizeX"], f["sizeZ"], f["views"]) == (2048, 200, 4)
    assert f["dtype"] == "uint16" and f["voxelSizeUm"] == [0.22, 0.22, 0.5]
    assert f["channels"] == ["GFP", "RFP"]
    json.dumps(d.to_dict(), allow_nan=False)


def test_zero_scaling_is_unknown_not_zero():
    xml = XML.replace(b"<Distance Id=\"Z\"><Value>5E-07</Value></Distance>",
                      b"<Distance Id=\"Z\"><Value>0</Value></Distance>")
    assert decode_czi(BytesSource(make_czi(xml))).fields["voxelSizeUm"] == [0.22, 0.22, None]


def test_json_safe_with_nan_in_xml():
    xml = XML.replace(b"<LensNA>1.0</LensNA>", b"<LensNA>NaN</LensNA>")
    d = decode_czi(BytesSource(make_czi(xml)))
    assert d.fields["numericalAperture"] is None
    json.dumps(d.to_dict(), allow_nan=False)


# ── real files, replayed from sparse recordings ───────────────────────────────

REAL = {
    # Lightsheet Z.1 (ZEN black / AIMApplication), whole 7 MB file; the writer
    # omits Size* for singleton dimensions and writes Z scaling as 0.
    "czi-image1": {
        "summary": "CZI image 1920x1920, 4 tiles, uint16, 1.822x1.822 um",
        "sizeX": 1920, "sizeY": 1920, "sizeZ": None, "sizeC": None, "tiles": 4,
        "voxelSizeUm": [1.822076, 1.822076, None], "channels": ["Cam1"], "exposureTimesMs": [30.0],
        "acquisitionMode": "SPIM", "microscope": "Lightsheet Z.1",
        "objective": "Fluar 2.5x/0.12 (User Medium : n=1.530)", "magnification": 5.0, "numericalAperture": 0.12,
        "software": "AIMApplication 9.3.6.393", "subBlockCount": 1, "attachmentCount": 3,
        "metadataPosition": 3296, "metadataBytes": 16489, "dtype": "uint16",
    },
    # Lightsheet Z.1 multi-view: metadata segment at byte 492,122,912 of a
    # 492 MB file (past the pixel data); the directory was not recorded, so its
    # counts read as unknown here, as they would on a server refusing that Range.
    "czi-mousebrain-z1": {
        "summary": "CZI image 1500x1000x41, 2 views, 2 illuminations, uint16, 2.567x2.567x7 um",
        "sizeX": 1500, "sizeY": 1000, "sizeZ": 41, "views": 2, "illuminations": 2,
        "voxelSizeUm": [2.566763, 2.566763, 7.0], "channels": ["Cam2"], "exposureTimesMs": [99.867],
        "acquisitionMode": "SPIM", "microscope": "Lightsheet Z.1", "magnification": 5.0,
        "subBlockCount": None, "attachmentCount": None, "hasAttachments": True,
        "metadataPosition": 492122912, "dtype": "uint16",
    },
    # Lattice Lightsheet 7 (ZEN blue 3.12), a max projection: Z is 1 and unscaled.
    "czi-lls7max": {
        "summary": "CZI image 754x707, 3 channels, 60 timepoints, uint16, 0.145x0.145 um",
        "sizeX": 754, "sizeY": 707, "sizeZ": 1, "sizeC": 3, "sizeT": 60,
        "voxelSizeUm": [0.144992, 0.144992, None], "channels": ["LLS", "LLS_2", "LLS_3"],
        "exposureTimesMs": [15.0, 15.0, 15.0], "acquisitionMode": "LatticeLightsheet",
        "microscope": "Lattice Lightsheet", "objective": "Objective Reflected Light",
        "magnification": 13.3, "numericalAperture": 0.44, "immersion": "Water",
        "software": "ZEN 3.12.107.07000", "acquisitionDate": "2025-08-28T11:12:33.5594323Z",
        "subBlockCount": 180, "attachmentCount": 4, "metadataBytes": 455802, "dtype": "uint16",
    },
}


@pytest.mark.parametrize("name", sorted(REAL))
def test_real_file_recording(name):
    d = decode_czi(load_sparse(name))
    assert d is not None and d.format == "czi"
    want = dict(REAL[name])
    assert d.summary == want.pop("summary")
    for k, v in want.items():
        assert d.fields[k] == v, k
    assert "xmlTruncated" not in d.fields
    json.dumps(d.to_dict(), allow_nan=False)


def test_real_file_decode_is_bounded():
    # Five reads at most, a few hundred KB, whatever the file size.
    from conftest import RecordingSource
    r = RecordingSource(load_sparse("czi-lls7max"))
    decode_czi(r)
    assert len(r.ranges) <= 5
    assert sum(len(b) for b in r.ranges.values()) < 500 * 1024
