"""Leica LIF decoder.

Golden: tests/fixtures/lif-pr2729.sparse.json, recorded from the OME sample
Leica-LIF/michael/PR2729_frameOrderCombinedScanTypes.lif (CC-BY). The
reference values come from readlif 0.6.6 on the same file: dims x=64, y=64,
z=3, t=2, m=4 (tiles), 2 channels, bit depth 8, scale 0.252 px/um in x and y
(3.968 um/px), 1/0.16058 = 6.2276 um in z; settings ObjectiveName
"C PLAN  10.0x0.22 POL", Magnification 10, NumericalAperture 0.22, ScanMode
xyzt, Immersion DRY.
"""

import json
import os
import struct

import scigantic_headers.lif  # noqa: F401  registers the decoder
from scigantic_headers import decode_bytes, has_decoder_for, read_for
from scigantic_headers.decoders import BytesSource
from scigantic_headers.lif import decode_lif

from conftest import FIXTURES_DIR, load_sparse


def _lif(xml: str, n_chars=None) -> bytes:
    body = xml.encode("utf-16-le")
    n = len(xml) if n_chars is None else n_chars
    return struct.pack("<iiBi", 0x70, len(body) + 5, 0x2A, n) + body


_XML = (
    '<LMSDataContainerHeader Version="2"><Element Name="Exp"><Children>'
    '<Element Name="Img1"><Data><Image><ImageDescription>'
    '<Channels><ChannelDescription DataType="0" Resolution="16" LUTName="Green"/>'
    '<ChannelDescription DataType="0" Resolution="16" LUTName="Red"/></Channels>'
    '<Dimensions><DimensionDescription DimID="1" NumberOfElements="512" Length="1.022e-4" Unit="m"/>'
    '<DimensionDescription DimID="2" NumberOfElements="256" Length="5.1e-5" Unit="m"/>'
    '<DimensionDescription DimID="3" NumberOfElements="11" Length="5e-6" Unit="m"/>'
    '<DimensionDescription DimID="4" NumberOfElements="3" Length="8.2" Unit="s"/>'
    '</Dimensions></ImageDescription>'
    '<Attachment Name="HardwareSetting" Application="LAS X" Software="LAS X 3.5" SystemTypeName="SP8">'
    '<ATLConfocalSettingDefinition ObjectiveName="HC PL APO 63x/1.40 OIL" Magnification="63" '
    'NumericalAperture="1.4" Immersion="OIL" ScanMode="xyzt" MicroscopeModel="DMI8"/>'
    '</Attachment></Image></Data></Element>'
    '<Element Name="Img2"><Data><Image><ImageDescription><Channels>'
    '<ChannelDescription DataType="1" Resolution="32"/></Channels><Dimensions>'
    '<DimensionDescription DimID="1" NumberOfElements="8" Length="1e-6" Unit="m"/>'
    '</Dimensions></ImageDescription></Image></Data></Element>'
    '</Children></Element></LMSDataContainerHeader>'
)


def test_registered_and_seeking():
    assert has_decoder_for("exp.lif")
    assert read_for("exp.lif").seek == 8 * 1024 * 1024


def test_synthetic_two_images():
    d = decode_bytes("e.lif", _lif(_XML))
    assert d.format == "lif"
    f = d.fields
    assert f["lifVersion"] == 2 and f["imageCount"] == 2 and f["images"] == ["Img1", "Img2"]
    assert (f["sizeX"], f["sizeY"], f["sizeZ"], f["sizeT"], f["sizeC"]) == (512, 256, 11, 3, 2)
    assert f["dtype"] == "uint16" and f["channels"] == ["Green", "Red"]
    assert f["voxelSizeUm"] == [0.2, 0.2, 0.5]
    assert f["objective"] == "HC PL APO 63x/1.40 OIL" and f["magnification"] == 63.0
    assert f["numericalAperture"] == 1.4 and f["immersion"] == "OIL" and f["scanMode"] == "xyzt"
    assert f["microscope"] == "DMI8" and f["software"] == "LAS X 3.5"
    assert f["xmlTruncated"] is False
    assert d.summary == "Leica LIF 512x256x11, 2 channels, 3 timepoints, uint16, 0.2x0.2x0.5 um, 2 images"


def test_float_channel_dtype_and_missing_hardware():
    xml = _XML.replace('Name="Img1"', 'Name="Skip"').replace("<Data><Image><ImageDescription><Channels>"
                                                             '<ChannelDescription DataType="0"', "<Data><Image><ImageDescription><Channels>"
                                                             '<ChannelDescription DataType="1"', 1)
    # Reorder so the float image is first: decode the second element only.
    xml2 = ('<LMSDataContainerHeader Version="2"><Element Name="Only"><Data><Image><ImageDescription>'
            '<Channels><ChannelDescription DataType="1" Resolution="32"/></Channels><Dimensions>'
            '<DimensionDescription DimID="1" NumberOfElements="8" Length="1e-6" Unit="m"/>'
            '<DimensionDescription DimID="2" NumberOfElements="1" Length="0" Unit="m"/>'
            '</Dimensions></ImageDescription></Image></Data></Element></LMSDataContainerHeader>')
    d = decode_lif(BytesSource(_lif(xml2)))
    assert d.fields["dtype"] == "float32" and d.fields["channels"] == ["Ch1"]
    assert d.fields["voxelSizeUm"] == [round(1e-6 / 7 * 1e6, 6), None, None]
    assert d.fields["objective"] is None and d.fields["software"] is None
    assert d.summary == "Leica LIF 8x1, float32"


def test_truncated_xml_still_yields_the_images_that_fit():
    full = _lif(_XML)
    cut = full[: 13 + 2 * (_XML.find("<Element Name=\"Img2\"") + 5)]
    d = decode_lif(BytesSource(cut))
    assert d is not None and d.fields["xmlTruncated"] is True
    assert d.fields["imageName"] == "Img1" and d.fields["sizeX"] == 512


def test_rejections():
    assert decode_lif(BytesSource(b"")) is None
    assert decode_lif(BytesSource(b"\x00" * 64)) is None
    assert decode_lif(BytesSource(struct.pack("<iiBi", 0x71, 10, 0x2A, 4) + b"x" * 64)) is None   # bad marker
    assert decode_lif(BytesSource(struct.pack("<iiBi", 0x70, 10, 0x2B, 4) + b"x" * 64)) is None   # bad test byte
    assert decode_lif(BytesSource(_lif("<Other/>" * 8))) is None                                  # not LIF XML
    assert decode_lif(BytesSource(_lif("<LMSDataContainerHeader"))) is None                       # unparseable
    assert decode_bytes("e.lif", _lif(_XML, n_chars=-5)) is None


def test_no_images_is_still_a_lif():
    d = decode_lif(BytesSource(_lif('<LMSDataContainerHeader Version="2"><Element Name="E"/></LMSDataContainerHeader>')))
    assert d.fields["imageCount"] == 0 and d.fields["imageName"] is None
    json.dumps(d.to_dict(), allow_nan=False)


def test_golden_pr2729_matches_readlif():
    src = load_sparse("lif-pr2729")
    d = decode_bytes("PR2729_frameOrderCombinedScanTypes.lif", b"")  # extension dispatch sanity
    assert d is None
    from scigantic_headers.decoders import decode_source
    d = decode_source("x.lif", src)
    f = d.fields
    # readlif: dims x=64 y=64 z=3 t=2 m=4, 2 channels, 8 bit
    assert (f["sizeX"], f["sizeY"], f["sizeZ"], f["sizeT"], f["sizeC"], f["tileCount"]) == (64, 64, 3, 2, 2, 4)
    assert f["dtype"] == "uint8" and f["channels"] == ["Green", "Red"]
    # readlif scale 0.252 px/um -> 3.968 um/px; z 1/0.160575 = 6.2276 um
    assert abs(f["voxelSizeUm"][0] - 1 / 0.252) < 1e-3
    assert abs(f["voxelSizeUm"][2] - 1 / 0.16057511583487422) < 1e-3
    assert f["objective"] == "C PLAN  10.0x0.22 POL" and f["magnification"] == 10.0
    assert f["numericalAperture"] == 0.22 and f["scanMode"] == "xyzt" and f["immersion"] == "DRY"
    assert f["imageName"] == "TileScan_001" and f["imageCount"] == 1
    with open(os.path.join(FIXTURES_DIR, "lif-pr2729.sparse.json")) as fh:
        assert d.to_dict() == json.load(fh)["expected"]
    json.dumps(d.to_dict(), allow_nan=False)
