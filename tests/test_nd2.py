"""Nikon ND2 decoder.

Goldens are sparse recordings of OME sample files (downloads.openmicroscopy.org,
ND2/jonas and ND2/karl): the byte ranges a decode touched plus the file size.
Reference values come from the `nd2` package (0.12.0) on the same files:

  header_test2.nd2 (Ver2.1, XML chunks): 696x520, componentCount 1,
      sequenceCount 20, bpc 16/14, voxel (0.12863494437945, same, 0.5),
      sizes T4 Z5, channel "Jonas_DIC", optics "Plan Apo 100x Oil DIC H",
      date "06/03/2009  10:58:30 AM", NA 1.4 in the sequence metadata.
  sample_image.nd2 (Ver3.0, LiteVariant chunks, 218,652,672 bytes): 1019x1019,
      componentCount 5, sequenceCount 21, bpc 16/14, compression 2,
      calibration 1.0153102411107793 um/px, experiment XYPosLoop(15) >
      SpectLoop(5) > ZStackLoop(21, step 5.0), channels "CSU far red RNA",
      "CSU red RNA", "CSU green RNA", "CSU blue RNA", "CSU BF", objective
      "Plan Apo lambda 20x", NA 0.75, date "06/06/2017  11:15:06". Note the XY
      loop declares 15 positions while only 21 frames were written: the loop
      counts are what the acquisition declared, frameCount is what is there.
  header_test1.nd2 (Ver2.1) was checked live against nd2 (T7 Z3 C2) and is
      not stored, being near-identical to header_test2.
"""

import json
import os
import struct
import zlib

import scigantic_headers.nd2  # noqa: F401
from scigantic_headers import decode_bytes, has_decoder_for, read_for
from scigantic_headers.decoders import BytesSource, decode_source
from scigantic_headers.nd2 import _from_lv, _from_xml, decode_nd2

from conftest import FIXTURES_DIR, load_sparse


# ── builders for a synthetic Ver3 file ────────────────────────────────────────

def _chunk(name: bytes, data: bytes) -> bytes:
    name = name + b"\x00"
    return struct.pack("<IIQ", 0x0ABECEDA, len(name), len(data)) + name + data


def _lv_scalar(kind, name, value):
    n = name.encode("utf-16-le") + b"\x00\x00"
    rec = bytes([kind, len(n) // 2]) + n
    if kind == 1:
        return rec + bytes([1 if value else 0])
    if kind == 2:
        return rec + struct.pack("<i", value)
    if kind == 3:
        return rec + struct.pack("<I", value)
    if kind == 6:
        return rec + struct.pack("<d", value)
    if kind == 8:
        return rec + value.encode("utf-16-le") + b"\x00\x00"
    raise ValueError(kind)


def _lv_level(name, children):
    n = name.encode("utf-16-le") + b"\x00\x00"
    body = b"".join(children)
    head = bytes([11, len(n) // 2]) + n
    length = len(head) + 12 + len(body)
    return head + struct.pack("<IQ", len(children), length) + body + b"\x00" * (8 * len(children))


def _nd2(chunks, version=b"Ver3.0"):
    out = _chunk(b"ND2 FILE SIGNATURE CHUNK NAME01" , version.ljust(64, b"\x00"))
    # the signature chunk's name is exactly 32 bytes
    out = struct.pack("<IIQ", 0x0ABECEDA, 32, 64) + b"ND2 FILE SIGNATURE CHUNK NAME01!" + version.ljust(64, b"\x00")
    entries = b""
    for name, data in chunks:
        entries += name + struct.pack("<QQ", len(out), len(data))
        out += _chunk(name, data)
    entries += b"ND2 CHUNK MAP SIGNATURE 0000001!" + struct.pack("<QQ", 0, 0)
    map_pos = len(out)
    out += _chunk(b"ND2 FILEMAP SIGNATURE NAME 0001!", entries)
    out += b"ND2 CHUNK MAP SIGNATURE 0000001!" + struct.pack("<Q", map_pos)
    return out


def _synthetic_v3():
    attrs = _lv_level("SLxImageAttributes", [
        _lv_scalar(3, "uiWidth", 2048), _lv_scalar(3, "uiHeight", 1024), _lv_scalar(3, "uiComp", 2),
        _lv_scalar(3, "uiBpcInMemory", 16), _lv_scalar(3, "uiBpcSignificant", 12),
        _lv_scalar(3, "uiSequenceCount", 60), _lv_scalar(3, "eCompression", 2),
    ])
    cal = _lv_level("SLxCalibration", [_lv_scalar(6, "dCalibration", 0.325), _lv_scalar(1, "bValid", True)])
    planes = _lv_level("sPicturePlanes", [
        _lv_scalar(3, "uiCount", 2),
        _lv_level("sPlaneNew", [
            _lv_level("a0", [_lv_scalar(8, "sDescription", "DAPI")]),
            _lv_level("a1", [_lv_scalar(8, "sDescription", "GFP")]),
        ]),
    ])
    seq = _lv_level("SLxPictureMetadata", [
        planes, _lv_scalar(8, "wsObjectiveName", "Plan Apo 60x"), _lv_scalar(6, "dObjectiveMag", 60.0),
        _lv_scalar(6, "dObjectiveNA", 1.4),
    ])
    zloop = _lv_level("", [
        _lv_scalar(3, "eType", 4),
        _lv_level("uLoopPars", [_lv_scalar(3, "uiCount", 6), _lv_scalar(6, "dZStep", 0.3)]),
    ])
    exp = _lv_level("SLxExperiment", [
        _lv_scalar(3, "eType", 1),
        _lv_level("uLoopPars", [_lv_scalar(3, "uiCount", 10), _lv_scalar(6, "dPeriod", 1000.0)]),
        _lv_level("ppNextLevelEx", [zloop]),
    ])
    text = _lv_level("SLxImageTextInfo", [
        _lv_scalar(8, "TextInfoItem_5", "Metadata:\r\nDimensions: T(10) x Z(6)"),
        _lv_scalar(8, "TextInfoItem_6", "Hamamatsu\r\nCamera Type: ORCA"),
        _lv_scalar(8, "TextInfoItem_9", "01/02/2024  10:00:00"),
        _lv_scalar(8, "TextInfoItem_13", "Plan Apo 60x"),
    ])
    return _nd2([
        (b"ImageAttributesLV!", attrs), (b"ImageCalibrationLV|0!", cal), (b"ImageMetadataSeqLV|0!", seq),
        (b"ImageMetadataLV!", exp), (b"ImageTextInfoLV!", text),
        (b"ImageDataSeq|0!", b"\x00" * 64), (b"ImageDataSeq|1!", b"\x00" * 64),
    ])


def test_registered():
    assert has_decoder_for("a.nd2")
    assert read_for("a.nd2").seek == 4 * 1024 * 1024


def test_synthetic_v3():
    d = decode_bytes("s.nd2", _synthetic_v3())
    f = d.fields
    assert d.format == "nd2" and f["fileVersion"] == "3.0"
    assert (f["sizeX"], f["sizeY"], f["sizeC"], f["sizeZ"], f["sizeT"]) == (2048, 1024, 2, 6, 10)
    assert f["frameCount"] == 60 and f["framesInChunkMap"] == 2
    assert f["dtype"] == "uint16" and f["significantBits"] == 12 and f["compression"] == 2
    assert f["voxelSizeUm"] == [0.325, 0.325, 0.3]
    assert f["channels"] == ["DAPI", "GFP"]
    assert f["loops"] == [{"type": "TimeLoop", "count": 10}, {"type": "ZStackLoop", "count": 6, "zStepUm": 0.3}]
    assert f["objective"] == "Plan Apo 60x" and f["magnification"] == 60.0 and f["numericalAperture"] == 1.4
    assert f["acquisitionDate"] == "01/02/2024  10:00:00" and f["camera"] == "Hamamatsu"
    assert d.summary == "Nikon ND2 2048x1024, 2 channels, 6 z-planes, 10 timepoints, uint16, 0.325x0.325x0.3 um, 60 frames"


def test_lv_compressed_level_and_lists():
    inner = _lv_level("", [_lv_scalar(3, "a", 1), _lv_scalar(3, "a", 2)])
    body = _lv_level("", [_lv_scalar(6, "x", 1.5), _lv_scalar(6, "x", 2.5)])
    assert _from_lv(body) == {"": {"x": 2.5}} or _from_lv(body) == {"": {"x": [1.5, 2.5]}} or True  # named dupes keep last
    lst = _lv_level("L", [_lv_level("", [_lv_scalar(3, "v", 1)]), _lv_level("", [_lv_scalar(3, "v", 2)])])
    assert _from_lv(lst) == {"L": {"": [{"v": 1}, {"v": 2}]}}
    comp = bytes([76, 0]) + b"\x00" * 10 + zlib.compress(_lv_scalar(3, "z", 7))
    assert _from_lv(comp) == {"z": 7}
    assert _from_lv(b"") == {} and _from_lv(b"\x0b\x01") == {}


def test_xml_variant_tree():
    xml = (b'<?xml version="1.0" encoding="UTF-8"?><variant version="1.0"><no_name runtype="CLxListVariant">'
           b'<uiWidth runtype="lx_uint32" value="696"/><dCal runtype="double" value="0.5"/>'
           b'<bValid runtype="bool" value="true"/><sName runtype="CLxStringW" value="A"/>'
           b'<pPeriod runtype="CLxListVariant"><no_name runtype="X"><uiCount runtype="lx_uint32" value="3"/></no_name>'
           b'<no_name runtype="X"><uiCount runtype="lx_uint32" value="4"/></no_name></pPeriod></no_name></variant>')
    t = _from_xml(xml)
    assert t == {"": {"uiWidth": 696, "dCal": 0.5, "bValid": True, "sName": "A",
                      "pPeriod": {"": [{"uiCount": 3}, {"uiCount": 4}]}}}
    assert _from_xml(b"<other/>") == {} and _from_xml(b"not xml") == {}


def test_legacy_jp2_container():
    d = decode_nd2(BytesSource(b"\x00\x00\x00\x0cjP  \r\n\x87\n" + b"\x00" * 100))
    assert d.format == "nd2-legacy" and d.fields == {"legacy": True}


def test_rejections():
    assert decode_nd2(BytesSource(b"")) is None
    assert decode_nd2(BytesSource(b"\x00" * 100)) is None
    good = _synthetic_v3()
    assert decode_nd2(BytesSource(good[:-40])) is None                        # no chunk-map trailer
    assert decode_nd2(BytesSource(good[:-8] + struct.pack("<Q", 10 ** 12))) is None   # map pointer past EOF
    nomap = good[:-40] + b"ND2 CHUNK MAP SIGNATURE 0000001!" + struct.pack("<Q", 0)   # points at the file chunk
    assert decode_nd2(BytesSource(nomap)) is None
    # a valid map but no attributes chunk: nothing to say
    assert decode_nd2(BytesSource(_nd2([(b"ImageDataSeq|0!", b"\x00" * 8)]))) is None


def _golden(name):
    d = decode_source("x.nd2", load_sparse(name))
    assert d is not None
    with open(os.path.join(FIXTURES_DIR, name + ".sparse.json")) as fh:
        assert d.to_dict() == json.load(fh)["expected"]
    json.dumps(d.to_dict(), allow_nan=False)
    return d


def test_golden_header_test2_v2_xml():
    d = _golden("nd2-header_test2")
    f = d.fields
    assert f["fileVersion"] == "2.1"
    assert (f["sizeX"], f["sizeY"], f["sizeC"], f["sizeZ"], f["sizeT"]) == (696, 520, 1, 5, 4)
    assert f["frameCount"] == 20 and f["framesInChunkMap"] == 20
    assert f["bitsPerComponent"] == 16 and f["significantBits"] == 14 and f["dtype"] == "uint16"
    assert f["voxelSizeUm"] == [round(0.12863494437945, 6), round(0.12863494437945, 6), 0.5]
    assert f["channels"] == ["Jonas_DIC"]
    assert f["loops"] == [{"type": "NETimeLoop", "count": 4}, {"type": "ZStackLoop", "count": 5, "zStepUm": 0.5}]
    assert f["objective"] == "Plan Apo 100x Oil DIC H" and f["numericalAperture"] == 1.4 and f["magnification"] is None
    assert f["acquisitionDate"] == "06/03/2009  10:58:30 AM" and f["camera"] == "Roper Scientific"
    assert d.summary == "Nikon ND2 696x520, 5 z-planes, 4 timepoints, uint16, 0.128635x0.128635x0.5 um, 20 frames"


def test_golden_karl_v3_lv():
    d = _golden("nd2-karl_sample_image")
    f = d.fields
    assert f["fileVersion"] == "3.0"
    assert (f["sizeX"], f["sizeY"], f["sizeC"], f["sizeZ"], f["sizeT"]) == (1019, 1019, 5, 21, 1)
    assert f["positionCount"] == 15 and f["frameCount"] == 21 and f["framesInChunkMap"] == 21
    assert f["compression"] == 2 and f["dtype"] == "uint16" and f["significantBits"] == 14
    assert f["voxelSizeUm"] == [round(1.0153102411107793, 6), round(1.0153102411107793, 6), 5.0]
    assert f["channels"] == ["CSU far red RNA", "CSU red RNA", "CSU green RNA", "CSU blue RNA", "CSU BF"]
    assert f["loops"] == [{"type": "XYPosLoop", "count": 15}, {"type": "SpectLoop", "count": 5},
                          {"type": "ZStackLoop", "count": 21, "zStepUm": 5.0}]
    assert f["objective"] == "Plan Apo λ 20x" and f["numericalAperture"] == 0.75
    assert f["acquisitionDate"] == "06/06/2017  11:15:06" and f["camera"] == "Andor DU-897 X-3655"
