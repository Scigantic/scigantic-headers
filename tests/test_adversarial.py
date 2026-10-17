"""Hostile files for the seeking decoders: pointers that loop back on
themselves, offsets past the end of the file, counts in the billions, XML
entity expansion. Each decode must terminate well inside the read cap and
return None or a JSON-safe header. Found nothing when first written; kept so a
future shortcut in a pointer walk fails here rather than in a tree walk."""

import json
import struct

import pytest

from scigantic_headers.decoders import SEEK_MAX_READS, Source, decode_source


class _Counting(Source):
    def __init__(self, data, size=None):
        super().__init__(1 << 30)
        self.d = data
        self.size = size if size is not None else len(data)
        self.n = 0

    def _fetch(self, o, n):
        self.n += 1
        return self.d[o:o + n]


def _decode(key, data, size=None, max_reads=64):
    src = _Counting(data, size)
    src.key = key
    hdr = decode_source(key, src)
    if hdr is not None:
        json.loads(json.dumps(hdr.to_dict(), allow_nan=False))
    assert src.n <= max_reads, "%s took %d reads" % (key, src.n)
    return hdr


def _ifd(width, nxt):
    return (struct.pack("<H", 3) + struct.pack("<HHII", 256, 4, 1, width) + struct.pack("<HHII", 257, 4, 1, 4)
            + struct.pack("<HHIHH", 258, 3, 1, 8, 0) + struct.pack("<I", nxt))


def test_tiff_ifd_self_loop_is_one_page():
    h = _decode("x.tif", b"II*\x00" + struct.pack("<I", 8) + _ifd(16, 8))
    assert h is not None and h.fields["pageCount"] == 1


def test_tiff_ten_thousand_ifds_are_capped():
    chain = b"II*\x00" + struct.pack("<I", 8)
    off, step = 8, len(_ifd(1, 0))
    for i in range(10000):
        chain += _ifd(16, off + step if i < 9999 else 0)
        off += step
    h = _decode("x.tif", chain, max_reads=SEEK_MAX_READS)
    assert h.fields["pageCountTruncated"] is True and h.fields["pageCount"] <= 64


@pytest.mark.parametrize("data", [
    b"II*\x00" + struct.pack("<I", 0xFFFFFFF0) + b"\0" * 100,                       # IFD0 past EOF
    b"II*\x00" + struct.pack("<I", 8) + struct.pack("<H", 65535) + b"\x01" * 4000,  # 65535 entries
    b"II+\x00\x08\x00\x00\x00" + struct.pack("<Q", 16) + struct.pack("<Q", 1 << 60) + b"\0" * 64,
    b"II*\x00" + struct.pack("<I", 8) + struct.pack("<H", 2) + struct.pack("<HHII", 256, 4, 1, 5)
    + struct.pack("<HHII", 270, 2, 0xFFFFFFFF, 8) + struct.pack("<I", 0),           # 4 GB description
])
def test_tiff_hostile_pointers_terminate(data):
    _decode("x.tif", data)


def _sb_v2(root):
    return (b"\x89HDF\r\n\x1a\n" + bytes([2, 8, 8, 0]) + struct.pack("<QQQQ", 0, 0xFFFFFFFFFFFFFFFF, 4096, root)
            + b"\0\0\0\0")


def test_hdf5_continuation_pointing_at_itself():
    oh = b"OHDR" + bytes([2, 0]) + struct.pack("<B", 40)
    oh += bytes([0x10]) + struct.pack("<H", 16) + b"\0" + struct.pack("<QQ", 48, 64) + b"\0" * 40
    _decode("x.h5", _sb_v2(48) + oh + b"\0" * 4000)


def test_hdf5_root_header_is_the_superblock():
    assert _decode("x.h5", _sb_v2(0) + b"\0" * 4000) is None


def test_hdf5_btree_child_is_itself():
    sb0 = (b"\x89HDF\r\n\x1a\n" + bytes([0, 0, 0, 0, 0, 8, 8, 0]) + struct.pack("<HHI", 4, 16, 0)
           + struct.pack("<QQQQ", 0, 0xFFFFFFFFFFFFFFFF, 8192, 0xFFFFFFFFFFFFFFFF)
           + struct.pack("<QQII", 0, 96, 1, 0) + struct.pack("<QQ", 200, 400))
    ohv1 = struct.pack("<BBHII", 1, 0, 1, 1, 24) + b"\0" * 4 + struct.pack("<HHB3x", 0x11, 16, 0) + struct.pack("<QQ", 200, 400)
    tree = (b"TREE" + bytes([0, 1]) + struct.pack("<H", 2) + struct.pack("<QQ", 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
            + struct.pack("<QQQQQ", 0, 200, 0, 200, 0))
    heap = b"HEAP" + bytes([0, 0, 0, 0]) + struct.pack("<QQQ", 64, 0, 480) + b"a\0" * 32
    data = sb0 + b"\0" * (96 - len(sb0)) + ohv1
    data += b"\0" * (200 - len(data)) + tree
    data += b"\0" * (400 - len(data)) + heap + b"\0" * 2000
    _decode("x.h5", data, max_reads=16)  # the cycle is caught, not walked to the read cap


def test_czi_hostile_metadata_positions():
    czi = bytearray(b"ZISRAWFILE" + b"\0" * 6 + struct.pack("<qq", 512, 512) + b"\0" * 512)
    struct.pack_into("<q", czi, 32 + 60, 0)
    struct.pack_into("<q", czi, 32 + 52, 1 << 62)
    _decode("x.czi", bytes(czi))                                   # metadata -> own header
    meta = bytearray(b"ZISRAWMETADATA\0\0" + struct.pack("<qq", 1 << 40, 1 << 40) + struct.pack("<ii", -5, 0)
                     + b"\0" * 256 + b"<ImageDocument>" * 100)
    struct.pack_into("<q", czi, 32 + 60, len(czi))
    _decode("x.czi", bytes(czi) + bytes(meta))                     # negative XML size
    struct.pack_into("<i", meta, 32, 1 << 30)
    _decode("x.czi", bytes(czi) + bytes(meta), size=1 << 31)      # 1 GB XML size


def test_nd2_hostile_chunk_maps():
    sig = b"ND2 CHUNK MAP SIGNATURE 0000001!"

    def chunk(name, data):
        return struct.pack("<IIQ", 0x0ABECEDA, len(name), len(data)) + name + data

    body = chunk(b"ImageAttributesLV!\0", b"\0" * 64)
    pos = len(body)
    cm = chunk(sig + b"\0", b"ImageAttributesLV!\0" + struct.pack("<QQ", 1 << 50, 1 << 50)
               + b"ImageMetadataLV!\0" + struct.pack("<QQ", pos, 1 << 40) + sig + struct.pack("<QQ", pos, 0))
    _decode("x.nd2", body + cm + sig + struct.pack("<Q", pos))
    _decode("x.nd2", body + cm + sig + struct.pack("<Q", len(body + cm)))
    _decode("x.nd2", body + sig + struct.pack("<Q", 1 << 60))


def test_lif_obf_klb_hostile_sizes():
    _decode("x.lif", struct.pack("<ii", 0x70, 1 << 30) + b"\x2a" + struct.pack("<i", 1 << 30) + b"<\0a\0" * 1000)
    obf = b"OMAS_BF\n\xff\xff" + struct.pack("<IQI", 1, 32, 0) + struct.pack("<Q", 0)
    st = (b"OMAS_BF_STACK\n\xff\xff" + struct.pack("<II", 1, 3) + struct.pack("<15I", 64, 64, 8, *[1] * 12)
          + struct.pack("<15d", *[1e-6] * 15) + struct.pack("<15d", *[0.0] * 15) + struct.pack("<5I", 4, 0, 0, 0, 0)
          + struct.pack("<QQQ", 0, 1 << 60, 32))
    _decode("x.msr", obf + b"\0" * (32 - len(obf)) + st + b"\0" * 512)   # next stack -> itself
    _decode("x.klb", bytes([2]) + struct.pack("<5I", *[0xFFFFFFFF] * 5) + struct.pack("<5f", 1, 1, 1, 1, 1)
            + bytes([1, 1]) + b"m" * 256 + struct.pack("<5I", 1, 1, 1, 1, 1) + b"\0" * 64)


def test_json_sidecars_hostile_sizes():
    _decode(".zattrs", b'{"multiscales":[{"axes":' + b"[" * 100000 + b"]" * 100000 + b"}]}")
    _decode(".zarray", b'{"zarr_format":2,"shape":[' + b"1," * 1000000 + b'1],"chunks":[1],"dtype":"<u2"}')
    _decode("attributes.json", b'{"n5":"2.0.0","dimensions":[' + b"9" * 2000000 + b"]}")


def test_xml_entities_are_not_expanded():
    laughs = (b'<?xml version="1.0"?><!DOCTYPE lolz [<!ENTITY lol "lol"><!ENTITY lol2 "&lol;&lol;&lol;&lol;">]>'
              b'<SpimData version="0.2"><a>&lol2;</a></SpimData>')
    assert _decode("x.xml", laughs) is None
    xxe = (b'<?xml version="1.0"?><!DOCTYPE OME [<!ENTITY xxe SYSTEM "file:///etc/passwd">]>'
           b'<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"><Image Name="&xxe;">'
           b'<Pixels SizeX="1" SizeY="1" SizeZ="1" SizeC="1" SizeT="1" Type="uint8" DimensionOrder="XYZCT"/></Image></OME>')
    h = _decode("x.ome.xml", xxe)
    assert h is None or "root" not in json.dumps(h.to_dict())
