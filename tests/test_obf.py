"""Imspector OBF / MSR decoder.

Goldens are sparse recordings of real .msr measurements from the msr-reader
project's test data (github.com/hoerlteam/msr-reader, MIT): the ranges a
decode touched, with the file size. Reference values were read with obffile
0.x (Christoph Gohlke) on the same files and agree with the Imspector
file-format documentation:

  bad44489_field_0_sted_0.msr  file v2; 2 stacks "ATTO 594_STED {0}",
      "STAR 635P_STED {0}", stack v7, rank 4, res (120, 110, 17, 1),
      len (3.6e-6, 3.3e-6, 2.55e-6, 1.0) m, dt 0x8 int16, zlib; footer labels
      ExpControl X/Y/Z/None, metres exponent 1 on the first three axes.
  mesc_yoyo2_img_001.msr       file v2; 2 stacks "gfp {0}" 250x250x1x1
      len 5e-5 m, "gfp {1}" 270x235x50x1.
  conf01.msr                   file v1; 1 unnamed stack 1024x1x1x1, int16,
      label "ExpControl T", no metres axis (a time trace).
  groupa_cell01.msr            file v1; 7 unnamed stacks (obffile keys its
      stack table by name and so collapses these to one; the linked list in
      the file has seven, which this decoder follows).
"""

import json
import os
import struct

import scigantic_headers.obf  # noqa: F401
from scigantic_headers import decode_bytes, has_decoder_for, read_for
from scigantic_headers.decoders import BytesSource, decode_source
from scigantic_headers.obf import decode_obf

from conftest import FIXTURES_DIR, load_sparse


def _stack(res, lengths, dt=0x4, comp=1, name=b"", descr=b"", data=b"", next_pos=0, version=5, footer=True):
    rank = len(res)
    res15 = list(res) + [1] * (15 - rank)
    len15 = list(lengths) + [1.0] * (15 - rank)
    off15 = [0.0] * 15
    hdr = struct.pack("<16sII15I15d15d5I3Q", b"OMAS_BF_STACK\n\xff\xff", version, rank, *res15, *len15, *off15,
                      dt, comp, 1, len(name), len(descr), 0, len(data), next_pos)
    out = hdr + name + descr + data
    if footer:
        # footer: size covers 4 + 60 + 60 + 4 + SI value (80) + 15 SI dims (1200) = 1408
        si_unit_m = struct.pack("<18id", 1, 1, *([0, 1] * 8), 1.0)
        si_unit_none = struct.pack("<18id", 0, 1, *([0, 1] * 8), 1.0)
        dims = b"".join(si_unit_m if i < min(rank, 3) else si_unit_none for i in range(15))
        ft = struct.pack("<I", 1408) + b"\x00" * 120 + struct.pack("<I", 0) + si_unit_none + dims
        labels = [b"ExpControl X", b"ExpControl Y", b"ExpControl Z", b"None"][:rank]
        for lab in labels:
            ft += struct.pack("<I", len(lab)) + lab
        out += ft
    return out


def _file(stacks, version=2, descr=b"<root/>"):
    first = 26 + len(descr) + (8 if version >= 2 else 0)
    head = struct.pack("<10sIQI", b"OMAS_BF\n\xff\xff", version, first, len(descr)) + descr
    if version >= 2:
        head += struct.pack("<Q", 0)
    # lay the stacks out back to back, patching next_stack_pos
    blobs = []
    pos = first
    positions = []
    for st in stacks:
        positions.append(pos)
        pos += len(st)
    body = b""
    for i, st in enumerate(stacks):
        nxt = positions[i + 1] if i + 1 < len(stacks) else 0
        # next_stack_pos is the last Q of the 368-byte header
        st = st[:360] + struct.pack("<Q", nxt) + st[368:]
        body += st
    return head + body


def test_registered():
    assert has_decoder_for("a.msr") and has_decoder_for("a.obf")
    assert read_for("a.msr").seek == 1024 * 1024


def test_synthetic_two_stacks_with_footer():
    s1 = _stack((120, 110, 17, 1), (3.6e-6, 3.3e-6, 2.55e-6, 1.0), dt=0x8, name=b"ATTO 594_STED {0}", data=b"\x00" * 40)
    s2 = _stack((64, 64), (6.4e-6, 6.4e-6), dt=0x40 | 0x40000000, comp=0, name=b"fft", data=b"\x00" * 8)
    d = decode_bytes("m.msr", _file([s1, s2]))
    f = d.fields
    assert d.format == "obf" and f["formatVersion"] == 2 and f["stackCount"] == 2
    assert f["name"] == "ATTO 594_STED {0}" and f["rank"] == 4 and f["sizes"] == [120, 110, 17, 1]
    assert (f["sizeX"], f["sizeY"], f["sizeZ"], f["sizeT"]) == (120, 110, 17, 1)
    assert f["dtype"] == "int16" and f["compression"] == "zlib"
    assert f["voxelSizeUm"] == [0.03, 0.03, 0.15]
    assert f["axisLabels"] == ["ExpControl X", "ExpControl Y", "ExpControl Z", "None"]
    assert f["stacks"][1] == {"name": "fft", "sizes": [64, 64], "dtype": "complex64"}
    assert f["description"] == "<root/>" and f["hasMetadata"] is False
    assert d.summary == "Imspector OBF 120x110x17, int16, 0.03x0.03x0.15 um, 2 stacks"


def test_version1_header_has_no_metadata_position_and_no_footer_means_assumed_metres():
    s = _stack((256, 256, 1), (2.56e-5, 2.56e-5, 1.0), dt=0x1, footer=False)
    d = decode_obf(BytesSource(_file([s], version=1, descr=b"")))
    f = d.fields
    assert f["formatVersion"] == 1 and f["description"] is None and f["hasMetadata"] is False
    assert f["axisLabels"] is None
    assert f["dtype"] == "uint8" and f["voxelSizeUm"] == [0.1, 0.1, None]


def test_dtype_table():
    for flag, name in [(0x1, "uint8"), (0x2, "int8"), (0x4, "uint16"), (0x8, "int16"), (0x10, "uint32"),
                       (0x20, "int32"), (0x40, "float32"), (0x80, "float64"), (0x1000, "uint64"),
                       (0x2000, "int64"), (0x10000, "bool"), (0x400, "rgb-uint8"), (0x80 | 0x40000000, "complex128")]:
        d = decode_obf(BytesSource(_file([_stack((4, 4), (1e-6, 1e-6), dt=flag, footer=False)])))
        assert d.fields["dtype"] == name, hex(flag)
    d = decode_obf(BytesSource(_file([_stack((4, 4), (1e-6, 1e-6), dt=0x12345, footer=False)])))
    assert d.fields["dtype"] is None


def test_stack_walk_is_bounded_and_cycle_safe():
    s = _stack((4, 4), (1e-6, 1e-6), footer=False)
    blob = _file([s])
    # point the single stack's next pointer at itself
    first = 26 + len(b"<root/>") + 8
    blob = blob[:first + 360] + struct.pack("<Q", first) + blob[first + 368:]
    d = decode_obf(BytesSource(blob))
    assert d.fields["stackCount"] == 1 and d.fields["stackCountTruncated"] is False


def test_rejections():
    assert decode_obf(BytesSource(b"")) is None
    assert decode_obf(BytesSource(b"OMAS_BF\n\xff\xff" + b"\x00" * 100)) is None     # first_stack_pos 0 -> no stacks
    assert decode_obf(BytesSource(b"OMAS_BF\n\x00\x00" + b"\x00" * 100)) is None     # wrong padding
    bad = _file([_stack((4, 4), (1e-6, 1e-6), footer=False)])
    assert decode_obf(BytesSource(bad[:40])) is None                                  # truncated at the stack
    rank0 = bad.replace(struct.pack("<I", 2), struct.pack("<I", 0), 1)
    assert decode_obf(BytesSource(rank0)) is None or decode_obf(BytesSource(rank0)).fields["rank"] > 0


def _golden(name):
    d = decode_source("x.msr", load_sparse(name))
    assert d is not None
    with open(os.path.join(FIXTURES_DIR, name + ".sparse.json")) as fh:
        assert d.to_dict() == json.load(fh)["expected"]
    json.dumps(d.to_dict(), allow_nan=False)
    return d.fields


def test_golden_bad44489_sted():
    f = _golden("obf-bad44489_sted")
    assert f["formatVersion"] == 2 and f["stackCount"] == 2 and f["hasMetadata"] is True
    assert [s["name"] for s in f["stacks"]] == ["ATTO 594_STED {0}", "STAR 635P_STED {0}"]
    assert f["stackVersion"] == 7 and f["sizes"] == [120, 110, 17, 1] and f["dtype"] == "int16"
    # the file holds float noise (3.5999998999e-06); the values are passed through as written
    assert [round(x, 9) for x in f["lengths"]] == [3.6e-06, 3.3e-06, 2.55e-06, 1.0]
    assert f["compression"] == "zlib"
    assert f["voxelSizeUm"] == [0.03, 0.03, 0.15]
    assert f["axisLabels"] == ["ExpControl X", "ExpControl Y", "ExpControl Z", "None"]


def test_golden_mesc_yoyo2():
    f = _golden("obf-mesc_yoyo2")
    assert f["stackCount"] == 2 and f["name"] == "gfp {0}" and f["sizes"] == [250, 250, 1, 1]
    assert f["stacks"][1] == {"name": "gfp {1}", "sizes": [270, 235, 50, 1], "dtype": "int16"}
    assert f["voxelSizeUm"] == [0.2, 0.2, None]


def test_golden_conf01_time_trace():
    f = _golden("obf-conf01")
    assert f["formatVersion"] == 1 and f["stackCount"] == 1 and f["name"] is None
    assert f["sizes"] == [1024, 1, 1, 1] and f["axisLabels"] == ["ExpControl T", "None", "None", "None"]
    assert f["voxelSizeUm"] == [None, None, None] and f["dtype"] == "int16"


def test_golden_groupa_seven_unnamed_stacks():
    f = _golden("obf-groupa_cell01")
    assert f["stackCount"] == 7 and f["stackCountTruncated"] is False
    assert f["sizes"] == [450, 490, 1, 1]
    assert f["stacks"][-1] == {"name": None, "sizes": [665, 1, 1, 1], "dtype": "float64"}  # the one obffile reports
