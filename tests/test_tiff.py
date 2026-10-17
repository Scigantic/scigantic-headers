"""TIFF family decoding: classic and BigTIFF headers, OME-TIFF, ImageJ
hyperstacks, Micro-Manager, plain stacks.

Real-file goldens (tests/fixtures/tiff-*.sparse.json, listed in
tiff-cases.json) were recorded from OME sample images (tubhiswt_C0.ome.tif,
4D-series.ome.tif / .tf8, multi-channel-time-series.ome.btf), two libtiff
pictures (smallliz.tif, quad-tile.tif) and the first 256 KiB of an ImageJ
hyperstack on Zenodo (record 10471496). Every expected dimension, dtype,
byte order, BigTIFF flag, page count, compression, photometric, resolution and
OME Pixels/Instrument value was checked against tifffile 2025 (TiffFile.pages,
.is_bigtiff, .ome_metadata, .imagej_metadata) before being frozen here.

The ImageJ and Micro-Manager synthetic files are hand-built below. The ImageJ
case reproduces a file tifffile wrote with imagej=True, resolution 1/0.325 and
spacing 1.5 um, which tifffile read back as 24 pages, channels 2, slices 4,
frames 3; the Micro-Manager case follows the published multipage TIFF layout
(summary JSON at byte 40, index map) and is checked by tifffile's
is_micromanager / micromanager_metadata in the development session (no real
Micro-Manager-native file was available; the Zenodo T-cell files turned out to
be ImageJ re-saves of a Micro-Manager acquisition).
"""

import json
import os
import struct

import pytest

from scigantic_headers import decode_bytes, has_decoder_for
from scigantic_headers.decoders import BytesSource, decode_source
from scigantic_headers.tiff import IFD_CHAIN_CAP, decode_tiff  # noqa: F401, registers the decoder
from conftest import FIXTURES_DIR, load_sparse

with open(os.path.join(FIXTURES_DIR, "tiff-cases.json")) as _fh:
    _CASES = json.load(_fh)["cases"]


# ── Builders ─────────────────────────────────────────────────────────────────

def _entry(bo, big, tag, typ, values, blob, size_of):
    """One IFD entry; out-of-line values are appended to `blob` (a list of
    bytes that the caller lays out after the IFD) and patched later."""
    if typ == 2:
        raw = values if isinstance(values, bytes) else values.encode() + b"\x00"
        count = len(raw)
    elif typ == 5:
        raw = b"".join(struct.pack(bo + "II", n, d) for n, d in values)
        count = len(values)
    else:
        fmt = {1: "B", 3: "H", 4: "I", 16: "Q", 12: "d"}[typ]
        raw = b"".join(struct.pack(bo + fmt, v) for v in values)
        count = len(values)
    inline_cap = 8 if big else 4
    if len(raw) <= inline_cap:
        return (tag, typ, count, raw.ljust(inline_cap, b"\x00"), None)
    return (tag, typ, count, None, raw)


def make_tiff(*, bo="<", big=False, width=64, height=32, bits=16, spp=1, sample_format=None,
              compression=1, photometric=1, description=None, software=None, xres=None, yres=None,
              resunit=None, tiled=False, pages=1, extra_entries=(), ifd_after_data=True,
              mm_header=None, chain_loop=False):
    """A valid TIFF with `pages` IFDs. Data is a few bytes per page. Out-of-line
    tag values follow each IFD. `mm_header` = (summary_json_bytes, index_entries)
    writes Micro-Manager's private header words and summary block."""
    off_fmt = "Q" if big else "I"
    cnt_fmt = "Q" if big else "H"
    entry_size = 20 if big else 12
    header_len = 16 if big else 8
    out = bytearray()
    out += (b"II" if bo == "<" else b"MM")
    if big:
        out += struct.pack(bo + "HHH", 43, 8, 0) + struct.pack(bo + "Q", 0)
    else:
        out += struct.pack(bo + "H", 42) + struct.pack(bo + "I", 0)
    if mm_header is not None:
        summary, index_entries = mm_header
        index_map = struct.pack(bo + "II", 3453623, len(index_entries)) + b"".join(
            struct.pack(bo + "5I", *e) for e in index_entries)
        # words: index map header, index map offset, display header, display offset,
        # comments header, comments offset, summary header, summary length
        summary_off = 40
        index_off = summary_off + len(summary)
        out += struct.pack(bo + "8I", 54773648, index_off, 483765892, 0, 99384722, 0, 2355492, len(summary))
        out += summary
        out += index_map
    data_bytes = max(1, width * height * spp * bits // 8)
    ifd_offsets = []
    for p in range(pages):
        entries = [
            _entry(bo, big, 256, 4, [width], None, None),
            _entry(bo, big, 257, 4, [height], None, None),
            _entry(bo, big, 258, 3, [bits] * spp, None, None),
        ]
        if compression is not None:
            entries.append(_entry(bo, big, 259, 3, [compression], None, None))
        entries.append(_entry(bo, big, 262, 3, [photometric], None, None))
        if description is not None and p == 0:
            entries.append(_entry(bo, big, 270, 2, description, None, None))
        entries.append(_entry(bo, big, 277, 3, [spp], None, None))
        if xres is not None:
            entries.append(_entry(bo, big, 282, 5, [xres], None, None))
        if yres is not None:
            entries.append(_entry(bo, big, 283, 5, [yres], None, None))
        if resunit is not None:
            entries.append(_entry(bo, big, 296, 3, [resunit], None, None))
        if software is not None and p == 0:
            entries.append(_entry(bo, big, 305, 2, software, None, None))
        if tiled:
            entries.append(_entry(bo, big, 322, 4, [16], None, None))
        if sample_format is not None:
            entries.append(_entry(bo, big, 339, 3, [sample_format], None, None))
        if p == 0:
            entries.extend(extra_entries)
        entries.sort(key=lambda e: e[0])
        # strip offsets: tag 273 inline, pointing at the data we emit first
        if ifd_after_data:
            data_off = len(out)
            out += b"\x00" * data_bytes
        entries.append(_entry(bo, big, 273, 4 if not big else 16, [len(out) if not ifd_after_data else data_off], None, None))
        entries.sort(key=lambda e: e[0])
        ifd_off = len(out)
        ifd_offsets.append(ifd_off)
        ifd_len = struct.calcsize(cnt_fmt) + len(entries) * entry_size + struct.calcsize(off_fmt)
        values_off = ifd_off + ifd_len
        ifd = bytearray(struct.pack(bo + cnt_fmt, len(entries)))
        tail = bytearray()
        for tag, typ, count, inline, raw in entries:
            if inline is not None:
                ifd += struct.pack(bo + "HH" + off_fmt, tag, typ, count) + inline
            else:
                ifd += struct.pack(bo + "HH" + off_fmt + off_fmt, tag, typ, count, values_off + len(tail))
                tail += raw + (b"\x00" if len(raw) % 2 else b"")
        ifd += struct.pack(bo + off_fmt, 0)  # next pointer, patched below
        out += ifd + tail
        if not ifd_after_data:
            out += b"\x00" * data_bytes
    # Link the chain.
    for i, off in enumerate(ifd_offsets):
        cnt = struct.unpack_from(bo + cnt_fmt, out, off)[0]
        next_pos = off + struct.calcsize(cnt_fmt) + cnt * entry_size
        if i + 1 < len(ifd_offsets):
            struct.pack_into(bo + off_fmt, out, next_pos, ifd_offsets[i + 1])
        elif chain_loop:
            struct.pack_into(bo + off_fmt, out, next_pos, ifd_offsets[0])
    struct.pack_into(bo + off_fmt, out, 8 if big else 4, ifd_offsets[0])
    return bytes(out)


# ── Golden fixtures from real files ───────────────────────────────────────────

@pytest.mark.parametrize("case", _CASES, ids=[c["name"] for c in _CASES])
def test_real_files_match_golden(case):
    d = decode_source(case["key"], load_sparse(case["fixture"]))
    assert d is not None, case["name"]
    assert d.to_dict() == case["expected"]


def test_goldens_cover_the_flavors():
    formats = {c["expected"]["format"] for c in _CASES}
    assert {"ome-tiff", "tiff", "imagej-tiff"} <= formats
    assert any(c["expected"]["fields"]["bigTiff"] for c in _CASES)
    assert any(c["expected"]["fields"]["byteOrder"] == "little" for c in _CASES)
    assert any(c["expected"]["fields"]["byteOrder"] == "big" for c in _CASES)


def test_golden_ome_tiff_details():
    by = {c["name"]: c["expected"] for c in _CASES}
    tub = by["tubhiswt_C0-ome-tif"]["fields"]
    assert (tub["sizeX"], tub["sizeY"], tub["sizeC"]) == (512, 512, 2)
    assert tub["microscope"] == "Nikon Eclipse TE300"
    assert tub["magnification"] == 100.0 and tub["numericalAperture"] == 1.4
    assert tub["acquisitionDate"] == "2013-01-15T17:02:40"
    four = by["4D-series-ome-tf8"]["fields"]
    assert four["bigTiff"] is True and four["byteOrder"] == "big"
    assert (four["sizeZ"], four["sizeT"], four["dtype"]) == (5, 7, "int8")
    assert four["pageCount"] == 35 and four["pageCountTruncated"] is False and four["planeCount"] == 35
    liz = by["smallliz-tif"]["fields"]
    assert liz["voxelSizeUm"] == [254.0, 254.0, None]   # 100 dpi
    assert liz["compression"] == "jpeg-old" and liz["samplesPerPixel"] == 3
    quad = by["quad-tile-tif"]["fields"]
    assert quad["tiled"] is True and quad["compression"] == "lzw" and quad["voxelSizeUm"] == [None, None, None]
    ij = by["imagej-tcell-hilo"]
    assert ij["format"] == "imagej-tiff"
    assert (ij["fields"]["sizeZ"], ij["fields"]["sizeT"], ij["fields"]["planeCount"]) == (74, 22, 1628)
    assert ij["fields"]["pageCountTruncated"] is True   # only the leading 256 KiB were available


# ── Synthetic: header variants ───────────────────────────────────────────────

@pytest.mark.parametrize("bo", ["<", ">"])
@pytest.mark.parametrize("big", [False, True])
def test_classic_and_bigtiff_both_orders(bo, big):
    d = decode_bytes("a.tif", make_tiff(bo=bo, big=big, width=100, height=50, bits=16, pages=3))
    assert d.format == "tiff"
    f = d.fields
    assert (f["width"], f["height"], f["dtype"], f["bigTiff"]) == (100, 50, "uint16", big)
    assert f["byteOrder"] == ("little" if bo == "<" else "big")
    assert f["pageCount"] == 3 and f["pageCountTruncated"] is False and f["sizeZ"] == 3
    assert d.summary.startswith("BigTIFF 100x50, 3 pages" if big else "TIFF 100x50, 3 pages")


def test_ifd_before_data_also_works():
    d = decode_bytes("a.tif", make_tiff(ifd_after_data=False, pages=2))
    assert d.fields["pageCount"] == 2


@pytest.mark.parametrize("bits,fmt,dtype", [(8, 1, "uint8"), (16, 1, "uint16"), (16, 2, "int16"),
                                            (32, 3, "float32"), (64, 3, "float64"), (32, 1, "uint32"), (1, 1, "bit")])
def test_dtype_from_bits_and_sample_format(bits, fmt, dtype):
    d = decode_bytes("a.tif", make_tiff(bits=bits, sample_format=fmt))
    assert d.fields["dtype"] == dtype


def test_resolution_units():
    inch = decode_bytes("a.tif", make_tiff(xres=(100, 1), yres=(100, 1), resunit=2)).fields["voxelSizeUm"]
    assert inch == [254.0, 254.0, None]
    cm = decode_bytes("a.tif", make_tiff(xres=(1000, 1), yres=(1000, 1), resunit=3)).fields["voxelSizeUm"]
    assert cm == [10.0, 10.0, None]
    none = decode_bytes("a.tif", make_tiff(xres=(1000, 1), yres=(1000, 1), resunit=1)).fields["voxelSizeUm"]
    assert none == [None, None, None]          # unit "none" is not micrometres
    zero = decode_bytes("a.tif", make_tiff(xres=(0, 1000), yres=(0, 1000), resunit=3)).fields["voxelSizeUm"]
    assert zero == [None, None, None]          # Bio-Formats writes 0/1000 when unknown


def test_compression_and_photometric_names():
    f = decode_bytes("a.tif", make_tiff(compression=5, photometric=2, spp=3, bits=8)).fields
    assert f["compression"] == "lzw" and f["photometric"] == "rgb" and f["samplesPerPixel"] == 3
    f = decode_bytes("a.tif", make_tiff(compression=50000)).fields
    assert f["compression"] == "zstd"
    f = decode_bytes("a.tif", make_tiff(compression=None)).fields
    assert f["compression"] is None            # tag absent (ImageJ omits it)


def test_software_and_text_description():
    f = decode_bytes("a.tif", make_tiff(software="HP IL v1.1", description="a plain caption")).fields
    assert f["software"] == "HP IL v1.1" and f["imageDescriptionKind"] == "text"


# ── Synthetic: OME-TIFF ──────────────────────────────────────────────────────

_OME = ('<?xml version="1.0" encoding="UTF-8"?><OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" '
        'Creator="test"><Instrument ID="I:0"><Microscope Manufacturer="Zeiss" Model="Lightsheet 7"/>'
        '<Objective ID="O:0" Model="20x/1.0" NominalMagnification="20" LensNA="1.0"/></Instrument>'
        '<Image ID="Image:0" Name="stack"><AcquisitionDate>2026-01-02T03:04:05</AcquisitionDate>'
        '<InstrumentRef ID="I:0"/><ObjectiveSettings ID="O:0"/>'
        '<Pixels ID="P:0" DimensionOrder="XYZCT" Type="uint16" SizeX="64" SizeY="32" SizeZ="10" SizeC="2" SizeT="1" '
        'PhysicalSizeX="0.5" PhysicalSizeY="0.5" PhysicalSizeZ="2" PhysicalSizeZUnit="µm">'
        '<Channel ID="C:0" Name="GFP" EmissionWavelength="510"/><Channel ID="C:1" Name="RFP"/>'
        '<TiffData/></Pixels></Image></OME>')


def test_ome_tiff_fields_and_summary():
    d = decode_bytes("s.ome.tif", make_tiff(description=_OME, pages=2))
    assert d.format == "ome-tiff"
    f = d.fields
    assert (f["sizeX"], f["sizeY"], f["sizeZ"], f["sizeC"], f["sizeT"]) == (64, 32, 10, 2, 1)
    assert f["voxelSizeUm"] == [0.5, 0.5, 2.0]
    assert f["channels"] == ["GFP", "RFP"] and f["emissionWavelengthsNm"] == [510.0, None]
    assert f["microscope"] == "Zeiss Lightsheet 7" and f["magnification"] == 20.0 and f["numericalAperture"] == 1.0
    assert f["planeCount"] == 20          # Z*C*T from the OME block, not the 2 IFDs present
    assert f["pageCount"] == 2
    assert f["imageDescriptionKind"] == "ome-xml" and f["flavor"] == "ome"
    assert d.summary == "OME-TIFF 64x32x10, 2 channels, uint16, 0.5x0.5x2 um"


def test_ome_tiff_in_bigtiff_big_endian():
    d = decode_bytes("s.ome.btf", make_tiff(bo=">", big=True, description=_OME))
    assert d.format == "ome-tiff" and d.fields["bigTiff"] is True


def test_imspector_flavor_is_marked():
    d = decode_bytes("s.ome.tif", make_tiff(description=_OME, software="ImSpector Pro 7.5"))
    assert d.format == "ome-tiff" and d.fields["flavor"] == "imspector"
    assert d.summary.startswith("ImSpector OME-TIFF")
    d = decode_bytes("s.tif", make_tiff(software="LaVision BioTec ImSpector"))
    assert d.format == "imspector-tiff" and d.fields["flavor"] == "imspector"


# ── Synthetic: ImageJ hyperstack (values cross-checked with tifffile) ────────

_IJ = "ImageJ=1.11a\nimages=24\nchannels=2\nslices=4\nframes=3\nhyperstack=true\nmode=grayscale\nloop=false\nspacing=1.5\nunit=um\n"


def test_imagej_hyperstack():
    d = decode_bytes("ij.tif", make_tiff(width=20, height=16, description=_IJ, xres=(40, 13), yres=(40, 13),
                                         resunit=1, software="tifffile.py", pages=4))
    assert d.format == "imagej-tiff"
    f = d.fields
    assert (f["sizeZ"], f["sizeC"], f["sizeT"], f["planeCount"]) == (4, 2, 3, 24)
    assert f["voxelSizeUm"] == [0.325, 0.325, 1.5]      # 13/40 um per pixel, spacing 1.5 um
    assert f["hyperstack"] is True and f["imagejVersion"] == "1.11a" and f["unit"] == "um"
    assert f["pageCount"] == 4 and f["pageCountTruncated"] is False
    assert d.summary == "ImageJ TIFF 20x16x4, 2 channels, 3 timepoints, uint16, 0.325x0.325x1.5 um, 24 planes"


def test_imagej_pixel_unit_gives_no_physical_size():
    d = decode_bytes("ij.tif", make_tiff(description="ImageJ=1.53f\nimages=5\nslices=5\nunit=pixel\n",
                                         xres=(1, 1), yres=(1, 1), resunit=1))
    assert d.fields["voxelSizeUm"] == [None, None, None] and d.fields["sizeZ"] == 5


# ── Synthetic: Micro-Manager ─────────────────────────────────────────────────

_MM_SUMMARY = json.dumps({"Summary": {
    "Channels": 2, "Slices": 10, "Frames": 5, "Positions": 1, "Width": 64, "Height": 32,
    "PixelType": "GRAY16", "PixelSize_um": 0.65, "z-step_um": 1.0, "ChNames": ["DAPI", "GFP"],
    "MicroManagerVersion": "2.0.1", "StartTime": "2026-02-03 10:11:12 +0100",
}}).encode()
_MM_TAG = json.dumps({"Camera": "Prime95B", "Binning": 2, "Exposure-ms": 50.0, "PixelSizeUm": 0.65}).encode()


def _mm_file(**kw):
    entries = [(51123, 2, len(_MM_TAG) + 1, None, _MM_TAG + b"\x00")]
    index = [(c, z, 0, 0, 0) for c in range(2) for z in range(10)]   # 20 planes for 2C x 10Z x 1T (truncated run)
    return make_tiff(width=64, height=32, bits=16, extra_entries=entries, mm_header=(_MM_SUMMARY, index), pages=2, **kw)


def test_micromanager_header_and_tag():
    d = decode_bytes("mm.tif", _mm_file())
    assert d.format == "micromanager-tiff"
    f = d.fields
    assert (f["sizeC"], f["sizeZ"], f["sizeT"], f["positions"]) == (2, 10, 5, 1)
    assert f["channels"] == ["DAPI", "GFP"] and f["voxelSizeUm"] == [0.65, 0.65, 1.0]
    assert f["pixelType"] == "GRAY16" and f["dtype"] == "uint16"
    assert f["micromanagerVersion"] == "2.0.1" and f["acquisitionDate"] == "2026-02-03 10:11:12 +0100"
    assert f["camera"] == "Prime95B" and f["binning"] == 2 and f["exposureMs"] == 50.0
    assert f["indexMapEntries"] == 20 and f["planeCount"] == 20   # the index map is what was actually written
    assert f["pageCount"] == 2
    assert d.summary == "Micro-Manager TIFF 64x32x10, 2 channels, 5 timepoints, uint16, 0.65x0.65x1 um"


def test_micromanager_from_description_summary_only():
    # Some exports carry the Summary JSON in ImageDescription with no private header.
    d = decode_bytes("mm.tif", make_tiff(description=_MM_SUMMARY.decode()))
    assert d.format == "micromanager-tiff"
    assert d.fields["imageDescriptionKind"] == "micromanager"
    assert d.fields["sizeC"] == 2 and d.fields["indexMapEntries"] is None and d.fields["planeCount"] == 100


# ── Totality and bounds ──────────────────────────────────────────────────────

def test_dispatch_extensions():
    for k in ("a.tif", "a.TIFF", "a.ome.tif", "a.ome.btf", "a.tf2", "a.tf8"):
        assert has_decoder_for(k), k
    assert not has_decoder_for("a.tf")


@pytest.mark.parametrize("bad", [
    b"",
    b"II",
    b"II\x2b\x00" + b"\x00" * 12,            # magic 43 but offset size not 8
    b"XX\x2a\x00" + b"\x08\x00\x00\x00",      # bad byte order
    b"II\x2a\x00" + b"\x00\x00\x00\x00",      # IFD offset 0
    b"II\x2a\x00" + b"\x08\x00\x00\x00" + b"\x00\x00",   # IFD with zero entries
    b"MM\x00\x2a" + b"\xff\xff\xff\xff",      # IFD offset past EOF
])
def test_rejects_non_tiff(bad):
    assert decode_bytes("a.tif", bad) is None


def test_zero_width_is_rejected():
    assert decode_bytes("a.tif", make_tiff(width=0)) is None


def test_chain_cap_and_truncation_flag():
    many = make_tiff(width=4, height=4, bits=8, pages=IFD_CHAIN_CAP + 10)
    f = decode_bytes("a.tif", many).fields
    assert f["pageCount"] == IFD_CHAIN_CAP and f["pageCountTruncated"] is True
    assert f["sizeZ"] is None and f["planeCount"] is None
    assert ">=%d pages" % IFD_CHAIN_CAP in decode_bytes("a.tif", many).summary


def test_self_referencing_chain_terminates():
    f = decode_bytes("a.tif", make_tiff(pages=2, chain_loop=True)).fields
    assert f["pageCount"] == 2 and f["pageCountTruncated"] is False


def test_truncated_file_reports_what_it_can():
    whole = make_tiff(pages=3, ifd_after_data=False)
    # Cut the file after the first IFD: later pages are unreachable.
    cut = whole[: len(whole) // 2]
    d = decode_bytes("a.tif", cut)
    assert d is not None and d.fields["pageCountTruncated"] is True


def test_huge_description_is_capped():
    d = decode_bytes("a.tif", make_tiff(description="x" * 70000))
    assert d is not None and d.fields["imageDescriptionKind"] == "text"


def test_json_safe():
    for blob in (make_tiff(description=_OME), _mm_file(), make_tiff(description=_IJ, xres=(40, 13), yres=(40, 13), resunit=1)):
        json.dumps(decode_bytes("a.tif", blob).to_dict(), allow_nan=False)


def test_sparse_source_outside_recorded_ranges_is_none():
    assert decode_source("a.tif", BytesSource(ranges={100: b"garbage"}, size=200)) is None


def test_magic_seeded_fuzz_is_total():
    # Valid header + random IFD bytes: drives the entry parsing, not just the magic check.
    import random
    rng = random.Random(4242)
    for _ in range(1500):
        n = rng.randint(0, 600)
        body = bytes(rng.getrandbits(8) for _ in range(n))
        big = rng.random() < 0.5
        bo = rng.choice(["<", ">"])
        if big:
            head = (b"II" if bo == "<" else b"MM") + struct.pack(bo + "HHHQ", 43, 8, 0, rng.choice([16, 8, 100]))
        else:
            head = (b"II" if bo == "<" else b"MM") + struct.pack(bo + "HI", 42, rng.choice([8, 10, 100]))
        try:
            d = decode_bytes("f.tif", head + body)
        except Exception as e:  # noqa: BLE001
            pytest.fail("raised: %r" % e)
        if d is not None:
            json.dumps(d.to_dict(), allow_nan=False)
