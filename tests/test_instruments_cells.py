"""Cell counters and solution analyzers (instruments/cells.py): golden cases
from real exports, declines across instruments and generic files, synthetic
edge cases, and a totality fuzz over every sniffer."""

import io
import json
import os
import random
import zipfile

import pytest

from scigantic_headers import DecodedHeader, decode_bytes
from scigantic_headers.decoders import BytesSource, INSTRUMENT_SEEK_BYTES
from scigantic_headers.instruments import INSTRUMENTS, TEXT_HEAD_BYTES
from scigantic_headers.instruments import cells

HERE = os.path.dirname(__file__)
FIXTURES = os.path.join(HERE, "fixtures", "instruments")
with open(os.path.join(FIXTURES, "cells-cases.json"), encoding="utf-8") as _fh:
    CASES = json.load(_fh)["cases"]


def _data(case):
    with open(os.path.join(FIXTURES, case["file"]), "rb") as fh:
        return fh.read()


def _sniff_all(data, key):
    """Formats of every sniffer of this family that accepts the bytes."""
    out = []
    for fn in cells.SNIFFERS:
        src = BytesSource(data, budget=INSTRUMENT_SEEK_BYTES)
        src.key = key
        try:
            h = fn(src)
        except Exception as e:  # a sniffer must never raise
            pytest.fail("%s raised %r on %s" % (fn.__name__, e, key))
        if h is not None:
            out.append(h.format)
    return out


def _assert_valid(h):
    assert isinstance(h, DecodedHeader)
    json.loads(json.dumps(h.to_dict(), allow_nan=False))


def test_cases_cover_every_instrument():
    assert len(CASES) >= 20
    names = {c["expected"]["format"] for c in CASES}
    family = {n for n, v in INSTRUMENTS.items() if v[2] in {
        "beckman_vi_cell_xr", "beckman_vi_cell_blu", "chemometec_nucleoview", "chemometec_nc_view",
        "revvity_matrix", "roche_cedex_hires", "roche_cedex_bioht", "beckman_pharmspec", "novabio_flex2"}}
    assert len(family) == 9
    assert names == family


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_golden(case):
    h = decode_bytes(case["key"], _data(case))
    assert h is not None
    _assert_valid(h)
    assert h.to_dict() == case["expected"]


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_only_its_own_sniffer_accepts(case):
    """Each fixture is claimed by its own instrument only, so no sniffer of
    this family takes another instrument's export."""
    got = set(_sniff_all(_data(case), case["key"]))
    assert got == {case["expected"]["format"]}


def test_identity_fields_match_the_table():
    xr = INSTRUMENTS[cells.XR]
    assert xr == ("Cell Counting", "Beckman Coulter", "beckman_vi_cell_xr", ["Brightfield"])
    assert INSTRUMENTS[cells.NUCLEOVIEW][3] == ["Dark Field"]
    assert INSTRUMENTS[cells.NCVIEW][3] == ["Fluorescence"]
    assert INSTRUMENTS[cells.PHARMSPEC][:1] == ("Solution Analyzer",)
    assert INSTRUMENTS[cells.PHARMSPEC][3] == ["Light Obscuration"]
    assert INSTRUMENTS[cells.FLEX2][3] == [
        "Metabolite Detection", "Blood Gas", "pH", "Osmolality", "Cell Counting"]


GENERIC = {
    "plain.csv": b"a,b,c\n1,2,3\n4,5,6\n",
    "samples.csv": b"Sample ID,Value\n001,2\n002,3\n",
    "well.csv": b"Well Name,Count\nA1,3\n",
    "date.csv": b"Date & Time,Sample ID,Sample Type\n1/1/2024,1,x\n",
    "cedex.csv": b"Cedex ID,Comment\n1,x\n",
    "semi.csv": b"Image;x;\nViability (%);50;\n",
    "names.csv": b"NAME,SAMPLE ID,TAGS\nx,1,\n",
    "notes.txt": b"Vi-CELL XR is a cell counter. Sample ID : none\n",
    "zero.txt": b"0\t2023-01-01 00:00:00\tsomething\n",
    "readme.txt": b"hello world\n" * 50,
    "data.xml": b'<?xml version="1.0"?><root><Sample ID="1"/></root>',
    "data.json": b'{"Sample ID": "1", "Viability (%)": 90}',
    "page.xls": b"<html><body><table><tr><td>hello</td></tr></table></body></html>",
    "empty.csv": b"",
    "bom.csv": b"\xef\xbb\xbf",
    "binary.xls": b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1" + b"\x00" * 600,
}


@pytest.mark.parametrize("key", sorted(GENERIC))
def test_generic_files_decline(key):
    assert _sniff_all(GENERIC[key], key) == []


def test_generic_csv_still_reaches_csv_decoder():
    h = decode_bytes("plain.csv", GENERIC["plain.csv"])
    assert h is None or h.format not in INSTRUMENTS


# ------------------------------------------------------------ synthetic


def _case(name):
    return next(c for c in CASES if c["name"] == name)


def test_xr_txt_utf16_without_bom():
    data = _data(_case("Beckman_Vi-Cell-XR_example07_instrumentOutput.txt"))
    h = decode_bytes("x.txt", data.decode("latin-1").encode("utf-16-le"))
    assert h is not None and h.fields["softwareVersion"] == "2.06.3"
    assert h.fields["serialNumber"] == "AW123456789"


def test_numeric_looking_sample_ids_stay_strings():
    text = ("Sample ID,Cell count,Viable cells,Total (x10^6) cells/mL,Viable (x10^6) cells/mL,"
            "Viability (%),Average diameter (μm),Average viable diameter (μm)\n"
            "007,1,1,1,1,99,1,1\n0010,1,1,1,1,99,1,1\n1e3,1,1,1,1,99,1,1\n")
    h = decode_bytes("x.csv", text.encode())
    assert h.format == cells.BLU
    assert h.fields["sampleIds"] == ["007", "0010", "1e3"]
    assert h.fields["sampleCount"] == 3


def test_bioht_serial_keeps_leading_zero():
    line = "0\t2023-09-17 13:04:06\t#ARC-FILE#\t1.1a\t2021-05-01\t2023-09-17\tCEDEX BIO HT\t012345\t6.0.0.1905 (1905)\tADMIN\n"
    row = "40\t2023-09-15 16:55:51\t0042\t\tSAM\t        \tGLN2B\t        \tmmol/L\t \t  2.45\t  0.17138\tR\n"
    h = decode_bytes("x.txt", (line + row).encode())
    assert h.fields["serialNumber"] == "012345"
    assert h.fields["sampleIds"] == ["0042"]
    assert h.fields["measurementCount"] == 1


def test_counts_omitted_when_file_exceeds_leading_read():
    line = "0\t2023-09-17 13:04:06\t#ARC-FILE#\t1.1a\t2021-05-01\t2023-09-17\tCEDEX BIO HT\t1\t6.0.0.1905 (1905)\tADMIN\n"
    row = "40\t2023-09-15 16:55:51\tS%d\t\tSAM\t        \tGLN2B\t        \tmmol/L\t \t  2.45\t  0.17138\tR\n"
    body = "".join(row % i for i in range(5000))
    data = (line + body).encode()
    assert len(data) > TEXT_HEAD_BYTES
    h = decode_bytes("x.txt", data)
    assert h.format == cells.BIOHT
    assert "sampleCount" not in h.fields and "measurementCount" not in h.fields
    assert h.summary == "Roche Cedex BioHT txt export"
    assert h.fields["serialNumber"] == "1"

    csv_hdr = "Cedex ID,System name,Sample identifer,Viable Cell Conc.\n"
    csv = csv_hdr + "".join("1,sys,S%06d,1\n" % i for i in range(20000))
    h = decode_bytes("x.csv", csv.encode())
    assert h.format == cells.HIRES and "sampleCount" not in h.fields
    assert len(h.fields["sampleIds"]) == cells.MAX_IDS


def test_nucleoview_offset_and_bad_dates():
    base = ("Assay; ;\nImage;x-1-c-a;\nOperator;me;\nViability (%%);68.2;\nLive (cells/ml);1;\n"
            "Total (cells/ml);2;\nDate time;%s\nTime zone offset;%s\n")
    h = decode_bytes("x.csv", (base % ("20231120 163016", "-5")).encode())
    assert h.fields["runDate"] == "2023-11-20T16:30:16-05:00"
    h = decode_bytes("x.csv", (base % ("20231120 163016", "5.5")).encode())
    assert h.fields["runDate"] == "2023-11-20T16:30:16"
    h = decode_bytes("x.csv", (base % ("20231399 163016", "0")).encode())
    assert "runDate" not in h.fields


def _xlsx(rows, sheet="Sheet1"):
    shared = []
    xml_rows = []
    for r, row in enumerate(rows, 1):
        cells_xml = []
        for c, v in enumerate(row):
            ref = chr(65 + c) + str(r)
            if v is None:
                continue
            if isinstance(v, str):
                shared.append(v)
                cells_xml.append('<c r="%s" t="s"><v>%d</v></c>' % (ref, len(shared) - 1))
            else:
                cells_xml.append('<c r="%s"><v>%s</v></c>' % (ref, v))
        if cells_xml:
            xml_rows.append('<row r="%d">%s</row>' % (r, "".join(cells_xml)))
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("xl/workbook.xml", '<workbook %s %s><sheets><sheet name="%s" sheetId="1" r:id="rId1"/></sheets></workbook>' % (ns, rns, sheet))
        z.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr("xl/sharedStrings.xml", "<sst %s>%s</sst>" % (ns, "".join("<si><t>%s</t></si>" % s for s in shared)))
        z.writestr("xl/worksheets/sheet1.xml", "<worksheet %s><sheetData>%s</sheetData></worksheet>" % (ns, "".join(xml_rows)))
    return buf.getvalue()


def test_xr_204_table_in_xlsx_is_not_the_report():
    """A 2.04 multi-sample table starts at row 4, the same row as the report's
    "Sample ID"; the next row tells them apart."""
    rows = [
        ["Vi-CELL XR 2.04"], ["c:\\x.xls"], ["Unit S/N:  0123"],
        ["Sample ID", "File name", "Cell type", "Sample date/time"],
        [None, None, None, None],
        ["001", "001.smp", "CHO", "2 Mar 2023  8:31:38 AM"],
        [2, "2.smp", "CHO", "2 Mar 2023  8:41:38 PM"],
        ["bad", "x", "x", "not a date"],
    ]
    h = decode_bytes("x.xlsx", _xlsx(rows))
    assert h.format == cells.XR and h.fields["exportFormat"] == "xlsx"
    assert h.fields["serialNumber"] == "0123"
    assert h.fields["sampleIds"] == ["001", "2"]
    assert h.fields["sampleCount"] == 2


def test_flex2_units_row_is_not_a_sample():
    text = ('"Date & Time","Sample ID","Sample Type","Gluc","pH","Vessel ID"\n'
            ',,,g / L,-,\n2/9/2025  10:40:00 AM,0001,NO CDV,4.5,7.1,R1\n')
    h = decode_bytes("SampleResults2025-02-09_104000.csv", text.encode())
    assert h.format == cells.FLEX2
    assert h.fields["sampleIds"] == ["0001"] and h.fields["sampleCount"] == 1
    assert h.fields["analytes"] == ["Gluc", "pH"]


# ----------------------------------------------------------------- fuzz

SEED = 4242


def _mutations(rng, data):
    yield data[: rng.randint(0, len(data))]
    b = bytearray(data)
    for _ in range(rng.randint(1, 40)):
        if b:
            b[rng.randrange(len(b))] = rng.getrandbits(8)
    yield bytes(b)
    lines = data.split(b"\n")
    rng.shuffle(lines)
    yield b"\n".join(lines)
    yield data.replace(b",", b"\t") if rng.random() < 0.5 else data.replace(b"\t", b",")
    yield bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 3000)))


def test_fuzz_totality():
    rng = random.Random(SEED)
    keys = ["x.csv", "x.txt", "x.xlsx", "x.xls"]
    seeds = [_data(c) for c in CASES] + list(GENERIC.values())
    n = 0
    for data in seeds:
        for _ in range(12):
            for m in _mutations(rng, data):
                for fn in cells.SNIFFERS:
                    src = BytesSource(m, budget=INSTRUMENT_SEEK_BYTES)
                    h = fn(src)
                    if h is not None:
                        _assert_valid(h)
                    n += 1
                h = decode_bytes(rng.choice(keys), m)
                if h is not None:
                    _assert_valid(h)
    assert n > 10000
