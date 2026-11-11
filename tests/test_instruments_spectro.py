"""Spectrophotometer exports (VISIONlite, NanoDrop 8000 / Eight / One, Qubit 4 /
Flex, Genesys 30, Genesys On-Board): golden decodes of real exports from
allotropy's test data, a cross-check of each decode against the ASM allotropy
produces for the same file, declines, and a totality fuzz."""

import glob
import io
import json
import os
import random
import zipfile

import pytest

from scigantic_headers import DecodedHeader, decode_bytes
from scigantic_headers.decoders import BytesSource
from scigantic_headers.instruments import INSTRUMENTS, spectro

HERE = os.path.dirname(__file__)
FX = os.path.join(HERE, "fixtures", "instruments")
CASES = json.load(open(os.path.join(FX, "spectro-cases.json"), encoding="utf-8"))["cases"]


def _bytes(case):
    with open(os.path.join(FX, case["file"]), "rb") as fh:
        return fh.read()


def _valid_or_none(h):
    if h is None:
        return
    assert isinstance(h, DecodedHeader)
    assert isinstance(h.format, str) and h.format
    json.loads(json.dumps(h.to_dict(), allow_nan=False))


def test_declared_instruments():
    table = {
        "Thermo Fisher Scientific VISIONlite": ("thermo_fisher_visionlite", ["Absorbance"]),
        "Thermo Fisher Scientific NanoDrop 8000": ("thermo_fisher_nanodrop_8000", ["Absorbance"]),
        "Thermo Fisher Scientific NanoDrop Eight": ("thermo_fisher_nanodrop_eight", ["Absorbance"]),
        "Thermo Fisher Scientific NanoDrop One": ("thermo_fisher_nanodrop_one", ["Absorbance"]),
        "Thermo Fisher Scientific Qubit 4": ("thermo_fisher_qubit4", ["Fluorescence"]),
        "Thermo Fisher Scientific Qubit Flex": ("thermo_fisher_qubit_flex", ["Fluorescence"]),
        "Thermo Fisher Scientific Genesys30": ("thermo_fisher_genesys30", ["Absorbance"]),
        "Thermo Fisher Scientific Genesys On-Board": ("thermo_fisher_genesys_on_board", ["Absorbance"]),
    }
    for name, (parser, modes) in table.items():
        category, vendor, p, m = INSTRUMENTS[name]
        assert (category, vendor, p, m) == ("Spectrophotometry", "Thermo Fisher Scientific", parser, modes)


def test_cases_cover_every_instrument():
    seen = set()
    for c in CASES:
        f = c["expected"]["fields"]
        seen.update(f.get("instrumentCandidates", [f.get("instrumentSoftware")]))
    assert {n for n in INSTRUMENTS if INSTRUMENTS[n][0] == "Spectrophotometry"} <= seen


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_golden(case):
    h = decode_bytes(case["key"], _bytes(case))
    assert h is not None
    assert {"format": h.format, "summary": h.summary, "fields": h.fields} == case["expected"]
    _valid_or_none(h)


@pytest.mark.parametrize("case", [c for c in CASES if "asm" in c], ids=[c["name"] for c in CASES if "asm" in c])
def test_agrees_with_allotropy_asm(case):
    """Every value both sides state must match what allotropy wrote to ASM for
    the same input file."""
    f = case["expected"]["fields"]
    asm = case["asm"]
    checked = 0
    assert f.get("sampleCount") == asm["sampleCount"]
    checked += 1
    if "wavelengthsNm" in asm or "wavelengthsNm" in f:
        assert sorted(f.get("wavelengthsNm", [])) == asm.get("wavelengthsNm", [])
        checked += 1
    for k in ("serialNumber", "softwareVersion", "excitation", "emission"):
        if k in asm or k in f:
            assert f.get(k) == asm.get(k), k
            checked += 1
    # NanoDrop One CSV: allotropy reads the assay from the file name, not the bytes.
    if "application" in f and "experimentType" in asm and case["name"] != "nanodrop_one_dsdna.csv":
        assert f["application"] == asm["experimentType"]
        checked += 1
    if "dataCubePoints" in asm:
        rng = f.get("wavelengthRangeNm") or f.get("spectrumRangeNm")
        assert rng == asm["dataCubeRangeNm"]
        if "measurementCount" in f:
            assert f["measurementCount"] == asm["dataCubePoints"]
        checked += 1
    sw = asm.get("softwareName")
    if sw and "instrumentCandidates" in f:
        # The ambiguous layout: allotropy's answer must be one of the candidates.
        assert any(sw in name for name in f["instrumentCandidates"])
    assert checked >= 1


def test_nanodrop_one_csv_does_not_claim_a_specific_assay():
    # allotropy takes dsDNA from the file name; the bytes say only "Nucleic Acid".
    case = next(c for c in CASES if c["name"] == "nanodrop_one_dsdna.csv")
    assert case["expected"]["fields"]["application"] == "Nucleic Acid"


def test_ambiguous_layout_names_both_programs():
    for name in ("visionlite_fixed_tab_utf16.csv", "genesys_on_board_fixed.csv"):
        f = next(c for c in CASES if c["name"] == name)["expected"]["fields"]
        assert "instrumentSoftware" not in f and "allotropyParser" not in f
        assert f["instrumentCandidates"] == [
            "Thermo Fisher Scientific VISIONlite",
            "Thermo Fisher Scientific Genesys On-Board",
        ]
        assert f["allotropyParserCandidates"] == ["thermo_fisher_visionlite", "thermo_fisher_genesys_on_board"]


def test_visionlite_rate_is_recognized_but_not_convertible():
    f = next(c for c in CASES if c["name"] == "visionlite_rate.csv")["expected"]["fields"]
    assert f["instrumentSoftware"] == "Thermo Fisher Scientific VISIONlite"
    assert "allotropyParser" not in f and f["allotropySupported"] is False


# ---------------------------------------------------------------- truncation


def test_large_export_identified_without_counts():
    """A file larger than the leading read is still identified; sample counts
    are left out because the table runs past what was read."""
    case = next(c for c in CASES if c["name"] == "nanodrop_eight_dsdna.txt")
    data = _bytes(case)
    lines = data.split(b"\r\n")
    body = b"\r\n".join(lines[:4]) + b"\r\n" + (lines[4] + b"\r\n") * 200
    assert len(body) > 256 * 1024
    h = decode_bytes("big.txt", body)
    assert h.format == "Thermo Fisher Scientific NanoDrop Eight"
    assert "sampleCount" not in h.fields
    assert h.fields["serialNumber"] == "NDE8675309"


def test_second_table_after_blank_line_is_not_counted():
    case = next(c for c in CASES if c["name"] == "nanodrop_8000_dna.txt")
    data = _bytes(case)
    assert b"Curve Type" in data  # the standard-curve table after the samples
    assert decode_bytes("x.txt", data).fields["sampleCount"] == 6


# ---------------------------------------------------------------- xlsx


def _xlsx(rows, sheet="Sheet1", date1904=False):
    shared = []
    xml_rows = []
    for r, row in enumerate(rows, 1):
        cells = []
        for c, v in enumerate(row):
            ref = chr(65 + c % 26) if c < 26 else "A" + chr(65 + c - 26)
            ref += str(r)
            if isinstance(v, str):
                shared.append(v.replace("&", "&amp;"))
                cells.append('<c r="%s" t="s"><v>%d</v></c>' % (ref, len(shared) - 1))
            elif v is not None:
                cells.append('<c r="%s"><v>%s</v></c>' % (ref, v))
        xml_rows.append('<row r="%d">%s</row>' % (r, "".join(cells)))
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    pr = '<workbookPr date1904="1"/>' if date1904 else ""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("xl/workbook.xml", '<workbook %s %s>%s<sheets><sheet name="%s" sheetId="1" r:id="rId1"/></sheets></workbook>' % (ns, rns, pr, sheet))
        z.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr("xl/sharedStrings.xml", "<sst %s>%s</sst>" % (ns, "".join("<si><t>%s</t></si>" % s for s in shared)))
        z.writestr("xl/worksheets/sheet1.xml", "<worksheet %s><sheetData>%s</sheetData></worksheet>" % (ns, "".join(xml_rows)))
    return buf.getvalue()


_Q4_COLS = ["Run ID", "Assay Name", "Test Name", "Test Date", "Qubit® tube conc.", "Units", "Original sample conc.", "Units",
            "Sample Volume (µL)", "Dilution Factor", "Std 1 RFU", "Std 2 RFU", "Std 3 RFU", "Excitation", "Emission", "Green RFU", "Far Red RFU"]


def test_synthetic_qubit4_xlsx_and_1904_dates():
    row = [45272, "RNA HS", "S1", 45272.5, 30.4, "ng/mL", 6.08, "ng/µL", 1, 200, 49.4, 21133.1, None, "Blue", "Green", 1303.1, None]
    h = decode_bytes("q.xlsx", _xlsx([_Q4_COLS, row, row]))
    assert h.format == "Thermo Fisher Scientific Qubit 4"
    assert h.fields["runDate"] == "2023-12-12" and h.fields["sampleCount"] == 2
    assert h.fields["application"] == "RNA HS"
    h = decode_bytes("q.xlsx", _xlsx([_Q4_COLS, row], date1904=True))
    assert h.fields["runDate"] == "2027-12-13"


def test_synthetic_nanodrop_one_protein_sheet():
    cols = ["Date", "Sample Name", "Protein A280 (mg/mL)", "A280"]
    h = decode_bytes("p.xlsx", _xlsx([cols, ["1/2/2024 3:04:05 PM", "BSA", 1.2, 0.8]], sheet="Protein A280 1_2_2024"))
    assert h.format == "Thermo Fisher Scientific NanoDrop One"
    assert h.fields["application"] == "Protein A280" and h.fields["wavelengthsNm"] == [280]


def test_assays_that_disagree_are_not_reported():
    row1 = [45272, "RNA HS", "S1", 45272.5, 1, "ng/mL", 1, "ng/µL", 1, 200, 1, 1, None, "Blue", "Green", 1, None]
    row2 = list(row1)
    row2[1] = "dsDNA HS"
    h = decode_bytes("q.xlsx", _xlsx([_Q4_COLS, row1, row2]))
    assert h is not None and "application" not in h.fields and h.fields["sampleCount"] == 2


def test_unrelated_xlsx_declines():
    assert decode_bytes("x.xlsx", _xlsx([["Date", "Sample Name", "Weight"], ["1/1/2024", "a", 1]])) is None
    assert decode_bytes("x.xlsx", _xlsx([["Run ID", "Test Date", "Assay Name", "Sample Name"], ["1", "2", "3", "4"]])) is None


# ---------------------------------------------------------------- declines

GENERIC = {
    "plain.csv": b"a,b,c\n1,2,3\n4,5,6\n",
    "sample.csv": b"Sample Name,Value\nx,1\ny,2\n",
    "sample_tab.csv": "Sample name\tA260\tNotes\r\nx\t1\thello\r\n".encode("utf-16"),
    "wells.csv": b"Well Position,Absorbance,Result\nA1,0.1,0.1\n",
    "plain.tsv": b"x\ty\n1\t2\n",
    "kv.txt": b"Application:\tdsDNA\nNotes:\tnothing\n",
    "notes.txt": b"Date: today\nhello world\n",
    "scan.csv": b"Scan,abc\nMode,ABS\n,,\nwavelength(nm) , ABS\n",
    "nm.csv": b"file,x,y\nnm,A\n400,0.1\n",
    "doc.xml": b"<?xml version='1.0'?><root><Sample Name='x'/></root>",
    "doc.json": b'{"Application": "dsDNA", "Serial Number": "NDE1"}',
    "empty.csv": b"",
    "empty.txt": b"",
}


@pytest.mark.parametrize("key", sorted(GENERIC))
def test_generic_files_decline(key):
    for sniff in spectro.SNIFFERS:
        assert sniff(BytesSource(GENERIC[key])) is None, (key, sniff.__name__)
    h = decode_bytes(key, GENERIC[key])
    assert h is None or not h.format.startswith("Thermo Fisher"), key


# Which sniffer is expected to claim each fixture; all others must decline.
OWNER = {
    "visionlite": "sniff_visionlite",
    "genesys_on_board": "sniff_visionlite",
    "nanodrop_8000": "sniff_nanodrop_8000",
    "nanodrop_eight": "sniff_nanodrop_eight",
    "nanodrop_one_dsdna.csv": "sniff_nanodrop_one_csv",
    "nanodrop_one_dsdna.xlsx": "sniff_nanodrop_one_xlsx",
    "qubit4.csv": "sniff_qubit_csv",
    "qubit4.xlsx": "sniff_qubit_xlsx",
    "qubit_flex.xlsx": "sniff_qubit_xlsx",
    "qubit_flex": "sniff_qubit_csv",
    "genesys30": "sniff_genesys30",
}


def _owner(name):
    for prefix in sorted(OWNER, key=len, reverse=True):
        if name.startswith(prefix):
            return OWNER[prefix]
    raise KeyError(name)


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_only_the_owning_sniffer_claims(case):
    data = _bytes(case)
    owner = _owner(case["name"])
    for sniff in spectro.SNIFFERS:
        got = sniff(BytesSource(data))
        if sniff.__name__ == owner:
            assert got is not None
        else:
            assert got is None, (case["name"], sniff.__name__)


def _other_family_fixtures():
    out = []
    for path in sorted(glob.glob(os.path.join(FX, "*", "*"))):
        if os.sep + "spectro" + os.sep in path or path.endswith(".md"):
            continue
        out.append(path)
    for path in sorted(glob.glob(os.path.join(HERE, "fixtures", "*"))):
        if os.path.isfile(path) and os.path.getsize(path) < 4 * 1024 * 1024:
            out.append(path)
    return out


@pytest.mark.parametrize("path", _other_family_fixtures())
def test_other_fixtures_decline(path):
    with open(path, "rb") as fh:
        data = fh.read()
    for sniff in spectro.SNIFFERS:
        assert sniff(BytesSource(data)) is None, (path, sniff.__name__)


# ---------------------------------------------------------------- totality fuzz

SEED = 1729


def _mutations(rng, data):
    yield data[: rng.randint(0, len(data))]
    b = bytearray(data)
    for _ in range(rng.randint(1, 12)):
        if b:
            b[rng.randrange(len(b))] = rng.getrandbits(8)
    yield bytes(b)
    b = bytearray(data)
    for _ in range(rng.randint(1, 6)):
        if b:
            i = rng.randrange(len(b))
            del b[i : i + rng.randint(1, 40)]
    yield bytes(b)
    lines = data.splitlines(True)
    rng.shuffle(lines)
    yield b"".join(lines)
    yield data.replace(b",", b"\t") if rng.random() < 0.5 else data.replace(b"\t", b",")


def test_fuzz_mutated_fixtures():
    rng = random.Random(SEED)
    for _ in range(6):
        for case in CASES:
            for m in _mutations(rng, _bytes(case)):
                for sniff in spectro.SNIFFERS:
                    _valid_or_none(sniff(BytesSource(m)))
                _valid_or_none(decode_bytes(case["key"], m))


def test_fuzz_random_bytes_and_prefixes():
    rng = random.Random(SEED + 1)
    prefixes = [
        b"", b"\xff\xfe", b"PK\x03\x04", b"sep=,\n", b"Scan,x\nMode,ABS\nLower,a\nUpper,\nLength,1e9\n,,\nwavelength(nm)\n",
        b"Application:\tdsDNA\nSerial Number:\tNDE1\n", b"Plate ID\tWell\tSample ID\tUser ID\tDate\tTime\tx\ty\n",
        b"Run ID,Assay Name,Test Name,Test Date,Qubit tube conc,Green RFU\n", b"Sample Name,500 nm   [A]\n",
        b"Sample Name,Dilution factor,Ordinate [A],Concentration [\n", b"a.csv,b,c,1/1/2024 1:00:00 AM\nnm,A\n",
        b"Date,Sample Name,Nucleic Acid(ng/uL),A260\n",
    ]
    for _ in range(1500):
        p = rng.choice(prefixes)
        tail = bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 400)))
        if rng.random() < 0.5:
            tail = bytes(rng.choice(b"0123456789.,\t\r\n -:/AnmE") for _ in range(rng.randint(0, 400)))
        data = p + tail
        for sniff in spectro.SNIFFERS:
            _valid_or_none(sniff(BytesSource(data)))
        for key in ("f.csv", "f.txt", "f.tsv", "f.xlsx"):
            _valid_or_none(decode_bytes(key, data))
