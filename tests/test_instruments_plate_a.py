"""Plate readers, family A: Gen5, Gen5 Image, SoftMax Pro, MARS, SMART
Control, Magellan, SkanIt. Golden fixtures are real exports from allotropy's
test suite (see fixtures/instruments/plate_a/NOTICE.md); their expected
decodes were checked against allotropy's ASM output for the same files."""

import io
import json
import os
import random
import zipfile

import pytest

from scigantic_headers import decode_bytes
from scigantic_headers.decoders import BytesSource, DecodedHeader
from scigantic_headers.instruments import INSTRUMENTS, TEXT_HEAD_BYTES
from scigantic_headers.instruments import plate_a as P

HERE = os.path.dirname(__file__)
FIX = os.path.join(HERE, "fixtures", "instruments", "plate_a")
with open(os.path.join(HERE, "fixtures", "instruments", "plate_a-cases.json"), encoding="utf-8") as fh:
    CASES = json.load(fh)["cases"]

SNIFFERS = {
    P.GEN5: P.sniff_gen5,
    P.GEN5_IMAGE: P.sniff_gen5_image,
    P.SOFTMAX: P.sniff_softmax_pro,
    P.MARS: P.sniff_mars,
    P.SMART: P.sniff_smart_control,
    P.MAGELLAN: P.sniff_magellan,
    P.SKANIT: P.sniff_skanit,
}


def _bytes(name):
    with open(os.path.join(FIX, name), "rb") as fh:
        return fh.read()


def _src(data):
    return BytesSource(data, budget=16 << 20)


def _valid_or_none(result):
    if result is None:
        return
    assert isinstance(result, DecodedHeader)
    assert result.format and isinstance(result.summary, str)
    json.loads(json.dumps(result.to_dict(), allow_nan=False))


def _xlsx(sheets, date1904=False):
    """A minimal .xlsx: {sheet name: rows of str/number}."""
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    shared, parts, entries, rels = [], [], [], []
    for i, (name, rows) in enumerate(sheets.items(), 1):
        xml_rows = []
        for r, row in enumerate(rows, 1):
            cells = []
            for c, v in enumerate(row):
                ref = chr(65 + c) + str(r)
                if v is None:
                    continue
                if isinstance(v, str):
                    shared.append(v)
                    cells.append('<c r="%s" t="s"><v>%d</v></c>' % (ref, len(shared) - 1))
                else:
                    cells.append('<c r="%s"><v>%s</v></c>' % (ref, v))
            xml_rows.append('<row r="%d">%s</row>' % (r, "".join(cells)))
        parts.append(("xl/worksheets/sheet%d.xml" % i, "<worksheet %s><sheetData>%s</sheetData></worksheet>" % (ns, "".join(xml_rows))))
        entries.append('<sheet name="%s" sheetId="%d" r:id="rId%d"/>' % (name, i, i))
        rels.append('<Relationship Id="rId%d" Target="worksheets/sheet%d.xml"/>' % (i, i))
    pr = '<workbookPr date1904="1"/>' if date1904 else ""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("xl/workbook.xml", "<workbook %s %s>%s<sheets>%s</sheets></workbook>" % (ns, rns, pr, "".join(entries)))
        z.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">%s</Relationships>' % "".join(rels))
        z.writestr("xl/sharedStrings.xml", "<sst %s>%s</sst>" % (ns, "".join("<si><t>%s</t></si>" % s for s in shared)))
        for n, x in parts:
            z.writestr(n, x)
    return buf.getvalue()


# ── golden fixtures ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_golden_fixture(case):
    data = _bytes(case["file"])
    got = decode_bytes(case["file"], data)
    assert got is not None
    assert got.to_dict() == case["expected"]


def test_fixtures_cover_every_instrument_and_are_small():
    names = {c["expected"]["format"] for c in CASES}
    assert names == set(SNIFFERS)
    for c in CASES:
        size = len(_bytes(c["file"]))
        assert size <= (200 if c["file"].endswith(".xlsx") else 32) * 1024


def test_identity_matches_the_instrument_table():
    table = {
        P.GEN5: ("Plate Reader", ["Absorbance", "Fluorescence", "Luminescence"], "agilent_gen5"),
        P.GEN5_IMAGE: ("Plate Reader", ["Optical Imaging"], "agilent_gen5_image"),
        P.SOFTMAX: ("Plate Reader", ["Absorbance", "Fluorescence", "Luminescence"], "moldev_softmax_pro"),
        P.MARS: ("Plate Reader", ["Absorbance", "Fluorescence", "Luminescence"], "bmg_mars"),
        P.SMART: ("Plate Reader", ["Fluorescence"], "bmg_labtech_smart_control"),
        P.MAGELLAN: ("Plate Reader", ["Absorbance"], "tecan_magellan"),
        P.SKANIT: ("Plate Reader", ["Absorbance", "Fluorescence", "Luminescence"], "thermo_skanit"),
    }
    for name, (cat, modes, parser) in table.items():
        c, _vendor, p, m = INSTRUMENTS[name]
        assert (c, m, p) == (cat, modes, parser)


# ── telling the family apart ──────────────────────────────────────────────────

@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_only_the_right_sniffer_claims_each_fixture(case):
    data = _bytes(case["file"])
    for name, sniff in SNIFFERS.items():
        got = sniff(_src(data))
        if name == case["expected"]["format"]:
            assert got is not None
        else:
            assert got is None, name


def test_gen5_and_gen5_image_split_on_the_read_steps():
    text = _bytes("gen5-image-cell-count.txt").decode("cp1252")
    assert P.sniff_gen5(_src(text.encode("cp1252"))) is None
    # The same header with an absorbance step instead of imaging is Gen5.
    photometric = text.replace("Read\tImage Montage", "Read\tAbsorbance Endpoint\r\n\tWavelengths:  450")
    assert P.sniff_gen5_image(_src(photometric.encode("cp1252"))) is None
    g = P.sniff_gen5(_src(photometric.encode("cp1252")))
    assert g.fields["readModes"][0] == "Absorbance" and g.fields["wavelengthsNm"] == [450]
    # Both kinds of step: Gen5 claims it and says imaging is present too.
    mixed = text.replace("Read\tImage Montage", "Read\tAbsorbance Endpoint\r\n\tWavelengths:  450\r\nRead\tImage Montage")
    assert P.sniff_gen5_image(_src(mixed.encode("cp1252"))) is None
    assert P.sniff_gen5(_src(mixed.encode("cp1252"))).fields["readModes"] == ["Absorbance", "Imaging"]
    # A Gen5 header with no recognizable read step stays Gen5 (no modes).
    bare = text.replace("Read\tImage Montage", "Shake\tLinear")
    assert P.sniff_gen5_image(_src(bare.encode("cp1252"))) is None
    assert "readModes" not in P.sniff_gen5(_src(bare.encode("cp1252"))).fields


def test_softmax_utf16_without_bom_and_utf8():
    data = _bytes("softmax-fluorescence-utf16.txt")
    assert data[:2] == b"\xff\xfe"
    text = data.decode("utf-16")
    expected = P.sniff_softmax_pro(_src(data)).fields
    for enc in ("utf-16-le", "utf-16-be", "utf-8"):
        assert P.sniff_softmax_pro(_src(text.encode(enc))).fields == expected


def test_counts_are_dropped_when_the_file_is_longer_than_the_head():
    data = _bytes("gen5-fluorescence-multiplate.txt")
    assert P.sniff_gen5(_src(data)).fields["plateCount"] == 3
    long = data + b"\r\n" * (TEXT_HEAD_BYTES // 2)
    f = P.sniff_gen5(_src(long)).fields
    assert "plateCount" not in f and f["softwareVersion"] == "3.12.08"
    smp = _bytes("softmax-luminescence-two-plates.txt")
    assert P.sniff_softmax_pro(_src(smp)).fields["plateCount"] == 2
    f = P.sniff_softmax_pro(_src(smp + b"\x00\x00" * TEXT_HEAD_BYTES)).fields
    assert "plateCount" not in f and f["blockCount"] == 2


def test_placeholders_are_not_reported():
    f = P.sniff_gen5(_src(_bytes("gen5-fluorescence-multiplate.txt"))).fields
    assert "instrumentModel" not in f and "serialNumber" not in f  # file says n/a


def test_mars_absorbance_label_with_two_numbers_reports_no_wavelength():
    # "Raw Data (450 1)": allotropy reads detection wavelength 1 nm; the label
    # is not a wavelength pair, so nothing is reported.
    f = P.sniff_mars(_src(_bytes("mars-absorbance-96-two-numbers.csv"))).fields
    assert "wavelengthsNm" not in f and "excitationWavelengthsNm" not in f
    assert f["plateFormat"] == 96 and f["sampleCount"] == 96


def test_skanit_serial_date_honors_the_1904_system():
    rows = [["Measurement results"], ["run.skax"], [45835.458333333336], [" "], ["Absorbance 1"], ["Wavelength: 450 nm"]]
    f = P.sniff_skanit(_src(_xlsx({"Absorbance 1_01": rows}))).fields
    assert f["runDate"] == "2025-06-27T11:00:00" and f["wavelengthsNm"] == [450]
    f = P.sniff_skanit(_src(_xlsx({"Absorbance 1_01": rows}, date1904=True))).fields
    assert f["runDate"] == "2029-06-28T11:00:00"


def test_synthetic_xlsx_exports():
    smart = _xlsx({
        "Microplate End point": [["User: u"]],
        "Table End point": [["Well", "A01", "A02", "B01"]],
        "Protocol Information": [["Test Name: t"], ["Date: 1/2/2025"], ["Time: 1:00:00 PM"], ["Measurement type:", "Fluorescence (FI)"],
                                 ["Microplate name:", "Greiner 384"], ["Excitation:", "485-20"], ["Emission:", "528-20"]],
    })
    f = P.sniff_smart_control(_src(smart)).fields
    assert (f["plateFormat"], f["sampleCount"], f["excitationWavelengthsNm"], f["wavelengthsNm"]) == (384, 3, [485], [528])
    mag = _xlsx({"Magellan Sheet 1": [["Raw data", "Well positions"], [0.1, "A1"], [0.2, "A2"],
                                      ["Date of measurement: 2024-02-03/Time of measurement: 9:05:00"], ["m.mth"], ["w.wsp"],
                                      ["405nm"], ["user"], ["SUNRISE"], ["Instrument serial number: 42"],
                                      ["Measurement mode: Absorbance"], ["Measurement wavelength: 405 nm"]]})
    f = P.sniff_magellan(_src(mag)).fields
    assert f["runDate"] == "2024-02-03T09:05:00" and f["instrumentName"] == "SUNRISE" and f["sampleCount"] == 2
    # Without the Magellan metadata block, a sheet with a "Well positions"
    # column is not claimed.
    assert P.sniff_magellan(_src(_xlsx({"Sheet1": [["Raw data", "Well positions"], [0.1, "A1"]]}))) is None


# ── declining everything else ─────────────────────────────────────────────────

GENERIC = {
    "plain.txt": b"Software Version is 1.0\nProcedure Details below\n",
    "tsv.txt": b"Software Version\t3.12.08\nPlate Number\tPlate 1\n",  # no Procedure Details
    "notes.txt": "Experiment notes, 96-well plate, 450 nm.\r\n".encode("utf-16"),
    "blocks.txt": b"##BLOCKS=x\nPlate:\tPlate1\n",
    "data.csv": b"Well,Value\nA1,0.1\nA2,0.2\n",
    "user.csv": b"User: me,Path: C:\\data\nsomething else\n",
    "table.csv": b",1,2,3\nA,1,2,3\nB,4,5,6\n",
    "doc.xml": b"<?xml version='1.0'?><root><Software>Gen5</Software></root>",
    "doc.json": b'{"Software Version": "3.12.08", "Procedure Details": []}',
    "empty.txt": b"",
    "binary.txt": bytes(range(256)) * 4,
}


@pytest.mark.parametrize("name", sorted(GENERIC))
def test_generic_files_are_declined(name):
    for sniff in SNIFFERS.values():
        assert sniff(_src(GENERIC[name])) is None
    got = decode_bytes(name, GENERIC[name])
    assert got is None or got.format not in SNIFFERS


def test_generic_and_other_xlsx_are_declined():
    books = [
        _xlsx({"Sheet1": [["a", "b"], [1, 2]]}),
        _xlsx({"Results": [["Measurement results"], ["not a session"]]}),
        _xlsx({"Absorbance 1_01": [["Measurement results"], ["file.xlsx"]]}),
        _xlsx({"Protocol Information": [["Test Name: x"]]}),
        _xlsx({"General information": [["", "Report generated with SW version", "", "", "Other Software 1.0"]]}),
    ]
    for data in books:
        for sniff in SNIFFERS.values():
            assert sniff(_src(data)) is None


def test_other_fixture_files_in_the_suite_are_declined():
    root = os.path.join(HERE, "fixtures")
    seen = 0
    for name in sorted(os.listdir(root)):
        path = os.path.join(root, name)
        if not os.path.isfile(path):
            continue
        with open(path, "rb") as fh:
            data = fh.read(TEXT_HEAD_BYTES)
        for sniff in SNIFFERS.values():
            assert sniff(_src(data)) is None, (name, sniff.__name__)
        seen += 1
    assert seen >= 20


def test_bids_tsv_and_json_still_decode_normally():
    j = decode_bytes("sub-01_eeg.json", b'{"SamplingFrequency": 500, "EEGChannelCount": 8}')
    assert j is not None and j.format not in SNIFFERS


# ── totality ──────────────────────────────────────────────────────────────────

SEED = 4242


def _seeds():
    gen5 = _bytes("gen5-absorbance-pathlength.txt")
    image = _bytes("gen5-image-cell-count.txt")
    smp = _bytes("softmax-fluorescence-utf16.txt")
    mars = _bytes("mars-fluorescence-384.csv")
    return [
        ("x.txt", b"\r\n\r\nSoftware Version\t3.12.08\r\n"),
        ("x.txt", b"Software Version\t3.1\r\nPlate Number\tP\r\nDate\t1\r\nTime\t2\r\nReader Type:\tX\r\nProcedure Details\r\n\r\n"),
        ("x.txt", b"##BLOCKS= 3\r\nPlate:\t"),
        ("x.csv", b"User: U,Path: P,Test run no.: 1\r\nTest name: T,Date: D,Time: T\r\n\r\nID1: a\r\nAbsorbance\r\n\r\nRaw Data ("),
        ("x.txt", gen5),
        ("x.txt", image),
        ("x.txt", smp),
        ("x.csv", mars),
    ]


def test_sniffers_total_on_random_and_seeded_bytes():
    rng = random.Random(SEED)
    seeds = _seeds()
    for _ in range(1500):
        key, base = rng.choice(seeds)
        buf = bytearray(base[: rng.randint(0, len(base))])
        # Corrupt, insert or append random bytes / tabs / commas / digits.
        for _ in range(rng.randint(0, 12)):
            op = rng.random()
            pos = rng.randint(0, len(buf))
            if op < 0.4:
                buf[pos:pos] = bytes(rng.choice(b"\t,\r\n0123456789:/ .ABCDEFGHnm-") for _ in range(rng.randint(1, 8)))
            elif op < 0.7 and buf:
                buf[min(pos, len(buf) - 1)] = rng.getrandbits(8)
            else:
                buf[pos:pos] = bytes(rng.getrandbits(8) for _ in range(rng.randint(1, 16)))
        data = bytes(buf)
        for sniff in SNIFFERS.values():
            try:
                _valid_or_none(sniff(_src(data)))
            except Exception as e:  # noqa: BLE001
                pytest.fail("%s raised on %d bytes: %r" % (sniff.__name__, len(data), e))
        _valid_or_none(decode_bytes(key, data))


def test_every_prefix_of_text_fixtures_is_safe():
    for c in CASES:
        if c["file"].endswith(".xlsx"):
            continue
        data = _bytes(c["file"])
        for cut in list(range(0, min(len(data), 600), 7)) + list(range(600, len(data), 211)):
            _valid_or_none(decode_bytes(c["file"], data[:cut]))


def test_xlsx_sniffers_total_on_damaged_workbooks():
    rng = random.Random(SEED + 1)
    books = [_bytes(c["file"]) for c in CASES if c["file"].endswith(".xlsx")]
    for _ in range(120):
        data = bytearray(rng.choice(books))
        op = rng.random()
        if op < 0.4:
            data = data[: rng.randint(0, len(data))]
        else:
            for _ in range(rng.randint(1, 20)):
                data[rng.randrange(len(data))] = rng.getrandbits(8)
        for sniff in SNIFFERS.values():
            try:
                _valid_or_none(sniff(_src(bytes(data))))
            except Exception as e:  # noqa: BLE001
                pytest.fail("%s raised on a damaged workbook: %r" % (sniff.__name__, e))
    # Random sheet contents in a well-formed workbook.
    for _ in range(200):
        rows = [[rng.choice(["Well positions", "Measurement results", "x.skax", "Test Name: t", "Measurement type:",
                             "Date of measurement: 2024-01-01/Time of measurement: 1:00:00", "a.mth", "Instrument serial number:",
                             "Wavelength: 9 nm", "Abs", "1", "2", "A", "", "Name", "Serial number", 3.5])
                 for _ in range(rng.randint(0, 4))] for _ in range(rng.randint(0, 12))]
        names = rng.sample(["Absorbance 1_01", "Protocol Information", "Table End point", "General information",
                            "Instrument information", "Magellan Sheet 1", "Session information"], rng.randint(1, 4))
        data = _xlsx({n: rows for n in names})
        for sniff in SNIFFERS.values():
            _valid_or_none(sniff(_src(data)))
