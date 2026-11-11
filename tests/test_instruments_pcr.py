"""qPCR / dPCR instrument exports: QuantStudio Real-Time PCR Software and
Design & Analysis, Absolute Q, Qiacuity, Diomni and CFX Maestro.

Golden inputs are allotropy's own test files (tests/fixtures/instruments/pcr,
see NOTICE.md). Each case in pcr-cases.json also records what allotropy's ASM
output says about the same file, and the decode is checked against it."""

import glob
import io
import json
import os
import random
import zipfile

import pytest

from scigantic_headers import DecodedHeader, decode_bytes
from scigantic_headers.decoders import BytesSource
from scigantic_headers.instruments import TEXT_HEAD_BYTES, pcr

HERE = os.path.dirname(__file__)
FIX = os.path.join(HERE, "fixtures", "instruments")
CASES = json.load(open(os.path.join(FIX, "pcr-cases.json"), encoding="utf-8"))["cases"]
FORMATS = {pcr.QS_RTPCR, pcr.QS_DA, pcr.ABSOLUTE_Q, pcr.QIACUITY, pcr.DIOMNI, pcr.CFX}

SNIFFERS = {
    "txt": [pcr.sniff_quantstudio_txt],
    "xlsx": [pcr.sniff_quantstudio_xlsx],
    "csv": [pcr.sniff_qiacuity_csv, pcr.sniff_cfx_csv, pcr.sniff_absolute_q_csv],
    "zip": [pcr.sniff_absolute_q_zip],
}
ALL_SNIFFERS = [s for v in SNIFFERS.values() for s in v]


def _read(rel):
    with open(os.path.join(FIX, rel), "rb") as fh:
        return fh.read()


def _norm(d):
    return json.loads(json.dumps(d))


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_golden_case(case):
    got = decode_bytes(case["key"], _read(case["file"]))
    assert got is not None
    assert _norm(got.to_dict()) == case["expected"]


# Where the decode and allotropy's ASM legitimately differ, and why.
ASM_DIFFERENCES = {
    # The RT-PCR parser writes a fixed "1.0"; the file states no version.
    ("rtpcr-*", "softwareVersion"),
    # allotropy names genotyping targets "<SNP assay>-<allele>"; the file has
    # the SNP assay, reported as snpAssays.
    ("rtpcr-txt-genotyping", "targets"),
    # The multichannel export lists channel combinations ("VIC JUN"); the
    # decode reports the single dyes.
    ("absoluteq-multichannel", "reporterDyes"),
    # The multichannel export has no well column; allotropy counts samples.
    ("absoluteq-multichannel", "wells"),
    # A constant allotropy writes for every Qiacuity file.
    ("qiacuity-*", "deviceIdentifier"),
}


def _expected_diff(name, key):
    for pat, k in ASM_DIFFERENCES:
        if k == key and (pat == name or (pat.endswith("*") and name.startswith(pat[:-1]))):
            return True
    return False


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_agrees_with_allotropy_asm(case):
    f = case["expected"]["fields"]
    asm = case["allotropyAsm"]
    mine = {
        "softwareVersion": f.get("softwareVersion"),
        "instrumentModel": f.get("instrumentModel"),
        "serialNumber": f.get("serialNumber"),
        "deviceIdentifier": f.get("instrumentName"),
        "targets": sorted(f["targets"]) if "targets" in f else None,
        "reporterDyes": sorted(f["reporterDyes"]) if "reporterDyes" in f else None,
        "wells": f.get("sampleCount"),
    }
    checked = 0
    for k, v in asm.items():
        if _expected_diff(case["name"], k):
            continue
        assert mine[k] == v, (k, mine[k], v)
        checked += 1
    assert checked >= 1
    # Nothing allotropy states is missing from the decode, bar the known gaps.
    assert f["allotropyParser"] in {"appbio_quantstudio", "appbio_quantstudio_designandanalysis",
                                    "appbio_absolute_q", "qiacuity_dpcr", "thermo_fisher_diomni",
                                    "cfxmaestro"}


def test_case_formats_cover_every_instrument():
    assert {c["expected"]["format"] for c in CASES} == FORMATS


# ── telling QuantStudio exports apart ─────────────────────────────────────────


def _xlsx(sheets):
    """A minimal .xlsx: sheets is [(name, rows)], rows lists of str/int."""
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    shared = []
    parts = []
    for name, rows in sheets:
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
            if cells:
                xml_rows.append('<row r="%d">%s</row>' % (r, "".join(cells)))
        parts.append("<worksheet %s><sheetData>%s</sheetData></worksheet>" % (ns, "".join(xml_rows)))
    buf = io.BytesIO()
    esc = lambda s: s.replace("&", "&amp;").replace("<", "&lt;")  # noqa: E731
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("xl/workbook.xml", "<workbook %s %s><sheets>%s</sheets></workbook>" % (
            ns, rns, "".join('<sheet name="%s" sheetId="%d" r:id="rId%d"/>' % (n, i + 1, i + 1)
                             for i, (n, _) in enumerate(sheets))))
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">%s</Relationships>'
                   % "".join('<Relationship Id="rId%d" Target="worksheets/sheet%d.xml"/>' % (i + 1, i + 1)
                             for i in range(len(sheets))))
        z.writestr("xl/sharedStrings.xml", "<sst %s>%s</sst>" % (ns, "".join("<si><t>%s</t></si>" % esc(s) for s in shared)))
        for i, p in enumerate(parts):
            z.writestr("xl/worksheets/sheet%d.xml" % (i + 1), p)
    return buf.getvalue()


LEGACY_HEADER = [
    ["Block Type", "384-Well Block"],
    ["Chemistry", "TAQMAN"],
    ["Experiment Name", "run 7"],
    ["Experiment Run End Time", "2024-05-20 19:39:52 PM PDT"],
    ["Experiment Type", "Comparative C\u0442 (\u0394\u0394C\u0442)"],
    ["Instrument Name", "QS5-01"],
    ["Instrument Serial Number", 1234],
    ["Quantification Cycle Method", "Ct"],
]
WELL_TABLE = [
    [None],
    ["Well", "Well Position", "Omit", "Sample Name", "Target Name", "Task", "Reporter"],
    [1, "A1", "false", "S1", "GAPDH", "UNKNOWN", "VIC"],
    [1, "A1", "false", "S1", "TGFB", "UNKNOWN", "FAM"],
    [2, "A2", "false", "S2", "GAPDH", "UNKNOWN", "VIC"],
]


def _legacy(model, first_sheet="Results"):
    rows = LEGACY_HEADER + [["Instrument Type", model]] + WELL_TABLE
    return _xlsx([(first_sheet, rows)])


def test_legacy_layout_from_design_and_analysis_v1():
    # Design & Analysis 1.x writes the Real-Time PCR layout but a trademark
    # sign in the instrument type (QuantStudio 3/5 run on Design & Analysis).
    h = decode_bytes("x.xlsx", _legacy("QuantStudio\u2122 5 System"))
    assert h.format == pcr.QS_DA
    f = h.fields
    assert f["allotropyParser"] == "appbio_quantstudio_designandanalysis"
    assert f["allotropyAlternateParser"] == "appbio_quantstudio"
    assert f["exportFormat"] == "xlsx (Real-Time PCR layout)"
    assert f["plateFormat"] == 384 and f["serialNumber"] == "1234"
    assert f["targets"] == ["GAPDH", "TGFB"] and f["reporterDyes"] == ["VIC", "FAM"]
    assert f["sampleCount"] == 2 and f["measurementCount"] == 3 and f["sampleNameCount"] == 2
    assert "softwareVersion" not in f


def test_legacy_layout_from_realtime_pcr_software():
    h = decode_bytes("x.xlsx", _legacy("QuantStudio(TM) 6 Flex System", "Sample Setup"))
    assert h.format == pcr.QS_RTPCR
    assert h.fields["allotropyParser"] == "appbio_quantstudio"
    assert h.fields["allotropyAlternateParser"] == "appbio_quantstudio_designandanalysis"
    assert h.fields["instrumentModel"] == "QuantStudio(TM) 6 Flex System"


def test_design_and_analysis_names_itself():
    rows = [
        ["File Name", "C:\\runs\\r.eds"],
        ["Instrument Type", "QuantStudio\u2122 7 Pro System"],
        ["Block Type", "96-Well 0.2-mL Block"],
        ["Instrument Serial Number", "1.23456789E8"],
        ["Run End Date/Time", "2024-10-17 10:29:52 AM EDT"],
        ["Quantification Cycle Method", "CT"],
        ["Software Name and Version", "Design & Analysis Software v2.8.1"],
        ["Plugin Name and Version", "Primary Analysis v1.9.0"],
        [None],
        ["Well", "Well Position", "Omit", "Sample", "Target", "Task", "Reporter"],
        [1, "A1", "false", "S", "T1", "UNKNOWN", "FAM"],
    ]
    h = decode_bytes("x.xlsx", _xlsx([("Results", rows), ("Amplification Data", [["Well"]])]))
    assert h.format == pcr.QS_DA
    f = h.fields
    assert f["softwareVersion"] == "2.8.1"
    assert f["serialNumber"] == "123456789"
    assert f["plateFormat"] == 96
    assert f["analysisPlugins"] == "Primary Analysis v1.9.0"
    assert f["runDate"] == "2024-10-17 10:29:52 AM EDT"
    assert f["exportFormat"] == "xlsx"
    assert "allotropyAlternateParser" not in f


def test_diomni_without_target_call_sheet_points_at_design_and_analysis():
    rows = [
        ["Instrument Type", "QuantStudio\u2122 7 Pro System"],
        ["Block Type", "384-Well Block"],
        ["Run End Data/Time", "2026-02-13 02:39:53 PM EST"],
        ["Software name", "Diomni"],
        ["Software version", "4.4.0"],
        [None],
        ["Well", "Omit", "Sample", "Target", "Task", "Reporter"],
        ["A1", "false", "S", "T", "UNKNOWN", "FAM"],
    ]
    h = decode_bytes("x.xlsx", _xlsx([("Primary_result", rows), ("Raw Data", [["Well"]])]))
    assert h.format == pcr.DIOMNI
    assert h.fields["softwareVersion"] == "4.4.0"
    assert h.fields["allotropyParser"] == "appbio_quantstudio_designandanalysis"


def test_xlsx_that_only_looks_tabular_is_declined():
    assert decode_bytes("x.xlsx", _xlsx([("Results", [["Block Type", "96"], ["Well", "x"]])])) is None
    assert decode_bytes("x.xlsx", _xlsx([("Sheet1", [["a", "b"], [1, 2]])])) is None
    # Real-Time PCR keys but an instrument that is not a QuantStudio.
    rows = [r for r in LEGACY_HEADER] + [["Instrument Type", "LightCycler 480"]] + WELL_TABLE
    assert decode_bytes("x.xlsx", _xlsx([("Results", rows)])) is None


def test_quantstudio_txt_in_utf16():
    raw = _read("pcr/appbio_quantstudio_example05.txt")
    want = decode_bytes("x.txt", raw).to_dict()
    for enc in ("utf-16", "utf-16-le"):
        got = decode_bytes("x.txt", raw.decode("utf-8").encode(enc))
        assert got is not None and got.to_dict() == want


def test_counts_withheld_when_the_well_table_is_cut():
    # A file longer than the leading read: header fields stay, counts that
    # would be partial are left out.
    raw = _read("pcr/appbio_quantstudio_example01.txt").decode("utf-8")
    top, rest = raw.split("[Sample Setup]", 1)
    lines = rest.split("\n")[1:]
    setup_cols = lines[0]
    body = []
    for line in lines[1:]:
        if not line.strip() or line.startswith("["):
            break
        body.append(line)
    setup_rows = "\n".join(body)
    filler = "\n".join([setup_rows] * (TEXT_HEAD_BYTES // max(1, len(setup_rows)) + 2))
    big = (top + "[Sample Setup]\n" + setup_cols + "\n" + filler).encode("utf-8")
    assert len(big) > TEXT_HEAD_BYTES
    h = pcr.sniff_quantstudio_txt(BytesSource(big, budget=1 << 26))
    assert h is not None
    assert h.fields["serialNumber"] == "INST_1"
    assert "sampleCount" not in h.fields and "targets" not in h.fields

    cfx = _read("pcr/cfxmaestro_example01.csv").decode("utf-8")
    head_line, body = cfx.split("\n", 1)
    big = (head_line + "\n" + body * (TEXT_HEAD_BYTES // len(body) + 2)).encode("utf-8")
    h = pcr.sniff_cfx_csv(BytesSource(big, budget=1 << 26))
    assert h is not None and "sampleCount" not in h.fields and h.fields["reporterDyes"] == ["FAM", "VIC"]


# ── declining everything else ─────────────────────────────────────────────────

GENERIC = {
    "txt": [b"", b"hello world\n", b"* a list item\n* another\n", b"* Key = Value\n* Other = 2\n\n[Section]\na\tb\n",
            b"* Block Type = 96-Well Block\n* Instrument Type = Some Reader\n* Chemistry = X\n* Experiment Name = n\n\n[Results]\nWell\n"],
    "csv": [b"", b"a,b,c\n1,2,3\n", b"Run,Date,Instrument,Plate\n1,2,3,4\n", b"sep=,\na,b\n1,2\n",
            b",Well,Fluor,Target\n,A01,FAM,x\n", b",,,FAM\nRun name,Group,Sample,Well\n"],
    "zip": [b"", b"PK\x03\x04garbage"],
    "xlsx": [b"", b"PK\x03\x04garbage", b"not a zip"],
}


def _zip(members):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for n, d in members.items():
            z.writestr(n, d)
    return buf.getvalue()


def test_generic_files_are_declined():
    GENERIC["zip"].append(_zip({"a_summary.csv": b"a,b\n1,2\n"}))
    GENERIC["zip"].append(_zip({"notes.txt": b"x"}))
    for ext, blobs in GENERIC.items():
        for b in blobs:
            for s in SNIFFERS[ext]:
                assert s(BytesSource(b, budget=1 << 20)) is None, (s.__name__, b[:40])
    for key, b in [("x.json", b'{"a": 1}'), ("x.xml", b"<a><b/></a>"), ("x.tsv", b"a\tb\n1\t2\n")]:
        h = decode_bytes(key, b)
        assert h is None or h.format not in FORMATS


def _accepting(data):
    return [s.__name__ for s in ALL_SNIFFERS if s(BytesSource(data, budget=1 << 26)) is not None]


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_exactly_one_sniffer_claims_each_fixture(case):
    assert len(_accepting(_read(case["file"]))) == 1


def _other_family_fixtures():
    out = []
    for d in sorted(glob.glob(os.path.join(FIX, "*"))):
        if os.path.isdir(d) and os.path.basename(d) != "pcr":
            out += [p for p in glob.glob(os.path.join(d, "**", "*"), recursive=True)
                    if os.path.isfile(p) and not p.endswith(".md")]
    return out


def test_other_families_fixtures_are_declined():
    # Other instrument families keep their fixtures beside these; none may be
    # claimed as a PCR export.
    for p in _other_family_fixtures():
        with open(p, "rb") as fh:
            data = fh.read()
        assert _accepting(data) == [], p
        h = decode_bytes(os.path.basename(p), data)
        assert h is None or h.format not in FORMATS, p


# ── totality ──────────────────────────────────────────────────────────────────

SEED = 1729


def _valid_or_none(h):
    if h is None:
        return
    assert isinstance(h, DecodedHeader) and h.format in FORMATS
    json.loads(json.dumps(h.to_dict(), allow_nan=False))


def _mutate(rng, data):
    b = bytearray(data)
    op = rng.randrange(6)
    if op == 0 and b:
        b = b[: rng.randrange(len(b))]
    elif op == 1:
        for _ in range(rng.randint(1, 20)):
            if b:
                b[rng.randrange(len(b))] = rng.getrandbits(8)
    elif op == 2:
        lines = bytes(b).split(b"\n")
        if lines:
            del lines[rng.randrange(len(lines))]
        b = bytearray(b"\n".join(lines))
    elif op == 3:
        lines = bytes(b).split(b"\n")
        if lines:
            i = rng.randrange(len(lines))
            lines.insert(i, bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 80))))
        b = bytearray(b"\n".join(lines))
    elif op == 4:
        b = bytearray(bytes(b).replace(b",", b",,", rng.randint(1, 5)))
    else:
        b = bytearray(bytes(b).replace(b"\t", b"", rng.randint(1, 50)))
    return bytes(b)


def test_sniffers_total_over_random_bytes():
    rng = random.Random(SEED)
    for _ in range(1500):
        buf = bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 2048)))
        for prefix in (b"", b"* Instrument Type = QuantStudio\n", b"sep=,\n", b",Well,Fluor,", b"Run,Date,Instrument,Plate,",
                       b"PK\x03\x04"):
            data = prefix + buf
            for s in ALL_SNIFFERS:
                try:
                    _valid_or_none(s(BytesSource(data, budget=1 << 20)))
                except Exception as e:  # noqa: BLE001
                    pytest.fail("%s raised on random input: %r" % (s.__name__, e))


@pytest.mark.parametrize("case", [c for c in CASES if c["go"]], ids=[c["name"] for c in CASES if c["go"]])
def test_sniffers_total_over_mutated_text_fixtures(case):
    rng = random.Random(SEED + len(case["name"]))
    data = _read(case["file"])
    for _ in range(300):
        m = _mutate(rng, data)
        for s in ALL_SNIFFERS:
            try:
                _valid_or_none(s(BytesSource(m, budget=1 << 20)))
            except Exception as e:  # noqa: BLE001
                pytest.fail("%s raised on a mutated %s: %r" % (s.__name__, case["name"], e))


def test_sniffers_total_over_mutated_xlsx():
    rng = random.Random(SEED)
    small = [c for c in CASES if c["file"].endswith((".xlsx", ".zip")) and len(_read(c["file"])) < 60_000]
    for c in small:
        data = _read(c["file"])
        for _ in range(25):
            m = _mutate(rng, data)
            for s in SNIFFERS["xlsx"] + SNIFFERS["zip"]:
                try:
                    _valid_or_none(s(BytesSource(m, budget=1 << 26)))
                except Exception as e:  # noqa: BLE001
                    pytest.fail("%s raised on a mutated %s: %r" % (s.__name__, c["name"], e))
    # Valid zips whose sheets are noise.
    for _ in range(40):
        rows = [[rng.choice(["Block Type", "Instrument Type", "Well", "Software name", "", "x"]),
                 rng.choice(["Diomni", "QuantStudio\u2122 5", "96", "", "1e400", "nan"])] for _ in range(rng.randint(0, 30))]
        name = rng.choice(["Results", "Primary_result", "Target Call", "Sample Setup"])
        _valid_or_none(pcr.sniff_quantstudio_xlsx(BytesSource(_xlsx([(name, rows)]), budget=1 << 26)))


@pytest.mark.parametrize("line,want", [
    ("", []), ("a,b", ["a", "b"]), (",a,", ["", "a", ""]), ('"a,b",c', ["a,b", "c"]),
    ('"a""b"', ['a"b']), ('"ab"c,d', ["abc", "d"]), ('a"b,c', ['a"b', "c"]), ('"unterminated', ["unterminated"]),
])
def test_csv_line_matches_the_csv_module(line, want):
    # Same cases as the Go port's TestPCRCSVLine; each also matches what the
    # csv module returns for a complete line.
    import csv
    assert pcr._csv_line(line) == want
    if line and line.count('"') % 2 == 0:
        assert next(csv.reader([line])) == want
