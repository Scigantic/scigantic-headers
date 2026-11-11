"""Flow cytometry, multi-analyte profiling and liquid handler exports: FlowJo,
FACSDiva, Luminex xPONENT and INTELLIFLEX, Bio-Plex Manager, Biomek and Echo.

Golden cases (fixtures/instruments/flow_map_lh-cases.json) are real exports
from allotropy's test data; scigantic-edge runs the same cases."""

import glob
import io
import json
import os
import random
import zipfile

import pytest

from scigantic_headers import decode_bytes
from scigantic_headers.decoders import BytesSource
from scigantic_headers.instruments import TEXT_HEAD_BYTES
from scigantic_headers.instruments import flow_map_lh as fam

HERE = os.path.dirname(__file__)
FIXTURES = os.path.join(HERE, "fixtures", "instruments")
CASES = json.load(open(os.path.join(FIXTURES, "flow_map_lh-cases.json"), encoding="utf-8"))["cases"]

MINE = {fam.FLOWJO, fam.FACSDIVA, fam.INTELLIFLEX, fam.XPONENT, fam.BIOPLEX,
        fam.BIOMEK, fam.ECHO_REFORMAT, fam.ECHO_CHERRY}

SNIFFERS = [
    fam.sniff_flowjo,
    fam.sniff_facsdiva,
    fam.sniff_bioplex,
    fam.sniff_luminex_csv,
    fam.sniff_luminex_xlsx,
    fam.sniff_biomek,
    fam.sniff_echo,
]


def _data(case):
    with open(os.path.join(FIXTURES, "flow_map_lh", case["file"]), "rb") as fh:
        return fh.read()


def test_cases_present():
    assert len(CASES) >= 20
    names = {c["expected"]["format"] for c in CASES}
    assert names == MINE


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_golden(case):
    got = decode_bytes(case["key"], _data(case))
    assert got is not None
    assert got.to_dict() == case["expected"]
    json.dumps(got.to_dict(), allow_nan=False)


# Values allotropy itself writes into the ASM for the same inputs (the .json
# next to each file in allotropy's testdata): software version, model, serial
# number and the number of measurement documents.
ALLOTROPY = {
    "flowjo_example_1.head.wsp": {"softwareVersion": "10.10.0", "instrumentModel": "FACSCelesta"},
    "flowjo_example_3.head.wsp": {"softwareVersion": "10.9.0", "instrumentModel": "ID7000"},
    "bio-rad_bio-plex_manager_example_01.head.xml": {
        "softwareVersion": "6.2.0.175", "serialNumber": "LX12345678912", "firmwareVersion": "2.6.1",
        "plateFormat": 96},
    "luminex_xPONENT_example02.csv": {
        "softwareVersion": "4.3.229.0", "instrumentModel": "FlexMAP 3D", "serialNumber": "FM3DD20137021",
        "plateFormat": 96, "measurementCount": 4, "batchName": "PCD00360_184221_20230517",
        "protocolName": "V116_pnUAD"},
    "luminex_xPONENT_example03.csv": {
        "softwareVersion": "4.3.229.0", "serialNumber": "FM3DD12345678", "measurementCount": 4},
    "luminex_intelliflex_example_01.csv": {
        "softwareVersion": "2.1.1015", "instrumentModel": "Intelliflex", "serialNumber": "SERIAL_12345",
        "measurementCount": 12, "runDate": "2024-10-18T11:36:00"},
    "luminex_intelliflex_single_dataset.csv": {"serialNumber": "IFLEXP24099001", "measurementCount": 3},
    "luminex_intelliflex_v2_1.csv": {
        "instrumentModel": "xMAP INTELLIFLEX\u00ae DR-SE", "serialNumber": "IFLEXP24226001",
        "measurementCount": 2, "protocolName": "HPV14_LIA_cLIA [v5]"},
    "luminex_intelliflex_v2_2.csv": {
        "instrumentModel": "xMAP INTELLIFLEX DR-SE", "serialNumber": "IFLXTEST001", "measurementCount": 10,
        "runDate": "2026-02-25T16:42:08"},
    "luminex_intelliflex_v2_2_excel.xlsx": {
        "instrumentModel": "xMAP INTELLIFLEX\u00ae DR-SE", "serialNumber": "IFLEXP24137001",
        "measurementCount": 2, "runDate": "2026-07-21T15:24:35"},
    "biomek_example_01.csv": {"serialNumber": "B8758519K13", "podHeadSerialNumbers": ["B8758919L19"],
                              "methodName": "!Strain Engineering\\Pfast Feas\\FF2_Induction"},
    "biomek_UnifiedTransfer09-14-2011.log": {"serialNumber": "1234", "measurementCount": 8},
    "echo_reformat_TransferOnly.csv": {
        "softwareVersion": "1.7.2", "instrumentModel": "Echo 525", "serialNumber": "E5XX-19119",
        "instrumentSoftwareVersion": "2.6.3", "transferCount": 4},
    "echo_reformat_SurveyVolumeExamples.csv": {"exceptionCount": 30},
    "echo_reformat_ExceptionsOnly.csv": {"exceptionCount": 4},
}


@pytest.mark.parametrize("name", sorted(ALLOTROPY))
def test_matches_allotropy(name):
    case = next(c for c in CASES if c["name"] == name)
    fields = decode_bytes(case["key"], _data(case)).fields
    for k, v in ALLOTROPY[name].items():
        assert fields.get(k) == v, k


def test_numeric_sample_names_stay_strings():
    case = next(c for c in CASES if c["name"] == "luminex_intelliflex_v2_1.csv")
    names = decode_bytes(case["key"], _data(case)).fields["sampleNames"]
    assert names == ["4074421429"]
    hdr = "INSTRUMENT TYPE,SERIAL NUMBER,SOFTWARE VERSION,PLATE NAME,PLATE START,WELL LOCATION,SAMPLE ID,R12: RP1 MEDIAN\n"
    rows = "".join("INTELLIFLEX,S1,Luminex INTELLIFLEX Bundle - 2.0.1,P,2024-01-02 03:04,A%d,%s,1\n" % (i, s)
                   for i, s in enumerate(["007", "1e3", "12.50", "-0", "NaN"]))
    f = decode_bytes("r.csv", (hdr + rows).encode()).fields
    assert f["sampleNames"] == ["007", "1e3", "12.50", "-0", "NaN"]
    assert f["measurementCount"] == 5 and f["sampleCount"] == 5
    assert f["softwareVersion"] == "2.0.1" and f["runDate"] == "2024-01-02T03:04:00"
    assert f["regionCount"] == 1


# --- per-instrument behavior on constructed inputs --------------------------

FLOWJO = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<!-- saved -->\n'
    '<Workspace version="20.0" flowJoVersion="10.8.1" xmlns:gating="http://www.isac-net.org/std/Gating-ML/v2.0/gating">\n'
    ' <Cytometers><Cytometer name="X" cyt="LSRFortessa &amp; HTS" serialnumber="R123" /></Cytometers>\n'
    ' <SampleList>\n'
    '  <Sample><SampleNode name="a.fcs" /></Sample>\n'
    '  <Sample><SampleNode name="b.fcs" /></Sample>\n'
    ' </SampleList>\n'
    '</Workspace>\n'
)


def test_flowjo_counts_samples_when_list_is_closed():
    f = decode_bytes("w.wsp", FLOWJO.encode()).fields
    assert f["softwareVersion"] == "10.8.1"
    assert f["instrumentModel"] == "LSRFortessa & HTS"
    assert f["serialNumber"] == "R123"
    assert f["sampleCount"] == 2
    # Same file as .xml goes to FlowJo too, not the generic XML decoder.
    assert decode_bytes("w.xml", FLOWJO.encode()).format == fam.FLOWJO
    cut = FLOWJO[: FLOWJO.index("</SampleList>")]
    assert "sampleCount" not in decode_bytes("w.wsp", cut.encode()).fields


def test_flowjo_utf16_and_declines():
    assert decode_bytes("w.wsp", FLOWJO.encode("utf-16")).fields["sampleCount"] == 2
    assert decode_bytes("w.wsp", b'<?xml version="1.0"?><Workspace version="2"/>') is None
    assert decode_bytes("w.wsp", b"<Other flowJoVersion='10'/>") is None


FACSDIVA = (
    "<?xml version='1.0' encoding='UTF-8'?>\n"
    '<bdfacs version="Version 9.0" release_version="Version 9.0">\n'
    '  <experiment name="Exp &lt;1&gt;">\n'
    '    <date>2020-01-02T03:04:05</date>\n'
    '    <export_time>2020-01-02T04:00:00</export_time>\n'
    '    <specimen name="S"><tube name="001"><date>2020-01-02T03:10:00</date>\n'
    '      <data_instrument_serial_number>V0001</data_instrument_serial_number>\n'
    '      <data_instrument_name>FACSymphony</data_instrument_name></tube>\n'
    '      <tube name="002"></tube></specimen>\n'
    '  </experiment>\n'
    '</bdfacs>\n'
)


def test_facsdiva():
    h = decode_bytes("e.xml", FACSDIVA.encode())
    f = h.fields
    assert h.format == fam.FACSDIVA
    assert f["softwareVersion"] == "9.0"
    assert f["experimentName"] == "Exp <1>"
    assert f["experimentDate"] == "2020-01-02T03:04:05"
    assert f["exportDate"] == "2020-01-02T04:00:00"
    assert f["instrumentModel"] == "FACSymphony" and f["serialNumber"] == "V0001"
    assert f["sampleCount"] == 2
    assert "sampleCount" not in decode_bytes("e.xml", FACSDIVA[:-12].encode()).fields


def test_bioplex_counts_wells_when_closed():
    doc = (
        '<?xml version="1.0"?><bpm:BioPlexManager xmlns:bpm="x" BioPlexManagerVersion="6.1" '
        'XmlExportDateTime="2024-01-01T00:00:00Z"><PlateDimensions Columns="24" Rows="16" TotalWells="384" />'
        '<Wells><Well RowNo="1" ColNo="1"><RunSettings><RegionsOfInterest><Region RegionNumber="1" />'
        '<Region RegionNumber="2" /></RegionsOfInterest></RunSettings><MachineInfo>'
        '<MicroControllerVersion>2.0</MicroControllerVersion><SerialNumber>LX1</SerialNumber></MachineInfo></Well>'
        '<Well RowNo="1" ColNo="2"></Well></Wells></bpm:BioPlexManager>'
    )
    f = decode_bytes("b.xml", doc.encode()).fields
    assert f["plateFormat"] == 384 and f["analyteCount"] == 2 and f["measurementCount"] == 2
    assert f["serialNumber"] == "LX1" and f["firmwareVersion"] == "2.0"


def _big_csv(head_lines, row, n):
    body = "".join(row % i for i in range(n))
    return ("".join(head_lines) + body).encode()


def test_xponent_truncated_read_reports_no_counts():
    case = next(c for c in CASES if c["name"] == "luminex_xPONENT_example02.csv")
    text = _data(case).decode("latin-1")
    pre, post = text.split('"Location","Sample"', 1)
    header_row, rest = post.split("\n", 1)
    big = (pre + '"Location","Sample"' + header_row + "\n"
           + "".join('"%d(1,A1)","Unknown%d","1","2","3"\r\n' % (i, i) for i in range(12000)))
    data = big.encode("latin-1")
    assert len(data) > TEXT_HEAD_BYTES
    f = decode_bytes("x.csv", data).fields
    assert f["softwareVersion"] == "4.3.229.0"
    assert "measurementCount" not in f and "sampleNames" not in f
    assert f["analyteCount"] == 2


def test_intelliflex_routing():
    case = next(c for c in CASES if c["name"] == "luminex_xPONENT_example02.csv")
    data = _data(case)
    assert decode_bytes("x.csv", data).format == fam.XPONENT
    swapped = data.replace(b'"FlexMAP 3D"', b'"Intelliflex"', 1)
    assert decode_bytes("x.csv", swapped).format == fam.INTELLIFLEX
    # A flat table that does not say INTELLIFLEX anywhere is not claimed.
    flat = b"INSTRUMENT TYPE,SERIAL NUMBER,SOFTWARE VERSION,WELL LOCATION,SAMPLE ID\nOTHER,1,2,A1,s\n"
    got = decode_bytes("x.csv", flat)
    assert got is None or got.format not in MINE


def _xlsx(rows, sheet="Results", date1904=False, numeric=()):
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    xml_rows = []
    for r, row in enumerate(rows, 1):
        cells = []
        for c, v in enumerate(row):
            ref = ""
            n = c + 1
            while n:
                n, rem = divmod(n - 1, 26)
                ref = chr(65 + rem) + ref
            ref += str(r)
            if (r, c) in numeric:
                cells.append('<c r="%s"><v>%s</v></c>' % (ref, v))
            else:
                cells.append('<c r="%s" t="inlineStr"><is><t>%s</t></is></c>' % (ref, v))
        xml_rows.append('<row r="%d">%s</row>' % (r, "".join(cells)))
    pr = '<workbookPr date1904="1"/>' if date1904 else ""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("xl/workbook.xml", '<workbook %s %s>%s<sheets><sheet name="%s" sheetId="1" r:id="rId1"/></sheets></workbook>' % (ns, rns, pr, sheet))
        z.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr("xl/worksheets/sheet1.xml", "<worksheet %s><sheetData>%s</sheetData></worksheet>" % (ns, "".join(xml_rows)))
    return buf.getvalue()


def test_intelliflex_xlsx_dates_and_numeric_ids():
    hdr = ["MODEL NAME", "INSTRUMENT TYPE", "SERIAL NUMBER", "PLATE START", "PLATE GEOMETRY", "WELL LOCATION", "SAMPLE ID"]
    rows = [hdr, ["xMAP INTELLIFLEX", "INTELLIFLEX", "S9", "45000.5", "Plate384Well", "A1", "0042"],
            ["xMAP INTELLIFLEX", "INTELLIFLEX", "S9", "45000.5", "Plate384Well", "B1", "17"]]
    data = _xlsx(rows, numeric={(2, 3), (3, 3), (3, 6)})
    f = decode_bytes("r.xlsx", data).fields
    assert f["runDate"] == "2023-03-15T12:00:00"
    assert f["sampleNames"] == ["0042", "17"]
    assert f["plateFormat"] == 384 and f["measurementCount"] == 2
    f = decode_bytes("r.xlsx", _xlsx(rows, date1904=True, numeric={(2, 3)})).fields
    assert f["runDate"] == "2027-03-16T12:00:00"
    assert decode_bytes("r.xlsx", _xlsx(rows, sheet="Sheet1")) is None


def test_biomek_layouts_and_declines():
    head = "Method = M\nLogged in user = u\nStarted 01/02/2020 10:00:00\nUnit serial number = None\nPod1 head serial number = None\n\n"
    f = decode_bytes("a.log", (head + "01/02/2020 10:00:01,Pod1,Aspirate,P1,L,,1,5,T\n").encode()).fields
    assert f["logLayout"] == "pipetting" and f["measurementCount"] == 1
    assert "serialNumber" not in f and "podHeadSerialNumbers" not in f and f["podCount"] == 1
    assert decode_bytes("a.log", b"Method = M\nsomething else\n") is None
    assert decode_bytes("a.log", b"2020-01-01 INFO started\n") is None


def test_echo_routing():
    base = "Run ID,1\nRun Date/Time,2024-01-01 01:02:03\nApplication Name,%s\nApplication Version,1.0\n\n[DETAILS]\nA,B\n1,2\n\n"
    assert decode_bytes("e.csv", (base % "Labcyte Echo Plate Reformat").encode()).format == fam.ECHO_REFORMAT
    assert decode_bytes("e.csv", (base % "Echo Cherry Pick").encode()).format == fam.ECHO_CHERRY
    got = decode_bytes("e.csv", (base % "Echo Dose Response").encode())
    assert got is None or got.format not in MINE


def test_dates():
    assert fam._date("2024-02-03 4:05 PM") == "2024-02-03T16:05:00"
    assert fam._date("2024-02-03 12:05:06 AM") == "2024-02-03T00:05:06"
    assert fam._date("2024-02-03 13:05 PM") == "2024-02-03 13:05 PM"
    assert fam._date("2024-13-03 10:00") == "2024-13-03 10:00"
    assert fam._date("5/17/2023 6:42:29 PM") == "5/17/2023 6:42:29 PM"
    assert fam._date("  ") is None


def test_cells():
    assert fam._cells('"a","b,c",,"d""e"') == ["a", "b,c", "", 'd"e']
    assert fam._cells("x") == ["x"]
    assert fam._cells('"open') == ["open"]


# --- declines -----------------------------------------------------------------

GENERIC = {
    "g.csv": b"a,b,c\n1,2,3\n4,5,6\n",
    "g.txt": b"hello world\nline two\n",
    "g.log": b"2024-01-01 12:00:00 INFO service started\n",
    "g.xml": b'<?xml version="1.0"?><root><item a="1"/></root>',
    "g.wsp": b"<html><body>not a workspace</body></html>",
    "g.json": b'{"a": 1, "b": [1, 2]}',
    "g.xlsx": _xlsx([["Instrument", "Reader"], ["Wells", "96"]], sheet="Results"),
    "epu.xml": b'<?xml version="1.0"?><MicroscopeImage xmlns="http://schemas.datacontract.org/2004/07/Fei.SharedObjects"><uniqueID>1</uniqueID></MicroscopeImage>',
    "ome.xml": b'<?xml version="1.0"?><OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"><Image ID="Image:0"/></OME>',
}


@pytest.mark.parametrize("key", sorted(GENERIC))
def test_generic_files_declined(key):
    src = BytesSource(GENERIC[key])
    for sniff in SNIFFERS:
        assert sniff(src) is None


def _foreign_files():
    """Other families' fixtures and this repo's own XML/JSON fixtures."""
    out = []
    for path in glob.glob(os.path.join(FIXTURES, "*", "*")):
        if os.path.basename(os.path.dirname(path)) != "flow_map_lh" and os.path.isfile(path):
            out.append(path)
    for pat in ("*.xml", "*.json", "*.ics"):
        out.extend(glob.glob(os.path.join(HERE, "fixtures", pat)))
    return sorted(out)


def test_foreign_fixtures_declined():
    files = _foreign_files()
    assert files
    for path in files:
        with open(path, "rb") as fh:
            data = fh.read(8 * 1024 * 1024)
        if os.path.basename(path).endswith(".sparse.json"):
            continue
        src = BytesSource(data)
        for sniff in SNIFFERS:
            assert sniff(src) is None, (path, sniff.__name__)


def test_existing_xml_decoders_unchanged():
    for path in glob.glob(os.path.join(HERE, "fixtures", "*.xml")):
        with open(path, "rb") as fh:
            data = fh.read()
        got = decode_bytes(os.path.basename(path), data)
        assert got is None or got.format not in MINE


def test_fcs_not_intercepted():
    data = b"FCS3.1    " + b" " * 48 + b"/$TOT/10/$PAR/2/"
    got = decode_bytes("x.fcs", data)
    assert got is None or got.format not in MINE


# --- totality -----------------------------------------------------------------

def _valid_or_none(h):
    if h is None:
        return
    assert isinstance(h.format, str) and h.format in MINE
    json.loads(json.dumps(h.to_dict(), allow_nan=False))


SEEDS = [b"<Workspace flowJoVersion='1'>", b"<bdfacs version='Version 1'>",
         b"<bpm:BioPlexManager BioPlexManagerVersion='1'>", b'"Program","xPONENT","","x"\n"Build","1"\n',
         b"INSTRUMENT TYPE,WELL LOCATION,SAMPLE ID,SOFTWARE VERSION\nINTELLIFLEX,",
         b"Method = m\nLogged in user = u\nUnit serial number = 1\n", b"Run ID,1\nApplication Name,Echo Cherry Pick\n[DETAILS]\n"]


def test_fuzz_total():
    rng = random.Random(4242)
    corpus = [_data(c) for c in CASES if not c["pythonOnly"]]
    alphabet = b'<>"\',=\n\r:/[] -.0123456789AaPpMm'
    for i in range(3000):
        kind = i % 3
        if kind == 0:
            buf = bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 600)))
        elif kind == 1:
            buf = rng.choice(SEEDS) + bytes(rng.choice(alphabet) for _ in range(rng.randint(0, 400)))
        else:
            base = bytearray(rng.choice(corpus)[: rng.randint(0, 6000)])
            for _ in range(rng.randint(0, 20)):
                if base:
                    base[rng.randrange(len(base))] = rng.choice(alphabet)
            buf = bytes(base)
        for sniff in SNIFFERS[:4] + SNIFFERS[5:]:
            _valid_or_none(sniff(BytesSource(buf)))
        key = rng.choice(["f.wsp", "f.xml", "f.csv", "f.log", "f.xlsx"])
        h = decode_bytes(key, buf)
        if h is not None:
            json.dumps(h.to_dict(), allow_nan=False)


def test_fuzz_xlsx_total():
    rng = random.Random(7)
    case = next(c for c in CASES if c["pythonOnly"])
    good = _data(case)
    for _ in range(200):
        b = bytearray(good)
        for _ in range(rng.randint(1, 30)):
            b[rng.randrange(len(b))] = rng.getrandbits(8)
        _valid_or_none(fam.sniff_luminex_xlsx(BytesSource(bytes(b))))
        _valid_or_none(fam.sniff_luminex_xlsx(BytesSource(bytes(b[: rng.randint(0, len(b))]))))
