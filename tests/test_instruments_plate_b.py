"""Plate readers, family B: Envision, Kaleido, MSD Methodical Mind, MSD
Discovery Workbench, CTL ImmunoSpot, Mabtech Apex, Unchained Labs Lunatic &
Stunner. Golden cases are real allotropy test exports (see
fixtures/instruments/plate_b/NOTICE.md); their values were checked against the
ASM allotropy writes for the same files."""

import io
import json
import os
import random
import zipfile

import pytest

from scigantic_headers import decode_bytes
from scigantic_headers.decoders import BytesSource, DecodedHeader
from scigantic_headers.instruments import TEXT_HEAD_BYTES, plate_b

HERE = os.path.dirname(__file__)
FIX = os.path.join(HERE, "fixtures", "instruments")
CASES = json.load(open(os.path.join(FIX, "plate_b-cases.json"), encoding="utf-8"))["cases"]


def _bytes(rel):
    with open(os.path.join(FIX, rel), "rb") as fh:
        return fh.read()


def _src(data):
    return BytesSource(data, budget=16 * 1024 * 1024)


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_golden(case):
    got = decode_bytes(case["key"], _bytes(case["file"]))
    assert got is not None
    assert got.to_dict() == case["expected"]


SNIFFERS = {
    "envision": (plate_b.sniff_envision, plate_b.ENVISION),
    "kaleido": (plate_b.sniff_kaleido, plate_b.KALEIDO),
    "methodical_mind": (plate_b.sniff_methodical_mind, plate_b.METHODICAL_MIND),
    "workbench_csv": (plate_b.sniff_workbench_csv, plate_b.WORKBENCH),
    "workbench_xlsx": (plate_b.sniff_workbench_xlsx, plate_b.WORKBENCH),
    "immunospot": (plate_b.sniff_immunospot, plate_b.IMMUNOSPOT),
    "apex": (plate_b.sniff_apex, plate_b.APEX),
    "lunatic_csv": (plate_b.sniff_lunatic_csv, plate_b.LUNATIC),
    "lunatic_xlsx": (plate_b.sniff_lunatic_xlsx, plate_b.LUNATIC),
}


def _xlsx(sheets):
    """A minimal .xlsx: {sheet name: rows}. Strings go through sharedStrings."""
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    shared = []
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        entries, rels = [], []
        for i, (name, rows) in enumerate(sheets.items(), 1):
            entries.append('<sheet name="%s" sheetId="%d" r:id="rId%d"/>' % (name, i, i))
            rels.append('<Relationship Id="rId%d" Target="worksheets/sheet%d.xml"/>' % (i, i))
            xml_rows = []
            for r, row in enumerate(rows, 1):
                cells = []
                for c, v in enumerate(row):
                    ref = chr(65 + c) + str(r)
                    if isinstance(v, str):
                        shared.append(v)
                        cells.append('<c r="%s" t="s"><v>%d</v></c>' % (ref, len(shared) - 1))
                    else:
                        cells.append('<c r="%s"><v>%s</v></c>' % (ref, v))
                xml_rows.append('<row r="%d">%s</row>' % (r, "".join(cells)))
            z.writestr("xl/worksheets/sheet%d.xml" % i,
                       "<worksheet %s><sheetData>%s</sheetData></worksheet>" % (ns, "".join(xml_rows)))
        z.writestr("xl/workbook.xml", "<workbook %s %s><sheets>%s</sheets></workbook>" % (ns, rns, "".join(entries)))
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">%s</Relationships>' % "".join(rels))
        z.writestr("xl/sharedStrings.xml", "<sst %s>%s</sst>" % (ns, "".join("<si><t>%s</t></si>" % s for s in shared)))
    return buf.getvalue()


GENERIC = {
    "generic.csv": b"id,name,value\n1,alpha,0.5\n2,beta,0.7\n",
    "generic-plate.csv": b",1,2,3\nA,1,2,3\nB,4,5,6\n",
    "generic.txt": b"Plate information\nThis is a note about a plate.\nSoftware version: 1.0\n",
    "generic.tsv.txt": b"Sample\tAssay\tWell\tSpot\nS1\tIL6\tA01\t1\n",
    "generic.xml": b'<?xml version="1.0"?><root><Plate id="1"/></root>',
    "generic.json": b'{"Plate ID": "P1", "Sample name": "x", "A260": 1.0}',
    "generic.xlsx": _xlsx({"Sheet1": [["Plate ID", "Sample name", "Value"], ["P1", "s", 1]]}),
    "generic-two.xlsx": _xlsx({"Data": [["a", "b"]], "Plate Info": [["Name:", "x"]]}),
    "empty": b"",
    "binary": bytes(range(256)) * 8,
}


def _all_inputs():
    """(label, owner sniffer names, bytes) for every fixture and generic input."""
    owners = {}
    for c in CASES:
        fmt = c["expected"]["format"]
        ext = c["key"].rsplit(".", 1)[-1]
        owners[c["file"]] = {k for k, (_, name) in SNIFFERS.items()
                             if name == fmt and (k.endswith(ext) or not k.endswith(("csv", "xlsx")))}
    out = [(rel, own, _bytes(rel)) for rel, own in owners.items()]
    other = os.path.join(FIX, "plate_b", "other")
    for n in sorted(os.listdir(other)):
        out.append(("plate_b/other/" + n, set(), _bytes("plate_b/other/" + n)))
    for n, data in GENERIC.items():
        out.append((n, set(), data))
    return out


INPUTS = _all_inputs()


@pytest.mark.parametrize("sniffer", sorted(SNIFFERS))
def test_sniffer_declines_what_it_does_not_own(sniffer):
    fn, name = SNIFFERS[sniffer]
    for label, owners, data in INPUTS:
        got = fn(_src(data))
        if sniffer in owners:
            assert got is not None and got.format == name, label
        else:
            assert got is None, (sniffer, label)


def test_msd_products_are_told_apart():
    mm = _bytes("plate_b/mm-two-plates.txt")
    wb = _bytes("plate_b/workbench-1.csv")
    # The same bytes under the other product's extension still name the right
    # product, or nothing.
    assert decode_bytes("x.txt", mm).format == plate_b.METHODICAL_MIND
    assert decode_bytes("x.csv", wb).format == plate_b.WORKBENCH
    assert plate_b.sniff_workbench_csv(_src(mm)) is None
    assert plate_b.sniff_methodical_mind(_src(wb)) is None
    got = decode_bytes("x.csv", mm)
    assert got is None or got.format not in (plate_b.METHODICAL_MIND, plate_b.WORKBENCH)
    assert decode_bytes("x.xlsx", _bytes("plate_b/workbench.xlsx")).format == plate_b.WORKBENCH


def test_counts_need_the_whole_file():
    data = _bytes("plate_b/envision-luminescence01.csv")
    assert decode_bytes("x.csv", data).fields["plateCount"] == 2
    padded = data + b"\r\n" * (TEXT_HEAD_BYTES // 2)
    got = decode_bytes("x.csv", padded)
    assert got.format == plate_b.ENVISION and "plateCount" not in got.fields
    wb = _bytes("plate_b/workbench-2.csv") + b"S999,Test Assay 1.0,H12,1\r\n" * 20000
    got = decode_bytes("x.csv", wb)
    assert got.format == plate_b.WORKBENCH
    assert "sampleCount" not in got.fields and "measurementCount" not in got.fields


def test_dates():
    d = plate_b._date
    assert d("2020-10-26T16:16:18.2284559-04:00") == "2020-10-26T16:16:18.2284559-04:00"
    assert d("2024-07-12 10:38:09") == "2024-07-12T10:38:09"
    assert d("Tue Oct 03 14:44:50 2023") == "2023-10-03T14:44:50"
    # Numeric day/month dates are kept as written: 01/08/2024 is January 8 in
    # one Lunatic export and could be August 1 in another.
    assert d("01/08/2024 10:15:59") == "01/08/2024 10:15:59"
    assert d("Xyz Oct 99 14:44:50 2023") == "Xyz Oct 99 14:44:50 2023"
    assert d("  ") is None


def test_xlsx_serial_dates():
    xl = _xlsx({"Plate Info": [["Plate Read Date:", 44105.46609953704], ["Machine ID:", "IRIS 1234"]],
                "Plate Database": [["Plate", "Well"]]})
    got = plate_b.sniff_apex(_src(xl))
    assert got.fields["runDate"] == "2020-10-01T11:11:11"
    assert got.fields["instrumentModel"] == "IRIS"
    assert got.fields["machineId"] == "IRIS 1234"


def test_synthetic_workbench_xlsx_with_empty_first_row_cells():
    xl = _xlsx({"MSD Data": [["x"]],
                "Workbench data": [["Plate_ABC123", "", ""],
                                   ["Sample", "Assay", "Well", "Spot", "Signal"],
                                   ["S1", "IL-6", "A01", 1, 100], ["S1", "IL-6", "A02", 1, 110],
                                   ["S2", "IL-6", "B01", 1, 90]]})
    got = decode_bytes("x.xlsx", xl)
    assert got.format == plate_b.WORKBENCH
    assert got.fields["sampleCount"] == 2 and got.fields["measurementCount"] == 3
    assert got.fields["assays"] == ["IL-6"]


# ---- totality ------------------------------------------------------------

SEED = 1729


def _json_safe(result):
    if result is None:
        return
    assert isinstance(result, DecodedHeader)
    json.loads(json.dumps(result.to_dict(), allow_nan=False))


def _mutations(rng, data):
    yield data[: rng.randint(0, len(data))]
    b = bytearray(data)
    for _ in range(rng.randint(1, 40)):
        if b:
            b[rng.randrange(len(b))] = rng.getrandbits(8)
    yield bytes(b)
    lines = data.split(b"\n")
    rng.shuffle(lines)
    yield b"\n".join(lines[: rng.randint(1, len(lines))])
    yield data.replace(b",", b"\t") if rng.random() < 0.5 else data.replace(b"\t", b",")
    yield data[: rng.randint(0, 400)] + bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 400)))


def test_sniffers_are_total_over_mutated_fixtures():
    rng = random.Random(SEED)
    text_inputs = [d for label, _, d in INPUTS if not d.startswith(b"PK")]
    for _ in range(40):
        for data in text_inputs:
            for m in _mutations(rng, data):
                for fn, _ in SNIFFERS.values():
                    _json_safe(fn(_src(m)))


def test_sniffers_are_total_over_random_bytes():
    rng = random.Random(SEED + 1)
    for _ in range(1500):
        buf = bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 2048)))
        for fn, _ in SNIFFERS.values():
            _json_safe(fn(_src(buf)))


def test_xlsx_sniffers_are_total_over_damaged_zips():
    rng = random.Random(SEED + 2)
    xlsx = [d for _, _, d in INPUTS if d.startswith(b"PK")]
    for _ in range(25):
        for data in xlsx:
            b = bytearray(data)
            for _ in range(rng.randint(1, 30)):
                b[rng.randrange(len(b))] = rng.getrandbits(8)
            for m in (bytes(b), data[: rng.randint(0, len(data))]):
                for fn, _ in SNIFFERS.values():
                    _json_safe(fn(_src(m)))
