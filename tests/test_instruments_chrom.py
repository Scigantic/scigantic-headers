"""Chromatography, electrophoresis and SPR exports (instruments/chrom.py):
golden cases from allotropy's test data, synthetic containers for the
branches those files do not reach, declines on everything else, and fuzzing."""

import io
import json
import os
import random
import struct
import zipfile

import pytest

from conftest import FIXTURES_DIR
from scigantic_headers import decode_bytes
from scigantic_headers.decoders import INSTRUMENT_SEEK_BYTES, BytesSource, DecodedHeader, decode_source
from scigantic_headers.instruments import INSTRUMENTS
from scigantic_headers.instruments import chrom
from scigantic_headers.instruments.chrom import partial_json

INST_DIR = os.path.join(FIXTURES_DIR, "instruments")
CASES = json.load(open(os.path.join(INST_DIR, "chrom-cases.json")))["cases"]

SNIFFERS = {
    chrom.OPENLAB: chrom.sniff_openlab,
    chrom.UNICORN: chrom.sniff_unicorn,
    chrom.EMPOWER: chrom.sniff_empower,
    chrom.CHROMELEON: chrom.sniff_chromeleon,
    chrom.TAPESTATION: chrom.sniff_tapestation,
    chrom.T200_CONTROL: chrom.sniff_t200_control,
    chrom.T200_EVALUATION: chrom.sniff_t200_evaluation,
    chrom.INSIGHT: chrom.sniff_insight,
}


def _source(path: str) -> BytesSource:
    if path.endswith(".sparse.json"):
        rec = json.load(open(path))
        return BytesSource(size=rec["size"], budget=INSTRUMENT_SEEK_BYTES,
                           ranges={int(o): bytes.fromhex(h) for o, h in rec["ranges"].items()})
    return BytesSource(open(path, "rb").read(), budget=INSTRUMENT_SEEK_BYTES)


def _src(data: bytes) -> BytesSource:
    return BytesSource(data, budget=INSTRUMENT_SEEK_BYTES)


def _valid_or_none(h):
    if h is None:
        return
    assert isinstance(h, DecodedHeader) and h.format and isinstance(h.summary, str)
    json.loads(json.dumps(h.to_dict(), allow_nan=False))


# ── identity table ───────────────────────────────────────────────────────────


def test_identities_match_the_instrument_table():
    want = {
        chrom.OPENLAB: ("Liquid Chromatography", ["Absorbance", "Fluorescence"], "agilent_openlab_cds"),
        chrom.UNICORN: ("Liquid Chromatography", ["Absorbance"], "cytiva_unicorn"),
        chrom.EMPOWER: ("Liquid Chromatography", ["Absorbance"], "benchling_empower"),
        chrom.CHROMELEON: ("Liquid Chromatography", ["Absorbance", "Fluorescence", "Conductivity"],
                           "benchling_chromeleon"),
        chrom.TAPESTATION: ("Electrophoresis", ["Fluorescence"], "agilent_tapestation_analysis"),
        chrom.T200_CONTROL: ("Binding Affinity Analyzer", ["Surface Plasmon Resonance"],
                             "cytiva_biacore_t200_control"),
        chrom.T200_EVALUATION: ("Binding Affinity Analyzer", ["Surface Plasmon Resonance"],
                                "cytiva_biacore_t200_evaluation"),
        chrom.INSIGHT: ("Binding Affinity Analyzer", ["Surface Plasmon Resonance"], "cytiva_biacore_insight"),
    }
    for name, (cat, modes, parser) in want.items():
        got = INSTRUMENTS[name]
        assert (got[0], got[3], got[2]) == (cat, modes, parser)


# ── golden cases ─────────────────────────────────────────────────────────────


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_golden_case(case):
    src = _source(os.path.join(INST_DIR, case["file"]))
    got = decode_source(case["key"], src)
    assert got is not None
    assert got.to_dict() == case["expected"]


def test_golden_cases_cover_every_instrument_but_t200_evaluation():
    # The T200 Evaluation files in allotropy's tests are 1-231 MB; that
    # instrument is covered by the synthetic compound file below.
    covered = {c["expected"]["format"] for c in CASES}
    assert covered == set(SNIFFERS) - {chrom.T200_EVALUATION}


def test_golden_values_read_from_the_files():
    by = {c["name"]: c["expected"]["fields"] for c in CASES}
    ol = by["openlab-sirius.rslt"]
    assert ol["softwareVersion"] == "2.5 - Build  2.5.0-REL_2.5.0.842"  # allotropy: same
    assert ol["instrumentName"] == "Sirius" and ol["injectionCount"] == 2
    assert ol["wavelengthsNm"] == [210, 228]
    assert ol["modules"][0] == {"type": "Pump", "name": "Quat. Pump", "model": "G7111B",
                                "serialNumber": "DEAEW08816", "firmwareVersion": "D.07.33 [0003]"}
    un = by["unicorn-1-trimmed.zip"]
    assert (un["instrumentModel"], un["systemName"]) == ("AKTA avant", "Avant1")
    assert un["wavelengthsNm"] == [280, 260]  # "UV 3_0" is a switched-off monitor
    em = by["empower-blanks-and-stds.json"]
    assert em["injectionCount"] == 7 and em["wavelengthsNm"] == [220]
    assert em["softwareVersion"].startswith("Empower 3 Software Build 3471")
    head = by["empower-example-01.head.json"]
    assert "channelCount" not in head and "injectionCount" not in head  # cut by the read
    ch = by["chromeleon-multi-signal.json"]
    assert ch["signals"] == ["UV_VIS_1", "UV_VIS_2"] and ch["wavelengthsNm"] == [220, 280]
    assert [m["model"] for m in ch["modules"]] == ["VF-P10-A", "VF-D40-A", "VF-A10-A"]
    assert "injectionCount" not in by["chromeleon-example.head.json"]
    ts = by["tapestation-example-01.xml"]
    assert (ts["serialNumber"], ts["sampleCount"], ts["screenTapeCount"]) == ("DEDAA00054", 4, 5)
    assert "sampleCount" not in by["tapestation-example-02.head.xml"]
    t2 = by["t200-control-immobilization.blr.sparse.json"]
    assert (t2["softwareVersion"], t2["instrumentModel"], t2["instrumentId"]) == ("2.0.1", "BiacoreT200", "1950172")
    assert t2["cycleCount"] == 4 and t2["runDate"] == "2020-01-16T16:59:38"
    ins = by["insight-transposed-flowcells.xlsx"]
    assert (ins["instrumentModel"], ins["serialNumber"], ins["cycleCount"]) == ("Biacore 8K+", "1234567", 3)


# ── truncation-tolerant JSON ─────────────────────────────────────────────────


def test_partial_json_complete_and_cut():
    v = partial_json('{"a": [1, 2.5, "x"], "b": {"c": true}}')
    assert v == {"a": [1, 2.5, "x"], "b": {"c": True}} and v.complete and v["a"].complete
    v = partial_json('{"a": [1, 2, 3], "b": {"c": "abc')
    assert not v.complete and v["a"].complete and v["b"] == {} and not v["b"].complete
    v = partial_json('{"a": [1, 2, 3')  # the trailing number may be cut: dropped
    assert v["a"] == [1, 2] and not v["a"].complete
    v = partial_json('{"a": tru')
    assert v == {} and not v.complete
    v = partial_json('{"s": "\\u00e9\\ud83d\\ude00\\n", "t": "\\ud800"}')
    assert v == {"s": "\u00e9\U0001F600\n", "t": "\ufffd"}
    assert partial_json("nope") is None and partial_json("") is None
    deep = partial_json("[" * 200)
    assert isinstance(deep, list)


# ── synthetic zips: OpenLab and Unicorn branches ─────────────────────────────


def _zip(members) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in members:
            z.writestr(name, data)
    return buf.getvalue()


def _openlab_bytes():
    return open(os.path.join(INST_DIR, "chrom", "openlab-sirius.rslt"), "rb").read()


def test_openlab_nested_zip():
    rslt = _openlab_bytes()
    outer = _zip([("Sirius.rslt", rslt), ("__MACOSX/._Sirius.rslt", b"\x00")])
    h = chrom.sniff_openlab(_src(outer))
    assert h is not None and h.fields["exportFormat"] == "zip holding an rslt result set"
    assert h.fields["injectionCount"] == 2
    assert decode_bytes("results.zip", outer).format == chrom.OPENLAB


def test_openlab_truncated_acaml_drops_counts():
    z = zipfile.ZipFile(io.BytesIO(_openlab_bytes()))
    name = next(n for n in z.namelist() if n.endswith(".acaml"))
    acaml = z.read(name)
    cut = acaml[: acaml.index(b"<Injections>") + 2000]
    h = chrom.sniff_openlab(_src(_zip([(name, cut)])))
    assert h is not None
    assert "injectionCount" not in h.fields and h.fields["sampleCount"] == 2
    assert h.fields["instrumentName"] == "Sirius"


def test_openlab_declines_other_zips():
    assert chrom.sniff_openlab(_src(_zip([("a.txt", b"hello")]))) is None
    assert chrom.sniff_openlab(_src(_zip([("x.acaml", b"<ACAML xmlns='urn:other'/>")]))) is None
    assert chrom.sniff_openlab(_src(_zip([("x.acaml", b"<!DOCTYPE x [<!ENTITY a 'b'>]><ACAML/>")]))) is None


def test_unicorn_requires_a_unicorn_result():
    uv = open(os.path.join(INST_DIR, "chrom", "unicorn-single-uv.zip"), "rb").read()
    z = zipfile.ZipFile(io.BytesIO(uv))
    members = [(n, z.read(n)) for n in z.namelist()]
    plain = [(n, d.replace(b"UNICORNVersion", b"Version") if n == "Result.xml" else d) for n, d in members]
    assert chrom.sniff_unicorn(_src(_zip(plain))) is None
    alone = [(n, d) for n, d in members if n == "Result.xml"]
    assert chrom.sniff_unicorn(_src(_zip(alone))) is None  # no chromatogram, no SystemData
    # With a UV curve named for its wavelength.
    renamed = [(n, d.replace(b"<Name>UV</Name>", b"<Name>UV 1_280</Name>")) for n, d in members]
    h = chrom.sniff_unicorn(_src(_zip(renamed)))
    assert h.fields["wavelengthsNm"] == [280] and h.fields["curves"][0] == "UV 1_280"


# ── synthetic OLE2 compound files: Biacore T200 ──────────────────────────────


def _cfb(streams) -> bytes:
    """A minimal version 3 compound file. `streams` maps a path tuple to bytes
    (a stream) or None (a storage). Streams under 4096 bytes go in the mini
    stream, as the format requires."""
    SS, MS = 512, 64
    nodes = {(): {"name": "Root Entry", "type": 5, "kids": []}}
    for path in sorted(streams, key=len):
        nodes[path] = {"name": path[-1], "type": 2 if streams[path] is not None else 1,
                       "kids": [], "data": streams[path]}
        nodes[path[:-1]]["kids"].append(path)
    order = [()] + [p for p in nodes if p != ()]
    index = {p: i for i, p in enumerate(order)}
    mini, minifat, big = b"", [], []
    for p in order[1:]:
        d = nodes[p].get("data")
        if d is None:
            continue
        if len(d) < 4096:
            nodes[p]["start"] = len(mini) // MS
            n = max(1, (len(d) + MS - 1) // MS)
            minifat += [len(minifat) + k + 1 for k in range(n - 1)] + [0xFFFFFFFE]
            mini += d.ljust(n * MS, b"\x00")
        else:
            big.append(p)
    n_dir = (len(order) * 128 + SS - 1) // SS
    n_mf = max(1, (len(minifat) * 4 + SS - 1) // SS)
    n_ms = max(1, (len(mini) + SS - 1) // SS)
    big_secs = {p: (len(nodes[p]["data"]) + SS - 1) // SS for p in big}
    n_data = n_dir + n_mf + n_ms + sum(big_secs.values())
    n_fat = 1
    while n_fat * 128 < n_data + n_fat:
        n_fat += 1
    fat = [0xFFFFFFFD] * n_fat
    def chain(n):
        s = len(fat)
        fat.extend([s + k + 1 for k in range(n - 1)] + [0xFFFFFFFE])
        return s
    dir0, mf0, ms0 = chain(n_dir), chain(n_mf), chain(n_ms)
    for p in big:
        nodes[p]["start"] = chain(big_secs[p])
    nodes[()]["start"], nodes[()]["size"] = ms0, len(mini)
    def key(p):
        return (len(nodes[p]["name"]), nodes[p]["name"].upper())
    dir_bytes = b""
    for p in order:
        nd = nodes[p]
        kids = sorted(nd["kids"], key=key)
        sib = {}
        if p:
            parent = sorted(nodes[p[:-1]]["kids"], key=key)
            i = parent.index(p)
            sib["right"] = index[parent[i + 1]] if i + 1 < len(parent) else 0xFFFFFFFF
        name = nd["name"].encode("utf-16-le") + b"\x00\x00"
        e = bytearray(128)
        e[:len(name)] = name
        struct.pack_into("<HBB", e, 64, len(name), nd["type"], 1)
        struct.pack_into("<III", e, 68, 0xFFFFFFFF, sib.get("right", 0xFFFFFFFF),
                         index[kids[0]] if kids else 0xFFFFFFFF)
        size = nd.get("size", len(nd["data"]) if nd.get("data") is not None else 0)
        struct.pack_into("<II", e, 116, nd.get("start", 0) if nd["type"] != 1 else 0, size)
        dir_bytes += bytes(e)
    h = bytearray(512)
    h[:8] = chrom.OLE_MAGIC
    struct.pack_into("<HHHHH", h, 24, 0x3E, 3, 0xFFFE, 9, 6)
    struct.pack_into("<IIIIIIIII", h, 40, 0, n_fat, dir0, 0, 4096, mf0, n_mf, 0xFFFFFFFE, 0)
    difat = list(range(n_fat)) + [0xFFFFFFFF] * (109 - n_fat)
    struct.pack_into("<109I", h, 76, *difat)
    fat += [0xFFFFFFFF] * (n_fat * 128 - len(fat))
    body = struct.pack("<%dI" % len(fat), *fat)
    body += dir_bytes.ljust(n_dir * SS, b"\x00")
    body += struct.pack("<%dI" % len(minifat), *minifat).ljust(n_mf * SS, b"\xff")
    body += mini.ljust(n_ms * SS, b"\x00")
    for p in big:
        body += nodes[p]["data"].ljust(big_secs[p] * SS, b"\x00")
    return bytes(h) + body


_ENV_CONTROL = (b"Application=Biacore T200 Control Software\nVersion=2.0.1\nRunTypeId=kineticsWizard\n"
                b"RunType=Kinetics/Affinity\nTemplateFile=C:\\Methods\\run 1.bwKin\nProcessingUnit=BiacoreT200\n"
                b"InstrumentId=1950172\nTimestamp=43082.6857180903\n")
_CHIP = b"Name=CM5\nId=171208-1131:1950172\nLotNo=10361662\nNoFcs=4\nNoSpots=1\n"


def _t200_control(app=_ENV_CONTROL, cycles=3, file_type=b"FileType=Result File\nFileTypeVersion=12\n"):
    s = {("\x03BIA filehandler info",): b"BDATA.DLL 1.0\x00",
         ("\x03BIA application info",): b"BIACOR~1.EX 0.0\x00",
         ("\x03BIA compability info",): file_type,
         ("Environment",): app, ("Chip",): _CHIP, ("AppData",): None,
         ("AppData", "Dip"): b"x" * 5000, ("RPoint Table",): b""}
    for c in range(1, cycles + 1):
        s[("_Cycle %d" % c,)] = None
        s[("_Cycle %d" % c, "_Window 1")] = None
    return _cfb(s)


def test_t200_control_synthetic():
    h = decode_bytes("run.blr", _t200_control())
    assert h is not None and h.format == chrom.T200_CONTROL
    f = h.fields
    assert (f["softwareVersion"], f["instrumentModel"], f["instrumentId"]) == ("2.0.1", "BiacoreT200", "1950172")
    assert (f["runType"], f["methodFile"], f["runDate"]) == ("Kinetics/Affinity", "run 1.bwKin", "2017-12-13T16:27:26")
    assert (f["chipType"], f["chipLotNumber"], f["flowCellCount"], f["cycleCount"]) == ("CM5", "10361662", 4, 3)
    assert h.summary == "Cytiva Biacore T200 Control result, Kinetics/Affinity, 3 cycles"


def test_t200_control_declines():
    other = _ENV_CONTROL.replace(b"Biacore T200 Control", b"Biacore X100 Control")
    assert chrom.sniff_t200_control(_src(_t200_control(app=other))) is None
    assert chrom.sniff_t200_control(_src(_t200_control(file_type=b"FileType=Something\n"))) is None
    assert chrom.sniff_t200_control(_src(_cfb({("Workbook",): b"\x09\x08" * 100}))) is None
    assert chrom.sniff_t200_evaluation(_src(_t200_control())) is None


def _t200_eval(total=b"7"):
    s = {("\x03BIA filehandler info",): b"BDATA.DLL 1.0\x00",
         ("\x03BIA application info",): b"T200Eval 3.1\x00",
         ("\x03BIA compability info",): b"FileType=T200 Evaluation File\nFileTypeVersion=5\n",
         ("Environment",): b"Application=Biacore T200 Evaluation Software\nVersion=3.2.1\n",
         ("DataManager",): b"Version=1\nTotalNumberOfCycles=" + total,
         ("_DataManager 1",): None,
         ("_DataManager 1", "Environment"): _ENV_CONTROL.replace(b"TemplateFile", b"MethodFile"),
         ("_DataManager 1", "Chip"): _CHIP,
         ("Evaluation",): None}
    for c in range(1, 6):
        s[("_DataManager 1", "_Cycle %d" % c)] = None
    for k in range(3):
        s[("Evaluation", "EvaluationItem%d" % k)] = None
        s[("Evaluation", "EvaluationItem%dBinary" % k)] = b"\x00" * 10
    return _cfb(s)


def test_t200_evaluation_synthetic():
    h = decode_bytes("eval.bme", _t200_eval())
    assert h is not None and h.format == chrom.T200_EVALUATION
    f = h.fields
    assert (f["softwareVersion"], f["applicationInfo"]) == ("3.2.1", "T200Eval 3.1")
    assert (f["acquisitionSoftware"], f["acquisitionSoftwareVersion"]) == ("Biacore T200 Control Software", "2.0.1")
    assert (f["instrumentModel"], f["instrumentId"], f["chipType"]) == ("BiacoreT200", "1950172", "CM5")
    assert f["cycleCount"] == 5  # counted from the _Cycle storages
    assert f["evaluationItemCount"] == 3
    assert chrom.sniff_t200_control(_src(_t200_eval())) is None


def test_t200_evaluation_identifies_when_the_directory_is_out_of_reach():
    # A large evaluation file spreads its directory over the whole file; past
    # the read budget only the file-type streams in the first directory
    # sector are reachable, and the sniffer reports identity alone.
    data = bytearray(_t200_eval())
    n_dir = (24 * 128 + 511) // 512  # 24 directory entries
    data[1536:1024 + n_dir * 512] = b"\x00" * ((n_dir - 1) * 512)  # all but the first directory sector
    h = chrom.sniff_t200_evaluation(_src(bytes(data)))
    assert h is not None
    assert "cycleCount" not in h.fields and "softwareVersion" not in h.fields


# ── synthetic .xlsx: Biacore Insight ─────────────────────────────────────────


def _xlsx(sheets) -> bytes:
    """sheets: [(name, rows)], rows of str/number cells."""
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    shared = []
    parts = []
    for name, rows in sheets:
        xr = []
        for r, row in enumerate(rows, 1):
            cells = []
            for c, v in enumerate(row):
                if v is None:
                    continue
                ref = chr(65 + c) + str(r)
                if isinstance(v, str):
                    shared.append(v)
                    cells.append('<c r="%s" t="s"><v>%d</v></c>' % (ref, len(shared) - 1))
                else:
                    cells.append('<c r="%s"><v>%s</v></c>' % (ref, v))
            xr.append('<row r="%d">%s</row>' % (r, "".join(cells)))
        parts.append("<worksheet %s><sheetData>%s</sheetData></worksheet>" % (ns, "".join(xr)))
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("xl/workbook.xml", "<workbook %s %s><sheets>%s</sheets></workbook>" % (ns, rns, "".join(
            '<sheet name="%s" sheetId="%d" r:id="rId%d"/>' % (n, i + 1, i + 1) for i, (n, _) in enumerate(sheets))))
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">%s</Relationships>'
                   % "".join('<Relationship Id="rId%d" Target="worksheets/sheet%d.xml"/>' % (i + 1, i + 1)
                             for i in range(len(sheets))))
        z.writestr("xl/sharedStrings.xml", "<sst %s>%s</sst>" % (ns, "".join("<si><t>%s</t></si>" % s for s in shared)))
        for i, p in enumerate(parts):
            z.writestr("xl/worksheets/sheet%d.xml" % (i + 1), p)
    return buf.getvalue()


_PROPS = [
    ["Properties"],
    ["Evaluation", None, None, "Run(s) in this evaluation"],
    ["Name", "Kinetics 1", None, "Name", "Run A"],
    ["Software"],
    ["Name", "Biacore Insight Evaluation"],
    ["Version", "6.0.7.1750"],
    ["Run 1"],
    ["Chip information", None, None, "Run information"],
    ["Chip type", "CM5", None, "Name", "Run A", None, "Cycles", 12.0],
    ["Lot number", 10361662.0, None, "Method", "Kinetics method", None, "Data collection rate", 10.0],
    ["First dock date", "x", None, "Instrument type", "Biacore 8K", None, "Instrument id", 2054321.0],
    [None, None, None, None, None, None, "Running buffer", "HBS-EP+"],
    ["Run 2"],
]


def test_insight_synthetic():
    data = _xlsx([("Properties", _PROPS), ("Report point table", [["Run", "Cycle"]]),
                  ("Evaluation - Kinetics", [["Group"]])])
    h = decode_bytes("export.xlsx", data)
    assert h is not None and h.format == chrom.INSIGHT
    f = h.fields
    assert (f["softwareVersion"], f["evaluationName"], f["runCount"]) == ("6.0.7.1750", "Kinetics 1", 2)
    assert (f["instrumentModel"], f["serialNumber"], f["methodName"]) == ("Biacore 8K", "2054321", "Kinetics method")
    assert (f["cycleCount"], f["dataCollectionRateHz"], f["chipLotNumber"]) == (12, 10, "10361662")
    assert f["runningBuffer"] == "HBS-EP+" and f["evaluationSheets"] == ["Evaluation - Kinetics"]


def test_insight_declines():
    no_rpt = _xlsx([("Properties", _PROPS)])
    assert chrom.sniff_insight(_src(no_rpt)) is None
    other = [r[:] for r in _PROPS]
    other[4] = ["Name", "Other Software"]
    assert chrom.sniff_insight(_src(_xlsx([("Properties", other), ("Report point table", [[]])]))) is None


# ── declines ─────────────────────────────────────────────────────────────────

GENERIC = {
    "plain.csv": b"a,b,c\n1,2,3\n",
    "plain.txt": b"Software Version\t3.1\nWavelength\t450\n",
    "plain.xml": b'<?xml version="1.0"?><File><Name>x</Name></File>',
    "plain.json": b'{"values": {"fields": {"a": 1}}, "injections": []}',
    "array.json": b"[1, 2, 3]",
    "bids.json": b'{"SamplingFrequency": 500, "EEGChannelCount": 8}',
    "empty.zip": _zip([("a.txt", b"x")]),
    "xlsx.xlsx": _xlsx([("Sheet1", [["a", 1]])]),
    "ole.bme": _cfb({("Workbook",): b"\x00" * 100}),
}


def _fixture_files():
    out = []
    for fam in sorted(os.listdir(INST_DIR)):
        d = os.path.join(INST_DIR, fam)
        if os.path.isdir(d):
            for n in sorted(os.listdir(d)):
                # Skip notices and subfolders (plate_b/other holds copies of
                # exports that belong to other families).
                if n != "NOTICE.md" and os.path.isfile(os.path.join(d, n)):
                    out.append((fam, os.path.join(d, n)))
    for n in ("instrument.ome.xml", "bdv-drosophila.xml", "bids-cases.json"):
        out.append(("repo", os.path.join(FIXTURES_DIR, n)))
    return out


def test_each_sniffer_declines_everything_else():
    owner = {c["file"].split("/", 1)[1]: c["expected"]["format"] for c in CASES}
    for fam, path in _fixture_files():
        own = owner.get(os.path.basename(path)) if fam == "chrom" else None
        for name, sniff in SNIFFERS.items():
            if name == own:
                continue
            assert sniff(_source(path)) is None, (name, path)
    for gname, data in GENERIC.items():
        for name, sniff in SNIFFERS.items():
            assert sniff(_src(data)) is None, (name, gname)
        h = decode_bytes(gname, data)
        assert h is None or h.format not in SNIFFERS


# ── totality ─────────────────────────────────────────────────────────────────

SEED = 1729


def _seeds():
    ol = _openlab_bytes()
    return [
        b"PK\x03\x04", ol, _zip([("x.acaml", b"<ACAML xmlns='urn:schemas-agilent-com:acaml21'><Doc><Content>")]),
        _zip([("Result.xml", b'<Result UNICORNVersion="7"><SystemName>a'), ("Chrom.1.Xml", b"<Chromatogram")]),
        b'{"values": {"fields": {"SampleSetId": 1, "SampleSetName": "a", "SystemName": "b"}, "channels": [{"fields": {"ChannelDescription": "220 nm"',
        b'{"sequence": {"sequence name": "s"}, "injections": [{"injection number": 1, "signals": [{"wavelength setting": 1e308',
        b"<File><FileInformation><Assay>RNA</Assay></FileInformation><ScreenTapes><ScreenTape>",
        _t200_control(), _t200_eval(),
        _xlsx([("Properties", _PROPS), ("Report point table", [[]])]),
    ]


def test_sniffers_total_on_mutated_inputs():
    rng = random.Random(SEED)
    seeds = _seeds() + [open(p, "rb").read() for f, p in _fixture_files()
                        if f == "chrom" and not p.endswith(".sparse.json")]
    for _ in range(600):
        buf = bytearray(rng.choice(seeds))
        op = rng.randrange(4)
        if op == 0 and buf:
            buf = buf[: rng.randrange(len(buf))]
        elif op == 1:
            for _ in range(rng.randint(1, 16)):
                if buf:
                    buf[rng.randrange(len(buf))] = rng.getrandbits(8)
        elif op == 2 and len(buf) > 8:
            a = rng.randrange(len(buf))
            buf[a:a + rng.randint(1, 64)] = bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 64)))
        data = bytes(buf)
        for name, sniff in SNIFFERS.items():
            try:
                _valid_or_none(sniff(_src(data)))
            except Exception as e:  # noqa: BLE001
                pytest.fail("%s raised on a mutated input (len %d): %r" % (name, len(data), e))


def test_sniffers_total_on_random_bytes():
    rng = random.Random(SEED + 1)
    prefixes = [b"", b"{", b"[", b"<", b"PK\x03\x04", chrom.OLE_MAGIC, b"\xef\xbb\xbf{", b"\xff\xfe"]
    for _ in range(1500):
        data = rng.choice(prefixes) + bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 3000)))
        for name, sniff in SNIFFERS.items():
            try:
                _valid_or_none(sniff(_src(data)))
            except Exception as e:  # noqa: BLE001
                pytest.fail("%s raised on random bytes (len %d): %r" % (name, len(data), e))
        for key in ("x.json", "x.xml", "x.zip", "x.rslt", "x.blr", "x.bme", "x.xlsx"):
            _valid_or_none(decode_bytes(key, data))


# ── optional cross-check against allotropy's own test data ───────────────────
#
# Set SCIGANTIC_ALLOTROPY_TESTDATA to an allotropy checkout's tests/parsers to
# decode every input there and compare with the ASM allotropy expects.

ALLOTROPY = os.environ.get("SCIGANTIC_ALLOTROPY_TESTDATA", "")


def _asm(path):
    d = json.load(open(path))
    return next(v for k, v in d.items() if k != "$asm.manifest")


@pytest.mark.skipif(not ALLOTROPY, reason="SCIGANTIC_ALLOTROPY_TESTDATA not set")
def test_against_allotropy_testdata():
    from scigantic_headers.sources import FileSource

    plan = [
        ("agilent_openlab_cds/testdata", ".rslt", ".json"),
        ("cytiva_unicorn/testdata", ".zip", ".json"),
        ("benchling_empower/testdata/input", ".json", None),
        ("benchling_chromeleon/testdata/input", ".json", None),
        ("agilent_tapestation_analysis/testdata", ".xml", ".json"),
        ("cytiva_biacore_t200_control/testdata", ".blr", ".json"),
        ("cytiva_biacore_t200_evaluation/testdata", ".bme", ".json"),
        ("cytiva_biacore_insight/testdata", ".xlsx", ".json"),
    ]
    n = 0
    for sub, ext, out_ext in plan:
        d = os.path.join(ALLOTROPY, sub)
        for name in sorted(os.listdir(d)):
            if not name.endswith(ext):
                continue
            path = os.path.join(d, name)
            h = decode_source(path, FileSource(path, INSTRUMENT_SEEK_BYTES))
            assert h is not None and h.format in SNIFFERS, path
            out = (os.path.join(d.replace("/input", "/output"), name) if out_ext is None
                   else path[: -len(ext)] + out_ext)
            if not os.path.exists(out):
                continue
            a = _asm(out)
            sw = a.get("data system document", {}).get("software version")
            if "softwareVersion" in h.fields:
                assert h.fields["softwareVersion"] == sw, path
            n += 1
    assert n >= 25
