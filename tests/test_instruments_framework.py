"""The instrument-export dispatch: sniffers run before an extension's own
decoder, decline cleanly, and an .xlsx is read through a bounded Source."""

import io
import zipfile

from scigantic_headers import decode_bytes, has_decoder_for, read_for
from scigantic_headers import decoders
from scigantic_headers.decoders import BytesSource, DecodedHeader
from scigantic_headers.instruments import head_text, text_of, xlsx_rows, xlsx_sheet_names


def _with_sniffer(ext, sniffer):
    saved = dict(decoders._INSTRUMENT_SNIFFERS)
    decoders.register_instrument_sniffer(ext, sniffer)
    return saved


def test_sniffer_runs_first_and_falls_back(monkeypatch):
    monkeypatch.setattr(decoders, "_INSTRUMENT_SNIFFERS", {})

    def sniff(src):
        t = head_text(src)
        if t and t.startswith("MAGIC"):
            return DecodedHeader("Test", "test", {})
        return None

    def boom(src):
        raise RuntimeError("a sniffer that raises is a decline")

    decoders.register_instrument_sniffer("json", boom)
    decoders.register_instrument_sniffer(["json", "zzz"], sniff)
    assert has_decoder_for("x.zzz")
    assert read_for("x.json").seek == decoders.INSTRUMENT_SEEK_BYTES
    assert decode_bytes("x.json", b"MAGIC 1").format == "Test"
    # Declined: the BIDS JSON decoder still gets the file.
    bids = decode_bytes("sub-01_eeg.json", b'{"SamplingFrequency": 500, "EEGChannelCount": 8}')
    assert bids is not None and bids.format != "Test"
    assert decode_bytes("x.zzz", b"nope") is None


def test_text_of_encodings():
    s = "Software Version\t3.11\r\n"
    assert text_of(s.encode("utf-16")) == s
    assert text_of(s.encode("utf-16-le")) == s
    assert text_of(b"\xef\xbb\xbf" + s.encode()) == s
    assert text_of("café au lait".encode("cp1252")) == "café au lait"
    assert text_of("café".encode("utf-8")[:-1]) == "caf"
    assert text_of(b"\x00\x01\x02\x03binary\x00") is None
    assert text_of(b"") is None


def _xlsx(rows, sheet="Results"):
    shared = []
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
        if cells:  # Excel omits empty rows
            xml_rows.append('<row r="%d">%s</row>' % (r, "".join(cells)))
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("xl/workbook.xml", '<workbook %s %s><sheets><sheet name="%s" sheetId="1" r:id="rId1"/></sheets></workbook>' % (ns, rns, sheet))
        z.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr("xl/sharedStrings.xml", "<sst %s>%s</sst>" % (ns, "".join("<si><t>%s</t></si>" % s for s in shared)))
        z.writestr("xl/worksheets/sheet1.xml", "<worksheet %s><sheetData>%s</sheetData></worksheet>" % (ns, "".join(xml_rows)))
    return buf.getvalue()


def test_xlsx_rows_through_a_source():
    data = _xlsx([["Instrument", "Reader"], ["Wells", 96], [], ["", "skip", "x"]])
    src = BytesSource(data, budget=1 << 20)
    assert xlsx_sheet_names(src) == ["Results"]
    rows = xlsx_rows(src)
    assert rows[0] == ["Instrument", "Reader"]
    assert rows[1] == ["Wells", "96"]
    assert rows[2] == []
    assert rows[3] == ["", "skip", "x"]
    assert xlsx_rows(src, sheet="Missing") is None
    assert xlsx_rows(BytesSource(b"not a zip")) is None
