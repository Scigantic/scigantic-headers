"""Illumina run-metadata decoding (RunInfo.xml / RunParameters.xml)."""

import json
import xml.etree.ElementTree as ET

from scigantic_headers import decode_illumina_run
from scigantic_headers.decoders import decode_bytes, has_decoder_for

RUNINFO = b"""<?xml version="1.0"?>
<RunInfo Version="2">
  <Run Id="200101_M00123_0007_000000000-ABCDE" Number="7">
    <Flowcell>000000000-ABCDE</Flowcell>
    <Instrument>M00123</Instrument>
    <Date>200101</Date>
    <Reads>
      <Read Number="1" NumCycles="151" IsIndexedRead="N"/>
      <Read Number="2" NumCycles="8" IsIndexedRead="Y"/>
      <Read Number="3" NumCycles="8" IsIndexedRead="Y"/>
      <Read Number="4" NumCycles="151" IsIndexedRead="N"/>
    </Reads>
    <FlowcellLayout LaneCount="1" SurfaceCount="2" SwathCount="1" TileCount="14"/>
  </Run>
</RunInfo>"""


def test_runinfo_core_fields():
    h = decode_illumina_run(RUNINFO)
    assert h.format == "illumina-run"
    assert h.fields["runId"] == "200101_M00123_0007_000000000-ABCDE"
    assert h.fields["instrument"] == "M00123"
    assert h.fields["flowcell"] == "000000000-ABCDE"
    assert h.fields["date"] == "200101"
    assert h.fields["laneCount"] == 1
    assert len(h.fields["reads"]) == 4
    assert h.fields["reads"][0] == {"numCycles": 151, "isIndex": False}
    assert h.fields["reads"][1] == {"numCycles": 8, "isIndex": True}


def test_dispatched_by_filename_not_extension():
    assert has_decoder_for("run/RunInfo.xml")
    assert has_decoder_for("run/runinfo.xml")  # case-insensitive
    assert not has_decoder_for("run/something.xml")  # .xml alone is not enough
    assert decode_bytes("RunInfo.xml", RUNINFO).format == "illumina-run"
    assert decode_bytes("other.xml", RUNINFO) is None


def test_matches_elementtree():
    h = decode_illumina_run(RUNINFO)
    root = ET.fromstring(RUNINFO)
    assert h.fields["runId"] == root.find(".//Run").get("Id")
    assert h.fields["instrument"] == root.findtext(".//Instrument")
    assert h.fields["flowcell"] == root.findtext(".//Flowcell")
    assert len(h.fields["reads"]) == len(root.findall(".//Read"))


def test_runparameters_best_effort():
    rp = b'<?xml version="1.0"?><RunParameters><RunId>xyz</RunId><ScannerID>NB500</ScannerID></RunParameters>'
    h = decode_illumina_run(rp)
    assert h.fields["runId"] == "xyz"
    assert h.fields["instrument"] == "NB500"
    assert h.fields["reads"] == []


def test_rejects_non_illumina():
    assert decode_illumina_run(b"<somethingElse/>") is None
    assert decode_illumina_run(b"") is None


def test_json_safe():
    json.dumps(decode_illumina_run(RUNINFO).to_dict(), allow_nan=False)
