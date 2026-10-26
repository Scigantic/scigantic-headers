"""Open Ephys binary-format manifest decoding.

Golden: tests/fixtures/openephys-cases.json, three real `structure.oebin`
files: the GUI's own demo manifest (0.4.5, four streams), a Neuropixels
session written by GUI 0.6.0, and a Rhythm FPGA session from GUI 0.4.4.1.
"""

import json
import os

from scigantic_headers import decode_bytes, decode_openephys_oebin, has_decoder_for
from scigantic_headers.decoders import read_for
from scigantic_headers.openephys import OPENEPHYS_READ_BYTES

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "openephys-cases.json")

with open(FIX) as _fh:
    CASES = json.load(_fh)["cases"]


def _case(name):
    return bytes.fromhex(next(c for c in CASES if name in c["name"])["hex"])


def test_extension_registered():
    assert has_decoder_for("recording1/structure.oebin")
    assert read_for("structure.oebin").leading == OPENEPHYS_READ_BYTES


def test_golden_fixture_matches():
    assert len(CASES) >= 3
    for case in CASES:
        d = decode_bytes(case["key"], bytes.fromhex(case["hex"]))
        assert d is not None and d.to_dict() == case["expected"], case["name"]


def test_four_demo_streams():
    d = decode_openephys_oebin(_case("plugin-GUI"))
    assert d.format == "openephys-oebin"
    assert d.summary == "Open Ephys recording, 4 continuous streams, 64 channels, 40000 Hz"
    f = d.fields
    assert f["guiVersion"] == "0.4.5" and f["continuousStreams"] == 4
    assert f["totalContinuousChannels"] == 64 and f["samplingRateHz"] == 40000.0
    first = f["continuous"][0]
    assert first["folderName"] == "data_stream_sine_wave/"
    assert first["channelCount"] == 16 and first["channelsListed"] == 16
    assert first["bitVolts"] == 0.05 and first["units"] == "uV"
    assert first["sourceProcessorName"] == "Demo source"


def test_neuropixels_session():
    f = decode_openephys_oebin(_case("spikewrap")).fields
    assert f["guiVersion"] == "0.6.0" and f["continuousStreams"] == 1
    assert f["totalContinuousChannels"] == 385
    assert f["streamNames"] == ["ProbeA-AP"]
    assert f["sourceProcessors"] == ["Neuropix-PXI"]
    # The sync channel scales differently from the 384 probe channels, so no
    # single bitVolts is claimed for the stream.
    assert "bitVolts" not in f["continuous"][0]
    assert f["eventStreams"] == 2
    assert f["events"][0]["channelName"] == "Neuropixels PXI Sync"


def test_rhythm_fpga_session():
    f = decode_openephys_oebin(_case("3DDD")).fields
    assert f["guiVersion"] == "0.4.4.1" and f["totalContinuousChannels"] == 35
    assert f["continuous"][0]["sourceProcessorName"] == "Rhythm FPGA"


def test_rejects_other_json():
    assert decode_openephys_oebin(b"") is None
    assert decode_openephys_oebin(b"[]") is None
    assert decode_openephys_oebin(b'{"name": "package"}') is None
    # Record lists present but all empty, and no GUI version: nothing to describe.
    assert decode_openephys_oebin(b'{"continuous": [], "events": []}') is None
    assert decode_openephys_oebin(b'{"GUI version": "0.6.0", "continuous": []}') is not None
