"""Neuralynx Cheetah and Pegasus header decoding.

Golden: tests/fixtures/neuralynx-cases.json. No public Neuralynx file was
reachable to copy, so each case is a header built to the per-version grammar
documented in python-neo's NlxHeader, the reference reader for this format.
neo 0.14.6 parsed all four and agreed on application, version, file type,
sampling rate, microvolts per bit, channel name, record size and the open and
close times. The last case keeps the 16 KiB of NUL padding a real file has.
"""

import json
import os

from scigantic_headers import decode_bytes, decode_neuralynx, has_decoder_for
from scigantic_headers.decoders import read_for
from scigantic_headers.neuralynx import NEURALYNX_HEADER_BYTES

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "neuralynx-cases.json")

with open(FIX) as _fh:
    CASES = json.load(_fh)["cases"]


def _case(name):
    return bytes.fromhex(next(c for c in CASES if name in c["name"])["hex"])


def test_extensions_registered():
    for ext in ("ncs", "nev", "nse", "ntt", "nvt"):
        assert has_decoder_for(f"session/CSC1.{ext}")
    assert read_for("x.ncs").leading == NEURALYNX_HEADER_BYTES


def test_golden_fixture_matches():
    assert len(CASES) >= 5
    for case in CASES:
        d = decode_bytes(case["key"], bytes.fromhex(case["hex"]))
        assert d is not None and d.to_dict() == case["expected"], case["name"]


def test_comment_line_dialect():
    d = decode_neuralynx(_case("Cheetah 5.5.1"))
    assert d.format == "neuralynx"
    assert d.summary == "Neuralynx CSC, CSC1, 32000 Hz, Cheetah 5.5.1"
    f = d.fields
    assert f["applicationName"] == "Cheetah" and f["applicationVersion"] == "5.5.1"
    assert f["fileType"] == "CSC" and f["recordKind"] == "continuous"
    assert f["channelName"] == "CSC1" and f["recordBytes"] == 1044
    assert f["samplingRateHz"] == 32000.0 and f["adMaxValue"] == 32767
    assert f["microvoltsPerBit"] == 0.030517578125
    assert f["timeOpened"] == "2013-11-29T17:05:16.793"
    assert f["timeClosed"] == "2013-11-29T18:03:13.603"
    assert f["dspLowCutFilterEnabled"] is True and f["inputInverted"] is True
    assert f["headerBytes"] == NEURALYNX_HEADER_BYTES


def test_property_dialect():
    f = decode_neuralynx(_case("Cheetah 5.7.4 CSC1.ncs (property")).fields
    # `-ApplicationName Cheetah "5.7.4 Development"` splits into both.
    assert f["applicationName"] == "Cheetah"
    assert f["applicationVersion"] == "5.7.4 Development"
    assert f["fileVersion"] == "3.4"
    assert f["timeOpened"] == "2017-02-16T17:56:04"
    assert f["originalFileName"] == "C:\\CheetahData\\2017-02-16_17-55-55\\CSC1.ncs"
    # A quoted value is stored without its quotes.
    assert f["referenceChannel"] == "Source 0 Reference 0"


def test_pegasus_event_file():
    f = decode_neuralynx(_case("Pegasus")).fields
    assert f["applicationName"] == "Pegasus" and f["applicationVersion"] == "2.1.1"
    assert f["fileType"] == "Event" and f["recordKind"] == "event"
    assert f["recordBytes"] == 184
    assert "samplingRateHz" not in f


def test_tetrode_repeats_keys_per_channel():
    f = decode_neuralynx(_case("tetrode")).fields
    assert f["fileType"] == "Spike" and f["recordKind"] == "spike"
    assert f["adChannelCount"] == 4 and f["waveformSamples"] == 32
    # Four ADBitVolts on one line; the first is the value a scalar field takes.
    assert f["microvoltsPerBit"] == 0.03662109375
    assert f["inputRange"] == "1200 1200 1200 1200"


def test_nul_padding_decodes_the_same():
    padded = _case("NUL padding")
    assert len(padded) == NEURALYNX_HEADER_BYTES
    assert decode_neuralynx(padded) == decode_neuralynx(_case("property"))


def test_rejects_other_files():
    assert decode_neuralynx(b"") is None
    assert decode_neuralynx(b"########" + b" not a neuralynx file\n") is None
    assert decode_neuralynx(b"-FileType CSC\n-SamplingFrequency 32000\n") is None
