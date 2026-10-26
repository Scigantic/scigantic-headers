"""SpikeGLX run metadata decoding.

Golden: tests/fixtures/spikeglx-cases.json, three real `.meta` files spanning
nine years of the format. Channel count, sampling rate, sample count and
microvolts per bit were checked against neo 0.14.6's extract_stream_info on
the two modern ones when the fixture was made; neo raises KeyError on the 2016
Phase 3A file, which this decoder still describes.
"""

import json
import os

from scigantic_headers import decode_bytes, decode_spikeglx_meta, has_decoder_for
from scigantic_headers.decoders import read_for
from scigantic_headers.spikeglx import SPIKEGLX_READ_BYTES

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "spikeglx-cases.json")

with open(FIX) as _fh:
    CASES = json.load(_fh)["cases"]


def _case(name):
    return bytes.fromhex(next(c for c in CASES if name in c["name"])["hex"])


def test_extension_registered():
    assert has_decoder_for("run_g0_t0.imec0.ap.meta")
    assert read_for("x.meta").leading == SPIKEGLX_READ_BYTES


def test_golden_fixture_matches():
    assert len(CASES) >= 3
    for case in CASES:
        d = decode_bytes(case["key"], bytes.fromhex(case["hex"]))
        assert d is not None and d.to_dict() == case["expected"], case["name"]


def test_neuropixels_2_four_shank():
    d = decode_spikeglx_meta(_case("1119617"))
    assert d.format == "spikeglx"
    assert d.summary == ("SpikeGLX imec.ap, 385 channels, 30000 Hz, 0.25 s, "
                         "Neuropixels 2.0 (4 shank)")
    f = d.fields
    assert f["channelCount"] == 385 and f["samplingRateHz"] == 30000.0
    assert f["apChannels"] == 384 and f["syncChannels"] == 1
    assert f["probeType"] == 24 and f["probePartNumber"] == "NP2010"
    assert f["probeSerial"] == "20403314871" and f["maxInt"] == 8192
    # int16 samples, no header: 5775000 / 2 / 385.
    assert f["sampleCount"] == 7500
    # imAiRangeMax / imMaxInt / fixed 2.0 gain, in microvolts; neo agrees.
    assert f["microvoltsPerBit"] == 0.762939453125
    assert f["runName"] == "1119617_LSE1_shank12"
    assert f["gateIndex"] == 0 and f["triggerIndex"] == 0
    assert f["device"] == "imec0" and f["stream"] == "ap"
    # The per-channel tables are counted, never copied.
    assert f["snsGeomMapEntries"] == 384 and "snsGeomMap" not in f


def test_neuropixels_1_gain_from_imro_table():
    f = decode_spikeglx_meta(_case("sample3B2")).fields
    assert f["probeType"] == 0 and f["probeModel"] == "Neuropixels 1.0"
    assert f["apGain"] == 500 and f["lfGain"] == 250
    # imAiRangeMax / 512 / apGain, in microvolts; neo agrees.
    assert f["microvoltsPerBit"] == 2.34375
    # A CatGT-concatenated run is `_g0_tcat`, with no trigger index.
    assert f["concatenated"] is True and "triggerIndex" not in f
    assert f["runName"] == "SC023_053119"


def test_phase_3a_file_without_type_key():
    f = decode_spikeglx_meta(_case("sample3.meta")).fields
    assert f["channelCount"] == 769 and f["samplingRateHz"] == 30000.0
    assert f["apChannels"] == 384 and f["lfChannels"] == 384
    # No imDatPrb_type and no imMaxInt, so no scaling is claimed.
    assert "probeType" not in f and "microvoltsPerBit" not in f
    assert f["durationSeconds"] == 120.0


def test_rejects_other_meta_files():
    assert decode_spikeglx_meta(b"") is None
    assert decode_spikeglx_meta(b"name=value\nother=1\n") is None
    # nSavedChans without a sample rate is not enough to claim the file.
    assert decode_spikeglx_meta(b"nSavedChans=385\n") is None
    assert decode_spikeglx_meta(b"nSavedChans=385\nimSampRate=30000\n") is not None
