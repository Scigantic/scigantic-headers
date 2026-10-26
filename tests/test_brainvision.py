"""BrainVision header and marker decoding.

Golden: tests/fixtures/brainvision-cases.json, real files from the
bids-standard/bids-examples corpus. Channel count, sampling rate and channel
names were checked against mne 1.13.2's read_raw_brainvision on all 64
non-empty `.vhdr` files in that corpus when the fixture was made; all 64
agreed.
"""

import json
import os

from scigantic_headers import (
    decode_brainvision_vhdr,
    decode_brainvision_vmrk,
    decode_bytes,
    has_decoder_for,
)
from scigantic_headers.brainvision import BRAINVISION_READ_BYTES
from scigantic_headers.decoders import read_for

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "brainvision-cases.json")

with open(FIX) as _fh:
    CASES = json.load(_fh)["cases"]


def _case(name):
    return bytes.fromhex(next(c for c in CASES if name in c["name"])["hex"])


def test_extensions_registered():
    assert has_decoder_for("eeg/sub-01_eeg.vhdr") and has_decoder_for("eeg/sub-01_eeg.vmrk")
    assert read_for("x.vhdr").leading == BRAINVISION_READ_BYTES


def test_golden_fixture_matches():
    assert len(CASES) >= 5
    for case in CASES:
        d = decode_bytes(case["key"], bytes.fromhex(case["hex"]))
        assert d is not None and d.to_dict() == case["expected"], case["name"]


def test_version_1_header():
    d = decode_brainvision_vhdr(_case("sub-cm4"))
    assert d.format == "brainvision"
    assert d.summary == "BrainVision, 61 channels, 1000 Hz, ieee_float_32 multiplexed"
    f = d.fields
    assert f["brainVisionVersion"] == "1.0" and f["channelCount"] == 61
    assert f["samplingIntervalMicroseconds"] == 1000.0 and f["samplingRateHz"] == 1000.0
    assert f["binaryFormat"] == "IEEE_FLOAT_32" and f["bytesPerSample"] == 4
    assert f["byteOrder"] == "little" and f["dataOrientation"] == "MULTIPLEXED"
    assert f["dataFile"] == "sub-cm4_task-FilteredSpeech_ieeg.eeg"
    assert f["channelNames"][0] == "G1" and len(f["channelNames"]) == 61
    assert f["channels"][0] == {"name": "G1", "resolution": 0.1}
    assert f["resolution"] == 0.1


def test_version_2_vectorized_header_with_coordinates():
    f = decode_brainvision_vhdr(_case("sub-35")).fields
    assert f["brainVisionVersion"] == "2.0" and f["dataOrientation"] == "VECTORIZED"
    assert f["dataType"] == "TIMEDOMAIN" and f["channelCount"] == 64
    assert f["samplingRateHz"] == 250.0 and f["dataPoints"] == 197235
    assert f["durationSeconds"] == 197235 / 250.0
    assert f["hasCoordinates"] is True


def test_markers():
    d = decode_brainvision_vmrk(_case("sub-05_task-matchingpennies_eeg.vmrk"))
    assert d.format == "brainvision-markers"
    f = d.fields
    assert f["markersRead"] == 301 and f["markersComplete"] is True
    assert f["markerTypes"] == {"Stimulus": 300, "New Segment": 1}
    assert len(f["markers"]) == 16  # detail is capped; the counts are not
    assert f["markers"][0]["type"] == "New Segment"


def test_marker_segment_timestamp():
    # `Mk1=New Segment,,1,1,0,20131129170516793000` carries the wall clock.
    text = (
        "Brain Vision Data Exchange Marker File, Version 1.0\n"
        "\n[Common Infos]\nDataFile=x.eeg\n"
        "\n[Marker Infos]\n"
        "Mk1=New Segment,,1,1,0,20131129170516793000\n"
        "Mk2=Stimulus,S  1,1024,1,0\n"
    ).encode()
    f = decode_brainvision_vmrk(text).fields
    assert f["segmentStartTime"] == "2013-11-29T17:05:16.793000"
    assert f["markers"][1] == {"type": "Stimulus", "description": "S  1", "position": 1024,
                               "size": 1, "channel": 0}


def test_escaped_comma_in_channel_name():
    text = (
        "Brain Vision Data Exchange Header File Version 1.0\n"
        "\n[Common Infos]\nNumberOfChannels=1\nSamplingInterval=1000\n"
        "\n[Channel Infos]\nCh1=Fp1\\1ref,,0.1,uV\n"
    ).encode()
    f = decode_brainvision_vhdr(text).fields
    assert f["channelNames"] == ["Fp1,ref"]
    assert f["channels"][0]["unit"] == "uV"


def test_rejects_other_ini_files():
    assert decode_brainvision_vhdr(b"[Common Infos]\nNumberOfChannels=1\n") is None
    assert decode_brainvision_vhdr(b"") is None
    # A marker file is not a header file and vice versa.
    assert decode_brainvision_vhdr(_case("sub-05_task-matchingpennies_eeg.vmrk")) is None
    assert decode_brainvision_vmrk(_case("sub-cm4_task-FilteredSpeech_ieeg.vhdr")) is None
