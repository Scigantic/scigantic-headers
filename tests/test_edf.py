"""EDF/EDF+ and BDF/BDF+ header decoding.

Golden: tests/fixtures/edf-cases.json. The first case is the real header of the
PhysioNet Sleep-EDF Expanded cassette recording SC4001E0-PSG.edf; the rest were
written by pyedflib 0.1.42. Every expected value was checked against pyedflib's
own reader on the same files when the fixture was made, with two known and
intended differences: pyedflib hides the annotations channel from its signal
list, and its start time carries the sub-second offset that lives in the
annotations record rather than the header this decoder reads.
"""

import json
import os

from scigantic_headers import decode_bytes, decode_edf, has_decoder_for
from scigantic_headers.decoders import read_for
from scigantic_headers.edf import EDF_READ_BYTES
from scigantic_headers.sources import is_decodable

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "edf-cases.json")

with open(FIX) as _fh:
    CASES = json.load(_fh)["cases"]


def _case(name):
    return next(c for c in CASES if name in c["name"])


def test_extensions_registered():
    assert has_decoder_for("sleep/SC4001E0-PSG.edf") and has_decoder_for("rec.bdf")
    assert is_decodable("sleep/SC4001E0-PSG.edf.gz")
    assert read_for("x.edf").leading == EDF_READ_BYTES


def test_golden_fixture_matches():
    assert len(CASES) >= 6
    for case in CASES:
        d = decode_bytes(case["key"], bytes.fromhex(case["hex"]))
        assert d is not None and d.to_dict() == case["expected"], case["name"]


def test_physionet_sleep_cassette_header():
    d = decode_edf(bytes.fromhex(_case("SC4001E0")["hex"]))
    assert d.format == "edf"
    assert d.summary == "EDF recording, 7 signals, 22:05:00"
    f = d.fields
    assert f["edfType"] == "EDF" and f["bitDepth"] == 16
    assert f["signalCount"] == 7 and f["headerBytes"] == 256 * 8
    assert f["numDataRecords"] == 2650 and f["dataRecordDurationSeconds"] == 30.0
    assert f["durationSeconds"] == 79500.0 and f["headerComplete"] is True
    # The header date is 24.04.89; only the EDF+ recording field says 1989.
    assert f["startDateTime"] == "1989-04-24T16:13:00"
    assert f["recording"] == "Startdate 24-APR-1989 X X X"
    assert f["subjectSex"] == "F" and f["subjectName"] == "Female_33yr"
    assert f["signalLabels"][0] == "EEG Fpz-Cz" and f["signalLabels"][-1] == "Event marker"
    eeg = f["signals"][0]
    assert eeg["physicalDimension"] == "uV" and eeg["samplesPerRecord"] == 3000
    assert eeg["samplingRateHz"] == 100.0 and eeg["digitalMin"] == -2048
    assert eeg["transducer"] == "Ag-AgCl electrodes"
    # Channels run at 100 Hz and 1 Hz, so no single rate is claimed.
    assert "samplingRateHz" not in f


def test_edf_plus_annotations_channel():
    f = decode_edf(bytes.fromhex(_case("test_generator.edf")["hex"])).fields
    assert f["edfType"] == "EDF+C" and f["annotationSignalCount"] == 1
    assert f["signals"][-1]["label"] == "EDF Annotations"
    assert f["signals"][-1]["annotations"] is True
    # The annotations channel is kept out of the recording's sampling rate.
    assert f["samplingRateHz"] == 200.0


def test_legacy_edf_with_annotations_channel():
    # Reserved field is blank, so this is plain EDF, but it still carries the
    # reserved annotations label; the label alone identifies the channel.
    f = decode_edf(bytes.fromhex(_case("test_legacy.edf")["hex"])).fields
    assert f["edfType"] == "EDF" and f["annotationSignalCount"] == 1
    assert f["samplingRateHz"] == 200.0


def test_bdf_is_24_bit():
    d = decode_edf(bytes.fromhex(_case("test_generator.bdf")["hex"]))
    f = d.fields
    assert f["edfType"] == "BDF+C" and f["bitDepth"] == 24
    assert f["signals"][0]["digitalMin"] == -8388608
    assert f["signals"][0]["digitalMax"] == 8388607
    assert f["annotationSignalCount"] == 1


def test_fractional_data_record_duration():
    f = decode_edf(bytes.fromhex(_case("datarec_generator_0_5")["hex"])).fields
    assert f["dataRecordDurationSeconds"] == 0.5
    assert f["durationSeconds"] == f["numDataRecords"] * 0.5


def test_rejects_non_edf():
    assert decode_edf(b"") is None
    assert decode_edf(b"\x00" * 1024) is None
    # Right magic, but the two redundant counts disagree.
    bad = bytearray(b" " * 256)
    bad[:8] = b"0       "
    bad[184:192] = b"9999    "
    bad[252:256] = b"2   "
    assert decode_edf(bytes(bad)) is None


def test_truncated_signal_header_is_reported():
    full = bytes.fromhex(_case("SC4001E0")["hex"])
    d = decode_edf(full[:256 + 100])
    assert d is not None and d.fields["headerComplete"] is False
    assert d.fields["signalCount"] == 7
    # The columns are laid out one after another, so a short read loses the
    # later ones entirely rather than truncating each signal.
    assert len(d.fields.get("signalLabels", [])) < 7
