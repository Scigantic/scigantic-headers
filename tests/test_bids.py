"""BIDS sidecar decoding.

Golden: tests/fixtures/bids-cases.json, real files from the
bids-standard/bids-examples corpus. The whole corpus was swept when the
fixture was made: 5853 `.json` and `.tsv` files decoded with no exception, and
the files that declined are the headerless derivative tables (motion,
timeseries, mixing, physio) and non-core derivative JSON this decoder is not
meant to claim.
"""

import json
import os

from scigantic_headers import decode_bids_json, decode_bids_tsv, decode_bytes, has_decoder_for
from scigantic_headers.bids import BIDS_READ_BYTES
from scigantic_headers.decoders import read_for

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "bids-cases.json")

with open(FIX) as _fh:
    CASES = json.load(_fh)["cases"]


def _case(name):
    return bytes.fromhex(next(c for c in CASES if name in c["name"])["hex"])


def test_extensions_registered():
    assert has_decoder_for("sub-01/func/sub-01_bold.json")
    assert has_decoder_for("participants.tsv")
    assert read_for("x.json").leading == BIDS_READ_BYTES
    # A name-dispatched decoder still wins over the .json extension.
    assert decode_bytes("zarr.json", b'{"zarr_format": 3, "node_type": "group"}') is not None


def test_golden_fixture_matches():
    assert len(CASES) >= 11
    for case in CASES:
        d = decode_bytes(case["key"], bytes.fromhex(case["hex"]))
        assert d is not None and d.to_dict() == case["expected"], case["name"]


def test_dataset_description():
    d = decode_bids_json(_case("dataset_description"))
    assert d.format == "bids-dataset"
    assert d.summary == "BIDS dataset, MNE sample dataset, BIDS 1.4.0, raw"
    f = d.fields
    assert f["BIDSVersion"] == "1.4.0" and f["License"] == "CC0"
    assert f["authorsCount"] == 2 and f["fundingCount"] == 11
    assert f["DatasetDOI"] == "10.18112/openneuro.ds000248.v1.2.4"


def test_mri_sidecar():
    f = decode_bids_json(_case("fMRI BOLD")).fields
    assert f["bidsModality"] == "mri" and f["TaskName"] == "Rest"
    assert f["RepetitionTime"] == 4.0 and f["EchoTime"] == 0.026
    # SliceTiming is 40 numbers; the count is reported, not the list.
    assert f["sliceCount"] == 40 and "SliceTiming" not in f


def test_ieeg_sidecar():
    d = decode_bids_json(_case("iEEG sidecar"))
    assert d.summary == "BIDS sidecar, ieeg, task FilteredSpeech, 1000 Hz"
    f = d.fields
    # iEEG keys win over the EEG keys the same file also carries.
    assert f["bidsModality"] == "ieeg"
    assert f["iEEGReference"] == "common average" and f["EEGChannelCount"] == 64
    assert f["RecordingDuration"] == 1204.041


def test_coordinate_system_sidecar():
    f = decode_bids_json(_case("coordinate system")).fields
    assert f["bidsModality"] == "coordsystem"
    assert f["MEGCoordinateSystem"] == "ElektaNeuromag"
    assert f["HeadCoilCoordinateUnits"] == "mm"


def test_data_dictionary():
    d = decode_bids_json(_case("data dictionary"))
    assert d.format == "bids-dictionary" and d.fields["columnCount"] == 4
    assert "Description" in next(iter(d.fields["columnDescriptions"].values()))


def test_events_table():
    d = decode_bids_tsv(_case("events with trial_type"))
    assert d.format == "bids-table"
    f = d.fields
    assert f["bidsFileType"] == "events" and f["rowsRead"] == 158
    assert f["rowsComplete"] is True and f["columns"][0] == "onset"
    assert f["onsetFirstSeconds"] == 0.061 and f["onsetLastSeconds"] == 600.409
    assert f["trialTypesDistinct"] == 4 and sum(f["trialTypes"].values()) == 158


def test_channels_and_electrodes_tables():
    channels = decode_bids_tsv(_case("channels table")).fields
    assert channels["bidsFileType"] == "channels" and channels["samplingRateHz"] == 1000.0
    assert channels["channelTypesDistinct"] == 1
    electrodes = decode_bids_tsv(_case("electrodes table")).fields
    assert electrodes["bidsFileType"] == "electrodes" and electrodes["electrodeCount"] == 64


def test_byte_order_mark_does_not_hide_the_columns():
    f = decode_bids_tsv(_case("UTF-8 BOM")).fields
    assert f["bidsFileType"] == "scans" and f["columns"] == ["filename", "acq_time"]


def test_truncated_table_says_so():
    data = _case("events with trial_type")
    assert len(data) < BIDS_READ_BYTES
    # Simulate a file longer than the bounded read.
    padded = data + b"\n" + b"0.0\t0.0\tpad\tn/a\tn/a\tn/a\t0\t0" * 20000
    f = decode_bids_tsv(padded[:BIDS_READ_BYTES]).fields
    assert f["rowsComplete"] is False
    assert "read" in decode_bids_tsv(padded[:BIDS_READ_BYTES]).summary


def test_declines_non_bids():
    assert decode_bids_json(b"") is None
    assert decode_bids_json(b"[1, 2, 3]") is None
    assert decode_bids_json(b'{"name": "pkg", "version": "1.0.0"}') is None
    assert decode_bids_tsv(b"a\tb\tc\n1\t2\t3\n") is None
    assert decode_bids_tsv(b"") is None
