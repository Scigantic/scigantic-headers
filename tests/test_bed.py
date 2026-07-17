"""BED interval-file decoding."""

import json

from scigantic_headers import decode_bed
from scigantic_headers.decoders import decode_bytes

BED6 = (
    b"track name=test\n"
    b"chr1\t1000\t2000\tfeat1\t500\t+\n"
    b"chr1\t3000\t4000\tfeat2\t900\t-\n"
)


def test_bed6_with_track_line():
    h = decode_bed(BED6)
    assert h.format == "bed"
    assert h.fields["numColumns"] == 6
    assert h.fields["bedType"] == "BED6"
    assert h.fields["hasTrackLine"] is True
    assert h.fields["firstChrom"] == "chr1"
    assert "strand" in h.fields["columns"]


def test_bed3():
    h = decode_bed(b"chr2\t0\t100\n")
    assert h.fields["numColumns"] == 3
    assert h.fields["hasTrackLine"] is False


def test_rejects_non_bed():
    assert decode_bed(b"chr1\tnotanumber\t2000\n") is None  # start/end not integers
    assert decode_bed(b"only one column here\n") is None
    assert decode_bed(b"") is None


def test_dispatch():
    assert decode_bytes("x.bed", BED6).format == "bed"


def test_json_safe():
    json.dumps(decode_bed(BED6).to_dict(), allow_nan=False)
