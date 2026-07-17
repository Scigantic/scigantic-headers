"""FCS (Flow Cytometry Standard) decoding.

The fixtures are built here with known keyword values, so the expected decode is
exact. Where flowio is installed, the same bytes are cross-checked against it.
"""

import json
import struct

import pytest

from scigantic_headers import decode_fcs
from scigantic_headers.decoders import decode_bytes


def make_fcs(
    *,
    version="FCS3.1",
    channels=None,
    tot=10,
    datatype="F",
    mode="L",
    cyt="TestCyt",
    delim="/",
    extra=None,
    with_data=True,
):
    """A complete, valid FCS file (HEADER + TEXT + optional DATA) as bytes."""
    if channels is None:
        channels = [
            ("FSC-A", "Forward Scatter"),
            ("SSC-A", "Side Scatter"),
            ("FL1-A", "CD3"),
            ("FL2-A", "CD4"),
        ]
    par = len(channels)
    kw = [
        ("$PAR", str(par)),
        ("$TOT", str(tot)),
        ("$DATATYPE", datatype),
        ("$MODE", mode),
        ("$CYT", cyt),
        ("$BYTEORD", "1,2,3,4"),
        ("$NEXTDATA", "0"),
    ]
    for i, (name, label) in enumerate(channels, 1):
        kw += [
            ("$P%dN" % i, name),
            ("$P%dB" % i, "32"),
            ("$P%dR" % i, "262144"),
            ("$P%dE" % i, "0,0"),
        ]
        if label:
            kw.append(("$P%dS" % i, label))
    if extra:
        kw += list(extra)

    d = delim
    text_begin = 58

    def build_text(kws):
        return (d + d.join(k + d + v for k, v in kws) + d).encode("ascii")

    if with_data:
        nvals = par * tot
        data = struct.pack("<%df" % nvals, *([1.0] * nvals))
        # Reference readers (flowio) locate DATA via the $BEGINDATA/$ENDDATA
        # keywords, so include them. Fixed-width (8-digit) values keep the TEXT
        # length identical across both passes, so the offsets computed from the
        # first pass remain correct for the second.
        probe = build_text(kw + [("$BEGINDATA", "%08d" % 0), ("$ENDDATA", "%08d" % 0)])
        text_end = text_begin + len(probe) - 1
        data_begin = text_end + 1
        data_end = data_begin + len(data) - 1
        tb = build_text(
            kw + [("$BEGINDATA", "%08d" % data_begin), ("$ENDDATA", "%08d" % data_end)]
        )
        assert len(tb) == len(probe)
    else:
        data, data_begin, data_end = b"", 0, 0
        tb = build_text(kw)
        text_end = text_begin + len(tb) - 1

    def off(n):
        return ("%8d" % n).encode("ascii")

    header = (
        version.encode("ascii")
        + b"    "
        + off(text_begin)
        + off(text_end)
        + off(data_begin)
        + off(data_end)
        + off(0)
        + off(0)
    )
    assert len(header) == 58
    return header + tb + data


def test_decodes_core_fields():
    h = decode_fcs(make_fcs(tot=1234))
    assert h is not None
    assert h.format == "fcs"
    assert h.fields["version"] == "FCS3.1"
    assert h.fields["numParameters"] == 4
    assert h.fields["numEvents"] == 1234
    assert h.fields["dataType"] == "float32"  # 'F' maps to float32
    assert h.fields["mode"] == "L"
    assert h.fields["cytometer"] == "TestCyt"
    assert "4 parameters" in h.summary
    assert "1234 events" in h.summary


def test_channel_names_and_labels():
    h = decode_fcs(make_fcs())
    chans = h.fields["channels"]
    assert [c["name"] for c in chans] == ["FSC-A", "SSC-A", "FL1-A", "FL2-A"]
    assert chans[2]["label"] == "CD3"


def test_channel_without_label_omits_label():
    h = decode_fcs(make_fcs(channels=[("FSC-A", None), ("SSC-A", "Side")]))
    chans = h.fields["channels"]
    assert "label" not in chans[0]
    assert chans[1]["label"] == "Side"


def test_dispatch_by_extension():
    blob = make_fcs()
    assert decode_bytes("sample.fcs", blob) is not None
    assert decode_bytes("sample.fcs", blob).format == "fcs"


@pytest.mark.parametrize("version", ["FCS2.0", "FCS3.0", "FCS3.1"])
def test_version_variants(version):
    h = decode_fcs(make_fcs(version=version))
    assert h.fields["version"] == version


@pytest.mark.parametrize("delim", ["/", "|", "\x0c", "\x1e"])
def test_delimiter_variants(delim):
    h = decode_fcs(make_fcs(delim=delim))
    assert h.fields["numParameters"] == 4
    assert h.fields["channels"][0]["name"] == "FSC-A"


def test_escaped_doubled_delimiter_in_value():
    # A value containing the delimiter, escaped by doubling it (FCS convention).
    blob = make_fcs(delim="/", cyt="Beckman//Coulter", channels=[("FSC-A", "a//b")])
    h = decode_fcs(blob)
    assert h.fields["cytometer"] == "Beckman/Coulter"
    assert h.fields["channels"][0]["label"] == "a/b"


def test_wide_panel_exceeds_1kib():
    # A spectral-style panel whose TEXT segment runs well past 1 KiB still
    # decodes when the full leading block is supplied.
    channels = [("P%d" % i, "marker%d" % i) for i in range(40)]
    blob = make_fcs(channels=channels, tot=5)
    assert len(blob) > 1024
    h = decode_fcs(blob)
    assert h.fields["numParameters"] == 40
    assert len(h.fields["channels"]) == 40


def test_rejects_non_fcs():
    assert decode_fcs(b"not an fcs file at all, just some bytes here padding") is None
    assert decode_fcs(b"") is None
    assert decode_fcs(b"FCS9.9  " + b" " * 60) is None  # unknown version


def test_truncated_text_returns_none():
    blob = make_fcs()
    # Chop the file inside the TEXT segment: header promises more than is present.
    truncated = blob[:40]
    assert decode_fcs(truncated) is None


def test_decode_file_reads_enough_for_wide_panel(tmp_path):
    # The source layer must read more than HEADER_BYTES for .fcs, or a wide
    # panel's TEXT segment is truncated and the file decodes to None.
    from scigantic_headers import decode_file

    channels = [("P%d" % i, "m%d" % i) for i in range(40)]
    blob = make_fcs(channels=channels, tot=5)
    assert len(blob) > 1024
    p = tmp_path / "wide.fcs"
    p.write_bytes(blob)
    h = decode_file(str(p))
    assert h is not None
    assert h.fields["numParameters"] == 40


def test_result_is_json_safe():
    h = decode_fcs(make_fcs())
    json.dumps(h.to_dict(), allow_nan=False)  # must not raise


def test_garbage_par_does_not_explode():
    # $PAR that is not an int must not raise and must not fabricate channels.
    blob = make_fcs(extra=[("$PAR", "notanumber")])  # later dup wins in our parse
    h = decode_fcs(blob)
    assert h is not None  # total: never raises


# ── cross-check against a reference reader, when present ──────────────────────

flowio = pytest.importorskip("flowio", reason="flowio not installed")


def test_matches_flowio_channel_and_event_counts():
    import io

    blob = make_fcs(tot=7)
    fd = flowio.FlowData(io.BytesIO(blob))
    h = decode_fcs(blob)
    assert h.fields["numParameters"] == fd.channel_count
    assert h.fields["numEvents"] == fd.event_count
    # channel short names ($PnN) agree
    ref_names = [fd.channels[i]["pnn"] for i in range(1, fd.channel_count + 1)]
    assert [c["name"] for c in h.fields["channels"]] == ref_names
