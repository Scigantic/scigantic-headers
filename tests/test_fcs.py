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


@pytest.mark.parametrize("version", ["FCS2.0", "FCS3.0", "FCS3.1", "FCS3.2"])
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


def test_truncated_header_returns_none():
    blob = make_fcs()
    # Chop the file inside the HEADER: no TEXT has been read at all.
    assert decode_fcs(blob[:40]) is None


def test_text_longer_than_read_decodes_what_was_read():
    # A TEXT segment that runs past the bytes given keeps the keywords that
    # were read, drops the pair the cut lands in, and says it is truncated.
    channels = [("P%d" % i, "marker%d" % i) for i in range(40)]
    blob = make_fcs(channels=channels, tot=5)
    cut = blob[:blob.index(b"marker20") + 3]
    h = decode_fcs(cut)
    assert h is not None and h.fields["textTruncated"] is True
    assert h.fields["numParameters"] == 40
    chans = h.fields["channels"]
    assert [c["name"] for c in chans] == ["P%d" % i for i in range(21)]
    # P20's label is the value the cut lands in: dropped, never half-read.
    assert [c.get("label") for c in chans] == ["marker%d" % i for i in range(20)] + [None]
    assert "textTruncated" not in decode_fcs(blob).fields


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

def test_matches_flowio_channel_and_event_counts():
    # Skip only this test without flowio. A module-level importorskip here
    # skipped every test in the file.
    flowio = pytest.importorskip("flowio", reason="flowio not installed")
    import io

    blob = make_fcs(tot=7)
    fd = flowio.FlowData(io.BytesIO(blob))
    h = decode_fcs(blob)
    assert h.fields["numParameters"] == fd.channel_count
    assert h.fields["numEvents"] == fd.event_count
    # channel short names ($PnN) agree
    ref_names = [fd.channels[i]["pnn"] for i in range(1, fd.channel_count + 1)]
    assert [c["name"] for c in h.fields["channels"]] == ref_names


# ── TEXT segment edge cases: legacy empty values, escapes, encodings, STEXT ───
#
# These build the TEXT bytes by hand so the exact delimiter runs are under test.

def raw_fcs(text, version=b"FCS3.1", tail=b""):
    """HEADER + the given TEXT bytes (+ tail), with no DATA segment."""
    text_begin = 58
    text_end = text_begin + len(text) - 1
    header = version + b"    " + b"%8d%8d" % (text_begin, text_end) + b"%8d" % 0 * 4
    assert len(header) == 58
    return header + text + tail


def test_legacy_empty_value_does_not_shift_later_pairs():
    # BD CellQuest (FCS 2.0) writes an empty value as a doubled delimiter:
    # `KEY\\NEXT\value\`. Read as an escaped literal, KEY's value would swallow
    # NEXT and every later keyword would become a value. Position decides: a
    # doubled delimiter straight after a keyword is an empty value.
    text = (
        b"\\$PAR\\2\\$TOT\\5\\$DATATYPE\\I\\$MODE\\L"
        b"\\&5Data File Prefix Part #1\\\\&6Data File Prefix Part #2\\\\"
        b"&8Acquisition Doc.\\LYMPH SUBSET ACQ\\$CYT\\FACScan\\$P1N\\FSC-H\\$P2N\\SSC-H"
        b"\\&13Analysis Doc.\\\\"
    )
    h = decode_fcs(raw_fcs(text, version=b"FCS2.0"))
    assert h.fields["cytometer"] == "FACScan"
    assert [c["name"] for c in h.fields["channels"]] == ["FSC-H", "SSC-H"]
    assert h.fields["numEvents"] == 5


def test_doubled_delimiter_inside_value_is_a_literal():
    # Spec escaping (3.2.7): `$P3F/488//10/` is the filter "488/10". A value may
    # also end in a literal (`C:/data/`, a run of three), and four in a row are
    # two literals (a Miltenyi writer's `561////10 nm`).
    text = (
        b"/$PAR/2/$TOT/1/$P1N/FSC/$P1S/488//10/$P2N/SSC/$P2S/561////10 nm/"
        b"$CYT/C://data///$MODE/L/"
    )
    h = decode_fcs(raw_fcs(text))
    assert h.fields["channels"][0]["label"] == "488/10"
    assert h.fields["channels"][1]["label"] == "561//10 nm"
    assert h.fields["cytometer"] == "C:/data/"
    assert h.fields["mode"] == "L"


def test_legacy_empty_value_and_escaped_literal_in_one_text():
    text = b"/$PAR/1/$TOT/3/NOTE//$P1N/FSC/$P1S/a//b/$CYT/Cyto/"
    h = decode_fcs(raw_fcs(text))
    assert h.fields["channels"] == [{"name": "FSC", "label": "a/b"}]
    assert h.fields["cytometer"] == "Cyto"
    assert h.fields["numEvents"] == 3


def test_split_round_trips_random_texts():
    # Property: any keywords without the delimiter and any values that do not
    # start with it (the spec's rule), empty values written the legacy way and
    # literals doubled, come back exactly.
    import random

    from scigantic_headers.fcs import _split

    rng = random.Random(2026)
    for delim in (b"/", b"\\", b"|", b"\x0c"):
        alphabet = [b"a", b"B", b"1", b" ", b"$", delim]
        for _ in range(300):
            pairs = []
            for _ in range(rng.randint(1, 12)):
                key = b"".join(rng.choice(alphabet[:-1]) for _ in range(rng.randint(1, 6)))
                if rng.random() < 0.3:
                    val = b""
                else:
                    val = rng.choice(alphabet[:-1]) + b"".join(
                        rng.choice(alphabet) for _ in range(rng.randint(0, 8))
                    )
                pairs.append((key, val))
            body = b"".join(k + delim + v.replace(delim, delim * 2) + delim for k, v in pairs)
            got = _split(body, delim[0])
            assert got == [x for kv in pairs for x in kv], (delim, pairs, got)


def test_latin1_value_is_kept_not_replaced():
    # FCS 2.0/3.0 TEXT from Windows instruments is not always UTF-8.
    text = b"/$PAR/1/$TOT/1/$P1N/FSC/$CYT/CELLQuest\xaa 3.3/"
    h = decode_fcs(raw_fcs(text, version=b"FCS2.0"))
    assert h.fields["cytometer"] == "CELLQuest\xaa 3.3"
    assert h.fields["cytometer"].encode("latin-1") == b"CELLQuest\xaa 3.3"


def test_missing_par_still_lists_channels_from_pnn():
    text = b"/$TOT/4/$P1N/FSC/$P2N/SSC/$P3N/FL1/$DATATYPE/F/"
    h = decode_fcs(raw_fcs(text))
    assert h.fields["numParameters"] is None
    assert [c["name"] for c in h.fields["channels"]] == ["FSC", "SSC", "FL1"]
    assert "? parameters" in h.summary


def test_channel_without_pnn_falls_back_to_pns():
    text = b"/$PAR/2/$TOT/1/$P1S/Forward Scatter/$P2N/SSC/$P2S/Side/"
    h = decode_fcs(raw_fcs(text))
    assert h.fields["channels"] == [{"name": "Forward Scatter"}, {"name": "SSC", "label": "Side"}]


def _with_stext(stext, leading_delim=True, offset_from_text_end=0, begin_override=None, end_override=None):
    """A file whose primary TEXT points at a supplemental TEXT placed after it."""
    primary = b"/$PAR/2/$TOT/1/$P1N/FSC/$P2N/SSC/$BEGINSTEXT/%08d/$ENDSTEXT/%08d/"
    length = len(primary % (0, 0))
    begin = 58 + length + offset_from_text_end
    seg = (b"/" if leading_delim else b"") + stext
    end = begin + len(seg) - 1
    text = primary % (begin_override if begin_override is not None else begin,
                      end_override if end_override is not None else end)
    return raw_fcs(text, tail=b"\x00" * offset_from_text_end + seg)


@pytest.mark.parametrize("leading_delim", [True, False])
def test_supplemental_text_is_read_and_primary_wins(leading_delim):
    # Optional keywords may live in the supplemental TEXT. $PAR is given twice;
    # the primary value (2) wins over the supplemental one (99).
    h = decode_fcs(_with_stext(b"$CYT/Aurora/$P1S/Forward/$PAR/99/", leading_delim=leading_delim, offset_from_text_end=16))
    assert h.fields["cytometer"] == "Aurora"
    assert h.fields["channels"][0] == {"name": "FSC", "label": "Forward"}
    assert h.fields["numParameters"] == 2


def test_supplemental_text_outside_the_read_is_ignored():
    blob = _with_stext(b"$CYT/Aurora/", offset_from_text_end=16)
    cut = blob[:blob.index(b"\x00" * 16)]  # everything up to the STEXT
    h = decode_fcs(cut)
    assert h.fields["numParameters"] == 2
    assert "cytometer" not in h.fields


def test_supplemental_text_with_bad_offsets_is_ignored():
    h = decode_fcs(_with_stext(b"$CYT/Aurora/", begin_override=500, end_override=400))
    assert h is not None and "cytometer" not in h.fields
    h = decode_fcs(_with_stext(b"$CYT/Aurora/", begin_override=99999990, end_override=99999999))
    assert h is not None and "cytometer" not in h.fields


def test_next_dataset_sets_multiple_datasets():
    # A Coulter .lmd carries an FCS 2.0 and an FCS 3.0 dataset; $NEXTDATA points
    # at the second. A single-dataset file has no flag at all.
    assert decode_fcs(make_fcs(extra=[("$NEXTDATA", "363602")])).fields["multipleDatasets"] is True
    assert "multipleDatasets" not in decode_fcs(make_fcs()).fields


def test_keyword_values_are_normalized():
    h = decode_fcs(make_fcs(extra=[("$DATATYPE", " f "), ("$CYT", "  Attune NxT  "), ("$TOT", " 42 ")]))
    assert h.fields["dataType"] == "float32"
    assert h.fields["cytometer"] == "Attune NxT"
    assert h.fields["numEvents"] == 42


def test_header_offset_with_trailing_byte():
    # Some LMD exports leave a stray byte after the digits of an offset field.
    blob = bytearray(make_fcs())
    assert blob[10:18] == b"      58"
    blob[10:18] = b"     58\x00"
    h = decode_fcs(bytes(blob))
    assert h is not None and h.fields["numParameters"] == 4


def test_leading_delimiter_run_in_body_is_skipped():
    # A TEXT that opens with two delimiters (an empty first keyword) still
    # yields the pairs after it.
    h = decode_fcs(raw_fcs(b"//$PAR/1/$TOT/9/$P1N/FSC/"))
    assert h.fields["numEvents"] == 9
    assert h.fields["channels"] == [{"name": "FSC"}]
