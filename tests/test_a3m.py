"""A3M multiple-sequence-alignment decoding."""

import json

import pytest

from scigantic_headers import decode_a3m, decode_bytes, has_decoder_for, is_decodable, read_for
from scigantic_headers.a3m import A3M_READ_BYTES

# A ColabFold complex MSA: the `#` chain line, numeric query name, MMseqs2
# tab-separated hit headers, a lower-case insertion and a gapped row.
COLABFOLD = (
    b"#5,4\t1,1\n"
    b">101\n"
    b"MKTAYWGHE\n"
    b">UniRef100_A0A001\t75\t0.590\t6.662E-14\t2\t9\t9\t1419\t1427\t1504\n"
    b"MKSAYWGHE\n"
    b">UniRef100_A0A002\t60\t0.400\t1.0E-10\t1\t9\t9\t1\t9\t20\n"
    b"MKTAYagWG--\n"
    b">MGYP000000000001\t50\t0.3\t1.0E-5\t1\t5\t9\t1\t5\t5\n"
    b"MKTAY----\n"
)


def test_colabfold_complex():
    h = decode_a3m(COLABFOLD)
    assert h.format == "a3m"
    f = h.fields
    assert f["queryName"] == "101"
    assert f["queryLength"] == 9
    assert f["chainLengths"] == [5, 4]
    assert f["chainCardinalities"] == [1, 1]
    assert f["sequenceCount"] == 4  # the query counts, as in every depth plot
    assert f["sequencesWithInsertions"] == 1
    assert f["lengthMismatches"] == 0
    # coverage: 9/9, 7/9, 5/9 -> mean 0.7778
    assert f["meanQueryCoverage"] == pytest.approx(0.7778)
    assert f["databases"] == {"UniRef100": 2, "MGnify": 1}
    assert f["firstHitHeader"].startswith("UniRef100_A0A001\t75\t0.590")
    assert f["queryCount"] == 1
    assert f["queries"] == [{"name": "101", "length": 9, "sequenceCount": 4}]
    assert f["truncated"] is False
    assert h.summary == "A3M MSA, query 9 aa (2 chains: 5+4), 4 sequences"


def test_homodimer_chain_line_expands_cardinality():
    h = decode_a3m(b"#9\t2\n>101\nMKTAYWGHE\n")
    assert h.fields["chainLengths"] == [9]
    assert h.fields["chainCardinalities"] == [2]
    assert "(2 chains: 9+9)" in h.summary


def test_monomer_chain_line_is_not_announced():
    h = decode_a3m(b"#9\t1\n>101\nMKTAYWGHE\n")
    assert h.fields["chainLengths"] == [9]
    assert h.summary == "A3M MSA, query 9 aa, 1 sequences"


def test_server_multi_query_nul_blocks():
    # The MMseqs2 / ColabFold server answers two queries in one file, blocks
    # separated by NUL, the file ending in one.
    data = (
        b">101\nMKTAY\n>UniRef100_X\t1\t1\t1\t1\t5\t5\t1\t5\t5\nMKSAY\n\x00"
        b">102\nWGHEQLN\n>UniRef100_Y\t1\t1\t1\t1\t7\t7\t1\t7\t7\nWGHE---\n"
        b">UniRef100_Z\t1\t1\t1\t1\t7\t7\t1\t7\t7\nWGHEQLN\n\x00"
    )
    h = decode_a3m(data)
    f = h.fields
    assert f["queryCount"] == 2
    assert f["queries"] == [
        {"name": "101", "length": 5, "sequenceCount": 2},
        {"name": "102", "length": 7, "sequenceCount": 3},
    ]
    assert f["queryName"] == "101" and f["sequenceCount"] == 2  # top-level fields describe the first block
    assert h.summary == "A3M MSA, 2 queries; first: query 5 aa, 2 sequences"


def test_hhsuite_annotation_rows_wrapped_lines_and_crlf():
    data = (
        b">ss_pred\r\nCCHHHHH\r\n"
        b">ss_conf\r\n9876543\r\n"
        b">d1abca_ b.1.1.1 (A:) some domain\r\n"
        b"MKTA\r\nYWG\r\n"
        b">tr|Q9XYZ1|Q9XYZ1_HUMAN Protein OS=Homo sapiens\r\n"
        b"MKsaTAY\r\nWG\r\n"
    )
    h = decode_a3m(data)
    f = h.fields
    assert f["annotationRows"] == ["ss_pred", "ss_conf"]
    assert f["queryName"] == "d1abca_"
    assert f["queryLength"] == 7
    assert f["sequenceCount"] == 2  # annotation rows are not sequences
    assert f["sequencesWithInsertions"] == 1
    assert f["lengthMismatches"] == 0
    assert f["meanQueryCoverage"] == 1.0
    assert f["databases"] == {"UniProtKB": 1}


def test_length_mismatch_is_counted_not_hidden():
    # An MSA whose rows do not fit the query: the realignment step was skipped,
    # or the file belongs to another sequence. Mismatched rows are counted and
    # left out of the coverage mean.
    data = b">query\nMKTAYWGHE\n>hit1\nMKTAYWGHE\n>hit2\nMKTAY\n>hit3\nMKTAYWGHEQQ\n"
    f = decode_a3m(data).fields
    assert f["sequenceCount"] == 4
    assert f["lengthMismatches"] == 2
    assert f["meanQueryCoverage"] == 1.0


def test_query_only_has_no_coverage():
    f = decode_a3m(b">query \nGSHMTGKSVDMY\n").fields
    assert f["queryName"] == "query"
    assert f["sequenceCount"] == 1
    assert f["meanQueryCoverage"] is None
    assert f["firstHitHeader"] is None
    assert f["databases"] is None


def test_dot_is_an_insertion_gap_and_star_a_match_state():
    f = decode_a3m(b">q\nMKTAY*\n>h\nMK..TAY*\n").fields
    assert f["queryLength"] == 6
    assert f["lengthMismatches"] == 0
    assert f["sequencesWithInsertions"] == 0


def test_truncated_read_drops_the_cut_record_and_says_so():
    rows = b">q\n" + b"A" * 60 + b"\n"
    hit = b">UniRef100_H\n" + b"A" * 60 + b"\n"
    n = (A3M_READ_BYTES - len(rows)) // len(hit) + 2
    data = (rows + hit * n)[:A3M_READ_BYTES]
    assert len(data) == A3M_READ_BYTES
    f = decode_a3m(data).fields
    assert f["truncated"] is True
    # Every complete record counts; the one the read cut through does not.
    whole = data.count(b"\n>")
    assert f["sequenceCount"] in (whole, whole + 1)
    assert f["lengthMismatches"] == 0
    assert decode_a3m(data).summary.endswith("+ sequences")


def test_rejects_non_alignment_text():
    assert decode_a3m(b"") is None
    assert decode_a3m(b"MKTAYWGHE\n>q\nMKTAY\n") is None  # sequence before any header
    assert decode_a3m(b">q\nMKT AY WGHE\n") is None  # not residue letters
    assert decode_a3m(b">q\n123456\n") is None
    assert decode_a3m(b"# just a comment\n") is None
    assert decode_a3m(b">\n\n") is None  # header with no residues
    assert decode_a3m(b"\x00\x00") is None


def test_odd_chain_lines_are_ignored_not_fatal():
    assert decode_a3m(b"#hello world\n>q\nMKT\n").fields["chainLengths"] is None
    assert decode_a3m(b"#1,2\t1\n>q\nMKT\n").fields["chainLengths"] is None  # unequal lists
    assert decode_a3m(b"#0\t1\n>q\nMKT\n").fields["chainLengths"] is None


def test_dispatch_and_read_size():
    assert has_decoder_for("msa/uniref.a3m")
    assert is_decodable("msa/uniref.a3m.gz")  # the reader inflates .gz before dispatch
    assert decode_bytes("x.a3m", COLABFOLD).format == "a3m"
    assert read_for("x.a3m").leading == A3M_READ_BYTES


def test_json_safe():
    json.dumps(decode_a3m(COLABFOLD).to_dict(), allow_nan=False)
