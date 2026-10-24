"""Decode an A3M multiple sequence alignment.

A3M is the FASTA-like alignment format HH-suite defined and the MSA format the
structure predictors read: ColabFold and the MMseqs2 server write it; AlphaFold,
Boltz, Chai, Protenix and OpenFold take it as input. The first record is the
query. In every later record an upper-case letter or `-` is a match state
aligned to a query column and a lower-case letter is an insertion relative to
the query, which carries no column; dropping the lower-case letters leaves
every row exactly the query's length.

Writers differ in the trimmings, and this reads each of them:

- ColabFold's complex MSA starts with `#L1,L2<TAB>C1,C2`: the chain lengths
  and copy counts of the concatenated query (`#231<TAB>1` for a monomer).
- The MMseqs2 / ColabFold server answers several queries in one file, blocks
  separated by a NUL byte, each with its own `>101`, `>102` query record.
- HH-suite puts `>ss_pred` / `>ss_conf` annotation rows before the query.
- A sequence may be wrapped over several lines (reformat.pl) or be one line.

Depth is the number that matters (prediction quality tracks it), so the record
count is reported for what the read covered, with `truncated` saying whether
the file went on past the read, together with the mean fraction of the query
each row covers and the rows whose aligned length does not match the query
(an MSA built for a different sequence, or a mis-staged file).
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

from .decoders import DecodedHeader, Read, register_decoder

# A ColabFold uniref MSA for a 300-residue query is 2 to 3 MB; the environmental
# one can pass 10 MB. Counting past the read is not possible, so `truncated`
# says when the count is a lower bound.
A3M_READ_BYTES = 4 * 1024 * 1024

# HH-suite annotation rows: not sequences, excluded from the depth, listed.
_ANNOTATION = {b"ss_pred", b"ss_conf", b"ss_dssp", b"aa_dssp", b"aa_pred", b"sa_dssp"}
_LOWER = bytes(range(ord("a"), ord("z") + 1))
_INSERT = _LOWER + b"."  # a `.` is an insertion-column gap in a2m-style rows
_SEQ_CHARS = frozenset(b"ABCDEFGHIJKLMNOPQRSTUVWXYZ" + _LOWER + b"-.*")
_QUERIES_MAX = 32
_INT_LIST = re.compile(rb"^\d+(,\d+)*$")
_HEADER_PREVIEW = 120

# Where a hit came from, read off its name. ColabFold's uniref.a3m is
# UniRef100_; its environmental MSA mixes MGnify proteins and SRA assemblies;
# jackhmmer / HHblits outputs carry UniProt `tr|` / `sp|` and UniRef90 names.
_DATABASES: List[Tuple[re.Pattern, str]] = [
    (re.compile(rb"^UniRef(100|90|50)_"), "UniRef%s"),
    (re.compile(rb"^(?:tr|sp)\|"), "UniProtKB"),
    (re.compile(rb"^UPI[0-9A-F]{10}"), "UniParc"),
    (re.compile(rb"^MGY[PS]\d"), "MGnify"),
    (re.compile(rb"^[SED]RR\d+"), "SRA"),
]


def _chain_spec(line: bytes) -> Tuple[Optional[List[int]], Optional[List[int]]]:
    """ColabFold's `#101,130<TAB>1,1` line: chain lengths and cardinalities, or
    (None, None) when the comment is something else."""
    parts = line[1:].strip().split()
    if len(parts) != 2 or not _INT_LIST.match(parts[0]) or not _INT_LIST.match(parts[1]):
        return None, None
    lengths = [int(x) for x in parts[0].split(b",")]
    cards = [int(x) for x in parts[1].split(b",")]
    if len(lengths) != len(cards) or min(lengths) <= 0 or min(cards) <= 0 or max(lengths + cards) >= 10 ** 9:
        return None, None
    return lengths, cards


def _records(block: bytes):
    """(header, sequence) pairs of one block, plus its leading `#` comment
    lines; None if the block is not alignment text."""
    comments: List[bytes] = []
    records: List[Tuple[bytes, bytearray]] = []
    for line in block.split(b"\n"):
        line = line.rstrip(b" \t")
        if not line:
            continue
        if line[:1] == b">":
            records.append((line[1:], bytearray()))
        elif not records:
            if line[:1] == b"#":
                comments.append(line)
                continue
            return None  # sequence text before any header: not A3M
        else:
            records[-1][1].extend(line)
    return comments, records


def _has_insertions(seq: bytes) -> bool:
    return len(seq.translate(None, _LOWER)) != len(seq)


def _database(name: bytes) -> Optional[str]:
    for pattern, tag in _DATABASES:
        m = pattern.match(name)
        if m:
            return tag % m.group(1).decode("ascii") if "%s" in tag else tag
    return None


def decode_a3m(data: bytes) -> Optional[DecodedHeader]:
    if not data:
        return None
    truncated = len(data) >= A3M_READ_BYTES
    blocks = [b for b in data.replace(b"\r", b"").split(b"\x00") if b.strip()]
    if not blocks:
        return None

    parsed = []
    for i, block in enumerate(blocks):
        r = _records(block)
        if r is None:
            return None
        comments, records = r
        if truncated and i == len(blocks) - 1 and records:
            records.pop()  # the read ended inside this record
        parsed.append((comments, records))

    chain_lengths = chain_cards = None
    for line in parsed[0][0]:
        chain_lengths, chain_cards = _chain_spec(line)
        if chain_lengths:
            break

    queries = []
    first = None
    for comments, records in parsed:
        annotations: List[str] = []
        query = None
        hits = []
        for header, seq in records:
            name = header.split(None, 1)[0] if header.strip() else b""
            if name in _ANNOTATION:
                annotations.append(name.decode("ascii"))
                continue
            if any(c not in _SEQ_CHARS for c in seq):
                return None
            if query is None:
                query = (name, header, bytes(seq))
            else:
                hits.append((name, header, bytes(seq)))
        if query is None:
            continue
        query_len = len(query[2].translate(None, _INSERT))
        if query_len == 0:
            continue
        if first is None:
            first = (query, hits, annotations, query_len)
        queries.append({
            "name": query[0].decode("utf-8", "replace") or None,
            "length": query_len,
            "sequenceCount": 1 + len(hits),
        })
    if first is None:
        return None

    query, hits, annotations, query_len = first
    with_insertions = 1 if _has_insertions(query[2]) else 0
    mismatches = 0
    coverage = 0.0
    covered = 0
    databases: Dict[str, int] = {}
    for name, header, seq in hits:
        match = seq.translate(None, _INSERT)
        if _has_insertions(seq):
            with_insertions += 1
        if len(match) != query_len:
            mismatches += 1
        else:
            coverage += 1.0 - match.count(b"-") / query_len
            covered += 1
        tag = _database(name)
        if tag:
            databases[tag] = databases.get(tag, 0) + 1

    first_hit = hits[0][1].decode("utf-8", "replace")[:_HEADER_PREVIEW] if hits else None
    fields = {
        "queryName": queries[0]["name"],
        "queryLength": query_len,
        "chainLengths": chain_lengths,
        "chainCardinalities": chain_cards,
        "sequenceCount": 1 + len(hits),
        "meanQueryCoverage": round(coverage / covered, 4) if covered else None,
        "sequencesWithInsertions": with_insertions,
        "lengthMismatches": mismatches,
        "annotationRows": annotations or None,
        "databases": databases or None,
        "firstHitHeader": first_hit,
        "queryCount": len(queries),
        "queries": queries[:_QUERIES_MAX],
        "truncated": truncated,
    }

    chains = ""
    if chain_lengths and sum(chain_cards) > 1:
        expanded = [str(n) for n, c in zip(chain_lengths, chain_cards) for _ in range(c)]
        chains = " (%d chains: %s)" % (len(expanded), "+".join(expanded))
    depth = "%d%s sequences" % (1 + len(hits), "+" if truncated else "")
    head = "A3M MSA, " if len(queries) == 1 else "A3M MSA, %d queries; first: " % len(queries)
    return DecodedHeader(
        format="a3m",
        summary="%squery %d aa%s, %s" % (head, query_len, chains, depth),
        fields=fields,
    )


register_decoder("a3m", decode_a3m, read=Read(leading=A3M_READ_BYTES))  # depth needs the rows
