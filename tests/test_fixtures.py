"""Golden-fixtures test — the drift guard.

Three decoders must agree on MRC: this library, the backend TS twin
(headerDecoders.ts), and the notebook (scigantic_empiar.parse_mrc_header).
fixtures/mrc-cases.json is the single source of truth: hex input -> expected
decode. This suite asserts the Python decoder matches it; the TS suite should
load the same file and assert the same, so a change to one decoder that the
other doesn't match fails a test rather than silently diverging.
"""

import json
import os

import pytest

from scigantic_headers import decode_bytes

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "mrc-cases.json")

with open(FIXTURES) as _fh:
    _CASES = json.load(_fh)["cases"]


@pytest.mark.parametrize("case", _CASES, ids=[c["name"] for c in _CASES])
def test_mrc_matches_golden_fixture(case):
    decoded = decode_bytes("fixture.mrc", bytes.fromhex(case["hex"]))
    assert decoded is not None, case["name"]
    assert decoded.to_dict() == case["expected"]


def test_fixture_file_is_nonempty():
    # Guard the guard: a truncated/empty fixtures file would make the parametrized
    # test vacuously pass, so assert it actually carries cases.
    assert len(_CASES) >= 4
