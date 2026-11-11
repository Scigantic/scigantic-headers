"""Every family's golden case decodes to its own instrument through the full
registry, with all families' sniffers registered together. Catches one
family's sniffer claiming another family's export."""

import glob
import json
import os

import pytest

from scigantic_headers import decode_bytes

FX = os.path.join(os.path.dirname(__file__), "fixtures", "instruments")


def _cases():
    out = []
    for f in sorted(glob.glob(os.path.join(FX, "*-cases.json"))):
        fam = os.path.basename(f)[: -len("-cases.json")]
        with open(f) as fh:
            for c in json.load(fh)["cases"]:
                path = os.path.join(FX, c["file"])
                if not os.path.exists(path):
                    path = os.path.join(FX, fam, c["file"])
                out.append(pytest.param(c.get("key") or os.path.basename(c["file"]), path,
                                        c["expected"]["format"], id="%s:%s" % (fam, c["name"])))
    return out


@pytest.mark.parametrize("key,path,fmt", _cases())
def test_case_decodes_to_its_own_instrument(key, path, fmt):
    with open(path, "rb") as fh:
        hdr = decode_bytes(key, fh.read())
    assert hdr is not None and hdr.format == fmt
