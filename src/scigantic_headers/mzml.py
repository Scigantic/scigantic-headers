"""Decode an mzML mass-spectrometry file's header metadata.

mzML is XML. The metadata worth surfacing sits in the document preamble: the
format version, the source file, the software, the instrument model (a cvParam
inside <instrumentConfiguration>), and the spectrum count and run start time
(attributes of <run> and <spectrumList>). The per-spectrum data comes far later
and is not read.

A truncated XML block will not go through a normal parser, so this pulls the few
attributes it needs with targeted patterns over the leading text, best-effort:
any field it cannot find is None. It reads a large leading block (see sources),
enough to reach <spectrumList>.
"""

from __future__ import annotations

import re
from typing import Optional

from .decoders import DecodedHeader, register_decoder

_NS = b"http://psi.hupo.org/ms/mzml"

_RE_VERSION = re.compile(rb'<mzML\b[^>]*\bversion="([^"]*)"')
_RE_SPECTRA = re.compile(rb'<spectrumList\b[^>]*\bcount="(\d+)"')
_RE_RUNTIME = re.compile(rb'<run\b[^>]*\bstartTimeStamp="([^"]*)"')
_RE_SOURCE = re.compile(rb'<sourceFile\b[^>]*\bname="([^"]*)"')
_RE_SOFTWARE = re.compile(rb'<software\b[^>]*\bid="([^"]*)"[^>]*\bversion="([^"]*)"')
_RE_CVNAME = re.compile(rb'<cvParam\b[^>]*\bname="([^"]*)"')


def _text(match, group=1):
    return match.group(group).decode("utf-8", "replace") if match else None


def _instrument_model(data: bytes) -> Optional[str]:
    """The instrument model is a cvParam inside <instrumentConfiguration>, ahead
    of the <componentList> (source/analyzer/detector). Take the first cvParam in
    that scope, falling back to any cvParam within the configuration."""
    start = data.find(b"<instrumentConfiguration")
    if start == -1:
        return None
    end = data.find(b"</instrumentConfiguration>", start)
    scope = data[start:end] if end != -1 else data[start:start + 8192]
    comp = scope.find(b"<componentList")
    head = scope[:comp] if comp != -1 else scope
    return _text(_RE_CVNAME.search(head)) or _text(_RE_CVNAME.search(scope))


def decode_mzml(data: bytes) -> Optional[DecodedHeader]:
    """Decode mzML header metadata from the file's leading bytes. Returns None if
    the bytes are not mzML (the <mzML> element and the mzML namespace must both
    be present)."""
    if b"<mzML" not in data or _NS not in data:
        return None

    version = _text(_RE_VERSION.search(data))
    spectra = _RE_SPECTRA.search(data)
    n_spectra = int(spectra.group(1)) if spectra else None
    model = _instrument_model(data)
    sm = _RE_SOFTWARE.search(data)
    software = "%s %s" % (_text(sm, 1), _text(sm, 2)) if sm else None

    fields = {
        "version": version,
        "numSpectra": n_spectra,
        "instrumentModel": model,
        "software": software,
        "sourceFile": _text(_RE_SOURCE.search(data)),
        "startTimeStamp": _text(_RE_RUNTIME.search(data)),
    }
    return DecodedHeader(
        format="mzml",
        summary="mzML %s, %s spectra, %s"
        % (version or "?", n_spectra if n_spectra is not None else "?", model or "unknown instrument"),
        fields=fields,
    )


register_decoder("mzml", decode_mzml)
