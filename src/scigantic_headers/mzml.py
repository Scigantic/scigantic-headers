"""Decode an mzML mass-spectrometry file's header metadata.

mzML is XML. The metadata worth surfacing sits in the document preamble: the
format version, the source file, the software, the instrument model (in
<instrumentConfiguration> or a param group it references), and the spectrum count and run start time
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

from .decoders import DecodedHeader, Read, register_decoder

_NS = b"http://psi.hupo.org/ms/mzml"

_RE_VERSION = re.compile(rb'<mzML\b[^>]*\bversion="([^"]*)"')
_RE_SPECTRA = re.compile(rb'<spectrumList\b[^>]*\bcount="(\d+)"')
_RE_RUNTIME = re.compile(rb'<run\b[^>]*\bstartTimeStamp="([^"]*)"')
_RE_SOURCE = re.compile(rb'<sourceFile\b[^>]*\bname="([^"]*)"')
_RE_SOFTWARE = re.compile(rb'<software\b[^>]*\bid="([^"]*)"[^>]*\bversion="([^"]*)"')
_RE_CVNAME = re.compile(rb'<cvParam\b[^>]*\bname="([^"]*)"')


def _text(match, group=1):
    return match.group(group).decode("utf-8", "replace") if match else None


# cvParams that can sit beside the model in an instrument configuration or param
# group but are not the model: the serial number and the customization note.
_NOT_MODEL = {b"MS:1000529", b"MS:1000032"}
_RE_CVTAG = re.compile(rb"<cvParam\b[^>]*>")
_RE_ACCESSION = re.compile(rb'\baccession="([^"]*)"')
_RE_GROUPREF = re.compile(rb'<referenceableParamGroupRef\b[^>]*\bref="([^"]*)"')


def _first_model_param(scope: bytes) -> Optional[str]:
    """The name of the first cvParam in `scope` that is not in _NOT_MODEL."""
    for tag in _RE_CVTAG.findall(scope):
        acc = _RE_ACCESSION.search(tag)
        if acc and acc.group(1) in _NOT_MODEL:
            continue
        name = _RE_CVNAME.search(tag)
        if name:
            return _text(name)
    return None


def _param_group(data: bytes, group_id: bytes) -> bytes:
    """The body of the <referenceableParamGroup> with this id, or b"" if it is
    not in `data`."""
    m = re.search(rb'<referenceableParamGroup\b[^>]*\bid="' + re.escape(group_id) + rb'"[^>]*>', data)
    if not m:
        return b""
    end = data.find(b"</referenceableParamGroup>", m.end())
    return data[m.end():end] if end != -1 else b""


def _instrument_model(data: bytes) -> Optional[str]:
    """The instrument model of the first <instrumentConfiguration>. The schema
    puts it ahead of the <componentList>, either as a cvParam or inside a
    referenceableParamGroup that the configuration points to (msconvert writes
    it that way). The components (source, analyzer, detector) are never the
    model, so they are not searched."""
    start = data.find(b"<instrumentConfiguration")
    if start == -1:
        return None
    end = data.find(b"</instrumentConfiguration>", start)
    scope = data[start:end] if end != -1 else data[start:start + 8192]
    comp = scope.find(b"<componentList")
    head = scope[:comp] if comp != -1 else scope
    model = _first_model_param(head)
    if model:
        return model
    for ref in _RE_GROUPREF.findall(head):
        model = _first_model_param(_param_group(data, ref))
        if model:
            return model
    return None


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


register_decoder("mzml", decode_mzml, read=Read(leading=1024 * 1024))  # XML preamble runs to <spectrumList>
