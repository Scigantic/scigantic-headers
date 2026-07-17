"""Decode a DICOM medical-image file header.

A DICOM file is a 128-byte preamble, the magic 'DICM', a File Meta group (0002,
always explicit VR little-endian) that names the transfer syntax, then the
dataset in that syntax. This reads the technical description of the image and
deliberately does NOT read patient identifiers (name, id, birth date): it answers
"what is this file", not "whose data is it". PixelData is never read.
"""

from __future__ import annotations

from typing import Optional

from .decoders import DecodedHeader, register_decoder

# Explicit-VR value representations that use a 2-byte reserved field + 4-byte
# length instead of a 2-byte length.
_LONG_VR = {b"OB", b"OW", b"OF", b"OD", b"OL", b"SQ", b"UT", b"UN", b"UC", b"UR"}
_META = 0x0002
_PIXEL_DATA = (0x7FE0, 0x0010)

_IMPLICIT_LE = "1.2.840.10008.1.2"
_EXPLICIT_BE = "1.2.840.10008.1.2.2"

# Technical tags only, with how to interpret the value. No patient tags.
_WANTED = {
    (0x0008, 0x0060): ("modality", "str"),
    (0x0008, 0x0070): ("manufacturer", "str"),
    (0x0008, 0x1030): ("studyDescription", "str"),
    (0x0008, 0x103E): ("seriesDescription", "str"),
    (0x0008, 0x0016): ("sopClassUid", "str"),
    (0x0028, 0x0002): ("samplesPerPixel", "us"),
    (0x0028, 0x0010): ("rows", "us"),
    (0x0028, 0x0011): ("columns", "us"),
    (0x0028, 0x0100): ("bitsAllocated", "us"),
}


def _read_element(data: bytes, i: int, explicit: bool, little: bool):
    """(group, element, value, next_i) or None if truncated / undefined-length."""
    bo = "little" if little else "big"
    if i + 8 > len(data):
        return None
    group = int.from_bytes(data[i:i + 2], bo)
    elem = int.from_bytes(data[i + 2:i + 4], bo)
    i += 4
    if explicit:
        vr = data[i:i + 2]
        i += 2
        if vr in _LONG_VR:
            i += 2  # reserved
            if i + 4 > len(data):
                return None
            length = int.from_bytes(data[i:i + 4], bo)
            i += 4
        else:
            if i + 2 > len(data):
                return None
            length = int.from_bytes(data[i:i + 2], bo)
            i += 2
    else:
        if i + 4 > len(data):
            return None
        length = int.from_bytes(data[i:i + 4], bo)
        i += 4
    if length == 0xFFFFFFFF:  # undefined-length sequence: stop, out of scope
        return None
    if i + length > len(data):
        return None
    return group, elem, data[i:i + length], i + length


def decode_dicom(data: bytes) -> Optional[DecodedHeader]:
    if len(data) < 132 or data[128:132] != b"DICM":
        return None

    # File Meta group 0002 is always explicit VR little-endian.
    transfer = None
    i = 132
    while True:
        r = _read_element(data, i, explicit=True, little=True)
        if r is None:
            break
        group, elem, value, nxt = r
        if group != _META:
            break  # dataset begins at i (this element is not meta); do not advance
        if (group, elem) == (_META, 0x0010):
            transfer = value.rstrip(b"\x00 ").decode("ascii", "replace")
        i = nxt

    explicit, little = True, True
    if transfer == _IMPLICIT_LE:
        explicit = False
    elif transfer == _EXPLICIT_BE:
        little = False

    found = {}
    wanted = dict(_WANTED)
    guard = 0
    while wanted and guard < 200000:
        guard += 1
        r = _read_element(data, i, explicit, little)
        if r is None:
            break
        group, elem, value, nxt = r
        if (group, elem) == _PIXEL_DATA:
            break
        spec = wanted.pop((group, elem), None)
        if spec:
            name, kind = spec
            if kind == "us" and len(value) >= 2:
                found[name] = int.from_bytes(value[:2], "little" if little else "big")
            elif kind == "str":
                found[name] = value.rstrip(b"\x00 ").decode("utf-8", "replace").strip() or None
        i = nxt

    dims = ""
    if found.get("rows") and found.get("columns"):
        dims = ", %dx%d" % (found["columns"], found["rows"])
    fields = {
        "modality": found.get("modality"),
        "rows": found.get("rows"),
        "columns": found.get("columns"),
        "bitsAllocated": found.get("bitsAllocated"),
        "samplesPerPixel": found.get("samplesPerPixel"),
        "manufacturer": found.get("manufacturer"),
        "studyDescription": found.get("studyDescription"),
        "seriesDescription": found.get("seriesDescription"),
        "sopClassUid": found.get("sopClassUid"),
        "transferSyntax": transfer,
    }
    return DecodedHeader(
        format="dicom",
        summary="DICOM %s%s" % (found.get("modality") or "image", dims),
        fields=fields,
    )


register_decoder("dcm", decode_dicom)
