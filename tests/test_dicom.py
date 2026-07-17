"""DICOM header decoding, including the PHI-avoidance guarantee."""

import json
import struct

from scigantic_headers import decode_dicom
from scigantic_headers.decoders import decode_bytes


def _pad(v: bytes) -> bytes:
    return v + (b" " if len(v) % 2 else b"")


def _explicit(group, elem, vr, value):
    out = struct.pack("<HH", group, elem) + vr
    if vr in (b"OB", b"OW", b"SQ", b"UN", b"UT"):
        out += b"\x00\x00" + struct.pack("<I", len(value))
    else:
        out += struct.pack("<H", len(value))
    return out + value


def _implicit(group, elem, value):
    return struct.pack("<HH", group, elem) + struct.pack("<I", len(value)) + value


def make_dicom(*, modality=b"MR", rows=256, cols=320, bits=16,
               transfer=b"1.2.840.10008.1.2.1", implicit=False, include_phi=False):
    def el(g, e, vr, v):
        return _implicit(g, e, v) if implicit else _explicit(g, e, vr, v)

    meta = _explicit(0x0002, 0x0010, b"UI", _pad(transfer))
    ds = b""
    ds += el(0x0008, 0x0016, b"UI", _pad(b"1.2.840.10008.5.1.4.1.1.4"))
    ds += el(0x0008, 0x0060, b"CS", _pad(modality))
    ds += el(0x0008, 0x0070, b"LO", _pad(b"ACME"))
    if include_phi:
        ds += el(0x0010, 0x0010, b"PN", _pad(b"DOE^JOHN"))     # PatientName
        ds += el(0x0010, 0x0020, b"LO", _pad(b"MRN0001"))      # PatientID
    ds += el(0x0028, 0x0002, b"US", struct.pack("<H", 1))
    ds += el(0x0028, 0x0010, b"US", struct.pack("<H", rows))
    ds += el(0x0028, 0x0011, b"US", struct.pack("<H", cols))
    ds += el(0x0028, 0x0100, b"US", struct.pack("<H", bits))
    return b"\x00" * 128 + b"DICM" + meta + ds


def test_explicit_vr_little_endian():
    h = decode_dicom(make_dicom())
    assert h.format == "dicom"
    assert h.fields["modality"] == "MR"
    assert h.fields["rows"] == 256
    assert h.fields["columns"] == 320
    assert h.fields["bitsAllocated"] == 16
    assert h.fields["manufacturer"] == "ACME"
    assert h.fields["transferSyntax"] == "1.2.840.10008.1.2.1"
    assert "MR" in h.summary and "320x256" in h.summary


def test_implicit_vr_little_endian():
    h = decode_dicom(make_dicom(transfer=b"1.2.840.10008.1.2", implicit=True))
    assert h.fields["modality"] == "MR"
    assert h.fields["rows"] == 256
    assert h.fields["columns"] == 320
    assert h.fields["bitsAllocated"] == 16


def test_never_surfaces_patient_identifiers():
    h = decode_dicom(make_dicom(include_phi=True))
    blob = json.dumps(h.to_dict())
    assert "patient" not in "".join(h.fields.keys()).lower()
    assert "DOE" not in blob and "MRN0001" not in blob
    # the technical fields are still read past the patient elements
    assert h.fields["rows"] == 256


def test_rejects_non_dicom():
    assert decode_dicom(b"\x00" * 200) is None     # no DICM magic at 128
    assert decode_dicom(b"short") is None


def test_dispatch():
    assert decode_bytes("img.dcm", make_dicom()).format == "dicom"


def test_json_safe():
    json.dumps(decode_dicom(make_dicom()).to_dict(), allow_nan=False)
