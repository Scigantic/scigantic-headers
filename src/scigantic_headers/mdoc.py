"""Decode a SerialEM .mdoc, the acquisition record written next to an image stack or movie.

SerialEM (and Thermo Fisher Tomography, which writes the same format) saves a
plain-text `.mdoc` beside each stack, tilt series, montage or saved movie. It holds
the settings the data was recorded with, image by image: tilt angle, dose, exposure
time, pixel spacing, defocus, stage position, magnification, number of frames and
the time of each exposure. The MRC or TIFF it describes usually carries none of
this, so the mdoc is the only place to read it from.

Layout:

    PixelSpacing = 1.35          global lines describe the whole stack
    Voltage = 300
    ImageFile = TS_01.mrc
    [T = SerialEM: Titan Krios ...  05-Jun-22  16:13:00]    free-text titles
    [T =     Tilt axis angle = 85.3, binning = 1  spot = 9  camera = 0]
    [ZValue = 0]                 one section per image in the stack
    TiltAngle = -0.01
    ExposureDose = 3.1
    ...

Three section kinds exist. `[ZValue = n]` is one image of a stack (one tilt, or one
piece of a montage). `[FrameSet = n]` is the whole of a single saved movie; the
small mdoc SerialEM writes beside each movie has one and no ZValue at all.
`[MontSection = n]` holds the settings of a montage as a whole; it is not an image,
so it is reported separately and left out of the per-image values. IMOD also reads
an image-series form (`ImageSeries = 1`, then one `[Image = name.tif]` section per
separate file); those sections count as images, and `imageFiles` lists the names.

Units are SerialEM's own and are kept: pixel spacing in angstrom, voltage in kV,
dose in electrons per square angstrom, defocus in micrometers (negative is
underfocus), exposure in seconds, stage position in micrometers. The fields that
`compare_acquisition` reads (`voltageKv`, `pixelSizeA`, `doseEA2`, `fractions`,
`magnification`) use the same names as the EPU decoder, so an mdoc can be checked
against a collection sheet, a RELION optics table or an EPU file directly.

`pixelSizeA` is the pixel spacing SerialEM recorded for the saved image. Movies saved
in super-resolution or with different binning from the image can sit on another
grid; check the movie's own dimensions before relying on it for the frames.
"""

from __future__ import annotations

import math
import re
from datetime import datetime
from typing import Dict, List, Optional, Tuple

from .decoders import DecodedHeader, Read, register_decoder

# A 60-tilt series is about 60 KiB; a large montage runs to a few hundred. Past this
# the per-image values are reported for the sections read, and `truncated` says so.
MDOC_READ_BYTES = 1024 * 1024

_SECTION_RE = re.compile(r"^\[\s*(ZValue|FrameSet|MontSection)\s*=\s*(-?\d+)\s*\]$")
# IMOD's image-series form: one section per separate image file, named in the header.
_IMAGE_SECTION_RE = re.compile(r"^\[\s*Image\s*=\s*(.+?)\s*\]$")
_TITLE_RE = re.compile(r"^\[\s*T\s*=(.*)\]$")
# SerialEM's keys are identifiers; a few camera-specific ones add '-', '(' and ')'.
_KV_RE = re.compile(r"^([A-Za-z][A-Za-z0-9_\-()]*)\s*=\s?(.*)$")
# Two title forms carry the tilt axis: "Tilt axis angle = 85.3, binning = 1 ..."
# (SerialEM 3 and later) and "TiltAxisAngle = 85.3  Binning = 1 ..." (older builds).
_TILT_AXIS_RE = re.compile(r"(?:Tilt axis angle|TiltAxisAngle)\s*=\s*(-?\d+(?:\.\d+)?)", re.I)
_TITLE_DATE_RE = re.compile(r"(\d{1,2}-[A-Za-z]{3}-\d{2,4}\s+\d{1,2}:\d{2}:\d{2})")
_VERSION_RE = re.compile(r"SerialEM Version\s+([0-9][0-9A-Za-z.\-]*)")

_MAX_TITLES = 10


def _num(v: Optional[str]) -> Optional[float]:
    """First number of a value, or None. SerialEM writes NaN for a reading it does
    not have (CameraLength in TEM mode); that is absent, not a number."""
    if v is None:
        return None
    parts = v.split()
    if not parts:
        return None
    try:
        x = float(parts[0])
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def _nums(v: Optional[str]) -> Optional[List[float]]:
    """Every number of a value (StagePosition is two, MinMaxMean three). Any token
    that is not a finite number makes the whole value absent, not a partial list."""
    if v is None:
        return None
    parts = v.split()
    if not parts:
        return None
    out = []
    for p in parts:
        try:
            x = float(p)
        except ValueError:
            return None
        if not math.isfinite(x):
            return None
        out.append(x)
    return out


def _int(v: Optional[str]) -> Optional[int]:
    x = _num(v)
    return int(x) if x is not None and x == int(x) and abs(x) < 2**53 else None


def _ints(v: Optional[str]) -> Optional[List[int]]:
    xs = _nums(v)
    if xs is None or any(x != int(x) or abs(x) >= 2**53 for x in xs):
        return None
    return [int(x) for x in xs]


def _text(v: Optional[str]) -> Optional[str]:
    v = (v or "").strip()
    return v or None


def _positive(x: Optional[float]) -> Optional[float]:
    return x if x is not None and x > 0 else None


def _iso(stamp: Optional[str]) -> Optional[str]:
    """'05-Jun-22  16:13:00' -> '2022-06-05T16:13:00'. SerialEM writes local time
    with no zone, so none is added. None if the stamp does not parse."""
    if not stamp:
        return None
    s = " ".join(stamp.split())
    for fmt in ("%d-%b-%y %H:%M:%S", "%d-%b-%Y %H:%M:%S"):
        try:
            return datetime.strptime(s, fmt).isoformat()
        except ValueError:
            continue
    return None


def _parse(text: str) -> Tuple[Dict[str, str], List[str], List[Tuple[str, int, Dict[str, str]]]]:
    glob: Dict[str, str] = {}
    titles: List[str] = []
    sections: List[Tuple[str, int, Dict[str, str]]] = []
    cur: Optional[Dict[str, str]] = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        m = _SECTION_RE.match(line)
        if m:
            cur = {}
            sections.append((m.group(1), int(m.group(2)), cur))
            continue
        m = _IMAGE_SECTION_RE.match(line)
        if m and not line.startswith("[T"):
            cur = {"__file": m.group(1)}
            sections.append(("Image", len(sections), cur))
            continue
        m = _TITLE_RE.match(line)
        if m:
            titles.append(m.group(1).strip())
            continue
        m = _KV_RE.match(line)
        if not m:
            continue  # a line cut off by the read, or one SerialEM's own readers would skip
        key, value = m.group(1), m.group(2).strip()
        # The mdoc SerialEM writes beside a saved movie has its title without
        # brackets ("T = SerialEM: ..."); it is still a title, not a setting.
        if key == "T" and cur is None:
            titles.append(value)
            continue
        (cur if cur is not None else glob)[key] = value
    return glob, titles, sections


def decode_mdoc(data: bytes) -> Optional[DecodedHeader]:
    if not data:
        return None
    if data.startswith(b"\xef\xbb\xbf"):
        data = data[3:]
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        # Windows builds can write a title or note in the system code page.
        text = data.decode("latin-1")
    if "\x00" in text:
        return None

    glob, titles, sections = _parse(text)

    images = [s for s in sections if s[0] in ("ZValue", "FrameSet", "Image")]
    montage_sections = [s for s in sections if s[0] == "MontSection"]
    # A real file has a section, or (written before the first image lands) a global
    # PixelSpacing or ImageFile line. Anything else is not an mdoc.
    if not sections and "PixelSpacing" not in glob and "ImageFile" not in glob:
        return None

    first = images[0][2] if images else {}

    def pick(key: str) -> Optional[str]:
        """Global value if present, else the first image's."""
        return glob[key] if key in glob else first.get(key)

    tilt_axis = None
    for t in titles:
        m = _TILT_AXIS_RE.search(t)
        if m:
            tilt_axis = float(m.group(1))
            break

    software = None
    m = _VERSION_RE.search(glob.get("Version", ""))
    if m:
        software = "SerialEM " + m.group(1)
    else:
        for t in titles:
            head = t.split(":", 1)[0].strip()
            if ":" in t and head in ("SerialEM", "Tomography"):
                software = head
                break

    title_date = None
    for t in titles:
        m = _TITLE_DATE_RE.search(t)
        if m:
            title_date = m.group(1)
            break

    tilts = [x for x in (_num(s[2].get("TiltAngle")) for s in images) if x is not None]
    doses = [x for x in (_num(s[2].get("ExposureDose")) for s in images) if x is not None]
    stamps = [_text(s[2].get("DateTime")) for s in images]
    stamps = [s for s in stamps if s]
    pieces = sum(1 for s in images if "PieceCoordinates" in s[2])

    total_dose = round(sum(doses), 4) if doses and any(d > 0 for d in doses) else None
    first_dose = _positive(_num(first.get("ExposureDose")))

    fields: Dict[str, object] = {
        # The values compare_acquisition reads, in the names the EPU decoder uses.
        "voltageKv": _positive(_num(pick("Voltage"))),
        "pixelSizeA": _positive(_num(pick("PixelSpacing"))),
        "magnification": _positive(_int(first.get("Magnification"))),
        "doseEA2": first_dose,  # one image (one tilt, or one movie); see totalDoseEA2
        "fractions": _positive(_int(first.get("NumSubFrames"))),
        "exposureS": _positive(_num(first.get("ExposureTime"))),
        # The stack as a whole.
        "software": software,
        "imageFile": _text(glob.get("ImageFile")),
        "imageSize": _ints(glob.get("ImageSize")),
        "dataMode": _int(glob.get("DataMode")),
        "montage": (_int(glob.get("Montage")) == 1) if "Montage" in glob else bool(montage_sections) or None,
        "tiltAxisAngle": tilt_axis,
        "titles": titles[:_MAX_TITLES],
        "sectionKind": images[0][0] if images else None,
        "imageCount": len(images),
        "montageSectionCount": len(montage_sections),
        "pieceCount": pieces or None,
        "imageFiles": [s[2]["__file"] for s in images if s[0] == "Image"][:_MAX_TITLES] or None,
        "truncated": len(data) >= MDOC_READ_BYTES,
        # Over every image read.
        "tiltAngleMin": min(tilts) if tilts else None,
        "tiltAngleMax": max(tilts) if tilts else None,
        "totalDoseEA2": total_dose,
        "acquiredAt": _iso(stamps[0] if stamps else title_date),
        "lastAcquiredAt": _iso(stamps[-1]) if len(stamps) > 1 else None,
        # The first image, as SerialEM wrote it (the settings are the same for every
        # tilt of a series; this is also the first file to land in a live session).
        "firstImage": {
            "tiltAngle": _num(first.get("TiltAngle")),
            "stagePosition": _nums(first.get("StagePosition")),
            "stageZ": _num(first.get("StageZ")),
            "defocusUm": _num(first.get("Defocus")),
            "targetDefocusUm": _num(first.get("TargetDefocus")),
            "imageShift": _nums(first.get("ImageShift")),
            "binning": _num(first.get("Binning")),
            "spotSize": _int(first.get("SpotSize")),
            "magIndex": _int(first.get("MagIndex")),
            "intensity": _num(first.get("Intensity")),
            "doseRate": _num(first.get("DoseRate")),
            "countsPerElectron": _num(first.get("CountsPerElectron")),
            "cameraIndex": _int(first.get("CameraIndex")),
            "rotationAngle": _num(first.get("RotationAngle")),
            "filterSlitAndLoss": _nums(first.get("FilterSlitAndLoss")),
            "subFramePath": _text(first.get("SubFramePath")),
            "frameDosesAndNumber": _nums(first.get("FrameDosesAndNumber")),
            "pieceCoordinates": _ints(first.get("PieceCoordinates")),
            "navigatorLabel": _text(first.get("NavigatorLabel")),
            "dateTime": _text(first.get("DateTime")),
        } if images else None,
    }
    if montage_sections:
        ms = montage_sections[0][2]
        fields["montageSection"] = {
            "fullMontageSize": _ints(ms.get("FullMontSize")),
            "pixelSizeA": _positive(_num(ms.get("PixelSpacing"))),
            "magnification": _positive(_int(ms.get("Magnification"))),
            "binning": _num(ms.get("Binning")),
        }

    summary = _summary(fields)
    return DecodedHeader(format="serialem-mdoc", summary=summary, fields=fields)


def _summary(f: Dict[str, object]) -> str:
    kind, n = f["sectionKind"], f["imageCount"]
    if kind == "FrameSet":
        head = "mdoc for one movie" if n == 1 else "mdoc, %d frame sets" % n
    elif f["montage"]:
        head = "montage mdoc, %d pieces" % (f["pieceCount"] or n)
    elif n:
        head = "mdoc, %d image%s" % (n, "" if n == 1 else "s")
    else:
        head = "mdoc, no images yet"
    parts = [head]
    if f["tiltAngleMin"] is not None and n > 1 and f["tiltAngleMax"] != f["tiltAngleMin"]:
        parts.append("tilt %g to %g deg" % (f["tiltAngleMin"], f["tiltAngleMax"]))
    if f["pixelSizeA"]:
        parts.append("%.4g A/px" % f["pixelSizeA"])
    if f["voltageKv"]:
        parts.append("%g kV" % f["voltageKv"])
    if f["totalDoseEA2"] and n > 1:
        parts.append("%.4g e/A^2 total" % f["totalDoseEA2"])
    elif f["doseEA2"]:
        parts.append("%.4g e/A^2" % f["doseEA2"])
    if f["fractions"]:
        parts.append("%d frame%s" % (f["fractions"], "" if f["fractions"] == 1 else "s"))
    if f["truncated"]:
        parts.append("first %d KiB read" % (MDOC_READ_BYTES // 1024))
    return ", ".join(parts)


register_decoder("mdoc", decode_mdoc, read=Read(leading=MDOC_READ_BYTES))
