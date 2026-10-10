"""Decode the acquisition metadata of a Thermo Fisher EPU image XML.

EPU writes one XML next to each movie (`FoilHole_..._Data_....xml`, root element
`MicroscopeImage`). It carries the settings the movie was recorded with: voltage,
pixel size, camera, binning, number of fractions, exposure, dose, applied defocus.
These are the values a processing run needs and the ones most often copied in by
hand from somewhere else, so they are worth reading from the source.

This is NOT registered for automatic dispatch. `.xml` is far too generic an
extension to claim (mzML and Illumina use it too), and the EPU file name carries no
fixed marker. Call `decode_epu_xml` on bytes you already know are an EPU file; it
returns None for anything without a `MicroscopeImage` root.

Units in the file, and what is returned here:
  AccelerationVoltage  volts          -> voltageKv (kV)
  pixelSize            metres         -> pixelSizeA (angstrom, AS RECORDED, see below)
  AppliedDefocus       metres         -> appliedDefocusUm (micrometres, negative = underfocus)
  Dose                 electrons/m^2  -> doseEA2 (electrons/angstrom^2)

`pixelSizeA` is the pixel size of the grid EPU recorded on. A camera read out with
2x2 binning reports a pixel twice the physical one, and a super-resolution factor
means the saved movie can sit on a finer grid still. Which grid a saved movie uses
depends on the export setting, so the file alone cannot say. `epu_pixel_size_candidates`
lists the plausible values; check the movie's own dimensions to choose.
"""

from __future__ import annotations

import re
from typing import Dict, Optional

from .decoders import DecodedHeader

_NS = rb"(?:\w+:)?"  # EPU XML prefixes its collection elements (a:, b:); match with or without


def _kv(data: bytes, key: str) -> Optional[str]:
    """Value of a `CustomData` key/value pair."""
    pat = (_NS + rb"Key>" + re.escape(key.encode()) + rb"</" + _NS + rb"Key>\s*<" + _NS +
           rb"Value[^>]*>([^<]*)</")
    m = re.search(b"<" + pat, data)
    return m.group(1).decode("utf-8", "replace").strip() if m else None


def _num(text: Optional[str]) -> Optional[float]:
    try:
        return float(text) if text is not None else None
    except ValueError:
        return None


def _first_float(data: bytes, pattern: bytes) -> Optional[float]:
    m = re.search(pattern, data, re.S)
    return _num(m.group(1).decode()) if m else None


def decode_epu_xml(data: bytes) -> Optional[DecodedHeader]:
    if b"<MicroscopeImage" not in data:
        return None

    volts = _first_float(data, rb"<AccelerationVoltage>([^<]*)</AccelerationVoltage>")
    px_m = _first_float(data, rb"<pixelSize>\s*<x>\s*<numericValue>([^<]*)</numericValue>")
    mag = _first_float(data, rb"<NominalMagnification>([^<]*)</NominalMagnification>")
    fractions = _first_float(data, rb"<" + _NS + rb"NumberOffractions\b[^>]*>(\d+)<")

    # The first positive ExposureTime is the camera's; a later plateCamera one is 0.
    exposure = next((v for v in (_num(x.decode()) for x in re.findall(rb"<ExposureTime>([^<]*)</ExposureTime>", data))
                     if v and v > 0), None)

    bin_m = re.search(rb"<Binning\b[^>]*>\s*<" + _NS + rb"x>(\d+)</" + _NS + rb"x>\s*<" + _NS + rb"y>(\d+)<", data)
    binning = (int(bin_m.group(1)), int(bin_m.group(2))) if bin_m else None
    ro = re.search(rb"<ReadoutArea\b[^>]*>\s*<" + _NS + rb"height>(\d+)</" + _NS + rb"height>\s*<" +
                   _NS + rb"width>(\d+)<", data)
    readout = {"height": int(ro.group(1)), "width": int(ro.group(2))} if ro else None

    superres = _num(_kv(data, "SuperResolutionFactor"))
    detector = _kv(data, "DetectorCommercialName") or _kv(data, "Detectors[EF-CCD].CommercialName")
    counting = _kv(data, "ElectronCountingEnabled") or _kv(data, "Detectors[EF-CCD].ElectronCounted")
    defocus_m = _num(_kv(data, "AppliedDefocus"))
    dose_m2 = _num(_kv(data, "Dose"))
    slit = _first_float(data, rb"<EnergySelectionSlitWidth>([^<]*)</EnergySelectionSlitWidth>")
    sw = re.search(rb"<ApplicationSoftware>([^<]*)</ApplicationSoftware>", data)
    sw_v = re.search(rb"<ApplicationSoftwareVersion>([^<]*)</ApplicationSoftwareVersion>", data)

    fields: Dict[str, object] = {
        "voltageKv": volts / 1000.0 if volts else None,
        "pixelSizeA": px_m * 1e10 if px_m else None,
        "magnification": int(mag) if mag else None,
        "detector": detector,
        "electronCounting": (counting or "").lower() == "true" if counting is not None else None,
        "binning": list(binning) if binning else None,
        "superResolutionFactor": int(superres) if superres else None,
        "readoutAreaPx": readout,
        "fractions": int(fractions) if fractions else None,
        "exposureS": exposure,
        "doseEA2": dose_m2 / 1e20 if dose_m2 else None,
        "appliedDefocusUm": defocus_m * 1e6 if defocus_m is not None else None,
        "energyFilterSlitEv": slit,
        "phasePlate": (_kv(data, "PhasePlateUsed") or "").lower() == "true" if _kv(data, "PhasePlateUsed") is not None else None,
        "software": ((sw.group(1).decode().strip() + " " + (sw_v.group(1).decode().strip() if sw_v else "")).strip()
                     if sw else None),
    }
    summary = "EPU image: %s kV, %s A/px as recorded, %s, %s fractions" % (
        ("%g" % fields["voltageKv"]) if fields["voltageKv"] else "?",
        ("%.4g" % fields["pixelSizeA"]) if fields["pixelSizeA"] else "?",
        detector or "unknown detector",
        fields["fractions"] if fields["fractions"] else "?",
    )
    return DecodedHeader(format="epu-xml", summary=summary, fields=fields)


def epu_pixel_size_candidates(fields: Dict[str, object]) -> Dict[str, float]:
    """Plausible pixel sizes (angstrom) for the saved movie, from decoded EPU fields.

    `asRecorded` is what EPU stored. `perBinning` divides by the binning, the
    physical pixel. `perBinningAndSuperRes` also divides by the super-resolution
    factor, the grid of a super-resolution movie. Only candidates that differ from
    `asRecorded` are added, so a plain unbinned camera returns one entry."""
    px = fields.get("pixelSizeA")
    if not px:
        return {}
    out = {"asRecorded": float(px)}
    b = fields.get("binning")
    s = fields.get("superResolutionFactor")
    per_bin = px / b[0] if b and b[0] and b[0] > 1 else None
    if per_bin:
        out["perBinning"] = per_bin
    if s and s > 1:
        out["perBinningAndSuperRes"] = (per_bin or px) / s
    return out
