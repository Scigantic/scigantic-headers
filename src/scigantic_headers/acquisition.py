"""Cross-check acquisition settings that come from different places.

A cryo-EM processing run needs voltage, pixel size, dose, frame count and Cs. They
usually arrive from several sources that were each typed or exported by someone
else: a collection sheet, the EPU or SerialEM file, the movie header, the RELION
or CryoSPARC optics table, the EMDB deposit. When two of them disagree, one is
wrong, and the wrong one silently caps the resolution rather than failing: a
processing run that used 300 kV for data recorded at 200 kV ran to the end and
simply stalled several angstrom short of the published map.

`compare_acquisition` takes the values from any number of sources, finds the fields
where they disagree beyond a tolerance, and says how serious it is. It does not
decide which source is right; it makes the disagreement visible before compute is
spent. Sources are plain dicts using the key names this library already returns
(`pixelSizeA`, `voltageKv` from `read_star_optics` and the MRC decoder; the rest from
`decode_epu_xml`), so they can be passed straight in.

    from scigantic_headers import compare_acquisition, read_star_optics
    findings = compare_acquisition({
        "sheet": {"voltageKv": 300, "pixelSizeA": 0.84, "doseEA2": 42.4, "fractions": 40},
        "relion": read_star_optics("Import/movies.star"),
        "deposit": {"voltageKv": 200, "pixelSizeA": 0.88},
    }, target_resolution_a=3.3)
    for f in findings:
        print(f.severity, f.field, f.message)
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Dict, List, Mapping, Optional

from .epu import epu_pixel_size_candidates

_SEVERITY_ORDER = {"error": 0, "warn": 1, "info": 2}


@dataclass(frozen=True)
class Finding:
    """One disagreement (or risk) in the acquisition settings."""

    field: str
    severity: str  # "error" | "warn" | "info"
    message: str
    values: Dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


def _present(sources: Mapping[str, Mapping[str, object]], key: str) -> Dict[str, float]:
    out = {}
    for name, src in sources.items():
        v = src.get(key) if src else None
        if isinstance(v, (int, float)) and not isinstance(v, bool) and v == v and v > 0:
            out[name] = float(v)
    return out


def _fmt(values: Mapping[str, float]) -> str:
    return ", ".join("%s=%g" % (k, v) for k, v in values.items())


def _spread(values: Mapping[str, float]) -> float:
    lo, hi = min(values.values()), max(values.values())
    return hi / lo - 1.0


def compare_acquisition(
    sources: Mapping[str, Mapping[str, object]],
    *,
    target_resolution_a: Optional[float] = None,
) -> List[Finding]:
    """Compare acquisition values across sources. Returns findings, worst first.

    An empty list means every field present in two or more sources agrees within
    tolerance (and, if `target_resolution_a` is given, the pixel size can reach it).
    A field only one source reports cannot be checked and is not flagged."""
    findings: List[Finding] = []

    volts = _present(sources, "voltageKv")
    if len(volts) >= 2 and max(volts.values()) - min(volts.values()) > 0.5:
        findings.append(Finding(
            "voltageKv", "error",
            "voltage differs between sources (%s). A wrong voltage corrupts the CTF and "
            "the dose weighting at high resolution without any error being raised." % _fmt(volts),
            dict(volts)))

    pix = _present(sources, "pixelSizeA")
    if len(pix) >= 2:
        s = _spread(pix)
        if s > 0.01:
            hint = ""
            for name, src in sources.items():
                for label, cand in epu_pixel_size_candidates(dict(src or {})).items():
                    if label == "asRecorded":
                        continue
                    for other_name, other in pix.items():
                        # Calibrated and nominal pixel sizes routinely differ by a few
                        # percent, so match loosely; this is a pointer, not a verdict.
                        if other_name != name and abs(cand / other - 1) < 0.06:
                            hint = (" %s=%g is within %.1f%% of %g, the %s pixel size derived from '%s' "
                                    "(recorded %g). Check which grid the movies are saved on from "
                                    "their own dimensions."
                                    % (other_name, other, 100 * abs(cand / other - 1), cand, label, name,
                                       src.get("pixelSizeA")))
                            break
                    if hint:
                        break
                if hint:
                    break
            findings.append(Finding(
                "pixelSizeA", "error" if s > 0.05 else "warn",
                "pixel size differs by %.1f%% (%s).%s" % (100 * s, _fmt(pix), hint), dict(pix)))

    dose = _present(sources, "doseEA2")
    if len(dose) >= 2:
        s = _spread(dose)
        if s > 0.10:
            findings.append(Finding(
                "doseEA2", "error" if s > 0.25 else "warn",
                "total dose differs by %.0f%% (%s). Dose weighting uses the per-frame dose." % (100 * s, _fmt(dose)),
                dict(dose)))

    frames = _present(sources, "fractions")
    if len(frames) >= 2 and len({int(v) for v in frames.values()}) > 1:
        findings.append(Finding(
            "fractions", "error",
            "frame count differs between sources (%s)." % _fmt(frames), dict(frames)))

    cs = _present(sources, "csMm")
    if len(cs) >= 2 and max(cs.values()) - min(cs.values()) > 0.1:
        findings.append(Finding("csMm", "warn", "Cs differs between sources (%s)." % _fmt(cs), dict(cs)))

    mag = _present(sources, "magnification")
    if len(mag) >= 2 and _spread(mag) > 0.01:
        findings.append(Finding("magnification", "warn",
                                "nominal magnification differs (%s)." % _fmt(mag), dict(mag)))

    if target_resolution_a:
        for px in sorted({round(v, 4) for v in pix.values()}):
            nyq = 2 * px
            if nyq > target_resolution_a:
                findings.append(Finding(
                    "pixelSizeA", "error",
                    "at %g A/px the Nyquist limit is %.2f A, coarser than the %g A target; the "
                    "target cannot be reached on this grid." % (px, nyq, target_resolution_a),
                    {"pixelSizeA": px, "nyquistA": nyq}))
            elif nyq * 1.1 > target_resolution_a:
                findings.append(Finding(
                    "pixelSizeA", "warn",
                    "at %g A/px the Nyquist limit is %.2f A, within 10%% of the %g A target; "
                    "there is almost no room above it." % (px, nyq, target_resolution_a),
                    {"pixelSizeA": px, "nyquistA": nyq}))

    findings.sort(key=lambda f: _SEVERITY_ORDER[f.severity])
    return findings
