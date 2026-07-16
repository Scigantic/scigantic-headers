"""Read acquisition optics from a RELION STAR file.

A raw MRC movie header carries dimensions, frame count, and dtype, but usually
not the pixel size (CELLA is zero). That number lives in the acquisition
metadata, which for RELION/CryoSPARC workflows is a STAR file. This extracts it
(plus voltage and Cs) so a context record can carry a real pixel size instead of
null.

STAR is a text format; this is a small, focused parser for the optics fields,
not a general STAR library. It handles the RELION 3.1+ `data_optics` loop, the
older key-value form, and the legacy detector-pixel-size / magnification pair.
Pure stdlib, no dependencies.
"""

from __future__ import annotations

from typing import Dict, Optional

# The optics block sits at the top of a RELION star file (before the large
# data_particles loop), so a bounded read captures it without loading gigabytes.
STAR_SCAN_BYTES = 256 * 1024


def parse_relion_optics(text: str) -> Dict[str, float]:
    """Extract {pixelSizeA, voltageKv, sphericalAberrationMm} (those present)
    from RELION STAR text. Returns an empty dict if none are found."""
    kv: Dict[str, str] = {}
    loops = []  # (columns, first_data_row)
    lines = text.splitlines()
    i, n = 0, len(lines)

    while i < n:
        s = lines[i].strip()
        i += 1
        if not s or s.startswith("#"):
            continue
        if s == "loop_":
            cols = []
            while i < n:
                t = lines[i].strip()
                if t.startswith("_"):
                    cols.append(t.split()[0])  # tag, dropping any "#n"
                    i += 1
                elif not t or t.startswith("#"):
                    i += 1
                else:
                    break
            first = None
            while i < n:
                t = lines[i].strip()
                if not t or t.startswith("#"):
                    i += 1
                    continue
                if t.startswith("_") or t.startswith("data_") or t == "loop_":
                    break
                first = t.split()
                i += 1
                break
            if first:
                loops.append((cols, first))
        elif s.startswith("_"):
            parts = s.split()
            if len(parts) >= 2:
                kv[parts[0]] = parts[1]

    def get(tag: str) -> Optional[str]:
        if tag in kv:
            return kv[tag]
        for cols, row in loops:
            if tag in cols:
                idx = cols.index(tag)
                if idx < len(row):
                    return row[idx]
        return None

    def as_float(v: Optional[str]) -> Optional[float]:
        try:
            return float(v) if v is not None else None
        except ValueError:
            return None

    out: Dict[str, float] = {}

    px = as_float(get("_rlnImagePixelSize")) or as_float(get("_rlnMicrographPixelSize"))
    if px is None:
        # Legacy form: detector pixel size (micron) / magnification, in angstrom.
        dps = as_float(get("_rlnDetectorPixelSize"))
        mag = as_float(get("_rlnMagnification"))
        if dps is not None and mag:
            px = dps / mag * 1e4
    if px is not None:
        out["pixelSizeA"] = round(px, 4)

    volt = as_float(get("_rlnVoltage"))
    if volt is not None:
        out["voltageKv"] = round(volt, 1)

    cs = as_float(get("_rlnSphericalAberration"))
    if cs is not None:
        out["sphericalAberrationMm"] = round(cs, 3)

    return out


def read_star_optics(path: str, max_bytes: int = STAR_SCAN_BYTES) -> Dict[str, float]:
    """Read a STAR file's leading `max_bytes` and parse its optics. Bounded so a
    huge data_particles table is never fully loaded, the optics block is first.
    Returns {} if the file can't be read or has no optics."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read(max_bytes)
    except OSError:
        return {}
    return parse_relion_optics(text)
