"""RELION STAR optics parsing — the fields that complete a card when the raw
MRC header omits pixel size. Fixtures match the exact text RELION writes."""

import os
import tempfile

from scigantic_headers import parse_relion_optics, read_star_optics

# RELION 3.1+ data_optics loop, as written by relion_refine / motioncorr.
RELION_31 = """
# version 30001

data_optics

loop_
_rlnOpticsGroupName #1
_rlnOpticsGroup #2
_rlnMicrographOriginalPixelSize #3
_rlnVoltage #4
_rlnSphericalAberration #5
_rlnImagePixelSize #6
opticsGroup1            1     0.885000   300.000000     2.700000     1.060000

data_particles

loop_
_rlnCoordinateX #1
_rlnCoordinateY #2
1234.5  6789.0
"""

# Older key-value style (no loop for optics).
RELION_KV = """
data_

_rlnImagePixelSize   1.35
_rlnVoltage          200.0
_rlnSphericalAberration  2.7
"""

# Legacy detector-pixel-size + magnification (pre-ImagePixelSize).
RELION_LEGACY = """
data_

loop_
_rlnDetectorPixelSize #1
_rlnMagnification #2
_rlnVoltage #3
14.0   130435.0   300.0
"""


def test_relion_31_loop():
    o = parse_relion_optics(RELION_31)
    assert o["pixelSizeA"] == 1.06         # _rlnImagePixelSize wins over Original
    assert o["voltageKv"] == 300.0
    assert o["sphericalAberrationMm"] == 2.7


def test_relion_key_value():
    o = parse_relion_optics(RELION_KV)
    assert o["pixelSizeA"] == 1.35
    assert o["voltageKv"] == 200.0


def test_relion_legacy_detector_magnification():
    o = parse_relion_optics(RELION_LEGACY)
    # 14.0 micron / 130435 * 1e4 = 1.0733... A/px
    assert round(o["pixelSizeA"], 3) == 1.073
    assert o["voltageKv"] == 300.0


# The exact column layout of a real deposition (EMPIAR-10084/micrographs_ctf.star):
# optics live in the MAIN data loop, not a data_optics block, and pixel size comes
# via DetectorPixelSize / Magnification. Verified against the live file.
REAL_EMPIAR_10084 = """
data_

loop_
_rlnMicrographName #1
_rlnCtfImage #2
_rlnDefocusU #3
_rlnDefocusV #4
_rlnDefocusAngle #5
_rlnVoltage #6
_rlnSphericalAberration #7
_rlnAmplitudeContrast #8
_rlnMagnification #9
_rlnDetectorPixelSize #10
_rlnCtfFigureOfMerit #11
micrographs/Hb_087.mrc CtfFind/job004/Hb_087.ctf:mrc 4999.59 4890.79 -34.41 300.000000 2.600000 0.070000 10000.000000 1.050000 0.024571
"""


def test_real_empiar_10084_main_loop_ctf():
    o = parse_relion_optics(REAL_EMPIAR_10084)
    assert o["pixelSizeA"] == 1.05          # 1.05 micron / 10000 * 1e4
    assert o["voltageKv"] == 300.0
    assert o["sphericalAberrationMm"] == 2.6


def test_micrograph_pixel_size_fallback():
    text = "data_optics\nloop_\n_rlnMicrographPixelSize #1\n0.82\n"
    assert parse_relion_optics(text)["pixelSizeA"] == 0.82


def test_empty_and_irrelevant():
    assert parse_relion_optics("") == {}
    assert parse_relion_optics("data_junk\nloop_\n_rlnCoordinateX #1\n1.0\n") == {}


def test_read_star_optics_from_file():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "particles.star")
        with open(p, "w") as fh:
            fh.write(RELION_31)
        assert read_star_optics(p)["pixelSizeA"] == 1.06


def test_read_star_optics_missing_file():
    assert read_star_optics("/no/such/file.star") == {}
