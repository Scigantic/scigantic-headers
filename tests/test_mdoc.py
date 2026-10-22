"""SerialEM .mdoc decoding.

The mdoc documents below are synthetic: they follow the layout and key names of
real SerialEM and Thermo Fisher Tomography output, but every value is made up, so
no collection data lives in this public repo. During development the decoder was
also run over 116 public mdoc files from GitHub and agreed with mdocfile
(github.com/teamtomo/mdocfile) on every value compared, for all 111 that mdocfile
can read.
"""

import json

from scigantic_headers import compare_acquisition, decode_bytes, decode_mdoc, has_decoder_for

TILT_SERIES = b"""PixelSpacing = 1.35
Voltage = 300
ImageFile = TS_01.mrc
ImageSize = 4096 4096
DataMode = 1
Version = SerialEM Version 4.1.10 64-bit,  built Jan  1 2024  10:00:00

[T = SerialEM: Example Krios   01-Mar-24  09:15:00    ]

[T =     Tilt axis angle = 85.3, binning = 1  spot = 9  camera = 0]

[ZValue = 0]
TiltAngle = 0.01
StagePosition = 10.5 -20.25
StageZ = 3.5
Magnification = 81000
Intensity = 0.12
ExposureDose = 3.1
PixelSpacing = 1.35
SpotSize = 9
Defocus = -3.02
ImageShift = 0.1 -0.2
RotationAngle = 175.3
ExposureTime = 0.8
Binning = 1
MagIndex = 31
CountsPerElectron = 32
MinMaxMean = -12 400 30.5
TargetDefocus = -3
SubFramePath = X:\\frames\\TS_01_000_0.0.tif
NumSubFrames = 8
FrameDosesAndNumber = 0.3875 8
DateTime = 01-Mar-24  09:16:02
FilterSlitAndLoss = 20 0
CameraLength = NaN

[ZValue = 1]
TiltAngle = 3.0
ExposureDose = 3.1
ExposureTime = 0.8
NumSubFrames = 8
DateTime = 01-Mar-24  09:16:40

[ZValue = 2]
TiltAngle = -3.0
ExposureDose = 3.2
ExposureTime = 0.8
NumSubFrames = 8
DateTime = 01-Mar-24  09:17:15
"""

# The file SerialEM writes beside each saved movie: title without brackets, one FrameSet.
FRAME_MDOC = b"""T = SerialEM: Example Krios                                  02-Mar-24  11:00:00

Voltage = 300
[FrameSet = 0]
TiltAngle = 12.0
StagePosition = 1 2
Magnification = 105000
ExposureDose = 2.5
PixelSpacing = 0.83
ExposureTime = 2
Binning = 1
NumSubFrames = 40
DateTime = 02-Mar-24  11:00:05
"""

# A montage: the [MontSection] describes the montage as a whole and is not an image.
MONTAGE = b"""PixelSpacing = 120.5
ImageFile = grid.mrc
ImageSize = 2880 2046
Montage = 1
DataMode = 1

[T = SerialEM: Example Krios   03-Mar-24  08:00:00    ]

[ZValue = 0]
PieceCoordinates = 0 0 0
TiltAngle = 0.0
Magnification = 135
ExposureDose = 0
PixelSpacing = 120.5
DateTime = 03-Mar-24  08:01:00

[ZValue = 1]
PieceCoordinates = 2600 0 0
TiltAngle = 0.0
Magnification = 135
ExposureDose = 0
DateTime = 03-Mar-24  08:01:30

[MontSection = 0]
TiltAngle = 0.0
Magnification = 135
PixelSpacing = 120.5
Binning = 2
FullMontSize = 5480 2046
"""


def test_tilt_series_fields():
    h = decode_mdoc(TILT_SERIES)
    f = h.fields
    assert h.format == "serialem-mdoc"
    assert f["voltageKv"] == 300 and f["pixelSizeA"] == 1.35
    assert f["magnification"] == 81000 and f["doseEA2"] == 3.1
    assert f["fractions"] == 8 and f["exposureS"] == 0.8
    assert f["software"] == "SerialEM 4.1.10"
    assert f["imageFile"] == "TS_01.mrc" and f["imageSize"] == [4096, 4096] and f["dataMode"] == 1
    assert f["tiltAxisAngle"] == 85.3
    assert f["sectionKind"] == "ZValue" and f["imageCount"] == 3 and f["montageSectionCount"] == 0
    assert (f["tiltAngleMin"], f["tiltAngleMax"]) == (-3.0, 3.0)
    assert f["totalDoseEA2"] == 9.4
    assert f["acquiredAt"] == "2024-03-01T09:16:02" and f["lastAcquiredAt"] == "2024-03-01T09:17:15"
    first = f["firstImage"]
    assert first["stagePosition"] == [10.5, -20.25] and first["imageShift"] == [0.1, -0.2]
    assert first["defocusUm"] == -3.02 and first["targetDefocusUm"] == -3
    assert first["subFramePath"] == "X:\\frames\\TS_01_000_0.0.tif"
    assert first["frameDosesAndNumber"] == [0.3875, 8]
    assert h.summary == "mdoc, 3 images, tilt -3 to 3 deg, 1.35 A/px, 300 kV, 9.4 e/A^2 total, 8 frames"


def test_frame_mdoc_with_unbracketed_title():
    f = decode_mdoc(FRAME_MDOC).fields
    assert f["sectionKind"] == "FrameSet" and f["imageCount"] == 1
    assert f["software"] == "SerialEM"
    assert f["titles"][0].startswith("SerialEM: Example Krios")
    # global PixelSpacing is absent here, so the frame set's own is used
    assert f["pixelSizeA"] == 0.83 and f["voltageKv"] == 300
    assert f["fractions"] == 40 and f["doseEA2"] == 2.5
    assert decode_mdoc(FRAME_MDOC).summary.startswith("mdoc for one movie")


def test_montage_section_is_not_an_image():
    f = decode_mdoc(MONTAGE).fields
    assert f["montage"] is True
    assert f["imageCount"] == 2 and f["pieceCount"] == 2 and f["montageSectionCount"] == 1
    assert f["firstImage"]["pieceCoordinates"] == [0, 0, 0]
    assert f["montageSection"] == {"fullMontageSize": [5480, 2046], "pixelSizeA": 120.5,
                                   "magnification": 135, "binning": 2.0}
    # a zero dose is "not calibrated", not a measurement
    assert f["doseEA2"] is None and f["totalDoseEA2"] is None
    assert decode_mdoc(MONTAGE).summary.startswith("montage mdoc, 2 pieces")


def test_tomography_title_and_older_tilt_axis_form():
    data = (b"PixelSpacing = 2.0\nImageFile = p.mrc\n"
            b"[T = Tomography: EXAMPLE    04-Mar-24  10:00:00]\n"
            b"[T = TiltAxisAngle = -174.2  Binning = 1  SpotSize = 7]\n"
            b"[ZValue = 0]\nTiltAngle = 1\n")
    f = decode_mdoc(data).fields
    assert f["software"] == "Tomography" and f["tiltAxisAngle"] == -174.2
    assert f["acquiredAt"] == "2024-03-04T10:00:00"


def test_imod_image_series():
    data = b"ImageSeries = 1\n\n[Image = a.tif]\nTiltAngle = 3\n\n[Image = b.tif]\nTiltAngle = 4\n"
    f = decode_mdoc(data).fields
    assert f["imageCount"] == 2 and f["imageFiles"] == ["a.tif", "b.tif"]
    assert (f["tiltAngleMin"], f["tiltAngleMax"]) == (3, 4)


def test_global_only_file_written_before_the_first_image():
    f = decode_mdoc(b"PixelSpacing = 1.1\nImageFile = new.mrc\nImageSize = 10 10\n").fields
    assert f["imageCount"] == 0 and f["firstImage"] is None and f["pixelSizeA"] == 1.1


def test_crlf_bom_and_latin1():
    crlf = b"\xef\xbb\xbf" + TILT_SERIES.replace(b"\n", b"\r\n")
    assert decode_mdoc(crlf).fields == decode_mdoc(TILT_SERIES).fields
    latin = TILT_SERIES.replace(b"Example Krios", b"Example Krios \xb5")
    assert decode_mdoc(latin).fields["titles"][0].startswith("SerialEM: Example Krios \u00b5")


def test_bad_values_are_absent_not_errors():
    data = (b"PixelSpacing = abc\nVoltage = NaN\n[ZValue = 0]\nTiltAngle = inf\n"
            b"StagePosition = 1 x\nMagnification = 1.5\nExposureDose = -2\n")
    f = decode_mdoc(data).fields
    assert f["pixelSizeA"] is None and f["voltageKv"] is None and f["tiltAngleMin"] is None
    assert f["firstImage"]["stagePosition"] is None and f["magnification"] is None and f["doseEA2"] is None
    json.dumps(decode_mdoc(data).to_dict(), allow_nan=False)


def test_truncated_read_is_reported():
    from scigantic_headers.mdoc import MDOC_READ_BYTES
    big = TILT_SERIES + b"".join(b"[ZValue = %d]\nTiltAngle = 1\nExposureDose = 1\n" % i
                                 for i in range(3, 40000))
    h = decode_mdoc(big[:MDOC_READ_BYTES])
    assert h.fields["truncated"] is True and "first 1024 KiB read" in h.summary
    assert decode_mdoc(TILT_SERIES).fields["truncated"] is False


def test_rejects_what_is_not_an_mdoc():
    for data in (b"", b"\x00\x01\x02" * 50, b"hello world\n", b"# a comment\nkey: value\n",
                 b"<MicroscopeImage/>", b"[section]\nname = x\n"):
        assert decode_mdoc(data) is None


def test_dispatched_by_extension():
    assert has_decoder_for("TS_01.mrc.mdoc") and has_decoder_for("frames/x.tif.MDOC")
    assert decode_bytes("TS_01.mrc.mdoc", TILT_SERIES).format == "serialem-mdoc"


def test_feeds_compare_acquisition():
    mdoc = decode_mdoc(TILT_SERIES).fields
    out = compare_acquisition({"mdoc": mdoc, "sheet": {"voltageKv": 200, "pixelSizeA": 1.35, "fractions": 8}})
    assert [(x.field, x.severity) for x in out] == [("voltageKv", "error")]
