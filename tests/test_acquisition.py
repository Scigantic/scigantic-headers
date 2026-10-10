"""EPU XML decoding and cross-source acquisition checks.

The EPU document below is synthetic: it follows the structure of a real EPU image
XML but every value is made up, so no collection data lives in this public repo.
The voltage case is modelled on a real incident with public numbers: a run that used
300 kV for EMPIAR-10581 / EMD-0731, which the deposit records at 200 kV.
"""

from scigantic_headers import (
    compare_acquisition,
    decode_epu_xml,
    epu_pixel_size_candidates,
)

EPU = b"""<MicroscopeImage xmlns="http://schemas.datacontract.org/2004/07/Fei.SharedObjects" xmlns:i="http://www.w3.org/2001/XMLSchema-instance">
<CustomData xmlns:a="http://schemas.microsoft.com/2003/10/Serialization/Arrays">
<a:KeyValueOfstringanyType><a:Key>Dose</a:Key><a:Value xmlns:b="http://www.w3.org/2001/XMLSchema" i:type="b:double">4.0E+21</a:Value></a:KeyValueOfstringanyType>
<a:KeyValueOfstringanyType><a:Key>DetectorCommercialName</a:Key><a:Value xmlns:b="http://www.w3.org/2001/XMLSchema" i:type="b:string">Example Detector</a:Value></a:KeyValueOfstringanyType>
<a:KeyValueOfstringanyType><a:Key>AppliedDefocus</a:Key><a:Value xmlns:b="http://www.w3.org/2001/XMLSchema" i:type="b:double">-2.2E-06</a:Value></a:KeyValueOfstringanyType>
<a:KeyValueOfstringanyType><a:Key>SuperResolutionFactor</a:Key><a:Value xmlns:b="http://www.w3.org/2001/XMLSchema" i:type="b:int">2</a:Value></a:KeyValueOfstringanyType>
<a:KeyValueOfstringanyType><a:Key>ElectronCountingEnabled</a:Key><a:Value xmlns:b="http://www.w3.org/2001/XMLSchema" i:type="b:boolean">true</a:Value></a:KeyValueOfstringanyType>
<a:KeyValueOfstringanyType><a:Key>PhasePlateUsed</a:Key><a:Value xmlns:b="http://www.w3.org/2001/XMLSchema" i:type="b:boolean">false</a:Value></a:KeyValueOfstringanyType>
</CustomData>
<SpatialScale><pixelSize><x><numericValue>1.92E-10</numericValue></x><y><numericValue>1.92E-10</numericValue></y></pixelSize></SpatialScale>
<microscopeData>
<acquisition><camera>
<Binning xmlns:a="http://schemas.datacontract.org/2004/07/System.Drawing"><a:x>2</a:x><a:y>2</a:y></Binning>
<ExposureTime>2.0</ExposureTime>
<ReadoutArea xmlns:a="http://schemas.datacontract.org/2004/07/System.Drawing"><a:height>2046</a:height><a:width>2880</a:width></ReadoutArea>
<FractionationSettings><b:NumberOffractions xmlns:b="x">30</b:NumberOffractions></FractionationSettings>
</camera><plateCamera><ExposureTime>0</ExposureTime></plateCamera></acquisition>
<core><ApplicationSoftware>EPU</ApplicationSoftware><ApplicationSoftwareVersion>3.0.0.1</ApplicationSoftwareVersion></core>
<gun><AccelerationVoltage>300000</AccelerationVoltage></gun>
<optics><EnergyFilter><EnergySelectionSlitWidth>20</EnergySelectionSlitWidth></EnergyFilter>
<TemMagnification><NominalMagnification>130000</NominalMagnification></TemMagnification></optics>
</microscopeData></MicroscopeImage>"""


def test_epu_fields():
    h = decode_epu_xml(EPU)
    assert h.format == "epu-xml"
    f = h.fields
    assert f["voltageKv"] == 300.0
    assert abs(f["pixelSizeA"] - 1.92) < 1e-9
    assert f["binning"] == [2, 2] and f["superResolutionFactor"] == 2
    assert f["fractions"] == 30 and f["exposureS"] == 2.0
    assert abs(f["doseEA2"] - 40.0) < 1e-9          # 4.0e21 e/m^2
    assert abs(f["appliedDefocusUm"] + 2.2) < 1e-9
    assert f["electronCounting"] is True and f["phasePlate"] is False
    assert f["detector"] == "Example Detector" and f["magnification"] == 130000
    assert f["energyFilterSlitEv"] == 20.0 and f["readoutAreaPx"] == {"height": 2046, "width": 2880}
    assert f["software"] == "EPU 3.0.0.1"


def test_epu_rejects_other_xml():
    assert decode_epu_xml(b"<RunInfo><Run/></RunInfo>") is None


def test_epu_pixel_candidates():
    f = decode_epu_xml(EPU).fields
    c = epu_pixel_size_candidates(f)
    assert abs(c["asRecorded"] - 1.92) < 1e-9
    assert abs(c["perBinning"] - 0.96) < 1e-9
    assert abs(c["perBinningAndSuperRes"] - 0.48) < 1e-9
    # unbinned, no super-resolution: only the recorded value
    assert list(epu_pixel_size_candidates({"pixelSizeA": 1.0})) == ["asRecorded"]


def test_voltage_mismatch_is_an_error():
    # The rehearsal incident: sheet/pipeline said 300 kV, the deposit says 200 kV.
    out = compare_acquisition({"pipeline": {"voltageKv": 300.0, "pixelSizeA": 0.8844},
                               "deposit": {"voltageKv": 200.0, "pixelSizeA": 0.88}})
    assert out[0].field == "voltageKv" and out[0].severity == "error"
    assert out[0].values == {"pipeline": 300.0, "deposit": 200.0}
    # 0.8844 vs 0.88 is 0.5%, under the 1% tolerance, so nothing is said about it
    assert [f.field for f in out] == ["voltageKv"]


def test_agreeing_sources_give_no_findings():
    assert compare_acquisition({
        "sheet": {"voltageKv": 300, "pixelSizeA": 0.84, "doseEA2": 42.4, "fractions": 40},
        "star": {"voltageKv": 300.0, "pixelSizeA": 0.8401},
        "header": {"pixelSizeA": 0.84},
    }) == []


def test_single_source_is_not_flagged():
    assert compare_acquisition({"only": {"voltageKv": 300, "pixelSizeA": 1.0}}) == []


def test_pixel_size_mismatch_points_at_the_binning_candidate():
    src = {"epu": {"pixelSizeA": 1.92, "binning": [2, 2], "superResolutionFactor": 2},
           "sheet": {"pixelSizeA": 0.97}}
    out = compare_acquisition(src)
    assert out[0].field == "pixelSizeA" and out[0].severity == "error"
    assert "perBinning" in out[0].message and "dimensions" in out[0].message


def test_dose_and_frames():
    # 42.4 vs 34.0 is a 25% gap: worth saying, but dose weighting is forgiving
    out = compare_acquisition({"a": {"doseEA2": 42.4, "fractions": 40},
                               "b": {"doseEA2": 34.0, "fractions": 40}})
    assert [(f.field, f.severity) for f in out] == [("doseEA2", "warn")]
    out = compare_acquisition({"a": {"doseEA2": 42.4}, "b": {"doseEA2": 30.0}})
    assert [(f.field, f.severity) for f in out] == [("doseEA2", "error")]
    out = compare_acquisition({"a": {"fractions": 40}, "b": {"fractions": 49}})
    assert out[0].field == "fractions" and out[0].severity == "error"


def test_nyquist_against_target():
    out = compare_acquisition({"x": {"pixelSizeA": 1.77}}, target_resolution_a=3.3)
    assert out[0].severity == "error" and "Nyquist" in out[0].message
    out = compare_acquisition({"x": {"pixelSizeA": 1.55}}, target_resolution_a=3.3)
    assert out[0].severity == "warn"            # 3.10 A Nyquist, within 10% of 3.3
    assert compare_acquisition({"x": {"pixelSizeA": 0.88}}, target_resolution_a=3.3) == []


def test_findings_are_ordered_worst_first_and_serialisable():
    out = compare_acquisition({"a": {"voltageKv": 300, "csMm": 2.7, "doseEA2": 40.0},
                               "b": {"voltageKv": 200, "csMm": 2.0, "doseEA2": 40.0}})
    assert [f.severity for f in out] == ["error", "warn"]
    assert out[0].to_dict()["field"] == "voltageKv"


def test_hint_tolerates_a_calibrated_pixel_size():
    # recorded 1.7675 at 2x2 binning -> 0.8837 per binning, against a stated 0.84 (5.2% apart)
    src = {"epu": {"pixelSizeA": 1.7675, "binning": [2, 2], "superResolutionFactor": 2},
           "sheet": {"pixelSizeA": 0.84}}
    out = compare_acquisition(src)
    assert out[0].field == "pixelSizeA" and "perBinning" in out[0].message
    assert "5.2%" in out[0].message
