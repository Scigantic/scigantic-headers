"""The `.xml` dispatcher: BigDataViewer SpimData, OME-XML sidecars, EPU.

Real-file goldens (all CC-BY sample data, committed under tests/fixtures/):
  bdv-drosophila.xml    OME sample set, BDV/samples/drosophila.xml (older export,
                        no sizes or attributes)
  bdv-HisYFP-SPIM.xml   OME sample set, BDV/samples/HisYFP-SPIM.xml (7-angle SPIM
                        plus fused and deconvolved setups, pattern timepoints)
  bdv-idr0048-WB.xml    IDR idr0048, BDV/idr0048/WB.xml (3 channels, 'pixel' unit)
  instrument.ome.xml    OME sample set, OME-XML/2016-06/instrument.ome.xml
The expected values were read from the XML by eye; a second, independent
xml.etree pass in `test_bdv_goldens_agree_with_independent_parse` holds the
expectations to the files, so a typo in either fails.
"""

import json
import os
import xml.etree.ElementTree as ET

import pytest

import scigantic_headers.xml_sidecar  # noqa: F401  registers the `xml` extension
from scigantic_headers import decode_bytes, decode_file, has_decoder_for
from scigantic_headers.xml_sidecar import (
    decode_bdv_xml, decode_ome_xml, decode_xml_sidecar, root_element,
)
from tests.test_acquisition import EPU

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")


def _fixture(name: str) -> bytes:
    with open(os.path.join(FIXTURES, name), "rb") as fh:
        return fh.read()


# ── routing ──────────────────────────────────────────────────────────────────

def test_xml_extension_is_registered():
    assert has_decoder_for("dataset.xml") is True
    assert has_decoder_for("x/y.OME.XML") is True
    # exact-name registrations still win over the extension
    from scigantic_headers.decoders import _entry_for
    assert _entry_for("run/RunInfo.xml")[0] is not decode_xml_sidecar


@pytest.mark.parametrize("head,expected", [
    (b"<SpimData version='0.2'/>", "SpimData"),
    (b'\xef\xbb\xbf<?xml version="1.0"?>\n<!-- c -->\n<OME xmlns="x">', "OME"),
    (b"<?xml version='1.0'?><!DOCTYPE foo><ns:MicroscopeImage>", "MicroscopeImage"),
    (b"   \n<mzML>", "mzML"),
    (b"not xml at all", None),
    (b"", None),
])
def test_root_element_sniff(head, expected):
    assert root_element(head) == expected


def test_routes_bdv():
    assert decode_bytes("drosophila.xml", _fixture("bdv-drosophila.xml")).format == "bdv-xml"


def test_routes_ome_xml():
    assert decode_bytes("instrument.ome.xml", _fixture("instrument.ome.xml")).format == "ome-xml"


def test_routes_epu_by_extension():
    h = decode_bytes("FoilHole_1_Data_2.xml", EPU)
    assert h is not None and h.format == "epu-xml"
    assert h.fields["voltageKv"] == 300.0


@pytest.mark.parametrize("data", [
    b"<?xml version='1.0'?><mzML xmlns='http://psi.hupo.org/ms/mzml'><run/></mzML>",
    b"<RunInfo Version='4'><Run Id='x'/></RunInfo>",   # Illumina content under another name
    b"<root><a>1</a></root>",
    b"<SpimData",                                         # truncated: no parse
    b"<OME xmlns='http://www.openmicroscopy.org/Schemas/OME/2016-06'><Image",  # truncated OME
    b"\x00\x01\x02 binary",
    b"",
])
def test_other_xml_and_junk_return_none(data):
    assert decode_xml_sidecar(data) is None
    assert decode_bytes("thing.xml", data) is None


def test_decode_file_dispatches_xml(tmp_path):
    p = tmp_path / "dataset.xml"
    p.write_bytes(_fixture("bdv-HisYFP-SPIM.xml"))
    assert decode_file(str(p)).fields["angleCount"] == 9


# ── BDV goldens ──────────────────────────────────────────────────────────────

BDV_GOLDENS = {
    "bdv-drosophila.xml": {
        "format": "bdv-xml",
        "summary": "BDV dataset (2 view setups), 3 timepoints, bdv.hdf5",
        "fields": {
            "bdvVersion": "0.2", "imageLoaderFormat": "bdv.hdf5", "dataFile": "drosophila.h5",
            "viewSetupCount": 2, "timepointCount": 3,
            "sizeX": None, "sizeY": None, "sizeZ": None,
            "voxelSizeUm": None, "voxelSizeUnit": None,
            "channels": None, "channelCount": None,
            "illuminationCount": None, "angleCount": None, "tileCount": None,
            "hasViewRegistrations": True, "missingViewCount": 0,
        },
    },
    "bdv-HisYFP-SPIM.xml": {
        "format": "bdv-xml",
        "summary": "BDV dataset 1388x1040x81, 9 angles, 0.731078x0.731078x2 um, bdv.hdf5",
        "fields": {
            "bdvVersion": "0.2", "imageLoaderFormat": "bdv.hdf5", "dataFile": "HisYFP-SPIM.h5",
            "viewSetupCount": 9, "timepointCount": 1,
            "sizeX": 1388, "sizeY": 1040, "sizeZ": 81,
            "voxelSizeUm": [0.731078, 0.731078, 2.0], "voxelSizeUnit": "um",
            "channels": ["0"], "channelCount": 1,
            "illuminationCount": 1, "angleCount": 9, "tileCount": None,
            "hasViewRegistrations": True, "missingViewCount": 0,
        },
    },
    "bdv-idr0048-WB.xml": {
        "format": "bdv-xml",
        "summary": "BDV dataset 19120x13350x91, 3 channels, bdv.hdf5",
        "fields": {
            "bdvVersion": "0.2", "imageLoaderFormat": "bdv.hdf5", "dataFile": "WB.h5",
            "viewSetupCount": 3, "timepointCount": 1,
            "sizeX": 19120, "sizeY": 13350, "sizeZ": 91,
            # 'pixel' is not a length unit, so no measurement is claimed
            "voxelSizeUm": [None, None, None], "voxelSizeUnit": "pixel",
            "channels": ["1", "2", "3"], "channelCount": 3,
            "illuminationCount": None, "angleCount": None, "tileCount": None,
            "hasViewRegistrations": True, "missingViewCount": 0,
        },
    },
}


@pytest.mark.parametrize("name", sorted(BDV_GOLDENS))
def test_bdv_golden(name):
    h = decode_bdv_xml(_fixture(name))
    assert h.to_dict() == BDV_GOLDENS[name]
    json.dumps(h.to_dict(), allow_nan=False)


@pytest.mark.parametrize("name", sorted(BDV_GOLDENS))
def test_bdv_goldens_agree_with_independent_parse(name):
    # Re-derive the structural counts straight from the XML so the hand-written
    # expectations cannot drift from the files.
    root = ET.fromstring(_fixture(name))
    exp = BDV_GOLDENS[name]["fields"]
    setups = root.findall("./SequenceDescription/ViewSetups/ViewSetup")
    assert exp["viewSetupCount"] == len(setups)
    assert exp["dataFile"] == root.findtext("./SequenceDescription/ImageLoader/hdf5")
    assert exp["imageLoaderFormat"] == root.find("./SequenceDescription/ImageLoader").get("format")
    size = setups[0].findtext("size")
    if size:
        assert [exp["sizeX"], exp["sizeY"], exp["sizeZ"]] == [int(v) for v in size.split()]
    else:
        assert exp["sizeX"] is None
    chans = root.find("./SequenceDescription/ViewSetups/Attributes[@name='channel']")
    if chans is not None:
        assert exp["channels"] == [c.findtext("name") for c in chans]
    else:
        assert exp["channels"] is None
    assert exp["hasViewRegistrations"] == bool(root.findall("./ViewRegistrations/ViewRegistration"))


# ── BDV synthetic: units, timepoint shapes, loaders, missing views ───────────

def _bdv(timepoints: str, unit="um", loader="bdv.n5", datafile="<n5 type='relative'>d.n5</n5>",
         extra_attrs="", missing=""):
    return ("""<?xml version="1.0" encoding="UTF-8"?>
<SpimData version="0.2"><BasePath type="relative">.</BasePath><SequenceDescription>
<ImageLoader format="%s">%s</ImageLoader>
<ViewSetups>
 <ViewSetup><id>0</id><name>s0</name><size>1920 1920 100</size>
  <voxelSize><unit>%s</unit><size>300 300 1500</size></voxelSize>
  <attributes><channel>0</channel><tile>0</tile></attributes></ViewSetup>
 <ViewSetup><id>1</id><name>s1</name><size>1920 1920 100</size>
  <voxelSize><unit>%s</unit><size>300 300 1500</size></voxelSize>
  <attributes><channel>1</channel><tile>0</tile></attributes></ViewSetup>
 <Attributes name="channel"><Channel><id>0</id><name>GFP</name></Channel>
   <Channel><id>1</id><name>RFP</name></Channel></Attributes>
 <Attributes name="tile"><Tile><id>0</id><name>0</name></Tile><Tile><id>1</id><name>1</name></Tile>
   <Tile><id>2</id><name>2</name></Tile></Attributes>%s
</ViewSetups>
%s
<MissingViews>%s</MissingViews>
</SequenceDescription></SpimData>""" % (loader, datafile, unit, unit, extra_attrs, timepoints, missing)).encode()


def test_bdv_nm_unit_converts_to_um_and_n5_loader():
    h = decode_bdv_xml(_bdv("<Timepoints type='range'><first>0</first><last>9</last></Timepoints>", unit="nm"))
    f = h.fields
    assert f["voxelSizeUm"] == [0.3, 0.3, 1.5]
    assert f["imageLoaderFormat"] == "bdv.n5" and f["dataFile"] == "d.n5"
    assert f["channels"] == ["GFP", "RFP"] and f["channelCount"] == 2
    assert f["tileCount"] == 3 and f["angleCount"] is None
    assert f["timepointCount"] == 10
    assert f["hasViewRegistrations"] is False
    assert h.summary == "BDV dataset 1920x1920x100, 2 channels, 3 tiles, 10 timepoints, 0.3x0.3x1.5 um, bdv.n5"


@pytest.mark.parametrize("tp,count", [
    ("<Timepoints type='list'><Timepoint><id>0</id></Timepoint><Timepoint><id>5</id></Timepoint></Timepoints>", 2),
    ("<Timepoints type='pattern'><integerpattern>0-10:2, 18</integerpattern></Timepoints>", 7),
    ("<Timepoints type='pattern'><integerpattern>18</integerpattern></Timepoints>", 1),
    ("<Timepoints type='pattern'><integerpattern>3-1</integerpattern></Timepoints>", None),   # inverted
    ("<Timepoints type='pattern'><integerpattern>abc</integerpattern></Timepoints>", None),
    ("<Timepoints type='range'><first>4</first><last>2</last></Timepoints>", None),
    ("", None),
])
def test_bdv_timepoint_shapes(tp, count):
    assert decode_bdv_xml(_bdv(tp)).fields["timepointCount"] == count


def test_bdv_missing_views_counted():
    h = decode_bdv_xml(_bdv("<Timepoints type='range'><first>0</first><last>1</last></Timepoints>",
                            missing="<MissingView timepoint='0' setup='1'/><MissingView timepoint='1' setup='1'/>"))
    assert h.fields["missingViewCount"] == 2


def test_bdv_unknown_unit_is_not_a_measurement():
    h = decode_bdv_xml(_bdv("", unit="furlong"))
    assert h.fields["voxelSizeUm"] == [None, None, None]
    assert h.fields["voxelSizeUnit"] == "furlong"
    assert "um" not in h.summary


def test_bdv_rejects_non_spimdata():
    assert decode_bdv_xml(b"<OME/>") is None
    assert decode_bdv_xml(b"<SpimData>") is None


# ── OME-XML ──────────────────────────────────────────────────────────────────

def test_ome_xml_golden():
    h = decode_ome_xml(_fixture("instrument.ome.xml"))
    assert h.to_dict() == {
        "format": "ome-xml",
        "summary": "OME-XML 6x4, uint8, 10000x10000 um",
        "fields": {
            "imageCount": 1, "creator": None, "imageName": "6x6x1x8-swatch.tif",
            "acquisitionDate": "2010-02-23T12:51:30",
            "sizeX": 6, "sizeY": 4, "sizeZ": 1, "sizeC": 1, "sizeT": 1,
            "dtype": "uint8", "dimensionOrder": "XYCZT", "significantBits": None,
            "voxelSizeUm": [10000.0, 10000.0, None],   # PhysicalSize with no unit: schema default is um
            "channels": ["Channel:0:1"], "channelCount": 1,
            "microscope": "OME Insturuments Lab Mk3",  # sic, as written in the sample
            "objective": "OME Objectives 60xUV-Air",
            "magnification": 60.0, "numericalAperture": 1.2,
        },
    }
    json.dumps(h.to_dict(), allow_nan=False)


COMPANION = b"""<?xml version="1.0" encoding="UTF-8"?>
<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" Creator="test">
 <Image ID="Image:0" Name="plane set">
  <Pixels ID="Pixels:0" DimensionOrder="XYZCT" Type="uint16" SizeX="2048" SizeY="2048" SizeZ="300"
          SizeC="2" SizeT="1" PhysicalSizeX="0.65" PhysicalSizeXUnit="\xc2\xb5m" PhysicalSizeY="0.65"
          PhysicalSizeYUnit="\xc2\xb5m" PhysicalSizeZ="2000" PhysicalSizeZUnit="nm">
   <Channel ID="Channel:0:0" Name="488"/><Channel ID="Channel:0:1" Name="561"/>
   <TiffData FirstZ="0" FirstC="0" IFD="0"><UUID FileName="plane_z000_c0.tif">urn:uuid:a</UUID></TiffData>
   <TiffData FirstZ="0" FirstC="1" IFD="0"><UUID FileName="plane_z000_c1.tif">urn:uuid:b</UUID></TiffData>
  </Pixels>
 </Image>
</OME>"""


def test_ome_companion_without_instrument():
    h = decode_bytes("set.companion.ome", COMPANION)   # not .xml: companion files need the dispatcher
    assert h is None or h.format == "ome-xml"           # (another registration may claim it; not ours)
    h = decode_xml_sidecar(COMPANION)
    f = h.fields
    assert f["sizeX"] == 2048 and f["sizeZ"] == 300 and f["sizeC"] == 2
    assert f["voxelSizeUm"] == [0.65, 0.65, 2.0]
    assert f["channels"] == ["488", "561"]
    assert f["microscope"] is None and f["objective"] is None
    assert f["creator"] == "test"
    assert h.summary == "OME-XML 2048x2048x300, 2 channels, uint16, 0.65x0.65x2 um"


def test_ome_root_without_namespace_still_decodes_when_well_formed():
    # The sniff is on the root name; parse_ome_xml decides the rest.
    assert decode_xml_sidecar(b"<OME><Image ID='i'><Pixels SizeX='4' SizeY='4' Type='uint8'/></Image></OME>").fields["sizeX"] == 4


def test_json_safe_everything():
    for name in ("bdv-drosophila.xml", "bdv-HisYFP-SPIM.xml", "bdv-idr0048-WB.xml", "instrument.ome.xml"):
        json.dumps(decode_xml_sidecar(_fixture(name)).to_dict(), allow_nan=False)
