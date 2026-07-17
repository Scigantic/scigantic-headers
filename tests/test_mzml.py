"""mzML (mass spectrometry) header decoding.

The samples are built here as well-formed mzML with known values, so the
expected decode is exact and is cross-checked against stdlib ElementTree (a
compliant XML parser) reading the same bytes.
"""

import json
import xml.etree.ElementTree as ET

from scigantic_headers import decode_mzml
from scigantic_headers.decoders import decode_bytes

NS = {"ms": "http://psi.hupo.org/ms/mzml"}


def make_mzml(
    *,
    version="1.1.0",
    n_spectra=42,
    instrument="Q Exactive",
    software=("Xcalibur", "2.2"),
    source="mysample.raw",
    start="2020-05-01T10:00:00Z",
    indexed=True,
    spectra_elems=1,
    pad=0,
):
    """A well-formed mzML document as bytes."""
    sw_id, sw_ver = software
    start_attr = ' startTimeStamp="%s"' % start if start else ""
    filler = "<!-- %s -->" % ("x" * pad) if pad else ""
    spectra = "".join(
        '<spectrum index="%d" id="scan=%d" defaultArrayLength="0">'
        '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/></spectrum>'
        % (i, i + 1)
        for i in range(spectra_elems)
    )
    body = (
        '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="%s" id="s1">'
        '<cvList count="1"><cv id="MS" fullName="PSI-MS" URI="http://x"/></cvList>'
        "<fileDescription>"
        '<fileContent><cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum"/></fileContent>'
        '<sourceFileList count="1"><sourceFile id="RAW1" name="%s" location="file:///data">'
        '<cvParam cvRef="MS" accession="MS:1000768" name="Thermo nativeID format"/></sourceFile>'
        "</sourceFileList></fileDescription>"
        "%s"
        '<softwareList count="1"><software id="%s" version="%s">'
        '<cvParam cvRef="MS" accession="MS:1000532" name="%s"/></software></softwareList>'
        '<instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">'
        '<cvParam cvRef="MS" accession="MS:1001911" name="%s"/>'
        '<componentList count="1"><source order="1">'
        '<cvParam cvRef="MS" accession="MS:1000398" name="nanoelectrospray"/></source></componentList>'
        "</instrumentConfiguration></instrumentConfigurationList>"
        '<run id="run1" defaultInstrumentConfigurationRef="IC1"%s>'
        '<spectrumList count="%d" defaultDataProcessingRef="dp1">%s</spectrumList>'
        "</run></mzML>"
        % (version, source, filler, sw_id, sw_ver, sw_id, instrument, start_attr, n_spectra, spectra)
    )
    head = '<?xml version="1.0" encoding="utf-8"?>'
    if indexed:
        doc = '%s<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">%s</indexedmzML>' % (head, body)
    else:
        doc = head + body
    return doc.encode("utf-8")


def test_decodes_core_fields():
    h = decode_mzml(make_mzml(n_spectra=1234))
    assert h is not None
    assert h.format == "mzml"
    assert h.fields["version"] == "1.1.0"
    assert h.fields["numSpectra"] == 1234
    assert h.fields["instrumentModel"] == "Q Exactive"
    assert h.fields["software"] == "Xcalibur 2.2"
    assert h.fields["sourceFile"] == "mysample.raw"
    assert h.fields["startTimeStamp"] == "2020-05-01T10:00:00Z"


def test_summary_mentions_count_and_instrument():
    h = decode_mzml(make_mzml(n_spectra=1000, instrument="timsTOF"))
    assert "1000 spectra" in h.summary
    assert "timsTOF" in h.summary


def test_dispatch_by_extension():
    blob = make_mzml()
    assert decode_bytes("run.mzml", blob).format == "mzml"
    assert decode_bytes("run.mzML", blob).format == "mzml"


def test_non_indexed_document():
    h = decode_mzml(make_mzml(indexed=False, version="1.0.0"))
    assert h.fields["version"] == "1.0.0"
    assert h.fields["numSpectra"] == 42


def test_software_id_and_version():
    h = decode_mzml(make_mzml(software=("Proteome Discoverer", "2.5")))
    assert h.fields["software"] == "Proteome Discoverer 2.5"


def test_missing_start_timestamp_is_none():
    h = decode_mzml(make_mzml(start=None))
    assert h.fields["startTimeStamp"] is None
    assert h.fields["numSpectra"] == 42  # other fields still decode


def test_rejects_non_mzml():
    assert decode_mzml(b"<?xml version='1.0'?><root>not mzml</root>") is None
    assert decode_mzml(b"just some bytes, definitely not xml") is None
    assert decode_mzml(b"") is None
    # <mzML present but the namespace is not: reject.
    assert decode_mzml(b'<mzML version="1.1.0">no namespace here</mzML>') is None


def test_matches_elementtree():
    blob = make_mzml(n_spectra=99, instrument="Orbitrap Fusion", start="2021-06-06T00:00:00Z")
    h = decode_mzml(blob)
    root = ET.fromstring(blob)
    sl = root.find(".//ms:spectrumList", NS)
    run = root.find(".//ms:run", NS)
    ic = root.find(".//ms:instrumentConfiguration/ms:cvParam", NS)
    mzml_el = root.find(".//ms:mzML", NS)
    assert h.fields["numSpectra"] == int(sl.get("count"))
    assert h.fields["startTimeStamp"] == run.get("startTimeStamp")
    assert h.fields["instrumentModel"] == ic.get("name")
    assert h.fields["version"] == mzml_el.get("version")


def test_decode_file_reads_past_1kib_preamble(tmp_path):
    from scigantic_headers import decode_file

    blob = make_mzml(n_spectra=7, pad=3000)
    assert blob.index(b"<spectrumList") > 1024  # a 1 KiB read would miss it
    p = tmp_path / "big.mzML"
    p.write_bytes(blob)
    h = decode_file(str(p))
    assert h is not None
    assert h.fields["numSpectra"] == 7  # None if only HEADER_BYTES were read


def test_gz_transparent(tmp_path):
    import gzip

    from scigantic_headers import decode_file

    p = tmp_path / "s.mzML.gz"
    with gzip.open(p, "wb") as f:
        f.write(make_mzml(n_spectra=5))
    h = decode_file(str(p))
    assert h is not None
    assert h.fields["numSpectra"] == 5


def test_result_is_json_safe():
    json.dumps(decode_mzml(make_mzml()).to_dict(), allow_nan=False)
