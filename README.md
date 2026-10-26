# scigantic-headers

<p align="center">
    <a href="https://github.com/Scigantic/scigantic-headers/actions/workflows/ci.yml">
        <img alt="CI" src="https://github.com/Scigantic/scigantic-headers/actions/workflows/ci.yml/badge.svg" /></a>
    <a href="https://pypi.org/project/scigantic-headers/">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/scigantic-headers" /></a>
    <a href="https://pypi.org/project/scigantic-headers/">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/scigantic-headers" /></a>
    <a href="https://github.com/Scigantic/scigantic-headers/blob/main/LICENSE">
        <img alt="License" src="https://img.shields.io/github/license/Scigantic/scigantic-headers" /></a>
</p>

Read the metadata of a scientific file and return its fields (dimensions, data
type, pixel size, columns, row count) as a dict, without reading the rest of the
file. Decodes formats across cryo-EM, light microscopy and light-sheet imaging,
arrays, tabular data, flow cytometry, mass spec, sequencing, genome annotation,
structural biology, and electrophysiology (the full list is below). The decode functions have no
dependencies; separate reader functions fetch the bytes from a file or a URL:
the leading bytes for most formats, the trailing bytes for Parquet (whose schema
is in a footer), and a bounded handful of reads wherever the pointers lead for
TIFF, HDF5, CZI, LIF and ND2 (whose metadata sits behind internal offsets).

```python
from scigantic_headers import decode_file
hdr = decode_file("session/frame.mrc")
hdr.summary            # "MRC stack 4096x4096x16, float32"
hdr.fields["dtype"]    # "float32"
hdr.fields["nz"]       # 16   (frames)
```

Zero runtime dependencies and header-only reads, so it runs air-gapped, on
amd64 or arm64, and installs into a slim image with nothing else to pull.

## Install

    pip install scigantic-headers

Published to PyPI on a `v*` tag by `.github/workflows/publish.yml`. Pin the
version in anything that consumes it (`scigantic-headers==X.Y.Z`). With no
package index at build time, install from a wheel or a checkout instead:

    pip install ./scigantic-headers

## Why it exists

You often want to know what a file is, its shape, dtype, and key acquisition
fields, without reading the whole thing. A scientific file's first ~1 KiB is a
structured header that already says so. This reads that, and nothing more.

Decoders today:

- **MRC / MRCS**: cryo-EM micrographs, movie stacks, tilt series, EMDB maps.
- **NPY**: NumPy arrays. Cross-checked against numpy.
- **NIfTI-1 and NIfTI-2**: neuroimaging volumes, `.nii` and `.nii.gz`. Both
  header layouts (348 and 540 bytes), big- and little-endian; `niftiVersion`
  says which. Cross-checked against nibabel.
- **NRRD** (`.nrrd`, detached `.nhdr`): the volume format of teem, ITK, 3D
  Slicer and the Allen Institute's mouse brain atlases (CCFv3). Version,
  dimension, sizes, dtype, encoding, endianness, space, space directions and
  the per-axis spacing derived from them (or `spacings`), origin, units, kinds,
  labels, measurement frame, the data file of a detached header, and the
  `key:=value` pairs Slicer writes. Verified on the real CCFv3 annotation
  header and cross-checked against pynrrd.
- **CryoSPARC `.cs`**: a structured-array dataset. Reports record count and the
  field schema.
- **SerialEM `.mdoc`**: the acquisition record SerialEM (and Thermo Fisher
  Tomography) writes beside a tilt series, montage or saved movie. Voltage, pixel
  size, dose, frame count, exposure, magnification and tilt axis, the tilt range
  and total dose over every image, the first image's stage position and defocus,
  and the montage settings. Handles tilt series (`[ZValue]`), single-movie files
  (`[FrameSet]`), montages (`[MontSection]`) and IMOD image series (`[Image]`).
  Cross-checked against mdocfile on 111 public mdoc files.
- **Parquet**: reads the footer, not a leading header, and returns the column
  names, physical types, and row count. Cross-checked against pyarrow.
- **FCS**: flow cytometry (Flow Cytometry Standard). Returns parameter and event
  counts, data type, and per-channel names. Cross-checked against flowio.
- **mzML**: mass spectrometry. Returns the spectrum count, instrument model,
  software, source file, and run start time from the XML preamble. Cross-checked
  against a compliant XML parse.
- **FASTQ**: sequencing reads. Parses the Illumina read-name convention in the
  first record for instrument, run, flowcell, and lane, plus the read length.
- **Illumina run**: RunInfo.xml / RunParameters.xml. Instrument, flowcell, run
  id, date, and read structure. Dispatched by file name, not extension.
- **VCF**: variant calls. Version, sample names, reference, and meta-line counts
  from the header.
- **SAM**: alignment header. Version, sort order, reference sequences, read
  groups, and programs.
- **PDB**: structure header. PDB id, classification, deposition date, title,
  experiment method, and resolution.
- **mmCIF / CIF**: structures (the modern PDB format). Entry id, title,
  keywords, experiment method, deposition date, and resolution (including the
  cryo-EM resolution item).
- **A3M**: multiple sequence alignments, the MSA format ColabFold and the
  MMseqs2 server write and AlphaFold, Boltz, Chai, Protenix and OpenFold read.
  Query name and length, depth (sequences in the read, `truncated` when the
  file goes on past it), ColabFold's `#` chain lengths and copy counts, the
  mean fraction of the query each row covers, rows with insertions, rows whose
  aligned length does not fit the query (an un-realigned or mis-staged MSA),
  HH-suite `ss_pred`/`ss_conf` annotation rows, the databases the hits came
  from (UniRef, UniProtKB, UniParc, MGnify, SRA) and the first hit's header.
  The server's NUL-separated multi-query files report every query. Checked
  against an independent parse on 188 ColabFold, Protenix, Boltz, OpenFold3
  and jackhmmer files.
- **GenBank**: annotated sequence records. LOCUS (name, length, molecule type,
  topology), definition, accession, version, and organism.
- **GFF3 / GTF**: genome annotation. Version and a preview of the sources and
  feature types.
- **BED**: genome intervals. BED flavor (column count) and track-line presence.
- **DICOM**: medical images. Modality, dimensions, bit depth, manufacturer, and
  study/series description. Reads technical fields only, never patient
  identifiers.

Electrophysiology and BIDS (added in 0.10.0). What an EEG, iEEG or Neuropixels
rig writes beside its raw samples, and the sidecars that describe a dataset
organised to the Brain Imaging Data Structure:

- **EDF / EDF+ / BDF / BDF+** (`.edf`, `.bdf`): the clinical and research
  electrophysiology container. Flavour and continuity (`EDF+C` / `EDF+D`), bit
  depth, start date and time (resolved to a four-digit year from the EDF+
  recording field), the EDF+ patient and recording subfields, data record
  count and duration, recording length, and per signal the label, units,
  physical and digital range, samples per record, derived sampling rate,
  transducer and prefiltering. The annotations channel is named and kept out
  of the recording's sampling rate. Cross-checked against pyedflib.
- **BrainVision** (`.vhdr`, `.vmrk`): Brain Products' header and marker files,
  one of the three container formats BIDS allows for EEG and iEEG. Format
  version, data file and orientation, channel count, sampling interval and
  rate, binary format and byte order, data points and duration, per-channel
  name, reference, resolution and unit, and whether electrode coordinates are
  present; from the marker file, the marker count, the count per type, the
  first markers and the segment start timestamp. Cross-checked against mne on
  64 public files.
- **SpikeGLX** (`.meta`): the run metadata SpikeGLX writes beside every
  Neuropixels binary, which is the only description that headerless int16
  matrix has. Stream type and kind, sample rate, channel counts split into AP,
  LF and sync, run name with gate and trigger index, file size and duration,
  derived sample count, probe type, model, part number and serial, slot, port
  and dock, and the microvolts per bit for the two probe generations the
  metadata reference defines it for. The per-channel tables are counted, not
  copied. Cross-checked against neo.
- **Neuralynx** (`.ncs`, `.nev`, `.nse`, `.ntt`, `.nvt`): the 16 KiB ASCII
  header every Cheetah and Pegasus file opens with, in both dialects (the
  older `##` comment lines with `-CheetahRev`, and the newer properties with
  `-ApplicationName` and `-TimeCreated`). File type and record kind, record
  size, application and version, channel name, sampling rate, AD range and
  microvolts per bit, DSP filter settings, and the open and close times.
  Checked against python-neo's header grammar.
- **Open Ephys** (`structure.oebin`): the JSON manifest of a binary-format
  recording directory, the only description of its headerless `.dat`
  matrices. GUI version, continuous streams with folder, stream and processor
  names, sample rate, channel count and the scaling where it is uniform
  across the stream, plus the event and spike streams.
- **BIDS sidecars** (`.json`, `.tsv`): `dataset_description.json`, the JSON
  metadata sidecars, the column data dictionaries, and the tabular files
  (events, participants, sessions, samples, scans, channels, electrodes,
  optodes). Both extensions are generic, so a file is claimed only when it
  carries a key the standard defines or one of its column signatures; key
  names are kept exactly as the standard spells them. Swept over the whole
  bids-standard/bids-examples corpus.

Light microscopy and light-sheet (added in 0.6.0; most of these follow internal
pointers, see "Pointer-based formats" below):

- **TIFF / BigTIFF / OME-TIFF** (`.tif`, `.tiff`, `.btf`, `.tf2`, `.tf8`):
  dimensions, dtype, compression, photometric, tiling, page count (chain
  followed up to a cap, with the truncation stated). The ImageDescription is
  read where it lives and classified: **OME-XML** (sizes X/Y/Z/C/T, physical
  sizes in um, channels, instrument, objective, acquisition date, via the
  shared OME-XML reader), **ImageJ** hyperstacks (slices/channels/frames,
  spacing and unit), **Micro-Manager** (summary JSON, private header, or tag
  51123), ImSpector/LaVision stacks by software name. Cross-checked against
  tifffile.
- **Zeiss CZI**: file header, metadata segment (XML, wherever the writer put it,
  including 490 MB into a light-sheet file), directory and attachment counts.
  Sizes X/Y/Z/C/T plus scenes, tiles, illuminations, views, rotations, phases;
  pixel type, scaling in um, channels and exposure, acquisition mode (SPIM,
  lattice light-sheet, widefield), microscope, objective, software. Verified on
  LSM 980, Axioscan, Lightsheet Z.1 and Lattice Lightsheet 7 files; cross-checked
  against czifile.
- **HDF5** and the imaging formats built on it, without libhdf5: superblocks
  v0 to v3, object headers v1 and v2, symbol tables and dense (fractal heap plus
  v2 B-tree) links and attributes, dataspace, datatype, layout, filters.
  Attributes are read from the object header the links point to, never scanned
  for, so a relocated header's stale copy cannot leak in.
  - **Imaris `.ims`**: X/Y/Z from DataSetInfo/Image, voxel size from the extents,
    channels with emission wavelengths, resolution levels, timepoints, the
    level-0 block dtype/chunks/compression, objective, recording date. An `.ims`
    that is a TIFF (Imaris IMS3) is reported as such through the TIFF decoder.
  - **BigDataViewer `.h5`**: setup and timepoint counts, the per-setup
    resolutions and subdivisions tables, the level-0 cells block shape (z, y,
    x), dtype and compression. Two Range requests describe a 51 GB file.
  - **Luxendo `.lux.h5`**: the generic structure labelled, with the `Data`
    block shape. No public sample was available; shape only.
  - **NWB `.nwb`** (Neurodata Without Borders, the format of the DANDI
    archive): version, identifier, session description and start time,
    session id, institution, lab, experimenters, keywords, publications, the
    subject (id, species, sex, age, genotype, strain), devices, electrode
    groups and electrode count, imaging planes, every acquisition series with
    its type, shape, dtype, unit and rate, processing modules and their
    interfaces, unit count and columns, interval tables with row counts,
    stimuli, and the schema namespaces the file embeds, with any extension
    (`ndx-*`) named. Targeted lookups, not a walk, so a session file with
    thousands of objects still reads in a few hundred bounded reads.
    Verified against h5py on four DANDI files from four labs and on
    pynwb-written files.
  - **NWB on Zarr** (`.nwb.zarr`, a directory): the same NWB model written by
    `hdmf-zarr` into a Zarr store. Version, namespace, object id, devices,
    electrode groups and count, imaging planes, every acquisition series with
    its type, shape, dtype, unit and rate, processing modules and their
    interfaces, unit count and columns, interval tables with row counts and
    columns, stimuli, and the cached schema namespaces with the version each
    was cached at. Both container versions: Zarr v3 reads the root
    `zarr.json`, Zarr v2 the consolidated `.zmetadata`. Unlike the HDF5
    flavour this cannot report the scalar strings (identifier, session
    description, subject): in Zarr those are chunk data behind a compressor
    this library does not carry, so they are absent and `stringValuesUnread`
    names the paths that hold them. Verified against a real DANDI 000719
    store and an hdmf-zarr-written one.
  - **generic `.h5` / `.hdf5` / `.he5`**: root attributes (fixed and
    variable-length strings, numbers), groups, and datasets with shape, dtype,
    layout, chunks and filters, bounded and with the truncation stated.
- **KLB** (Keller Lab Block): the 319-byte header: X/Y/Z/C/T, dtype, bzip2 /
  zlib / none, block size, pixel size (KLB stores no unit; micrometres by
  convention), the 256-character metadata field.
- **ICS** (Image Cytometry Standard, `.ics` + `.ids`): version, layout order and
  sizes, dtype, byte order, compression, scale converted to um, units, labels,
  sensor/history entries (Huygens writes the objective NA there).
- **Leica LIF**: the UTF-16 XML header: image count and names, dimensions by
  DimID, dtype from channel resolution, voxel size from physical lengths,
  channels, objective, magnification, NA, software version. Cross-checked
  against readlif.
- **Nikon ND2**: the chunk map at the end of the file leads to the attribute,
  calibration, text-info and metadata chunks (both the Ver2 XML and the Ver3
  LiteVariant encodings are decoded for the needed values): width, height,
  components, bit depth, frame count, experiment loops (Z/T/positions), um per
  pixel, z step, channels, objective and NA, acquisition date, camera,
  compression. About 125 KB of a 218 MB file is read. Cross-checked against the
  nd2 package. Legacy JPEG 2000 ND2 is recognised only.
- **Imspector OBF / `.msr`** (Abberior, LaVision BioTec): file header, stack
  chain (count, names, rank, sizes, dtype, compression, physical lengths to um
  for the axes the footer marks as metres, axis labels). Verified on five real
  Abberior STED `.msr` files (hoerlteam/msr-reader test data) and cross-checked
  against obffile; no ImSpector-written LaVision `.msr` was available.
- **Zarr stores**: OME-Zarr / OME-NGFF 0.4 (`.zattrs`) and 0.5 (`zarr.json`,
  Zarr v3): axes, levels, level-0 scale to um, channels, plate/labels markers;
  plain `.zarray` / `.zgroup` arrays and groups (shape, chunks, dtype, codec);
  **N5** `attributes.json` (dimensions, block size, dtype, compression, pixel
  resolution). These are dispatched by file name. A store *directory* decodes
  through the sidecar it contains, and an OME-Zarr store also merges its level-0
  array's shape and dtype, so `decode_file("image.zarr")` reads like one image.
  Cross-checked against zarr on IDR data.
- **XML sidecars**: the generic `.xml` extension is now claimed by a dispatcher
  that sniffs the root element and returns None for anything it does not know:
  **BigDataViewer / BigStitcher `SpimData`** (loader format, data file, setups,
  timepoints, size, voxel size, channels, illuminations, angles, tiles,
  registrations), **OME-XML** sidecars and companion files (`.ome.xml`, same
  fields as OME-TIFF), and **EPU** `MicroscopeImage` (previously opt-in only).
  Illumina RunInfo.xml / RunParameters.xml keep their file-name dispatch, which
  wins over the extension.

`.gz` is transparent. The reader decompresses just the leading block, so a
gzipped `.nii.gz` or `.mrc.gz` decodes without inflating the whole file, and
dispatch sees the inner format.

Adding a format is a pure `bytes -> DecodedHeader | None` function plus one
`register_decoder` call. That call also declares how the bytes are fetched
(`read=Read(leading=...)`, `read=Read(footer=...)` for a schema in the tail like
Parquet, or `read=Read(seek=budget)` for a pointer-based format, whose decoder
then takes a `Source` instead of bytes), so the read strategy lives next to the
decoder, not in a table elsewhere. Nothing about the dispatch, the readers, or
the batch path is format-specific: the readers ask the registry how to read and
never name a format. NPY and NIfTI landed with zero plumbing changes.

## Pointer-based formats

TIFF keeps its image directories wherever the writer put them (libtiff and
Bio-Formats put them after the pixel data); HDF5 is a graph of object headers
reached from a superblock; a CZI's metadata segment can sit 490 MB into the
file; an ND2 starts from a chunk map at its end. No fixed leading or trailing
read decodes these, so a decoder for one of them takes a `Source`: a
random-access view with `read_at(offset, n)` and `read_tail(n)` that is total
(short reads instead of errors) and bounded three ways: a byte budget declared
with the decoder, a cap on reads, and, over HTTP, a cap on Range requests
(`URL_MAX_REQUESTS`). Reads from a URL are coalesced into 64 KiB cached chunks,
so a walk over neighbouring structures costs one round trip. The decode is
still a few KB of a file of any size; it is just not the first few KB. Measured
on the OME sample server:

    BDV drosophila.h5 (620 MB)                2 Range requests
    BDV WB.h5 (51 GB)                         2 Range requests
    Imaris Convallaria .ims (7 MB)            6 Range requests
    CZI MouseBrain light-sheet (492 MB)       5 Range requests (metadata at 492 MB)
    OME-TIFF tubhiswt (270 KB, 1 page)        3 Range requests
    OME-TIFF multi-channel .btf (21 pages)    3 Range requests (chain not chased)
    Leica LIF PR2729 (227 KB)                 2 Range requests
    Nikon ND2 sample_image (218 MB)           6 Range requests
    KLB img.klb                               2 Range requests

A plain TIFF whose directories trail each plane costs a request per page until
the chain cap stops it. An OME-TIFF states its plane count in the XML, so its
chain is not chased at all: `planeCount` comes from the OME block and
`pageCount` stops at the first directory, marked truncated.

`decode_bytes(key, data)` still works for these formats (the bytes become an
in-memory Source), and `BytesSource(ranges={offset: bytes})` replays a sparse
recording, which is how the golden fixtures hold a 500 MB file's few KB of
metadata (`tests/fixtures/*.sparse.json`).

## Pixel size for raw movies

### Cross-check the acquisition settings

Voltage, pixel size, dose and frame count usually arrive from several places typed
or exported by different people. A wrong one does not fail; it quietly caps the
resolution. `compare_acquisition` lines the sources up and reports disagreements:

```python
from scigantic_headers import compare_acquisition, decode_epu_xml, decode_file, read_star_optics

epu = decode_epu_xml(open("FoilHole_..._Data_....xml", "rb").read()).fields
mdoc = decode_file("TS_01.mrc.mdoc").fields      # SerialEM, if the session used it
for f in compare_acquisition({
    "sheet": {"voltageKv": 300, "pixelSizeA": 0.84, "doseEA2": 42.4, "fractions": 40},
    "epu": epu,
    "serialem": mdoc,
    "relion": read_star_optics("Import/movies.star"),
}, target_resolution_a=3.3):
    print(f.severity, f.field, f.message)
```

It does not decide which source is right. EPU records the pixel size of the binned
readout grid, so `epu_pixel_size_candidates(epu)` lists the likely values for the
saved movie; compare with the movie's own dimensions. `decode_epu_xml` is not
auto-dispatched (`.xml` is too generic), so call it on a file you know is EPU's.
An mdoc's `doseEA2` is the dose of one image (one tilt, or one movie), the value to
compare with a per-movie dose; `totalDoseEA2` is the sum over the series.

A raw MRC movie header has no pixel size: CELLA is 0, so `decode_*` returns
`pixelSizeA` as `None`. The value is in the RELION STAR or CryoSPARC `.cs` file
the workflow writes in the session, and the library reads it. `read_session_optics`
finds that file next to (or above) the data file and returns the optics:

```python
from scigantic_headers import decode_file, read_session_optics

hdr = decode_file("Movies/frame.mrc")       # header fields; pixelSizeA is None
opt = read_session_optics("Movies/frame.mrc")
# {'pixelSizeA': 1.05, 'voltageKv': 300.0, 'source': 'relion-star:particles.star'}
```

`read_star_optics(path)` and `read_cryosparc_optics(path)` read a specific optics
file when you already know which one. The RELION 3.1 `data_optics` loop, the
older key-value form, optics columns in the main data loop, and the legacy
detector-pixel-size / magnification pair are all handled.

## Layout

    src/scigantic_headers/
      decoders.py   pure core: registry + MRC / NPY / NIfTI / .cs decoders
      parquet.py    Parquet footer decode (small Thrift-compact reader)
      fcs.py        FCS (flow cytometry) header + TEXT-segment decode
      mzml.py       mzML (mass spec) XML-preamble metadata decode
      fastq.py      FASTQ first-record / read-name decode
      illumina.py   Illumina RunInfo.xml / RunParameters.xml decode
      vcf.py        VCF header decode
      sam.py        SAM alignment-header decode
      pdb.py        PDB structure-header decode
      mmcif.py      mmCIF / CIF structure-header decode
      a3m.py        A3M multiple-sequence-alignment decode
      genbank.py    GenBank flat-file header decode
      gff.py        GFF3 / GTF annotation-header decode
      bed.py        BED interval-file decode
      dicom.py      DICOM technical-metadata decode (no patient data)
      tiff.py       TIFF / BigTIFF / OME-TIFF / ImageJ / Micro-Manager (seeking)
      ome_xml.py    shared OME-XML reader (OME-TIFF, .ome.xml, companion files)
      czi.py        Zeiss CZI segments + metadata XML (seeking)
      hdf5.py       HDF5 walker; Imaris, BigDataViewer, Luxendo, generic (seeking)
      lif.py        Leica LIF XML header (seeking)
      nd2.py        Nikon ND2 chunk map + LiteVariant chunks (seeking)
      obf.py        Imspector OBF / .msr stack chain (seeking)
      klb.py        Keller Lab Block header
      ics.py        Image Cytometry Standard header
      zarr.py       OME-Zarr / Zarr v2+v3 / N5 sidecars (by file name)
      xml_sidecar.py  .xml dispatcher: BDV SpimData, OME-XML, EPU
      mdoc.py       SerialEM / Tomography .mdoc acquisition record
      sources.py    read bytes from file or URL (leading, footer, or seeking); bounded parallel batch
      star.py       RELION STAR optics reader
      cryosparc.py  CryoSPARC .cs optics reader
      optics.py     read_session_optics: find a data file's optics file
      cli.py        scigantic-headers <file|url|--dir>
      benchmark.py  scigantic-headers-bench, measures the speed levers
    tests/          pytest (507 tests, incl. fuzz, golden and sparse real-file fixtures)

## Robustness

A decoder is handed the leading bytes of arbitrary files, so it must be total:
for any input it returns None or a valid, JSON-safe DecodedHeader. It never
raises and never emits NaN or inf. This is enforced, not assumed.

- **Fuzz** (`test_fuzz.py`): thousands of random and magic-seeded byte buffers
  through every decoder. Each result must serialize with `allow_nan=False`.
- **Golden fixtures** (`fixtures/mrc-cases.json`): hex input to exact decode. A
  change to a decoder that the fixture does not expect fails a test rather than
  passing silently.
- Non-finite floats (a garbage NIfTI pixdim, an inf MRC cell) sanitize to null.

## Speed

The decode is microseconds: a 1 KiB read and a few `struct` unpacks. Optimizing
that is a rounding error. All the speed is in I/O, and two things carry it. Both
are measured, not asserted (run `scigantic-headers-bench`).

**Read the header, not the file.** The read is bounded (`HEADER_BYTES`, 1 KiB, for
most formats; a few hundred KB for text formats with long headers, declared per
decoder), so its cost is independent of file size and a multi-GB movie decodes as
fast as a small one.

    header read + decode (1 KiB):        15 us
    full-file read + decode (200 MB):  22.7 ms      (about 1,500x cheaper)

**Parallelize the I/O, bounded.** Header reads across files are independent and
I/O-bound, so a thread pool overlaps their latency (the read releases the GIL).
Measured over 8 headers pulled by HTTP Range from a remote archive:

    serial   (1 worker):   25.7 s
    parallel (8 workers):   0.7 s        (about 35x)

That 35x is large because the server's per-request latency is about 3 s, so
overlapping eight requests recovers a lot of dead wait. On low-latency local
storage the speedup is smaller. The win scales with per-read latency, which is
why the pool is bounded and tunable (default 8; past a point more connections
hit server throttling or disk-queue thrash and get slower).

Two more, by construction:

- **Filter before I/O.** `has_decoder_for` is a cheap extension check. The batch
  and walk paths use it to skip a file before any syscall.
- **Zero deps, fast import.** The dtype table is a plain dict, not numpy. There
  is nothing to install and nothing to import, and it runs air-gapped.

The biggest system-level lever lives in the caller, not here: decode once when a
file lands, cache the small result, and answer every later query from it without
re-decoding.

## Use

    from scigantic_headers import decode_file, decode_paths, decode_urls, iter_decodable_files

    decode_file("s/frame.mrc")                              # one local file
    decode_file("acq/image.zarr")                           # a Zarr or N5 store directory
    decode_paths(iter_decodable_files("/data/sessions"))    # a tree, parallel
    decode_urls([...], workers=8)                           # remote, by Range

    # CLI
    scigantic-headers session/frame.mrc
    scigantic-headers https://example.org/path/img.mrcs
    scigantic-headers https://example.org/path/sample.ims     # HDF5 over Range
    scigantic-headers --dir /data/sessions --workers 8

    pytest
    scigantic-headers-bench            # reproduce the numbers above

## Limitations

- **HDF5 is walked, not fully parsed.** Groups and attribute sets are read
  from compact storage, symbol tables, and dense storage with v2 B-trees of
  depth 0 or 1 and fractal heaps with direct blocks; a group with thousands of
  links (a BDV file with thousands of timepoints) is listed up to a cap, with
  the truncation stated. Huge heap objects, virtual datasets, external links
  and soft links are not followed. Attribute values are decoded for strings and
  small numeric arrays; others are reported by type.
- **Vendor formats without a public sample** (Micro-Manager-native TIFF,
  ImSpector TIFF, Luxendo `.lux.h5`) were implemented from the published
  layouts and tested on constructed files; their tests say so. Olympus
  OIB (an OLE compound document) and OIR, Hamamatsu DCIMG, and raw `.bin` /
  `.dat` / `.raw` exports (which have no header to read) are not decoded.
- **Neuralynx has no public sample** that was reachable, so its fixtures are
  headers built to the per-version grammar documented in python-neo's
  `NlxHeader` and parsed by neo, not copies of a rig's output. The binary
  records after the header are not read, so a file's record count is not
  reported.
- **Legacy ND2** (JPEG2000 based) is recognised but not described.
- BAM's BGZF framing and per-record structure are not handled; plain gzip is.
- **A Zarr store's chunk data is never opened**, so an NWB-on-Zarr record
  carries the structure, shapes, units and rates but not the scalar strings
  the HDF5 flavour reads out of its object headers. A store written without
  consolidated metadata is named and versioned but not described, since
  walking it would mean one read per node.
- **Output is not a standard.** The fields are an ad-hoc dict, not an
  interchange format such as Allotrope ASM. Map the dict to a standard if you
  need one.

## License

MIT. See LICENSE.
