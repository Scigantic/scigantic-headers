# scigantic-headers

Read the metadata of a scientific file and return its fields (dimensions, data
type, pixel size, columns, row count) as a dict, without reading the rest of the
file. Decodes formats across cryo-EM, imaging, arrays, tabular data, flow
cytometry, mass spec, sequencing, genome annotation, and structural biology (the
full list is below). The decode
functions have no dependencies; separate reader functions fetch the bytes from a
file or a URL (the leading bytes for most formats, the trailing bytes for
Parquet, whose schema is in a footer).

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
- **NIfTI-1**: neuroimaging volumes, `.nii` and `.nii.gz`. Cross-checked against
  nibabel. Big- and little-endian.
- **CryoSPARC `.cs`**: a structured-array dataset. Reports record count and the
  field schema.
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
- **GenBank**: annotated sequence records. LOCUS (name, length, molecule type,
  topology), definition, accession, version, and organism.
- **GFF3 / GTF**: genome annotation. Version and a preview of the sources and
  feature types.
- **BED**: genome intervals. BED flavor (column count) and track-line presence.
- **DICOM**: medical images. Modality, dimensions, bit depth, manufacturer, and
  study/series description. Reads technical fields only, never patient
  identifiers.

`.gz` is transparent. The reader decompresses just the leading block, so a
gzipped `.nii.gz` or `.mrc.gz` decodes without inflating the whole file, and
dispatch sees the inner format.

Adding a format is a pure `bytes -> DecodedHeader | None` function plus one
`register_decoder` call. That call also declares how many bytes the decoder needs
and from which end (`read=Read(leading=...)`, or `read=Read(footer=...)` for a
schema in the tail like Parquet), so the read strategy lives next to the decoder,
not in a table elsewhere. Nothing about the dispatch, the readers, or the batch
path is format-specific: the readers ask the registry how to read and never name
a format. NPY and NIfTI landed with zero plumbing changes.

## Pixel size for raw movies

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
      genbank.py    GenBank flat-file header decode
      gff.py        GFF3 / GTF annotation-header decode
      bed.py        BED interval-file decode
      dicom.py      DICOM technical-metadata decode (no patient data)
      sources.py    read leading bytes from file or URL; bounded parallel batch
      star.py       RELION STAR optics reader
      cryosparc.py  CryoSPARC .cs optics reader
      optics.py     read_session_optics: find a data file's optics file
      cli.py        scigantic-headers <file|url|--dir>
      benchmark.py  scigantic-headers-bench, measures the speed levers
    tests/          pytest (202 tests, incl. fuzz + golden fixtures)

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
    decode_paths(iter_decodable_files("/data/sessions"))    # a tree, parallel
    decode_urls([...], workers=8)                           # remote, by Range

    # CLI
    scigantic-headers session/frame.mrc
    scigantic-headers https://example.org/path/img.mrcs
    scigantic-headers --dir /data/sessions --workers 8

    pytest
    scigantic-headers-bench            # reproduce the numbers above

## Limitations

- **Pointer-based formats are not supported.** HDF5 and TIFF store their layout
  behind internal pointers, so reading the header or footer alone is not enough.
  Parquet (a footer) is supported; formats that need to chase offsets through the
  file are not. Plain gzip is handled (the reader inflates the leading block),
  but BAM's BGZF framing and per-record structure are not.
- **Output is not a standard.** The fields are an ad-hoc dict, not an
  interchange format such as Allotrope ASM. Map the dict to a standard if you
  need one.

## License

MIT. See LICENSE.
