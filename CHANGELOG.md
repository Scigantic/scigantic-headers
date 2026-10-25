# Changelog

## 0.9.0

- NIfTI-2: the 540-byte header (`n+2` / `ni2` magic, 64-bit dims, double
  pixdim) decodes through the same `nifti` decoder as NIfTI-1, so an axis past
  32767 voxels reads correctly. A `niftiVersion` field (1 or 2) is added to
  both. Verified against nibabel-written files, both byte orders.
- NRRD: a decoder for `.nrrd` and detached `.nhdr` headers (teem, ITK, 3D
  Slicer, the Allen CCFv3 atlases): sizes, dtype, encoding, endian, space,
  space directions and derived per-axis spacing, origin, units, kinds, labels,
  measurement frame, data file, `key:=value` pairs, and whether the header fit
  the read. Verified on the real CCFv3 `annotation_25.nrrd` header and against
  pynrrd.
- NWB: `.nwb` files (Neurodata Without Borders, the DANDI archive format)
  decode through the HDF5 reader as their own format: version, identifier,
  session description, start time and id, institution, lab, experimenters,
  keywords, publications, subject, devices, electrode groups and count,
  imaging planes, acquisition series (type, shape, dtype, unit, rate),
  processing modules and interfaces, unit count and columns, interval tables,
  stimuli, embedded schema namespaces and extensions. Targeted lookups rather
  than a walk. Verified against h5py on four DANDI files and pynwb output.
- HDF5: a soft link in a v1 symbol-table group (an undefined object address
  in the symbol node) no longer ends the listing of that node; the entries
  after it are read. An NWB TimeSeries keeps its `electrode` link in front of
  `starting_time`, so the sampling rate of every MIES series was unreadable.
- HDF5: the global heap scan widens in steps (8 KiB, 256 KiB, 1 MiB, never
  past the collection's own size) instead of a fixed 8 KiB window, so strings
  on the large collections a session file builds are read.
- HDF5: identical reads and group listings are memoized within one decode,
  so repeated lookups through the same group no longer count against the
  read cap; no byte is fetched that was not fetched before.

## 0.8.0

- A3M: a decoder for multiple sequence alignments (`.a3m`, `.a3m.gz`), the
  MSA format ColabFold and the MMseqs2 server write and AlphaFold, Boltz,
  Chai, Protenix and OpenFold read. Reports the query name and length, the
  depth within a 4 MiB read with `truncated` when the file goes on, ColabFold's
  `#` chain lengths and cardinalities, mean query coverage, rows with
  insertions, rows whose aligned length does not fit the query, HH-suite
  annotation rows, the source databases of the hits and the first hit's
  header. The server's NUL-separated multi-query files list every query.
  Agrees with an independent parse on 188 real files from ColabFold,
  Protenix, Boltz, OpenFold3 and jackhmmer.

## 0.7.2

- FCS: the TEXT splitter now tells a legacy empty value from an escaped
  delimiter by position. BD CellQuest and other FCS 2.0 / 3.0 writers put an
  empty value down as a doubled delimiter (`KEY\\NEXT\value\`); read as an
  escape, KEY swallowed NEXT and every later pair swapped keyword and value.
  A keyword or value may not begin with the delimiter (FCS 3.1 3.2.7), so a
  doubled delimiter straight after a keyword is an empty value, and inside a
  value it is still a literal (`$P3F/488//10/` is `488/10`). On flowio's own
  `data1.fcs` this reads the 4 empty CellQuest pairs that flowio 1.4.0 itself
  misreads; the other 31 corpus files give identical keywords. A keyword that
  contains the delimiter is the one layout this cannot read, and it costs
  only that pair.
- FCS: a supplemental TEXT segment that lies within the leading read is parsed
  too, so `$CYT` and `$PnS` placed there (optional keywords may be) are found.
  Primary TEXT keywords win.
- FCS: a channel without `$PnN` takes its `$PnS` as the name; when `$PAR` is
  missing or unusable, channels are listed from the `$PnN` keywords present.
- FCS: values that are not UTF-8 (CellQuest's Latin-1 `\xaa`) decode as
  Latin-1 instead of becoming replacement characters; `$DATATYPE` and `$CYT`
  are stripped and `$DATATYPE` upper-cased; a HEADER offset with a stray byte
  after the digits is read; `multipleDatasets` is set when `$NEXTDATA` points
  at a second dataset (Coulter LMD files hold an FCS 2.0 and an FCS 3.0 copy).
- FCS: the splitter works on `bytes.split` runs instead of a per-byte loop
  (1 MiB of TEXT in 13 ms instead of ~300 ms).

## 0.7.1

Fixes found by cross-checking these decoders against other readers.

- KLB: `dataOffset` was 319, the end of the fixed header. The blocks start after
  the offset table that follows it, one 8-byte entry per block, so it is now
  319 + 8 x the block count. On the OME sample `img.klb` that is 447, where a
  bzip2 stream starts, and 447 plus the table's last entry is the file size.
- FCS: FCS 3.2 files were rejected. The leading read is now 1 MiB instead of
  256 KiB, and a TEXT segment longer than that decodes as far as it goes, with
  `textTruncated` set, instead of returning None.
- mzML: the instrument model came from the first cvParam in the instrument
  configuration, which in msconvert output is the ion source (for example
  "electrospray ionization"), because msconvert keeps the model in a
  referenceableParamGroup. The decoder now follows that reference, skips the
  serial number, and never reads the model from the component list. 8 of 54
  real and test mzML files now report the right model; none change otherwise.
- BigDataViewer XML: a timepoint pattern range that ends below its start, such
  as `9-4`, counts down the way spimdata does instead of reading as unknown.
- Tests: the FCS tests were all skipped unless flowio was installed, because of
  a module-level `importorskip`. Only the flowio cross-check skips now.

## 0.7.0

SerialEM `.mdoc` decoding.

- New decoder for `.mdoc`, the acquisition record SerialEM and Thermo Fisher
  Tomography write beside a tilt series, montage or saved movie, dispatched by
  extension with a 1 MiB leading read. Also importable as `decode_mdoc`.
- Returns `voltageKv`, `pixelSizeA`, `doseEA2`, `fractions`, `magnification` and
  `exposureS` under the names the EPU decoder uses, so an mdoc goes straight into
  `compare_acquisition`; plus the software and version, tilt axis, tilt range,
  total dose, acquisition times, the first image's stage position, defocus and
  frame path, and montage settings.
- All four section kinds: `[ZValue]` (stack images), `[FrameSet]` (the file
  beside a single saved movie, including its unbracketed title line),
  `[MontSection]` (montage settings, kept out of the per-image values) and IMOD's
  `[Image = name]` series. NaN and unparseable values decode as absent; CRLF, a
  UTF-8 byte order mark and Latin-1 text are accepted.
- Checked on 116 public mdoc files from GitHub: every value compared agrees with
  mdocfile on the 111 it can read, and the other five (deliberately malformed test
  files) decode with the bad values left out. 46,400 mutations of those files
  raised nothing.

## 0.6.0

Light microscopy and light-sheet formats, and the reader model they needed.

- New `Source` abstraction for pointer-based formats: `Read(seek=budget)`
  registers a decoder that takes a bounded random-access `Source` instead of
  bytes. `FileSource` and `UrlSource` (HTTP Range, 64 KiB cached chunks, a cap
  on requests) implement it; `BytesSource` wraps bytes or a sparse set of ranges.
  `decode_source(key, src)` serves every registered format from a Source. A
  directory decodes through the sidecar file it contains (Zarr, N5 stores).
- TIFF / BigTIFF / OME-TIFF / ImageJ / Micro-Manager / ImSpector decoder
  (`.tif`, `.tiff`, `.btf`, `.tf2`, `.tf8`), with a shared OME-XML reader.
- Zeiss CZI decoder (`.czi`).
- HDF5 walker without libhdf5 (superblock v0-v3, object headers v1/v2, symbol
  tables, dense links and attributes via fractal heaps and v2 B-trees), with
  Imaris `.ims` (incl. the TIFF-based IMS3), BigDataViewer `.h5`, Luxendo
  `.lux.h5` and generic `.h5` / `.hdf5` / `.he5` flavors.
- Leica LIF (`.lif`), Nikon ND2 (`.nd2`), Imspector OBF (`.obf`, `.msr`),
  Keller Lab Block (`.klb`), ICS (`.ics`) decoders.
- Zarr sidecars by file name: `zarr.json`, `.zattrs`, `.zarray`, `.zgroup`
  (OME-Zarr / OME-NGFF 0.4 and 0.5, plain Zarr v2 and v3) and N5
  `attributes.json`. An OME-Zarr store directory merges its level-0 array's
  shape and dtype into the result.
- The `.xml` extension is now dispatched: a sniff of the root element routes
  BigDataViewer `SpimData`, OME-XML (`.ome.xml`, companion files) and EPU
  `MicroscopeImage`; any other XML decodes to None. File-name registrations
  (Illumina) still win. `decode_epu_xml` is unchanged and still importable.
- Behavior change: `has_decoder_for("anything.xml")` is now True, so a tree
  walk opens `.xml` files (1 MiB leading read) and discards the unknown ones.
- Stress-tested before release: a sweep of 390 files on the OME sample server
  by URL (one ND2 Ver3 variant with a 4128-byte chunk-map name was rejected and
  is now read), 200k random and magic-seeded buffers, 14k mutations of the
  real-file fixtures, and hand-built hostile files (self-referencing pointers,
  offsets past EOF, XML entity expansion), kept as `tests/test_adversarial.py`.
- Golden fixtures for the seeking decoders are sparse recordings of real
  sample files (`tests/fixtures/*.sparse.json`, replayed through
  `BytesSource(ranges=...)`); `tests/conftest.py` carries the recorder.

## 0.5.0

- Added `decode_epu_xml`: reads the acquisition settings out of a Thermo Fisher
  EPU image XML (voltage, pixel size as recorded, binning, super-resolution
  factor, fractions, exposure, dose, applied defocus, detector, energy-filter
  slit). Not registered for automatic dispatch, because `.xml` is too generic an
  extension to claim. `epu_pixel_size_candidates` lists the plausible pixel sizes
  for the saved movie, since EPU records the pixel size of the binned grid.
- Added `compare_acquisition` and `Finding`: cross-check voltage, pixel size, dose,
  frame count, Cs and magnification across any number of sources (a collection
  sheet, EPU, a movie header, a RELION or CryoSPARC optics table, a deposit) and
  report disagreements with a severity. With `target_resolution_a` it also flags a
  pixel size whose Nyquist limit cannot reach the target. Motivated by a processing
  run that used 300 kV for data recorded at 200 kV and silently stalled several
  angstrom short of the published resolution.

## 0.4.1

- Added `[project.urls]` (Homepage, Repository, Issues) — the prior release had
  none, despite being live on PyPI.
- Expanded `classifiers` from 2 entries to a full set: development status,
  intended audience, per-version Python classifiers (3.9-3.12, matching the
  CI matrix), scientific-computing/bio-informatics topics, and `Typing ::
  Typed`.
- Added a `py.typed` marker (the package is pure Python with type hints but
  shipped without one) and declared it in `[tool.setuptools.package-data]` so
  it lands in the wheel.
- `authors` now includes the `support@scigantic.com` contact email, matching
  the sibling `scigantic-wwpdb` package.
- Added README badges (CI, PyPI version, license, Python versions).
- GitHub repo topics set: `bioinformatics`, `cryo-em`, `file-formats`,
  `parser`, `scientific-computing`.
- `requires-python = ">=3.9"` was checked for an upper bound and deliberately
  left open: this package is zero-dependency, so there's no third-party
  ceiling to defend against, and Python itself doesn't need one the way a
  dependency does.

## 0.4.0

Internal simplification: the read strategy now lives with the decoder, not in a
table elsewhere. `register_decoder` takes a `read=` argument saying how many bytes
the decoder needs and from which end (`Read(leading=N)`, or `Read(footer=N)` for a
schema in the tail like Parquet). The readers consult the registry and no longer
name any format, so adding a format really is one pure function plus one call.

- New public API: `Read` (the read spec) and `read_for(key)` (the spec the readers
  will use for a key). Additive; existing code is unaffected. `register_decoder`
  gained an optional `read=` keyword.
- Removed the per-format read-size tables and the Parquet-specific footer re-read
  from `sources`, which is now format-agnostic.
- Behavior change: a Parquet footer larger than the 1 MiB trailing read now decodes
  to `None` (no context), the same as any other header that does not fit, instead
  of a precise re-read. This affects only pathologically large schemas and removes
  an inconsistency where local files got the re-read but remote (HTTP Range) did not.
- The HTTP User-Agent for URL reads now reports the installed version instead of a
  hardcoded `0.1`.

## 0.3.0

Broadens coverage from cryo-EM/imaging/tabular to the formats a wet lab produces
across mass spec, sequencing, and structural biology. New decoders:

- mzML (mass spectrometry): spectrum count, instrument model, software, source
  file, run start time.
- FASTQ: instrument/run/flowcell/lane from the Illumina read name, read length.
- Illumina run metadata: RunInfo.xml / RunParameters.xml, dispatched by file
  name (a new file-name registry, since `.xml` alone is too generic).
- VCF, SAM: header fields (version, samples/references, meta-line counts).
- PDB, mmCIF: structure headers (id, classification, method, resolution).
- GenBank: LOCUS line plus definition, accession, version, organism.
- GFF3/GTF, BED: annotation and interval headers.
- DICOM: modality, dimensions, bit depth, manufacturer. Technical fields only,
  never patient identifiers.

Now 17 formats. Fuzz covers every registered decoder, including the file-name
dispatch path.

## 0.2.0

- FCS (flow cytometry) decoder: parameter and event counts, data type, and
  per-channel names, read from the HEADER offset table and TEXT segment. The
  leading read grows for `.fcs`, whose keyword segment can exceed 1 KiB.
  Cross-checked against flowio.
- README and source docstrings describe only the library; references to a
  surrounding system have been removed.

## 0.1.3

- `__version__` now derives from the installed package metadata instead of a
  hardcoded string, so it always matches the release and cannot drift. (0.1.1
  and 0.1.2 reported `0.1.0` from `scigantic_headers.__version__`.)

## 0.1.2

- Relicensed under MIT (was proprietary evaluation-only).

## 0.1.1

- README describes installing from PyPI with a version pin, not vendoring a
  copy into the monorepo (which no longer happens).
- Package summary lists the actual formats and drops the marketing phrasing.

No code changes.

## 0.1.0

First release.

- Header decoders (pure, zero-dependency): MRC/MRCS (cryo-EM), NPY (numpy
  arrays), NIfTI-1 (neuroimaging, `.nii` and `.nii.gz`), CryoSPARC `.cs`, and
  Parquet (column names, physical types, row count, read from the footer).
- Acquisition-optics readers: RELION STAR and CryoSPARC `.cs` (pixel size,
  voltage, spherical aberration). These recover the pixel size a raw MRC movie
  header omits. `read_session_optics(path)` finds a data file's optics file in
  its session directory and returns the pixel size in one call.
- Bounded, `.gz`-transparent readers; bounded-parallel batch over files or HTTP
  Range.
- Verified against reference implementations (numpy, nibabel) and real EMPIAR
  data; MRC held to a golden fixture shared with the backend TypeScript twin.
