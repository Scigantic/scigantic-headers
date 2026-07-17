# Changelog

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
