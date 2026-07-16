# Changelog

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
