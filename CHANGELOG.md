# Changelog

## 0.1.0

First release.

- Header decoders (pure, zero-dependency): MRC/MRCS (cryo-EM), NPY (numpy
  arrays), NIfTI-1 (neuroimaging, `.nii` and `.nii.gz`), CryoSPARC `.cs`.
- Acquisition-optics readers: RELION STAR and CryoSPARC `.cs` (pixel size,
  voltage, spherical aberration). These recover the pixel size a raw MRC movie
  header omits. `read_session_optics(path)` finds a data file's optics file in
  its session directory and returns the pixel size in one call.
- Bounded, `.gz`-transparent readers; bounded-parallel batch over files or HTTP
  Range.
- Verified against reference implementations (numpy, nibabel) and real EMPIAR
  data; MRC held to a golden fixture shared with the backend TypeScript twin.
