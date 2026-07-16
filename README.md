# scigantic-headers

Read the metadata of a scientific file and return its fields (dimensions, data
type, pixel size, columns, row count) as a dict, without reading the rest of the
file. Decodes MRC, NPY, NIfTI, CryoSPARC `.cs`, and Parquet. The decode
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

One implementation, shared. The cloud probe's TypeScript twin
(`backend/src/services/headerDecoders.ts`), the notebook, and the edge
context-producer describe a file the same way. Same byte offsets, same
`mode -> dtype` table, held to the same test fixtures so they cannot drift.

Edge use: the context-producer (`infrastructure/edge/context-producer`) runs
this next to an instrument, with no network, on amd64 or arm64. The zero
dependencies and header-only reads are what let it run there.

## Install

    pip install scigantic-headers

Published to PyPI on a `v*` tag by `.github/workflows/publish.yml`. Pin the
version in anything that consumes it (`scigantic-headers==0.1.0`); the notebook
and edge-producer images in the monorepo install it that way.

It has zero runtime dependencies, so it runs air-gapped and installs into a slim
image without pulling anything else. Where there is no package index at build
time, install from a wheel or a checkout of this repo instead:

    pip install ./scigantic-headers

## Why it exists

A schema card or context record wants to say what a file is without opening it
over a FUSE mount. A scientific file's first ~1 KiB is a structured header that
already says so. This reads that, and nothing more.

Decoders today:

- **MRC / MRCS**: cryo-EM micrographs, movie stacks, tilt series, EMDB maps.
- **NPY**: NumPy arrays (ML, genomics, materials, not cryo-EM). Cross-checked
  against numpy.
- **NIfTI-1**: neuroimaging volumes, `.nii` and `.nii.gz` (not cryo-EM).
  Cross-checked against nibabel. Big- and little-endian.
- **CryoSPARC `.cs`**: a structured-array dataset. Reports record count and the
  field schema.
- **Parquet**: reads the footer, not a leading header, and returns the column
  names, physical types, and row count. Cross-checked against pyarrow.

`.gz` is transparent. The reader decompresses just the leading block, so a
gzipped `.nii.gz` or `.mrc.gz` decodes without inflating the whole file, and
dispatch sees the inner format.

Adding a format is a pure `bytes -> DecodedHeader | None` function plus one
`register_decoder` call. Nothing about the dispatch, the readers, or the batch
path is cryo-EM-specific. NPY and NIfTI landed with zero plumbing changes.

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
      sources.py    read leading bytes from file or URL; bounded parallel batch
      star.py       RELION STAR optics reader
      cryosparc.py  CryoSPARC .cs optics reader
      optics.py     read_session_optics: find a data file's optics file
      cli.py        scigantic-headers <file|url|--dir>
      benchmark.py  scigantic-headers-bench, measures the speed levers
    tests/          pytest (117 tests, incl. fuzz + golden fixtures)

## Robustness

A decoder is handed the leading bytes of arbitrary files, so it must be total:
for any input it returns None or a valid, JSON-safe DecodedHeader. It never
raises and never emits NaN or inf. This is enforced, not assumed.

- **Fuzz** (`test_fuzz.py`): thousands of random and magic-seeded byte buffers
  through every decoder. Each result must serialize with `allow_nan=False`.
- **Golden fixtures** (`fixtures/mrc-cases.json`): hex input to exact decode.
  This is the single source of truth the three MRC decoders (this library, the
  TypeScript twin, scigantic_empiar) are all held to, so none can silently
  diverge.
- Non-finite floats (a garbage NIfTI pixdim, an inf MRC cell) sanitize to null.

## Speed

The decode is microseconds: a 1 KiB read and a few `struct` unpacks. Optimizing
that is a rounding error. All the speed is in I/O, and two things carry it. Both
are measured, not asserted (run `scigantic-headers-bench`).

**Read the header, not the file.** The read is exactly `HEADER_BYTES` (1 KiB),
so its cost is independent of file size and a multi-GB movie decodes as fast as
a small one.

    header read + decode (1 KiB):        15 us
    full-file read + decode (200 MB):  22.7 ms      (about 1,500x cheaper)

**Parallelize the I/O, bounded.** Header reads across files are independent and
I/O-bound, so a thread pool overlaps their latency (the read releases the GIL).
Measured over 8 real EMPIAR headers pulled by HTTP Range:

    serial   (1 worker):   25.7 s
    parallel (8 workers):   0.7 s        (about 35x)

That 35x is large because EBI's per-request latency is about 3 s, so overlapping
eight requests recovers a lot of dead wait. On low-latency local storage the
speedup is smaller. The win scales with per-read latency, which is why the pool
is bounded and tunable (default 8, the EMPIAR sweet spot; past a point more
connections hit server throttling or disk-queue thrash and get slower).

Two more, by construction:

- **Filter before I/O.** `has_decoder_for` is a cheap extension check. The batch
  and walk paths use it to skip a file before any syscall.
- **Zero deps, fast import.** The dtype table is a plain dict, not numpy. There
  is nothing to install and nothing to import, and it runs air-gapped.

The biggest system-level lever lives in the caller, not here: decode once when a
file lands, write the context record, and answer every later query from the tiny
record without re-decoding. The edge context-producer and the cloud probe both
do this.

## Use

    from scigantic_headers import decode_file, decode_paths, decode_urls, iter_decodable_files

    decode_file("s/frame.mrc")                              # one local file
    decode_paths(iter_decodable_files("/mnt/sessions"))     # a tree, parallel
    decode_urls([...], workers=8)                           # remote, by Range

    # CLI
    scigantic-headers session/frame.mrc
    scigantic-headers https://ftp.ebi.ac.uk/empiar/.../img.mrcs
    scigantic-headers --dir /mnt/flashblade/sessions --workers 8

    pytest
    scigantic-headers-bench --empiar   # reproduce the numbers above

## Limitations

- **Pointer-based formats are not supported.** HDF5 and TIFF store their layout
  behind internal pointers, so reading the header or footer alone is not enough.
  Parquet (a footer) is supported; formats that need to chase offsets through the
  file are not. Plain gzip is handled (the reader inflates the leading block),
  but BAM's BGZF framing and per-record structure are not.
- **Output is not a standard.** The fields are an ad-hoc dict, not an
  interchange format such as Allotrope ASM. Map the dict to a standard if you
  need one.

## Repository

This is the library's home and the version published to PyPI. The Scigantic
monorepo installs it from PyPI; it is not vendored there. The paths above
(`backend/...`, `infrastructure/...`) point at the TypeScript twin, the
notebook, and the edge context-producer in that repo that use it.
