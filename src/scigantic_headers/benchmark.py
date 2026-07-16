"""Measure where the time actually goes, honestly, not by assertion.

Two things worth measuring, matching the two real speed levers:

  1. Bounded read vs full read. Decoding reads HEADER_BYTES; it does not read
     the file. This shows the header read is ~constant and tiny while a full
     read scales with file size, the single biggest win, and the reason a
     multi-GB movie costs the same to decode as a 1 KB one.

  2. Serial vs bounded-parallel over many *remote* files. Header reads are
     independent and I/O-bound, so a thread pool overlaps their latency. Run
     against real EMPIAR URLs so the numbers are real network latency, not a
     simulation, this reproduces the EMPIAR range-reader's finding at the
     file-set level.

    scigantic-headers-bench                 # local read-size bench (offline)
    scigantic-headers-bench --empiar        # + remote serial-vs-parallel (network)
"""

from __future__ import annotations

import argparse
import os
import struct
import tempfile
import time

from .decoders import decode_bytes
from .sources import DEFAULT_WORKERS, decode_file, decode_urls, read_leading_bytes

# A handful of real EMPIAR-10002 movie stacks (public, no auth). Header pulled
# via Range, nothing is downloaded.
_EMPIAR_BASE = "https://ftp.ebi.ac.uk/empiar/world_availability/10002/data"
_EMPIAR_FILES = [
    "15_movie_gc.mrcs", "16_movie_gc.mrcs", "17_movie_gc.mrcs", "18_movie_gc.mrcs",
    "20_movie_gc.mrcs", "21_movie_gc.mrcs", "22_movie_gc.mrcs", "23_movie_gc.mrcs",
]


def _make_mrc_file(path: str, size_mb: int) -> None:
    """Write a valid-header MRC file of `size_mb` megabytes so we can time a full
    read against a bounded read at a realistic size."""
    buf = bytearray(1024)
    struct.pack_into("<4i", buf, 0, 4096, 4096, 16, 2)
    struct.pack_into("<i", buf, 28, 4096)
    buf[208:212] = b"MAP "
    with open(path, "wb") as fh:
        fh.write(buf)
        fh.write(b"\0" * (size_mb * 1024 * 1024))


def _time(fn, repeat=5):
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def bench_read_size(size_mb: int = 200) -> None:
    print(f"\n[1] bounded read vs full read ({size_mb} MB file)")
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "big.mrc")
        _make_mrc_file(path, size_mb)

        def bounded():
            decode_bytes(path, read_leading_bytes(path))

        def full():
            with open(path, "rb") as fh:
                decode_bytes(path, fh.read())

        tb = _time(bounded)
        tf = _time(full)
        print(f"    header read + decode (1 KiB): {tb * 1e6:8.1f} us")
        print(f"    full file  read + decode ({size_mb} MB): {tf * 1e3:8.1f} ms")
        print(f"    bounded read is ~{tf / tb:,.0f}x cheaper, and flat in file size")


def bench_parallel_remote(workers: int = DEFAULT_WORKERS) -> None:
    urls = [f"{_EMPIAR_BASE}/{name}" for name in _EMPIAR_FILES]
    print(f"\n[2] serial vs {workers}-way parallel over {len(urls)} real EMPIAR headers (network)")
    t_serial = _time(lambda: decode_urls(urls, workers=1), repeat=1)
    t_par = _time(lambda: decode_urls(urls, workers=workers), repeat=1)
    got = len(decode_urls(urls, workers=workers))
    print(f"    serial (1 worker):        {t_serial:6.2f} s")
    print(f"    parallel ({workers} workers):     {t_par:6.2f} s")
    if t_par > 0:
        print(f"    speedup: {t_serial / t_par:.1f}x   ({got}/{len(urls)} decoded)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="scigantic-headers-bench", description=__doc__)
    ap.add_argument("--size-mb", type=int, default=200, help="file size for the read-size bench")
    ap.add_argument("--empiar", action="store_true", help="also run the remote parallel bench (needs network)")
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = ap.parse_args(argv)

    bench_read_size(args.size_mb)
    if args.empiar:
        bench_parallel_remote(args.workers)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
