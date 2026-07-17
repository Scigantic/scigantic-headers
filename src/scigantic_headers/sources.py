"""Read leading bytes from files or URLs and decode them, serially or in a
bounded parallel batch.

The decode (decoders.py) is microseconds; the cost here is I/O, so this is where
speed is won or lost. Two principles are baked in:

  * Bounded reads. We read a header (a few KiB, up to whatever a format's decoder
    declared it needs), never the file. For a multi-GB movie that is a millionfold
    saving over reading the whole thing, and the read cost is independent of file
    size. This module is format-agnostic: it asks the registry how many bytes a
    key needs and from which end, and never names a format itself.

  * Bounded parallelism. Header reads across files are independent, so a thread
    pool overlaps their I/O latency. Threads (not processes) are right because
    the work is I/O-bound, the GIL is released during the read. The win is real
    for high-latency storage (NFS, cold disk, remote HTTP) and marginal for warm
    local files; the pool is bounded and tunable because past a point more
    connections hurt (server throttling, disk queue thrash). Default 8.

URL reads use an HTTP Range request, so a header is pulled from a remote archive
without downloading the file.
"""

from __future__ import annotations

import os
import urllib.request
from concurrent.futures import ThreadPoolExecutor
import gzip
import zlib
from importlib.metadata import PackageNotFoundError, version as _pkg_version
from typing import Dict, Iterable, List, Optional, Tuple

from .decoders import HEADER_BYTES, DecodedHeader, decode_bytes, has_decoder_for, read_for

DEFAULT_WORKERS = 8

# Identify the client (and version) to the archives we range-read, with a contact.
# The version comes from the installed package metadata, so it tracks the release
# and cannot go stale in the string.
try:
    _VERSION = _pkg_version("scigantic-headers")
except PackageNotFoundError:  # a source checkout with no install metadata
    _VERSION = "0.0.0+unknown"
_UA = "scigantic-headers/%s (+https://scigantic.com; mailto:support@scigantic.com)" % _VERSION

# Enough compressed bytes to be sure the first HEADER_BYTES of a gzip stream
# decompress from them. Headers compress well, so this is generous.
_GZ_FETCH_BYTES = 128 * 1024


def _inner_key(key: str) -> str:
    """Drop a trailing '.gz' so dispatch sees the real format: 'scan.nii.gz' ->
    'scan.nii'. Compression is a transport concern the decoders never see."""
    return key[:-3] if key.lower().endswith(".gz") else key


def is_decodable(key: str) -> bool:
    """True if a decoder handles this key, seeing through a '.gz' wrapper."""
    return has_decoder_for(_inner_key(key))


def read_leading_bytes(path: str, n: int = HEADER_BYTES) -> bytes:
    """First `n` bytes of a local file, transparently decompressing a '.gz'.
    gzip.open reads incrementally, so a multi-GB .nii.gz is not fully inflated,
    only the leading block needed for `n` bytes."""
    if path.lower().endswith(".gz"):
        with gzip.open(path, "rb") as fh:
            return fh.read(n)
    with open(path, "rb") as fh:
        return fh.read(n)


def _gunzip_leading(compressed: bytes, n: int) -> bytes:
    """Inflate up to `n` bytes from the head of a gzip stream, tolerating a
    truncated tail (we only fetched a leading chunk)."""
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)  # 16 = gzip framing
    try:
        return d.decompress(compressed, n)
    except zlib.error:
        return b""


def read_leading_bytes_url(url: str, n: int = HEADER_BYTES, timeout: float = 30.0) -> bytes:
    """First `n` bytes of a remote file via an HTTP Range request, no full
    download. For a '.gz' URL, fetch a compressed chunk and inflate the head."""
    gz = url.lower().endswith(".gz")
    want = _GZ_FETCH_BYTES if gz else n
    req = urllib.request.Request(url, headers={"Range": f"bytes=0-{want - 1}", "User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read(want)
    return _gunzip_leading(raw, n) if gz else raw


def read_trailing_bytes(path: str, n: int) -> bytes:
    """Last `n` bytes of a local file (the whole file if smaller)."""
    with open(path, "rb") as fh:
        fh.seek(0, 2)
        size = fh.tell()
        fh.seek(max(0, size - n))
        return fh.read()


def read_trailing_bytes_url(url: str, n: int, timeout: float = 30.0) -> bytes:
    """Last `n` bytes of a remote file via an HTTP suffix Range request."""
    req = urllib.request.Request(url, headers={"Range": f"bytes=-{n}", "User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def decode_file(path: str) -> Optional[DecodedHeader]:
    """Decode one local file's header, or None if no decoder applies / it fails.
    The decoder's registered `Read` says how many bytes to fetch and from which
    end: a footer format (Parquet) is read from the tail, everything else from the
    start. Sees through a '.gz' wrapper for leading-header formats."""
    if not is_decodable(path):
        return None
    inner = _inner_key(path)
    spec = read_for(inner)
    try:
        if spec.footer is not None:
            data = read_trailing_bytes(path, spec.footer)
        else:
            data = read_leading_bytes(path, spec.leading)
    except (OSError, EOFError, gzip.BadGzipFile):
        return None
    return decode_bytes(inner, data)


def decode_url(url: str) -> Optional[DecodedHeader]:
    """Decode a remote file's header via a Range request, or None. The decoder's
    registered `Read` says how many bytes and from which end: a footer format uses
    a suffix Range, everything else a leading Range."""
    if not is_decodable(url):
        return None
    inner = _inner_key(url)
    spec = read_for(inner)
    try:
        if spec.footer is not None:
            data = read_trailing_bytes_url(url, spec.footer)
        else:
            data = read_leading_bytes_url(url, spec.leading)
    except Exception:
        return None
    return decode_bytes(inner, data)


def _batch(sources: Iterable[str], one, workers: int) -> List[Tuple[str, Optional[DecodedHeader]]]:
    # Pre-filter on extension (no I/O) so the pool only spends slots on files
    # that could decode ('.gz'-aware).
    candidates = [s for s in sources if is_decodable(s)]
    if not candidates:
        return []
    # No point spinning up more workers than items.
    workers = max(1, min(workers, len(candidates)))
    if workers == 1:
        return [(s, one(s)) for s in candidates]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(one, candidates))
    return list(zip(candidates, results))


def decode_paths(
    paths: Iterable[str], *, workers: int = DEFAULT_WORKERS
) -> Dict[str, DecodedHeader]:
    """Decode many local files concurrently. Returns {path: DecodedHeader} for
    those that decoded (files with no decoder or an invalid header are omitted)."""
    return {s: d for s, d in _batch(paths, decode_file, workers) if d is not None}


def decode_urls(
    urls: Iterable[str], *, workers: int = DEFAULT_WORKERS
) -> Dict[str, DecodedHeader]:
    """Decode many remote files concurrently via Range requests. This is where
    parallelism pays most, network latency dominates and overlaps cleanly."""
    return {s: d for s, d in _batch(urls, decode_url, workers) if d is not None}


def iter_decodable_files(root: str) -> Iterable[str]:
    """Walk `root` yielding only files a decoder handles ('.gz'-aware)."""
    for dirpath, _dirs, files in os.walk(root):
        for name in files:
            if is_decodable(name):
                yield os.path.join(dirpath, name)
