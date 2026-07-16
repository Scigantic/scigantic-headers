"""Read leading bytes from files or URLs and decode them, serially or in a
bounded parallel batch.

The decode (decoders.py) is microseconds; the cost here is I/O, so this is where
speed is won or lost. Two principles are baked in:

  * Bounded reads. We read exactly HEADER_BYTES, never the file. For a multi-GB
    movie that is a ~million-fold saving over reading the whole thing, and the
    read cost is independent of file size.

  * Bounded parallelism. Header reads across files are independent, so a thread
    pool overlaps their I/O latency. Threads (not processes) are right because
    the work is I/O-bound, the GIL is released during the read. The win is real
    for high-latency storage (NFS/FlashBlade, cold disk, remote HTTP) and
    marginal for warm local files; the pool is bounded and tunable because, as
    the EMPIAR range-reader found, past a point more connections hurt (server
    throttling, disk queue thrash). Default 8, the EMPIAR sweet spot.

URL reads use an HTTP Range request, so a header is pulled from a remote archive
without downloading the file, the same lazy access the FUSE mounts use.
"""

from __future__ import annotations

import os
import urllib.request
from concurrent.futures import ThreadPoolExecutor
import gzip
import zlib
from typing import Dict, Iterable, List, Optional, Tuple

from .decoders import HEADER_BYTES, DecodedHeader, decode_bytes, has_decoder_for

DEFAULT_WORKERS = 8
_UA = "scigantic-headers/0.1 (+https://scigantic.com; mailto:support@scigantic.com)"

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


def decode_file(path: str) -> Optional[DecodedHeader]:
    """Decode one local file's header, or None if no decoder applies / it fails.
    Sees through a '.gz' wrapper; skips the read entirely when nothing decodes."""
    if not is_decodable(path):
        return None
    try:
        data = read_leading_bytes(path)
    except (OSError, EOFError, gzip.BadGzipFile):
        return None
    return decode_bytes(_inner_key(path), data)


def decode_url(url: str) -> Optional[DecodedHeader]:
    """Decode a remote file's header via a Range request, or None. Sees through
    a '.gz' wrapper."""
    if not is_decodable(url):
        return None
    try:
        data = read_leading_bytes_url(url)
    except Exception:
        return None
    return decode_bytes(_inner_key(url), data)


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
