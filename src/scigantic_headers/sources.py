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

from .decoders import (
    HEADER_BYTES, DecodedHeader, Source, decode_bytes, decode_source, has_decoder_for,
    read_for, sidecar_names,
)

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


class FileSource(Source):
    """A `Source` over a local file: seek + bounded read per request. The file
    stays open for the decode and is closed by `close()` (or the `with` block
    in `decode_file`)."""

    def __init__(self, path: str, budget: int):
        super().__init__(budget)
        self._fh = open(path, "rb")
        self.size = os.fstat(self._fh.fileno()).st_size

    def _fetch(self, offset: int, n: int) -> bytes:
        self._fh.seek(offset)
        return self._fh.read(n)

    def close(self) -> None:
        self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


# A remote read is rounded up to this granularity and cached, so a decoder that
# walks a chain of small structures (HDF5 object headers, TIFF IFDs) near each
# other costs one round trip, not one per structure. 64 KiB is well under the
# per-request latency cost that dominates remote reads.
URL_CHUNK_BYTES = 64 * 1024

# Range requests one decode may send. Over HTTP each is a round trip, so this is
# the knob that keeps a remote decode to a handful of requests (a typical HDF5
# or TIFF decode needs 3 to 12; a pathological file is cut off, not chased).
URL_MAX_REQUESTS = 48


class UrlSource(Source):
    """A `Source` over an http(s) URL: one HTTP Range request per uncached read,
    rounded out to `URL_CHUNK_BYTES`. The file size comes from the first
    response's Content-Range (needed by footer-walking formats like ND2). A
    server that ignores Range (200 with the whole body) is detected and the
    decode gives up rather than downloading the file: a short read is returned."""

    def __init__(self, url: str, budget: int, timeout: float = 30.0):
        super().__init__(budget)
        self._url = url
        self._timeout = timeout
        self._cache: Dict[int, bytes] = {}  # chunk start -> bytes
        self._requests = 0
        self.size = self._probe_size()

    def _request(self, range_header: str) -> Tuple[bytes, Optional[int], int]:
        if self._requests >= URL_MAX_REQUESTS:
            return b"", None, 0
        self._requests += 1
        req = urllib.request.Request(self._url, headers={"Range": range_header, "User-Agent": _UA})
        with urllib.request.urlopen(req, timeout=self._timeout) as resp:
            status = resp.status
            total = None
            cr = resp.headers.get("Content-Range")
            if cr and "/" in cr:
                tail = cr.rsplit("/", 1)[1]
                total = int(tail) if tail.isdigit() else None
            if status != 206:
                return b"", total, status  # Range ignored: do not pull the body
            return resp.read(), total, status

    def _probe_size(self) -> Optional[int]:
        """Learn the file size from a 1-byte suffix Range; counts toward the
        read budget like any other request. Returns None if unlearnable."""
        try:
            _, total, _status = self._request("bytes=-1")
            return total
        except Exception:
            return None

    def _fetch(self, offset: int, n: int) -> bytes:
        end = offset + n
        if self.size is not None:
            end = min(end, self.size)
        out = bytearray()
        chunk = (offset // URL_CHUNK_BYTES) * URL_CHUNK_BYTES
        while chunk < end:
            buf = self._cache.get(chunk)
            if buf is None:
                last = chunk + URL_CHUNK_BYTES - 1
                if self.size is not None:
                    last = min(last, self.size - 1)
                buf, total, status = self._request("bytes=%d-%d" % (chunk, last))
                if status != 206:
                    break
                if self.size is None and total is not None:
                    self.size = total
                self._cache[chunk] = buf
            lo = max(offset, chunk) - chunk
            hi = min(end, chunk + URL_CHUNK_BYTES) - chunk
            piece = buf[lo:hi]
            out += piece
            if len(piece) < hi - lo:
                break  # EOF inside this chunk
            chunk += URL_CHUNK_BYTES
        return bytes(out)


def _decode_directory(path: str) -> Optional[DecodedHeader]:
    """A directory decodes through the first registered sidecar file it holds
    (zarr.json, .zattrs, attributes.json, ...), which is how a Zarr or N5 store
    gets a header. None if it holds none."""
    try:
        present = {n.lower() for n in os.listdir(path)}
    except OSError:
        return None
    for name in sidecar_names():
        if name in present:
            # Find the real-cased entry so the open works on case-sensitive filesystems.
            for entry in os.listdir(path):
                if entry.lower() == name:
                    hdr = decode_file(os.path.join(path, entry))
                    if hdr is not None:
                        return _complete_store(path, hdr)
                    break
    return None


def _complete_store(path: str, hdr: DecodedHeader) -> DecodedHeader:
    """A multiscale group names its level-0 array but not its shape; when that
    array is a child directory of the store, decode it too and merge."""
    level0 = hdr.fields.get("level0Path")
    if not isinstance(level0, str) or not level0.strip("/"):
        return hdr
    child = os.path.join(path, *level0.strip("/").split("/"))
    if not os.path.isdir(child):
        return hdr
    from .zarr import with_level0_array  # the one store format with this two-level shape
    return with_level0_array(hdr, _decode_directory(child))


def decode_file(path: str) -> Optional[DecodedHeader]:
    """Decode one local file's header, or None if no decoder applies / it fails.
    The decoder's registered `Read` says how to fetch: a footer format (Parquet)
    is read from the tail, a seeking format (TIFF, HDF5, CZI, ND2) follows its
    pointers through a bounded `FileSource`, everything else is read from the
    start. Sees through a '.gz' wrapper for leading-header formats. A directory
    (a Zarr or N5 store) decodes through the sidecar file it contains."""
    if os.path.isdir(path):
        return _decode_directory(path)
    if not is_decodable(path):
        return None
    inner = _inner_key(path)
    spec = read_for(inner)
    try:
        if spec.seek is not None:
            with FileSource(path, spec.seek) as src:
                return decode_source(inner, src)
        if spec.footer is not None:
            data = read_trailing_bytes(path, spec.footer)
        else:
            data = read_leading_bytes(path, spec.leading)
    except (OSError, EOFError, gzip.BadGzipFile):
        return None
    return decode_bytes(inner, data)


def decode_url(url: str) -> Optional[DecodedHeader]:
    """Decode a remote file's header via Range requests, or None. The decoder's
    registered `Read` says how: a footer format uses a suffix Range, a seeking
    format issues a bounded number of Range requests where its pointers lead,
    everything else a single leading Range."""
    if not is_decodable(url):
        return None
    inner = _inner_key(url)
    spec = read_for(inner)
    try:
        if spec.seek is not None:
            return decode_source(inner, UrlSource(url, spec.seek))
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
