"""scigantic-headers, decode a file or URL header into typed fields.

    scigantic-headers session/frame.mrc
    scigantic-headers https://example.org/path/img.mrcs
    scigantic-headers --dir /data/sessions --workers 8
    scigantic-headers image.zarr              # a Zarr or N5 store (directory)

A URL is read with HTTP Range requests, so a header is pulled from a remote
archive without downloading the file; pointer-based formats (TIFF, HDF5, CZI,
ND2) issue a handful of them.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

from .sources import decode_file, decode_paths, decode_url, is_decodable, iter_decodable_files


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="scigantic-headers", description=__doc__)
    ap.add_argument("source", nargs="?", help="a file path or http(s) URL")
    ap.add_argument("--dir", help="decode every recognized file under this directory")
    ap.add_argument("--workers", type=int, default=8, help="parallel workers for --dir")
    args = ap.parse_args(argv)

    if args.dir:
        results = decode_paths(iter_decodable_files(args.dir), workers=args.workers)
        print(json.dumps(
            {p: h.to_dict() for p, h in results.items()}, indent=2, default=str))
        print(f"# {len(results)} decoded", file=sys.stderr)
        return 0

    if not args.source:
        ap.error("give a file/URL, or --dir")
    is_store = os.path.isdir(args.source)  # a Zarr or N5 store decodes through its sidecar file
    if not is_store and not is_decodable(args.source):
        print(f"no decoder for '{args.source}' (extension not recognized)", file=sys.stderr)
        return 1

    hdr = decode_url(args.source) if args.source.startswith(("http://", "https://")) else decode_file(args.source)
    if hdr is None:
        print(f"header did not decode for {args.source} (invalid or unreadable)", file=sys.stderr)
        return 1
    print(json.dumps({"source": args.source, **hdr.to_dict()}, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
