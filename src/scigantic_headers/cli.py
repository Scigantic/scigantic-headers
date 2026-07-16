"""scigantic-headers, decode a file or URL header into typed fields.

    scigantic-headers session/frame.mrc
    scigantic-headers https://ftp.ebi.ac.uk/empiar/.../img.mrcs
    scigantic-headers --dir /mnt/flashblade/sessions --workers 8

A URL is read with an HTTP Range request, so a header is pulled from a remote
archive without downloading the file.
"""

from __future__ import annotations

import argparse
import json
import sys

from .decoders import has_decoder_for
from .sources import decode_file, decode_paths, decode_url, iter_decodable_files


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
    if not has_decoder_for(args.source):
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
