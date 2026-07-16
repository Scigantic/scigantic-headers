"""Recover acquisition optics for a data file from its session.

A raw MRC movie header has no pixel size (CELLA = 0). The value lives in a RELION
STAR or CryoSPARC .cs file the workflow writes in the same session. This finds
that file next to (or above) the data file and returns its optics, so a decoded
header can be completed in one call:

    hdr = decode_file(path)                 # header fields; pixelSizeA is None
    optics = read_session_optics(path)      # {'pixelSizeA': 1.05, 'voltageKv': 300, ...}
"""

from __future__ import annotations

import os
from typing import Dict

from .cryosparc import read_cryosparc_optics
from .star import read_star_optics


def read_session_optics(path: str, *, depth: int = 3) -> Dict[str, object]:
    """Find the acquisition optics for `path` from its session directory.

    Searches `path`'s directory and up to `depth` parent directories for a RELION
    `.star` or CryoSPARC `.cs`, and returns the first optics that carry a pixel
    size, with a `source` field naming the file. Returns {} if none is found.
    `path` may be a file or a directory. Pure stdlib; the reads are bounded (the
    optics block is at the top of a STAR, the first record of a .cs)."""
    d = path if os.path.isdir(path) else os.path.dirname(os.path.abspath(path))
    for _ in range(max(1, depth)):
        try:
            names = sorted(os.listdir(d))
        except OSError:
            names = []
        for name in names:
            if name.endswith(".star"):
                optics, kind = read_star_optics(os.path.join(d, name)), "relion-star"
            elif name.endswith(".cs"):
                optics, kind = read_cryosparc_optics(os.path.join(d, name)), "cryosparc-cs"
            else:
                continue
            if optics.get("pixelSizeA"):
                return {**optics, "source": f"{kind}:{name}"}
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return {}
