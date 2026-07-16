"""scigantic-headers, decode scientific file headers into typed context.

Bytes in, typed fields out. Pure, zero-dependency core (decoders) plus bounded
parallel I/O readers (sources). One implementation shared by the cloud probe's
twin, the notebook, and the edge box, so a file is described the same way
everywhere.

    from scigantic_headers import decode_file, decode_paths
    hdr = decode_file("session/frame.mrc")
    print(hdr.summary)  # "MRC stack 4096x4096x16, float32"
"""

from .cryosparc import parse_cryosparc_optics, read_cryosparc_optics
from .optics import read_session_optics
from .decoders import (
    HEADER_BYTES,
    DecodedHeader,
    decode_bytes,
    decode_cryosparc_header,
    decode_mrc_header,
    decode_nifti_header,
    decode_npy_header,
    extension_of,
    has_decoder_for,
    register_decoder,
)
from .sources import (
    DEFAULT_WORKERS,
    decode_file,
    decode_paths,
    decode_url,
    decode_urls,
    is_decodable,
    iter_decodable_files,
    read_leading_bytes,
    read_leading_bytes_url,
)
from .star import parse_relion_optics, read_star_optics

__version__ = "0.1.0"
__all__ = [
    "HEADER_BYTES", "DecodedHeader", "decode_bytes", "decode_mrc_header",
    "decode_nifti_header", "decode_npy_header", "decode_cryosparc_header",
    "parse_cryosparc_optics", "read_cryosparc_optics", "read_session_optics",
    "extension_of", "has_decoder_for", "register_decoder",
    "DEFAULT_WORKERS", "decode_file", "decode_paths", "decode_url",
    "decode_urls", "is_decodable", "iter_decodable_files", "read_leading_bytes",
    "read_leading_bytes_url", "parse_relion_optics", "read_star_optics",
    "__version__",
]
