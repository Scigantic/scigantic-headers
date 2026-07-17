"""scigantic-headers, decode scientific file headers into typed context.

Bytes in, typed fields out. Pure, zero-dependency core (decoders) plus bounded
parallel I/O readers (sources). One implementation, so a file is described the
same way everywhere it runs.

    from scigantic_headers import decode_file, decode_paths
    hdr = decode_file("session/frame.mrc")
    print(hdr.summary)  # "MRC stack 4096x4096x16, float32"
"""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

from .bed import decode_bed
from .cryosparc import parse_cryosparc_optics, read_cryosparc_optics
from .dicom import decode_dicom
from .fastq import decode_fastq
from .fcs import decode_fcs
from .genbank import decode_genbank
from .gff import decode_gff
from .illumina import decode_illumina_run
from .mmcif import decode_mmcif
from .mzml import decode_mzml
from .optics import read_session_optics
from .parquet import decode_parquet
from .pdb import decode_pdb
from .sam import decode_sam
from .vcf import decode_vcf
from .decoders import (
    HEADER_BYTES,
    DecodedHeader,
    Read,
    decode_bytes,
    decode_cryosparc_header,
    decode_mrc_header,
    decode_nifti_header,
    decode_npy_header,
    extension_of,
    has_decoder_for,
    read_for,
    register_decoder,
    register_decoder_for_name,
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

try:
    __version__ = _pkg_version("scigantic-headers")
except PackageNotFoundError:  # a source checkout with no install metadata
    __version__ = "0.0.0+unknown"
__all__ = [
    "HEADER_BYTES", "DecodedHeader", "Read", "decode_bytes", "decode_mrc_header",
    "decode_nifti_header", "decode_npy_header", "decode_cryosparc_header",
    "parse_cryosparc_optics", "read_cryosparc_optics", "read_session_optics",
    "decode_parquet", "decode_fcs", "decode_mzml",
    "decode_fastq", "decode_illumina_run", "decode_vcf", "decode_sam", "decode_pdb",
    "decode_mmcif", "decode_genbank", "decode_gff", "decode_bed", "decode_dicom",
    "extension_of", "has_decoder_for", "read_for", "register_decoder", "register_decoder_for_name",
    "DEFAULT_WORKERS", "decode_file", "decode_paths", "decode_url",
    "decode_urls", "is_decodable", "iter_decodable_files", "read_leading_bytes",
    "read_leading_bytes_url", "parse_relion_optics", "read_star_optics",
    "__version__",
]
