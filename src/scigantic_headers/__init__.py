"""scigantic-headers, decode scientific file headers into typed context.

Bytes in, typed fields out. Pure, zero-dependency core (decoders) plus bounded
parallel I/O readers (sources). One implementation, so a file is described the
same way everywhere it runs.

    from scigantic_headers import decode_file, decode_paths
    hdr = decode_file("session/frame.mrc")
    print(hdr.summary)  # "MRC stack 4096x4096x16, float32"
"""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

from .acquisition import Finding, compare_acquisition
from .bed import decode_bed
from .czi import decode_czi
from .hdf5 import decode_hdf5
from .ics import decode_ics
from .klb import decode_klb
from .lif import decode_lif
from .mdoc import decode_mdoc
from .nd2 import decode_nd2
from .obf import decode_obf
from .ome_xml import parse_ome_xml
from .tiff import decode_tiff
from .xml_sidecar import decode_bdv_xml, decode_ome_xml, decode_xml_sidecar
from .zarr import decode_n5_attributes, decode_zarr_json, decode_zarray, decode_zattrs, decode_zgroup
from .cryosparc import parse_cryosparc_optics, read_cryosparc_optics
from .dicom import decode_dicom
from .epu import decode_epu_xml, epu_pixel_size_candidates
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
    SEEK_MAX_READS,
    BytesSource,
    DecodedHeader,
    Read,
    Source,
    decode_bytes,
    decode_source,
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
    FileSource,
    UrlSource,
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
    "Finding", "compare_acquisition", "decode_epu_xml", "epu_pixel_size_candidates", "decode_mdoc",
    "HEADER_BYTES", "SEEK_MAX_READS", "DecodedHeader", "Read", "Source", "BytesSource",
    "FileSource", "UrlSource", "decode_bytes", "decode_source", "decode_mrc_header",
    "decode_nifti_header", "decode_npy_header", "decode_cryosparc_header",
    "parse_cryosparc_optics", "read_cryosparc_optics", "read_session_optics",
    "decode_parquet", "decode_fcs", "decode_mzml",
    "decode_fastq", "decode_illumina_run", "decode_vcf", "decode_sam", "decode_pdb",
    "decode_mmcif", "decode_genbank", "decode_gff", "decode_bed", "decode_dicom",
    "decode_czi", "decode_hdf5", "decode_ics", "decode_klb", "decode_lif", "decode_nd2", "decode_obf",
    "decode_tiff", "parse_ome_xml",
    "decode_bdv_xml", "decode_ome_xml", "decode_xml_sidecar",
    "decode_n5_attributes", "decode_zarr_json", "decode_zarray", "decode_zattrs", "decode_zgroup",
    "extension_of", "has_decoder_for", "read_for", "register_decoder", "register_decoder_for_name",
    "DEFAULT_WORKERS", "decode_file", "decode_paths", "decode_url",
    "decode_urls", "is_decodable", "iter_decodable_files", "read_leading_bytes",
    "read_leading_bytes_url", "parse_relion_optics", "read_star_optics",
    "__version__",
]
