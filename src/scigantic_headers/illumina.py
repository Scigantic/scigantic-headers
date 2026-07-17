"""Decode Illumina sequencing run metadata (RunInfo.xml / RunParameters.xml).

These sit at the root of a sequencer run folder and describe the run without the
reads: instrument, flowcell, run id, date, and read structure (cycle counts and
which reads are indices). RunInfo.xml is well standardized; RunParameters.xml
varies by platform, so its fields are extracted best-effort with alternates.

Dispatched by file name, since the .xml extension alone is too generic (see
decoders.register_decoder_for_name).
"""

from __future__ import annotations

import re
from typing import Optional

from .decoders import DecodedHeader, Read, register_decoder_for_name

_READ_EL = re.compile(rb"<Read\b[^>]*/?>")
_CYCLES = re.compile(rb'\bNumCycles="(\d+)"')
_INDEXED = re.compile(rb'\bIsIndexedRead="([YyNn])"')


def _first(data: bytes, patterns) -> Optional[str]:
    for p in patterns:
        m = re.search(p, data)
        if m:
            return m.group(1).decode("utf-8", "replace").strip() or None
    return None


def decode_illumina_run(data: bytes) -> Optional[DecodedHeader]:
    is_runinfo = b"<RunInfo" in data
    is_runparams = b"<RunParameters" in data
    if not (is_runinfo or is_runparams):
        return None

    run_id = _first(data, [
        rb'<Run\b[^>]*\bId="([^"]*)"',
        rb"<RunId>([^<]*)</RunId>",
        rb"<RunID>([^<]*)</RunID>",
    ])
    instrument = _first(data, [
        rb"<Instrument>([^<]*)</Instrument>",
        rb"<InstrumentName>([^<]*)</InstrumentName>",
        rb"<ScannerID>([^<]*)</ScannerID>",
        rb"<InstrumentID>([^<]*)</InstrumentID>",
    ])
    flowcell = _first(data, [
        rb"<Flowcell>([^<]*)</Flowcell>",
        rb"<FlowcellSerialBarcode>([^<]*)</FlowcellSerialBarcode>",
        rb"<Barcode>([^<]*)</Barcode>",
    ])
    date = _first(data, [rb"<Date>([^<]*)</Date>"])
    lanes = _first(data, [rb'<FlowcellLayout\b[^>]*\bLaneCount="(\d+)"'])

    reads = []
    for el in _READ_EL.findall(data):
        c = _CYCLES.search(el)
        if c:
            ix = _INDEXED.search(el)
            reads.append({
                "numCycles": int(c.group(1)),
                "isIndex": bool(ix and ix.group(1).upper() == b"Y"),
            })

    fields = {
        "runId": run_id,
        "instrument": instrument,
        "flowcell": flowcell,
        "date": date,
        "reads": reads,
        "laneCount": int(lanes) if lanes else None,
    }
    kind = "RunInfo" if is_runinfo else "RunParameters"
    return DecodedHeader(
        format="illumina-run",
        summary="Illumina run %s (%s), instrument %s, %d reads"
        % (run_id or "?", kind, instrument or "?", len(reads)),
        fields=fields,
    )


register_decoder_for_name(
    ["RunInfo.xml", "RunParameters.xml"], decode_illumina_run, read=Read(leading=256 * 1024)
)
