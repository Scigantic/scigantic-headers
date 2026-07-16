"""Decode a Parquet file's schema from its footer.

Parquet keeps its schema and row count in a footer, not a leading header:

    PAR1 ... row groups ... [FileMetaData][footer_length: 4 bytes LE][PAR1]

So this reads the trailing bytes (see sources.read_trailing_bytes), takes the
footer length from the last 8 bytes, and parses the FileMetaData that precedes
them. That metadata is Thrift-compact-encoded, so this module includes a small
reader for the parts we need: version, num_rows, and the schema (column names
and physical types). Pure stdlib, no pyarrow.

`decode_parquet` takes the *trailing* bytes of the file (unlike the other
decoders, which take the leading bytes), because that is where a Parquet schema
lives.
"""

from __future__ import annotations

from typing import Optional

from .decoders import DecodedHeader, register_decoder

# Parquet physical types (parquet.thrift Type enum).
_PARQUET_TYPE = {
    0: "BOOLEAN", 1: "INT32", 2: "INT64", 3: "INT96",
    4: "FLOAT", 5: "DOUBLE", 6: "BYTE_ARRAY", 7: "FIXED_LEN_BYTE_ARRAY",
}

_MAGIC = b"PAR1"


class _Compact:
    """Just enough of the Thrift compact protocol to walk FileMetaData."""

    def __init__(self, buf: bytes):
        self.b = buf
        self.i = 0

    def _u8(self) -> int:
        v = self.b[self.i]
        self.i += 1
        return v

    def _varint(self) -> int:
        r = shift = 0
        while True:
            x = self._u8()
            r |= (x & 0x7F) << shift
            if not (x & 0x80):
                return r
            shift += 7

    def _zigzag(self) -> int:
        n = self._varint()
        return (n >> 1) ^ -(n & 1)

    def _binary(self) -> bytes:
        n = self._varint()
        s = self.b[self.i:self.i + n]
        self.i += n
        return s

    def skip(self, t: int) -> None:
        """Consume a value of compact type `t` without keeping it."""
        if t in (1, 2):        # bool true/false: value is the type, no bytes
            return
        if t == 3:             # i8
            self.i += 1
        elif t in (4, 5, 6):   # i16 / i32 / i64
            self._varint()
        elif t == 7:           # double
            self.i += 8
        elif t == 8:           # binary / string
            self._binary()
        elif t in (9, 10):     # list / set
            h = self._u8()
            size = h >> 4
            et = h & 0x0F
            if size == 15:
                size = self._varint()
            for _ in range(size):
                self.skip(et)
        elif t == 11:          # map
            size = self._varint()
            if size:
                kv = self._u8()
                kt, vt = kv >> 4, kv & 0x0F
                for _ in range(size):
                    self.skip(kt)
                    self.skip(vt)
        elif t == 12:          # struct
            self.walk(lambda fid, ft: False)

    def walk(self, handler) -> None:
        """Read struct fields until STOP. `handler(field_id, type)` reads the
        value and returns True, or returns False to have the field skipped."""
        fid = 0
        while True:
            h = self._u8()
            if h == 0:                 # STOP
                return
            t = h & 0x0F
            delta = h >> 4
            fid = self._zigzag() if delta == 0 else fid + delta
            if not handler(fid, t):
                self.skip(t)


def _read_schema_element(c: _Compact) -> dict:
    el: dict = {}

    def handler(fid: int, t: int) -> bool:
        if fid == 1 and t == 5:        # type (i32)
            el["type"] = c._zigzag()
            return True
        if fid == 4 and t == 8:        # name (binary)
            el["name"] = c._binary().decode("utf-8", "replace")
            return True
        if fid == 5 and t == 5:        # num_children (i32)
            el["num_children"] = c._zigzag()
            return True
        return False

    c.walk(handler)
    return el


def _parse_footer(footer: bytes) -> dict:
    c = _Compact(footer)
    out: dict = {"version": None, "num_rows": None, "columns": []}

    def handler(fid: int, t: int) -> bool:
        if fid == 1 and t == 5:        # version (i32)
            out["version"] = c._zigzag()
            return True
        if fid == 3 and t == 6:        # num_rows (i64)
            out["num_rows"] = c._zigzag()
            return True
        if fid == 2 and t == 9:        # schema: list<SchemaElement>
            h = c._u8()
            size = h >> 4
            if size == 15:
                size = c._varint()
            for _ in range(size):
                el = _read_schema_element(c)
                if el.get("type") is not None:   # a leaf column has a physical type
                    out["columns"].append({
                        "name": el.get("name", ""),
                        "type": _PARQUET_TYPE.get(el["type"], f"type{el['type']}"),
                    })
            return True
        return False

    c.walk(handler)
    return out


def decode_parquet(tail: bytes) -> Optional[DecodedHeader]:
    """Decode a Parquet schema from the file's *trailing* bytes. `tail` must end
    with the file's last bytes (the footer, its 4-byte length, and PAR1).
    Returns None if the magic is wrong or the footer is not fully present in
    `tail` (read more trailing bytes and retry)."""
    if len(tail) < 8 or tail[-4:] != _MAGIC:
        return None
    footer_len = int.from_bytes(tail[-8:-4], "little")
    if footer_len <= 0 or footer_len + 8 > len(tail):
        return None
    footer = tail[-8 - footer_len:-8]
    try:
        meta = _parse_footer(footer)
    except Exception:
        return None
    cols = meta["columns"]
    nrows = meta["num_rows"]
    rows = nrows if nrows is not None else "?"
    return DecodedHeader(
        format="parquet",
        summary=f"Parquet, {rows} rows, {len(cols)} columns",
        fields={"numRows": nrows, "numColumns": len(cols), "columns": cols},
    )


register_decoder("parquet", decode_parquet)
