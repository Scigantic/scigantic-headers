"""Decode an HDF5 container's structure, and the imaging formats built on it
(Imaris `.ims`, BigDataViewer `.h5`, Luxendo `.lux.h5`, generic `.h5`/`.hdf5`),
without libhdf5.

HDF5 is pointer based: a superblock names the root group's object header, whose
messages point at the B-tree and heap (or hold the links) that name its
children, whose object headers hold their dataspace, datatype, layout, filters
and attributes. None of this sits at a fixed offset past the superblock, so
this is a seeking decoder: it follows those pointers through a bounded `Source`
and stops when the budget, the read cap, or its own walk limits are reached.
What it reads is a few KB, wherever the file put it.

Structures handled (HDF5 file format specification v3):
  superblock v0/v1 (symbol-table root) and v2/v3 (root object header address)
  object header v1 (messages padded to 8) and v2 ("OHDR" + "OCHK" chunks)
  messages: dataspace (1), datatype (3), link (6), layout (8), filters (11),
            attribute (12), continuation (16), symbol table (17), link info (2)
  group storage: v1 B-tree ("TREE" type 0) + "SNOD" symbol nodes + "HEAP" names;
                 v2 compact link messages. Dense (fractal heap) link storage is
                 not walked: such a group reports its children as unknown.
  attribute values: fixed strings (incl. Imaris' one-byte-string arrays), small
                 fixed-point and floating-point scalars/arrays. Larger or exotic
                 values are reported by type only.

Attributes are read from the object header the links point to, never from a
byte scan: an HDF5 writer that relocates an object header leaves the stale copy
in place, so a scan can surface values the file no longer holds.

Flavors, detected from the structure:
  imaris    root attribute ImarisDataSet; voxel size from DataSetInfo/Image
            ExtMin/ExtMax over X/Y/Z; channels from DataSetInfo/Channel N; the
            pyramid shape from DataSet/ResolutionLevel i/TimePoint j/Channel k.
            An `.ims` that is a TIFF (Imaris "IMS3" files start with II*) is
            reported as imaris-ims3 with the TIFF decode when that is available.
  bdv-h5    root children sNN (setups) and tNNNNN (timepoints) plus
            __DATA_TYPES__; resolutions/subdivisions per setup; the level-0
            cells dataset of the first timepoint and setup.
  luxendo   by extension (.lux.h5): the generic structure plus the first
            dataset's shape, labelled as Luxendo.
  hdf5      anything else: the root attributes, groups and datasets found.
"""

from __future__ import annotations

import math
import struct
from typing import Dict, List, Optional, Tuple

from .decoders import DecodedHeader, Read, Source, register_decoder

HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"
_TIFF_MAGICS = (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")

# Walk limits. These bound the decode independent of the byte budget: a file
# with ten thousand timepoint groups is described after this many, with the
# truncation stated, rather than walked to the end.
MAX_OBJECTS = 160        # object headers parsed per decode
MAX_CHILDREN = 64        # links listed per group
MAX_ATTR_BYTES = 4096    # attribute value bytes decoded
MAX_DEPTH = 5

_FILTER_NAMES = {
    1: "deflate", 2: "shuffle", 3: "fletcher32", 4: "szip", 5: "nbit", 6: "scaleoffset",
    32000: "lzf", 32001: "blosc", 32004: "lz4", 32008: "bitshuffle", 32015: "zstd",
    32017: "jpeg", 32018: "jpeg-xl", 32026: "blosc2",
}


def _finite(x):
    return x if x is not None and math.isfinite(x) else None


class _Object:
    """A parsed object header: its messages by type, in order."""

    __slots__ = ("address", "messages", "version", "truncated")

    def __init__(self, address: int):
        self.address = address
        self.messages: List[Tuple[int, bytes]] = []
        self.version = 0
        self.truncated = False

    def first(self, mtype: int) -> Optional[bytes]:
        for t, d in self.messages:
            if t == mtype:
                return d
        return None

    def all(self, mtype: int) -> List[bytes]:
        return [d for t, d in self.messages if t == mtype]


class _H5:
    """Bounded reader over a Source. Every method is total: a short read or a
    malformed structure yields None/[] rather than an exception, and the caller
    reports what it got."""

    def __init__(self, src: Source):
        self.src = src
        self.O = 8  # size of offsets
        self.L = 8  # size of lengths
        self.base = 0
        self.root: Optional[int] = None
        self.root_symtab: Optional[Tuple[int, int]] = None  # (btree, heap) for a v0/v1 root
        self.objects_parsed = 0
        self._cache: Dict[int, _Object] = {}

    # ── primitives ────────────────────────────────────────────────────────
    def read(self, off: int, n: int) -> bytes:
        return self.src.read_at(self.base + off, n)

    def _u(self, b: bytes, i: int, n: int) -> Optional[int]:
        if i < 0 or i + n > len(b):
            return None
        return int.from_bytes(b[i:i + n], "little")

    def _addr(self, b: bytes, i: int) -> Optional[int]:
        v = self._u(b, i, self.O)
        # The "undefined address" is all ones.
        return None if v is None or v == (1 << (8 * self.O)) - 1 else v

    # ── superblock ────────────────────────────────────────────────────────
    def open(self) -> bool:
        head = self.read(0, 2048)
        if not head:
            return False
        sb_off = None
        # The superblock may follow a user block at 512, 1024, 2048, ... bytes.
        for cand in (0, 512, 1024):
            if head[cand:cand + 8] == HDF5_SIGNATURE:
                sb_off = cand
                break
        if sb_off is None:
            for cand in (2048, 4096, 8192, 16384):
                chunk = self.read(cand, 8)
                if chunk == HDF5_SIGNATURE:
                    sb_off = cand
                    head = self.read(cand, 2048)
                    break
            else:
                return False
        else:
            head = head[sb_off:]
        if len(head) < 96:
            return False  # a superblock with its root entry needs at least this much
        version = head[8]
        if version in (0, 1):
            self.O, self.L = head[13], head[14]
            if self.O not in (2, 4, 8) or self.L not in (2, 4, 8):
                return False
            p = 24 if version == 0 else 28
            base = self._addr(head, p)
            # Root group symbol table entry follows the four addresses.
            ste = p + 4 * self.O
            oh = self._addr(head, ste + self.O)
            cache_type = self._u(head, ste + 2 * self.O, 4)
            if oh is None:
                return False
            self.base = (base or 0) + sb_off if sb_off and (base or 0) == 0 else (base or 0)
            self.root = oh
            if cache_type == 1:
                bt = self._addr(head, ste + 2 * self.O + 8)
                hp = self._addr(head, ste + 2 * self.O + 8 + self.O)
                if bt is not None and hp is not None:
                    self.root_symtab = (bt, hp)
            return True
        if version in (2, 3):
            self.O, self.L = head[9], head[10]
            if self.O not in (2, 4, 8) or self.L not in (2, 4, 8):
                return False
            base = self._addr(head, 12)
            root = self._addr(head, 12 + 3 * self.O)
            if root is None:
                return False
            self.base = (base or 0) + sb_off if sb_off and (base or 0) == 0 else (base or 0)
            self.root = root
            return True
        return False

    # ── object headers ────────────────────────────────────────────────────
    def object(self, addr: int) -> Optional[_Object]:
        if addr in self._cache:
            return self._cache[addr]
        if self.objects_parsed >= MAX_OBJECTS:
            return None
        self.objects_parsed += 1
        head = self.read(addr, 512)
        if len(head) < 16:
            return None
        obj = _Object(addr)
        if head[:4] == b"OHDR":
            self._parse_v2(obj, head)
        elif head[0] == 1:
            self._parse_v1(obj, head)
        else:
            return None
        self._cache[addr] = obj
        return obj

    def _parse_v1(self, obj: _Object, head: bytes) -> None:
        obj.version = 1
        nmsgs = self._u(head, 2, 2) or 0
        size = self._u(head, 8, 4) or 0
        if size > 1 << 20:
            obj.truncated = True
            return
        body = head[16:16 + size]
        if len(body) < size:
            body = self.read(obj.address + 16, size)
        self._messages_v1(obj, body, nmsgs, depth=0)

    def _messages_v1(self, obj: _Object, body: bytes, nmsgs: int, depth: int) -> None:
        i = 0
        while i + 8 <= len(body) and len(obj.messages) < 512:
            mtype = self._u(body, i, 2)
            msize = self._u(body, i + 2, 2)
            if mtype is None or msize is None:
                break
            data = body[i + 8:i + 8 + msize]
            if len(data) < msize:
                obj.truncated = True
                break
            i += 8 + msize
            if mtype == 0x10:  # continuation: offset, length
                off = self._addr(data, 0)
                ln = self._u(data, self.O, self.L)
                if off is not None and ln and ln < 1 << 20 and depth < 8:
                    more = self.read(off, ln)
                    if len(more) < ln:
                        obj.truncated = True
                    self._messages_v1(obj, more, nmsgs, depth + 1)
                continue
            if mtype != 0:
                obj.messages.append((mtype, data))

    def _parse_v2(self, obj: _Object, head: bytes) -> None:
        obj.version = 2
        flags = head[5]
        i = 6
        if flags & 0x20:
            i += 16
        if flags & 0x10:
            i += 4
        size_len = 1 << (flags & 3)
        chunk_size = self._u(head, i, size_len)
        i += size_len
        if chunk_size is None or chunk_size > 1 << 20:
            obj.truncated = True
            return
        body = head[i:i + chunk_size]
        if len(body) < chunk_size:
            body = self.read(obj.address + i, chunk_size)
            if len(body) < chunk_size:
                obj.truncated = True
        self._messages_v2(obj, body, flags, depth=0)

    def _messages_v2(self, obj: _Object, body: bytes, flags: int, depth: int) -> None:
        i = 0
        # The chunk ends with a 4-byte checksum; stop before it.
        end = len(body) - 4
        while i + 4 <= end and len(obj.messages) < 512:
            mtype = body[i]
            msize = self._u(body, i + 1, 2) or 0
            i += 4
            if flags & 4:
                i += 2
            data = body[i:i + msize]
            if len(data) < msize:
                obj.truncated = True
                break
            i += msize
            if mtype == 0x10:
                off = self._addr(data, 0)
                ln = self._u(data, self.O, self.L)
                if off is not None and ln and ln < 1 << 20 and depth < 8:
                    more = self.read(off, ln)
                    if more[:4] == b"OCHK":
                        self._messages_v2(obj, more[4:], flags, depth + 1)
                    else:
                        obj.truncated = True
                continue
            if mtype != 0:
                obj.messages.append((mtype, data))

    # ── groups ────────────────────────────────────────────────────────────
    def children(self, obj: _Object, symtab: Optional[Tuple[int, int]] = None):
        """[(name, address)] of a group's links, capped at MAX_CHILDREN, plus
        a flag saying whether the listing is complete. Non-groups give ([], True)."""
        out: List[Tuple[str, int]] = []
        complete = True
        st = obj.first(0x11)
        if st is not None and len(st) >= 2 * self.O:
            symtab = (self._addr(st, 0), self._addr(st, self.O))
        if symtab and symtab[0] is not None and symtab[1] is not None:
            heap = self._heap(symtab[1])
            complete = self._btree_group(symtab[0], heap, out, 0)
            return out, complete
        links = obj.all(0x6)
        if links:
            for ln in links:
                parsed = self._link(ln)
                if parsed:
                    out.append(parsed)
                if len(out) >= MAX_CHILDREN:
                    complete = len(links) <= MAX_CHILDREN
                    break
            return out, complete
        li = obj.first(0x2)
        if li is not None and len(li) >= 2:
            # Link Info: dense storage, the links live as messages in a fractal
            # heap indexed by a v2 B-tree on name.
            i = 2 + (8 if li[1] & 1 else 0)
            heap_addr = self._addr(li, i)
            btree_addr = self._addr(li, i + self.O)
            if heap_addr is None:
                return out, True  # an empty group: no heap was ever allocated
            if btree_addr is None:
                return [], False
            heap = _FractalHeap(self, heap_addr)
            recs = self._btree_v2_records(btree_addr)
            if recs is None:
                return [], False
            for rec in recs:
                # Type 5 record: name hash u32, then the heap ID.
                body = heap.get(rec[4:])
                parsed = self._link(body) if body else None
                if parsed:
                    out.append(parsed)
                if len(out) >= MAX_CHILDREN:
                    out.sort()
                    return out, len(recs) <= MAX_CHILDREN
            out.sort()  # the B-tree orders by name hash; names are what a reader expects
            return out, True
        return out, True

    def _heap(self, addr: int) -> Optional[Tuple[int, int]]:
        h = self.read(addr, 8 + 2 * self.L + self.O)
        if h[:4] != b"HEAP":
            return None
        size = self._u(h, 8, self.L)
        data = self._addr(h, 8 + 2 * self.L)
        if size is None or data is None:
            return None
        return data, size

    def _heap_name(self, heap, off: int) -> str:
        if heap is None or off < 0 or off >= heap[1]:
            return ""
        raw = self.read(heap[0] + off, min(256, heap[1] - off))
        return raw.split(b"\x00", 1)[0].decode("utf-8", "replace")

    def _btree_group(self, addr: int, heap, out: List[Tuple[str, int]], depth: int) -> bool:
        if depth > 6 or len(out) >= MAX_CHILDREN:
            return False
        node = self.read(addr, 8 + 2 * self.O + 2)
        if node[:4] != b"TREE" or node[4] != 0:
            return False
        level = node[5]
        used = self._u(node, 6, 2) or 0
        entry = self.L + self.O
        body = self.read(addr + 8 + 2 * self.O, entry * used + self.L)
        complete = True
        for k in range(used):
            child = self._addr(body, self.L + k * entry)
            if child is None:
                return False
            if level == 0:
                complete = self._snod(child, heap, out) and complete
            else:
                complete = self._btree_group(child, heap, out, depth + 1) and complete
            if len(out) >= MAX_CHILDREN:
                return False
        return complete

    def _snod(self, addr: int, heap, out: List[Tuple[str, int]]) -> bool:
        h = self.read(addr, 8)
        if h[:4] != b"SNOD":
            return False
        n = self._u(h, 6, 2) or 0
        ent = 2 * self.O + 24
        body = self.read(addr + 8, ent * n)
        for k in range(n):
            if len(out) >= MAX_CHILDREN:
                return False
            name_off = self._u(body, k * ent, self.O)
            oh = self._addr(body, k * ent + self.O)
            if name_off is None or oh is None:
                return False
            out.append((self._heap_name(heap, name_off), oh))
        return True

    def _link(self, m: bytes) -> Optional[Tuple[str, int]]:
        if len(m) < 3 or m[0] != 1:
            return None
        flags = m[1]
        i = 2
        ltype = 0
        if flags & 8:
            ltype = m[i]
            i += 1
        if flags & 4:
            i += 8
        if flags & 0x10:
            i += 1
        nlen = self._u(m, i, 1 << (flags & 3))
        i += 1 << (flags & 3)
        if nlen is None or i + nlen > len(m):
            return None
        name = m[i:i + nlen].decode("utf-8", "replace")
        i += nlen
        if ltype != 0:
            return None  # soft/external links do not lead to an object here
        addr = self._addr(m, i)
        return (name, addr) if addr is not None else None

    # ── dataspace / datatype / layout / filters ───────────────────────────
    def dataspace(self, m: Optional[bytes]) -> Optional[List[int]]:
        if m is None or len(m) < 2:
            return None
        ver, rank = m[0], m[1]
        if ver == 1:
            i = 8
        elif ver == 2:
            if len(m) > 3 and m[3] == 0:  # scalar
                return []
            if len(m) > 3 and m[3] == 2:  # null
                return []
            i = 4
        else:
            return None
        if rank > 32:
            return None
        dims = []
        for _ in range(rank):
            d = self._u(m, i, self.L)
            if d is None:
                return None
            dims.append(d)
            i += self.L
        return dims

    def datatype(self, m: Optional[bytes]) -> Optional[Dict[str, object]]:
        """{'dtype': name, 'size': bytes, 'class': n}"""
        if m is None or len(m) < 8:
            return None
        cls = m[0] & 0x0F
        bits = m[1] | (m[2] << 8) | (m[3] << 16)
        size = self._u(m, 4, 4) or 0
        if cls == 0:
            name = ("int" if bits & 8 else "uint") + str(size * 8)
        elif cls == 1:
            name = {2: "float16", 4: "float32", 8: "float64", 16: "float128"}.get(size, "float%d" % (size * 8))
        elif cls == 2:
            name = "bitfield%d" % (size * 8)
        elif cls == 3:
            name = "string(%d)" % size
        elif cls == 4:
            name = "opaque"
        elif cls == 5:
            name = "reference"
        elif cls == 6:
            name = "compound(%d)" % (bits & 0xFFFF)
        elif cls == 7:
            name = "enum"
        elif cls == 8:
            name = "array"
        elif cls == 9:
            name = "vlen-string" if (bits & 0xF) == 1 else "vlen"
        elif cls == 10:
            name = "time"
        else:
            name = "class%d" % cls
        return {"dtype": name, "size": size, "class": cls, "bits": bits}

    def layout(self, m: Optional[bytes]) -> Dict[str, object]:
        """{'layout': 'compact|contiguous|chunked|virtual', 'chunks': [...], 'address': int, 'size': int}"""
        out: Dict[str, object] = {}
        if m is None or len(m) < 2:
            return out
        ver = m[0]
        names = {0: "compact", 1: "contiguous", 2: "chunked", 3: "virtual"}
        if ver in (1, 2):
            ndims = m[1]
            cls = m[2]
            out["layout"] = names.get(cls, "class%d" % cls)
            i = 8
            if cls != 0:
                out["address"] = self._addr(m, i)
                i += self.O
            dims = []
            for _ in range(ndims):
                v = self._u(m, i, 4)
                if v is None:
                    break
                dims.append(v)
                i += 4
            if cls == 2 and dims:
                out["chunks"] = dims[:-1]  # last entry is the element size
            return out
        if ver == 3:
            cls = m[1]
            out["layout"] = names.get(cls, "class%d" % cls)
            if cls == 0:
                out["size"] = self._u(m, 2, 2)
            elif cls == 1:
                out["address"] = self._addr(m, 2)
                out["size"] = self._u(m, 2 + self.O, self.L)
            elif cls == 2:
                ndims = m[2]
                out["address"] = self._addr(m, 3)
                i = 3 + self.O
                dims = []
                for _ in range(ndims):
                    v = self._u(m, i, 4)
                    if v is None:
                        break
                    dims.append(v)
                    i += 4
                out["chunks"] = dims[:-1]
            return out
        if ver == 4:
            cls = m[1]
            out["layout"] = names.get(cls, "class%d" % cls)
            if cls == 1:
                out["address"] = self._addr(m, 2)
                out["size"] = self._u(m, 2 + self.O, self.L)
            elif cls == 2:
                ndims = m[3]
                enc = m[4]
                i = 5
                dims = []
                for _ in range(ndims):
                    v = self._u(m, i, enc)
                    if v is None:
                        break
                    dims.append(v)
                    i += enc
                out["chunks"] = dims[:-1]
            return out
        return out

    def filters(self, m: Optional[bytes]) -> List[str]:
        if m is None or len(m) < 2:
            return []
        ver, n = m[0], m[1]
        i = 8 if ver == 1 else 2
        out = []
        for _ in range(min(n, 8)):
            fid = self._u(m, i, 2)
            if fid is None:
                break
            i += 2
            name_len = 0
            if ver == 1 or fid >= 256:
                name_len = self._u(m, i, 2) or 0
                i += 2
            i += 2  # flags
            ncd = self._u(m, i, 2) or 0
            i += 2
            name = ""
            if name_len:
                name = m[i:i + name_len].split(b"\x00", 1)[0].decode("utf-8", "replace")
                if ver == 1:
                    name_len = (name_len + 7) // 8 * 8
                i += name_len
            i += 4 * ncd
            if ver == 1 and ncd % 2:
                i += 4
            out.append(_FILTER_NAMES.get(fid) or name or "filter%d" % fid)
        return out

    # ── attributes ────────────────────────────────────────────────────────
    def attributes(self, obj: _Object) -> Dict[str, object]:
        out: Dict[str, object] = {}
        for m in obj.all(0x0C)[:128]:
            parsed = self._attribute(m)
            if parsed:
                out[parsed[0]] = parsed[1]
        ai = obj.first(0x15)
        if ai is not None and len(ai) >= 2:
            # Attribute Info: dense storage, each attribute message in a fractal
            # heap, indexed by a v2 B-tree on name (type 8 records: heap ID,
            # message flags u8, creation order u32, name hash u32).
            i = 2 + (2 if ai[1] & 1 else 0)
            heap_addr = self._addr(ai, i)
            btree_addr = self._addr(ai, i + self.O)
            if heap_addr is not None and btree_addr is not None:
                heap = _FractalHeap(self, heap_addr)
                dense: Dict[str, object] = {}
                for rec in (self._btree_v2_records(btree_addr) or [])[:128]:
                    body = heap.get(rec[:heap.id_len])
                    parsed = self._attribute(body) if body else None
                    if parsed:
                        dense[parsed[0]] = parsed[1]
                for k in sorted(dense):  # B-tree order is by name hash
                    out[k] = dense[k]
        return out

    # ── v2 B-tree (name indexes of dense link/attribute storage) ──────────
    def _btree_v2_records(self, addr: int) -> Optional[List[bytes]]:
        """All records of a v2 B-tree, depth 0 or 1 (dense groups and attribute
        sets of imaging files are small). None if the tree is unreadable."""
        h = self.read(addr, 16 + self.O + 2 + self.L + 4)
        if h[:4] != b"BTHD":
            return None
        rec_size = self._u(h, 10, 2) or 0
        depth = self._u(h, 12, 2) or 0
        node_size = self._u(h, 6, 4) or 0
        root = self._addr(h, 16)
        root_nrec = self._u(h, 16 + self.O, 2) or 0
        if root is None or rec_size == 0 or rec_size > 64 or node_size > 1 << 20:
            return None
        if depth == 0:
            return self._btree_v2_leaf(root, root_nrec, rec_size)
        if depth == 1:
            node = self.read(root, 6 + root_nrec * rec_size + (root_nrec + 1) * (self.O + 8))
            if node[:4] != b"BTIN":
                return None
            out: List[bytes] = list(self._records(node, 6, root_nrec, rec_size))
            max_leaf = max(1, (node_size - 10) // rec_size)
            nrec_size = (max_leaf.bit_length() + 7) // 8
            p = 6 + root_nrec * rec_size
            for _ in range(root_nrec + 1):
                child = self._addr(node, p)
                n = self._u(node, p + self.O, nrec_size)
                p += self.O + nrec_size
                if child is None or n is None:
                    return None
                leaf = self._btree_v2_leaf(child, n, rec_size)
                if leaf is None:
                    return None
                out.extend(leaf)
                if len(out) > 4096:
                    break
            return out
        return None  # deeper trees: thousands of entries, out of scope

    def _btree_v2_leaf(self, addr: int, nrec: int, rec_size: int) -> Optional[List[bytes]]:
        if nrec > 4096:
            return None
        node = self.read(addr, 6 + nrec * rec_size)
        if node[:4] != b"BTLF":
            return None
        return list(self._records(node, 6, nrec, rec_size))

    @staticmethod
    def _records(node: bytes, start: int, n: int, size: int):
        for k in range(n):
            rec = node[start + k * size:start + (k + 1) * size]
            if len(rec) < size:
                return
            yield rec

    def _attribute(self, m: bytes) -> Optional[Tuple[str, object]]:
        if len(m) < 8:
            return None
        ver = m[0]
        name_len = self._u(m, 2, 2) or 0
        dt_len = self._u(m, 4, 2) or 0
        ds_len = self._u(m, 6, 2) or 0
        i = 8
        if ver == 3:
            i += 1  # name character set
        if ver == 1:
            pad = lambda n: (n + 7) // 8 * 8  # noqa: E731
        else:
            pad = lambda n: n  # noqa: E731
        name = m[i:i + name_len].split(b"\x00", 1)[0].decode("utf-8", "replace")
        i += pad(name_len)
        dt = self.datatype(m[i:i + dt_len])
        i += pad(dt_len)
        dims = self.dataspace(m[i:i + ds_len])
        i += pad(ds_len)
        data = m[i:]
        return name, self._value(dt, dims, data)

    def _value(self, dt, dims, data: bytes):
        """Decode a small attribute value, else describe it by type."""
        if dt is None:
            return None
        count = 1
        for d in (dims or []):
            count *= d
        size = dt["size"]
        total = count * size
        if total > MAX_ATTR_BYTES or total > len(data):
            return "<%s x%d>" % (dt["dtype"], count)
        cls = dt["class"]
        little = not (dt["bits"] & 1)
        if cls == 3:
            # Fixed string, or an array of 1-byte strings (Imaris writes its text this way).
            raw = data[:total]
            return raw.split(b"\x00", 1)[0].decode("utf-8", "replace") if raw else ""
        if cls == 0:
            vals = []
            for k in range(count):
                vals.append(int.from_bytes(data[k * size:(k + 1) * size], "little" if little else "big",
                                           signed=bool(dt["bits"] & 8)))
            return vals[0] if dims == [] or (dims is None and count == 1) else vals
        if cls == 1 and size in (4, 8):
            fmt = ("<" if little else ">") + ("f" if size == 4 else "d")
            vals = [_finite(struct.unpack_from(fmt, data, k * size)[0]) for k in range(count)]
            return vals[0] if dims == [] or (dims is None and count == 1) else vals
        if cls == 9 and (dt["bits"] & 0xF) == 1 and count == 1 and size == 16:
            # A variable-length string is a global heap reference: u32 length,
            # collection address O, object index u32.
            text = self._global_heap_string(data[:16])
            return text if text is not None else "<vlen-string>"
        return "<%s>" % dt["dtype"]

    def _global_heap_string(self, ref: bytes) -> Optional[str]:
        addr = self._addr(ref, 4)
        idx = self._u(ref, 4 + self.O, 4)
        if addr is None or not idx:
            return None
        col = self.read(addr, 8 + self.L + 8192)
        if col[:4] != b"GCOL":
            return None
        i = 8 + self.L
        for _ in range(256):
            oid = self._u(col, i, 2)
            size = self._u(col, i + 8, self.L)
            if oid is None or size is None:
                return None
            if oid == 0:
                return None
            if oid == idx:
                raw = col[i + 8 + self.L:i + 8 + self.L + size]
                if len(raw) < size:
                    raw = self.read(addr + i + 8 + self.L, size)
                return raw.split(b"\x00", 1)[0].decode("utf-8", "replace")
            i += 8 + self.L + (size + 7) // 8 * 8
        return None

    def _single_chunk(self, btree: int, ndims: int) -> Optional[int]:
        """Address of the only chunk of a one-chunk dataset, via its v1 B-tree
        (type 1) leaf: key (chunk size u32, filter mask u32, offsets L x ndims+1),
        child address O."""
        key = 8 + self.L * (ndims + 1)
        node = self.read(btree, 8 + 2 * self.O + key + self.O)
        if node[:4] != b"TREE" or node[4] != 1 or node[5] != 0 or (self._u(node, 6, 2) or 0) < 1:
            return None
        return self._addr(node, 8 + 2 * self.O + key)

    # ── datasets ──────────────────────────────────────────────────────────
    def dataset_info(self, obj: _Object) -> Optional[Dict[str, object]]:
        ds = obj.first(0x1)
        dt = obj.first(0x3)
        if ds is None or dt is None:
            return None
        shape = self.dataspace(ds)
        typ = self.datatype(dt)
        lay = self.layout(obj.first(0x8))
        info: Dict[str, object] = {
            "shape": shape,
            "dtype": typ["dtype"] if typ else None,
            "layout": lay.get("layout"),
            "chunks": lay.get("chunks"),
            "filters": self.filters(obj.first(0xB)),
        }
        return info

    def read_small_dataset(self, obj: _Object, max_bytes: int = 4096):
        """The values of a small contiguous numeric dataset (BDV's resolutions),
        or None."""
        info = self.dataset_info(obj)
        lay = self.layout(obj.first(0x8))
        typ = self.datatype(obj.first(0x3))
        if not info or typ is None or lay.get("address") is None:
            return None
        shape = info["shape"] or []
        count = 1
        for d in shape:
            count *= d
        size = typ["size"]
        if count * size > max_bytes:
            return None
        if lay.get("layout") == "contiguous":
            addr = lay["address"]
        elif lay.get("layout") == "chunked" and list(lay.get("chunks") or []) == list(shape) and not info["filters"]:
            # One unfiltered chunk holding the whole array (how BDV writes its
            # resolutions table): the v1 chunk B-tree at the layout address has
            # a single leaf entry whose child is the chunk.
            addr = self._single_chunk(lay["address"], len(shape))
            if addr is None:
                return None
        else:
            return None
        raw = self.read(addr, count * size)
        if len(raw) < count * size:
            return None
        little = not (typ["bits"] & 1)
        if typ["class"] == 1 and size in (4, 8):
            fmt = ("<" if little else ">") + ("f" if size == 4 else "d")
            vals = [_finite(struct.unpack_from(fmt, raw, k * size)[0]) for k in range(count)]
        elif typ["class"] == 0:
            vals = [int.from_bytes(raw[k * size:(k + 1) * size], "little" if little else "big",
                                   signed=bool(typ["bits"] & 8)) for k in range(count)]
        else:
            return None
        if len(shape) == 2:
            return [vals[r * shape[1]:(r + 1) * shape[1]] for r in range(shape[0])]
        return vals


class _FractalHeap:
    """Managed objects of a fractal heap ("FRHP"), as used by dense link and
    attribute storage. Supports a root direct block and a root indirect block
    whose children are direct blocks, which covers the heaps imaging writers
    produce; anything larger returns None for its objects."""

    def __init__(self, h: _H5, addr: int):
        self.h = h
        self.ok = False
        self._dblocks: Optional[List[Tuple[int, int, int]]] = None  # (heap offset, size, address)
        O, L = h.O, h.L
        hdr = h.read(addr, 16 + 12 * L + 2 + 2 * L + 2 + 2 + O + 2 + 16)
        if hdr[:4] != b"FRHP":
            return
        self.id_len = h._u(hdr, 5, 2) or 0
        filter_len = h._u(hdr, 7, 2) or 0
        self.flags = hdr[9]
        # Skip the accounting fields: next huge id (L), huge B-tree (O), free
        # space (L), free-space manager (O), then eight more lengths.
        p = 14 + L + O + L + O + 8 * L
        self.width = h._u(hdr, p, 2) or 0
        p += 2
        self.start_size = h._u(hdr, p, L) or 0
        p += L
        self.max_direct = h._u(hdr, p, L) or 0
        p += L
        self.max_heap_bits = h._u(hdr, p, 2) or 0
        p += 2 + 2  # starting number of rows
        self.root = h._addr(hdr, p)
        p += O
        self.cur_rows = h._u(hdr, p, 2) or 0
        self.filtered = filter_len > 0
        self.off_size = (self.max_heap_bits + 7) // 8
        self.len_size = self.id_len - 1 - self.off_size
        if self.root is None or self.width == 0 or self.start_size == 0 or self.len_size <= 0:
            return
        self.ok = True

    def _direct_blocks(self) -> List[Tuple[int, int, int]]:
        if self._dblocks is not None:
            return self._dblocks
        blocks: List[Tuple[int, int, int]] = []
        if self.cur_rows == 0:
            blocks.append((0, self.start_size, self.root))
        else:
            h = self.h
            entry = h.O + (h.L + 4 if self.filtered else 0)
            n = self.width * self.cur_rows
            ib = h.read(self.root, 5 + h.O + self.off_size + n * entry)
            if ib[:4] == b"FHIB":
                p = 5 + h.O + self.off_size
                off = 0
                for row in range(self.cur_rows):
                    size = self.start_size if row == 0 else self.start_size << (row - 1)
                    for _col in range(self.width):
                        addr = h._addr(ib, p)
                        p += entry
                        if addr is not None and size <= self.max_direct:
                            blocks.append((off, size, addr))
                        off += size
                    if len(blocks) > 256:
                        break
        self._dblocks = blocks
        return blocks

    def get(self, heap_id: bytes) -> Optional[bytes]:
        if not self.ok or len(heap_id) < 1:
            return None
        kind = (heap_id[0] >> 4) & 3
        if kind == 2:  # tiny: the object is in the ID itself
            n = (heap_id[0] & 0x0F) + 1
            return heap_id[1:1 + n]
        if kind != 0 or len(heap_id) < 1 + self.off_size + self.len_size:
            return None
        off = int.from_bytes(heap_id[1:1 + self.off_size], "little")
        length = int.from_bytes(heap_id[1 + self.off_size:1 + self.off_size + self.len_size], "little")
        if length == 0 or length > MAX_ATTR_BYTES * 4:
            return None
        # Heap offsets count from the start of a direct block's header, so an
        # object's file address is simply the block address plus the offset
        # into that block (verified on Imaris files: the first object sits at
        # heap offset 22, right after the 22-byte block header).
        for boff, bsize, baddr in self._direct_blocks():
            if boff <= off < boff + bsize:
                data = self.h.read(baddr + (off - boff), length)
                return data if len(data) == length else None
        return None


# ── structure walk ──────────────────────────────────────────────────────────

def _is_dataset(obj: _Object) -> bool:
    return obj.first(0x1) is not None and obj.first(0x3) is not None and obj.first(0x11) is None


def _walk(h: _H5, obj: _Object, path: str, depth: int, groups: List[str], datasets: List[Dict[str, object]],
          symtab=None) -> bool:
    """Collect group paths and dataset descriptions under `obj`. Returns False
    when the walk was cut short by a limit."""
    complete = True
    kids, kids_complete = h.children(obj, symtab)
    complete = complete and kids_complete
    for name, addr in kids:
        child = h.object(addr)
        if child is None:
            complete = False
            continue
        cpath = path + "/" + name if path else name
        if _is_dataset(child):
            if len(datasets) < 32:
                info = h.dataset_info(child) or {}
                info["path"] = cpath
                datasets.append(info)
            else:
                complete = False
        else:
            if len(groups) < 128:
                groups.append(cpath)
            else:
                complete = False
            if depth + 1 < MAX_DEPTH:
                complete = _walk(h, child, cpath, depth + 1, groups, datasets) and complete
            else:
                complete = False
    return complete


def _lookup(h: _H5, root: _Object, path: str, symtab=None) -> Optional[_Object]:
    """Follow a slash path from `root`, or None."""
    obj = root
    for part in path.split("/"):
        if not part:
            continue
        kids, _complete = h.children(obj, symtab if obj is root else None)
        addr = next((a for n, a in kids if n == part), None)
        if addr is None:
            return None
        obj = h.object(addr)
        if obj is None:
            return None
    return obj


def _num(text) -> Optional[float]:
    try:
        v = float(str(text).strip().split()[0])
    except (TypeError, ValueError, IndexError):
        return None
    return v if math.isfinite(v) else None


def _int(text) -> Optional[int]:
    v = _num(text)
    return int(v) if v is not None else None


# ── flavors ─────────────────────────────────────────────────────────────────

def _imaris(h: _H5, root: _Object, root_attrs: Dict[str, object], generic: Dict[str, object]) -> DecodedHeader:
    f: Dict[str, object] = {"imarisVersion": root_attrs.get("ImarisVersion")}
    info = _lookup(h, root, root_attrs.get("DataSetInfoDirectoryName") or "DataSetInfo", h.root_symtab)
    image_attrs: Dict[str, object] = {}
    channels: List[Dict[str, object]] = []
    if info is not None:
        kids, _ = h.children(info)
        for name, addr in kids:
            o = h.object(addr)
            if o is None:
                continue
            if name == "Image":
                image_attrs = h.attributes(o)
            elif name.startswith("Channel "):
                a = h.attributes(o)
                channels.append({"index": _int(name[8:]), "name": a.get("Name"),
                                 "emission": a.get("LSMEmissionWavelength"),
                                 "excitation": a.get("LSMExcitationWavelength")})
            elif name == "Imaris":
                a = h.attributes(o)
                f["application"] = " ".join(str(x) for x in (a.get("ManufactorString"), a.get("Version")) if x) or None
            elif name == "ImarisDataSet":
                a = h.attributes(o)
                f["creator"] = a.get("Creator")
            elif name == "CustomData":
                a = h.attributes(o)
                f["customDate"] = a.get("DateAndTime")
    channels.sort(key=lambda c: c["index"] if c["index"] is not None else 1 << 30)
    sx, sy, sz = (_int(image_attrs.get(k)) for k in ("X", "Y", "Z"))
    voxel: List[Optional[float]] = []
    for k, n in ((0, sx), (1, sy), (2, sz)):
        lo, hi = _num(image_attrs.get("ExtMin%d" % k)), _num(image_attrs.get("ExtMax%d" % k))
        voxel.append(_finite(round((hi - lo) / n, 6)) if lo is not None and hi is not None and n else None)
    unit = image_attrs.get("Unit")
    if unit and str(unit).strip().lower() in ("mm",):
        voxel = [v * 1000.0 if v is not None else None for v in voxel]
    elif unit and str(unit).strip().lower() in ("nm",):
        voxel = [v / 1000.0 if v is not None else None for v in voxel]

    # Pyramid shape and the level-0 data block.
    levels = timepoints = nchan = None
    data = None
    dataset = _lookup(h, root, root_attrs.get("DataSetDirectoryName") or "DataSet", h.root_symtab)
    if dataset is not None:
        rl, complete = h.children(dataset)
        levels = len(rl) if complete else None
        first = h.object(rl[0][1]) if rl else None
        if first is not None:
            tp, complete = h.children(first)
            timepoints = len(tp) if complete else None
            t0 = h.object(tp[0][1]) if tp else None
            if t0 is not None:
                ch, complete = h.children(t0)
                nchan = len(ch) if complete else None
                c0 = h.object(ch[0][1]) if ch else None
                if c0 is not None:
                    kids, _ = h.children(c0)
                    d = next((a for n, a in kids if n == "Data"), None)
                    dobj = h.object(d) if d is not None else None
                    if dobj is not None:
                        data = h.dataset_info(dobj)
    f.update({
        "sizeX": sx, "sizeY": sy, "sizeZ": sz,
        "sizeC": nchan if nchan is not None else (len(channels) or None),
        "sizeT": timepoints if timepoints is not None else _int(image_attrs.get("TimePoints")),
        "dtype": data["dtype"] if data else None,
        "voxelSizeUm": voxel,
        "unit": (str(unit).strip() or None) if unit is not None else None,
        "resolutionLevels": levels,
        "chunks": data["chunks"] if data else None,
        "compression": (data["filters"] or None) if data else None,
        "channels": [c["name"] for c in channels[:32]],
        "channelCount": nchan if nchan is not None else (len(channels) or None),
        "emissionWavelengths": [c["emission"] for c in channels[:32]] if any(c["emission"] for c in channels) else None,
        "recordingDate": image_attrs.get("RecordingDate"),
        "microscopeMode": image_attrs.get("MicroscopeMode"),
        "objective": image_attrs.get("LensPower") or image_attrs.get("Lenspower"),
        "numericalAperture": _num(image_attrs.get("NumericalAperture")),
        "imageName": image_attrs.get("Name"),
        "description": (image_attrs.get("Description") if image_attrs.get("Description") not in ("(description not specified)", "", None) else None),
    })
    dims = "x".join(str(d) for d in (sx, sy, sz) if d) if sx and sy else "?"
    parts = ["Imaris IMS %s" % dims]
    if f["sizeC"] and f["sizeC"] > 1:
        parts.append("%d channels" % f["sizeC"])
    if f["sizeT"] and f["sizeT"] > 1:
        parts.append("%d timepoints" % f["sizeT"])
    if f["dtype"]:
        parts.append(str(f["dtype"]))
    if all(v is not None for v in voxel):
        parts.append("x".join("%g" % v for v in voxel) + " " + (f["unit"] or "um"))
    if levels:
        parts.append("%d resolution levels" % levels)
    return DecodedHeader(format="imaris", summary=", ".join(parts), fields=f)


def _bdv(h: _H5, root: _Object, kids: List[Tuple[str, int]], complete: bool) -> DecodedHeader:
    setups = sorted(n for n, _ in kids if len(n) == 3 and n[0] == "s" and n[1:].isdigit())
    tps = sorted(n for n, _ in kids if len(n) == 6 and n[0] == "t" and n[1:].isdigit())
    f: Dict[str, object] = {
        "setupCount": len(setups) if complete else None,
        "timepointCount": len(tps) if complete else None,
        "setupCountAtLeast": len(setups) if not complete else None,
        "timepointCountAtLeast": len(tps) if not complete else None,
        "firstTimepoint": tps[0] if tps else None,
    }
    resolutions = subdivisions = None
    if setups:
        s0 = _lookup(h, root, setups[0], h.root_symtab)
        if s0 is not None:
            for name, addr in h.children(s0)[0]:
                o = h.object(addr)
                if o is None:
                    continue
                if name == "resolutions":
                    resolutions = h.read_small_dataset(o)
                elif name == "subdivisions":
                    subdivisions = h.read_small_dataset(o)
    f["resolutions"] = resolutions
    f["subdivisions"] = subdivisions
    f["resolutionLevels"] = len(resolutions) if resolutions else None
    cells = None
    if tps and setups:
        c = _lookup(h, root, "%s/%s/0/cells" % (tps[0], setups[0]), h.root_symtab)
        if c is not None:
            cells = h.dataset_info(c)
    shape = (cells or {}).get("shape") or []
    # BDV stores cells as [z, y, x].
    sz, sy, sx = (list(shape) + [None, None, None])[:3] if len(shape) == 3 else (None, None, None)
    f.update({
        "sizeX": sx, "sizeY": sy, "sizeZ": sz,
        "sizeC": None, "sizeT": f["timepointCount"],
        "dtype": (cells or {}).get("dtype"),
        "chunks": (cells or {}).get("chunks"),
        "compression": ((cells or {}).get("filters") or None) if cells else None,
        "voxelSizeUm": [None, None, None],  # the XML beside the file carries it
    })
    dims = "x".join(str(d) for d in (sx, sy, sz) if d) if sx else "?"
    parts = ["BigDataViewer HDF5 %s" % dims]
    if f["setupCount"]:
        parts.append("%d setups" % f["setupCount"])
    elif f["setupCountAtLeast"]:
        parts.append(">=%d setups" % f["setupCountAtLeast"])
    if f["timepointCount"] and f["timepointCount"] > 1:
        parts.append("%d timepoints" % f["timepointCount"])
    elif f["timepointCountAtLeast"]:
        parts.append(">=%d timepoints" % f["timepointCountAtLeast"])
    if f["dtype"]:
        parts.append(str(f["dtype"]))
    if f["resolutionLevels"]:
        parts.append("%d resolution levels" % f["resolutionLevels"])
    return DecodedHeader(format="bdv-h5", summary=", ".join(parts), fields=f)


def _generic(kind: str, label: str, root_attrs, groups, datasets, complete, ext_fields=None) -> DecodedHeader:
    first = datasets[0] if datasets else None
    f: Dict[str, object] = {
        "rootAttributes": {k: v for k, v in list(root_attrs.items())[:32]},
        "groups": groups[:64],
        "groupCount": len(groups) if complete else None,
        "datasets": datasets[:16],
        "datasetCount": len(datasets) if complete else None,
        "walkComplete": complete,
    }
    if first:
        shape = first.get("shape") or []
        f["firstDataset"] = first.get("path")
        f["shape"] = shape
        f["dtype"] = first.get("dtype")
        f["chunks"] = first.get("chunks")
        f["compression"] = first.get("filters") or None
    if ext_fields:
        f.update(ext_fields)
    parts = [label]
    if first and first.get("shape"):
        parts.append("%s %s" % (first["path"], "x".join(str(d) for d in first["shape"])))
        if first.get("dtype"):
            parts.append(str(first["dtype"]))
    parts.append("%d groups" % len(groups) + ("" if complete else "+"))
    parts.append("%d datasets" % len(datasets) + ("" if complete else "+"))
    return DecodedHeader(format=kind, summary=", ".join(parts), fields=f)


def _luxendo(h: _H5, root: _Object, root_attrs, groups, datasets, complete) -> DecodedHeader:
    # Luxendo (Bruker) .lux.h5: one 3-D block per view in a `Data` dataset
    # (z, y, x), the acquisition metadata in a JSON sidecar next to the file.
    first = next((d for d in datasets if d.get("path", "").split("/")[-1] == "Data"), datasets[0] if datasets else None)
    ext: Dict[str, object] = {}
    if first and len(first.get("shape") or []) == 3:
        z, y, x = first["shape"]
        ext.update({"sizeX": x, "sizeY": y, "sizeZ": z})
    hdr = _generic("luxendo-h5", "Luxendo HDF5", root_attrs, groups, datasets, complete, ext)
    return hdr


def decode_hdf5(src: Source) -> Optional[DecodedHeader]:
    """Decode an HDF5 container (or an Imaris IMS3 TIFF) from a Source. Total:
    None for anything that is not one, or that is cut off before the root
    group can be read."""
    try:
        return _decode(src)
    except Exception:
        return None


def _decode(src: Source) -> Optional[DecodedHeader]:
    head = src.read_at(0, 8)
    if head[:4] in _TIFF_MAGICS:
        return _tiff_container(src)
    h = _H5(src)
    if not h.open() or h.root is None:
        return None
    root = h.object(h.root)
    if root is None:
        return None
    root_attrs = h.attributes(root)
    kids, complete = h.children(root, h.root_symtab)
    names = {n for n, _ in kids}

    if root_attrs.get("ImarisDataSet") == "ImarisDataSet" or "DataSetInfo" in names and "DataSet" in names:
        return _imaris(h, root, root_attrs, {})
    if "__DATA_TYPES__" in names and any(n.startswith("s") and n[1:].isdigit() for n in names):
        return _bdv(h, root, kids, complete)

    groups: List[str] = []
    datasets: List[Dict[str, object]] = []
    walk_complete = _walk(h, root, "", 0, groups, datasets, h.root_symtab) and complete
    if (src.key or "").lower().endswith(".lux.h5"):
        return _luxendo(h, root, root_attrs, groups, datasets, walk_complete)
    return _generic("hdf5", "HDF5", root_attrs, groups, datasets, walk_complete)


def _tiff_container(src: Source) -> Optional[DecodedHeader]:
    """An `.ims` that is a TIFF (Imaris IMS3). Decoded by the TIFF decoder when
    present, relabelled; otherwise reported minimally."""
    try:
        from .tiff import decode_tiff  # type: ignore
    except Exception:
        decode_tiff = None
    if decode_tiff is not None:
        hdr = decode_tiff(src)
        if hdr is not None:
            fields = dict(hdr.fields)
            fields["container"] = hdr.format
            return DecodedHeader(format="imaris-ims3", summary="Imaris IMS3 (TIFF): " + hdr.summary, fields=fields)
    return DecodedHeader(format="imaris-ims3", summary="Imaris IMS3 (TIFF container)", fields={"container": "tiff"})


# 2 MiB covers a root group with hundreds of links, the pyramid walk of an
# Imaris file, and a BDV file's setup datasets; the walk limits above stop a
# pathological file well before that.
HDF5_BUDGET = 2 * 1024 * 1024

register_decoder(["h5", "hdf5", "he5", "hdf", "ims"], decode_hdf5, read=Read(seek=HDF5_BUDGET))
