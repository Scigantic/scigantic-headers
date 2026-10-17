"""The Source abstraction behind seeking decoders: budgets, request caps, sparse
replay, file and URL sources."""

import threading
from http.server import BaseHTTPRequestHandler, HTTPServer

import pytest

from scigantic_headers.decoders import (
    SEEK_MAX_READS, BytesSource, DecodedHeader, Read, decode_bytes, decode_source,
    register_decoder,
)
from scigantic_headers.sources import FileSource, UrlSource, URL_CHUNK_BYTES, decode_file


def test_bytes_source_reads_and_eof():
    s = BytesSource(b"0123456789")
    assert s.size == 10
    assert s.read_at(2, 3) == b"234"
    assert s.read_at(8, 5) == b"89"       # short at EOF
    assert s.read_at(20, 5) == b""        # past EOF
    assert s.read_at(-1, 5) == b""        # negative offset
    assert s.read_tail(4) == b"6789"
    assert BytesSource(b"ab").read_tail(10) == b"ab"


def test_byte_budget_is_enforced():
    s = BytesSource(b"x" * 100, budget=10)
    assert len(s.read_at(0, 8)) == 8
    assert len(s.read_at(0, 8)) == 2      # only 2 left
    assert s.read_at(0, 8) == b""         # spent


def test_read_count_cap_is_enforced():
    s = BytesSource(b"x" * 100000)
    for _ in range(SEEK_MAX_READS):
        assert s.read_at(0, 1) == b"x"
    assert s.read_at(0, 1) == b""


def test_sparse_ranges_replay_like_a_file():
    s = BytesSource(ranges={0: b"HEAD", 1000: b"TAIL"})
    assert s.size == 1004
    assert s.read_at(0, 4) == b"HEAD"
    assert s.read_at(1000, 2) == b"TA"
    assert s.read_at(500, 4) == b""        # unrecorded gap reads as nothing
    assert s.read_at(2, 10) == b"AD"       # bounded to the recorded range
    assert s.read_tail(4) == b"TAIL"


def _seeking(src):
    head = src.read_at(0, 4)
    if head != b"SEEK":
        return None
    tail = src.read_tail(4)
    return DecodedHeader(format="seektest", summary="ok", fields={"tail": tail.decode(), "size": src.size})


register_decoder("seektest", _seeking, read=Read(seek=4096))


def test_decode_bytes_hands_seeking_decoder_a_source():
    d = decode_bytes("a.seektest", b"SEEK" + b"." * 100 + b"END!")
    assert d.fields == {"tail": "END!", "size": 108}
    assert decode_bytes("a.seektest", b"NOPE") is None


def test_decode_source_serves_leading_and_footer_decoders_too(mrc):
    hdr = mrc(64, 64, 1, 2, mx=64, cella_x=64.0)
    assert decode_source("f.mrc", BytesSource(hdr + b"\0" * 5000)).fields["nx"] == 64
    assert decode_source("f.nothing", BytesSource(hdr)) is None


def test_file_source_and_decode_file(tmp_path):
    p = tmp_path / "x.seektest"
    p.write_bytes(b"SEEK" + b"\0" * 10000 + b"last")
    with FileSource(str(p), budget=1 << 20) as s:
        assert s.size == 10008
        assert s.read_at(10004, 4) == b"last"
    d = decode_file(str(p))
    assert d.format == "seektest" and d.fields["tail"] == "last"


def test_seeking_decoder_that_raises_is_contained():
    def boom(src):
        raise ValueError("bad")
    register_decoder("boomtest", boom, read=Read(seek=16))
    assert decode_bytes("a.boomtest", b"anything") is None


# ── UrlSource against a local Range-capable server ────────────────────────────

_BODY = b"SEEK" + bytes(range(256)) * 600 + b"TAIL"   # 153,608 bytes, > 2 chunks


class _RangeHandler(BaseHTTPRequestHandler):
    honor_range = True
    requests = []

    def log_message(self, *a):
        pass

    def do_GET(self):
        rng = self.headers.get("Range")
        type(self).requests.append(rng)
        if not self.honor_range or not rng:
            self.send_response(200)
            self.send_header("Content-Length", str(len(_BODY)))
            self.end_headers()
            self.wfile.write(_BODY)
            return
        spec = rng.split("=", 1)[1]
        if spec.startswith("-"):
            n = int(spec[1:])
            start, end = max(0, len(_BODY) - n), len(_BODY) - 1
        else:
            a, b = spec.split("-")
            start = int(a)
            end = min(int(b), len(_BODY) - 1) if b else len(_BODY) - 1
        chunk = _BODY[start:end + 1]
        self.send_response(206)
        self.send_header("Content-Range", "bytes %d-%d/%d" % (start, end, len(_BODY)))
        self.send_header("Content-Length", str(len(chunk)))
        self.end_headers()
        self.wfile.write(chunk)


@pytest.fixture
def range_server():
    srv = HTTPServer(("127.0.0.1", 0), _RangeHandler)
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    _RangeHandler.requests.clear()
    _RangeHandler.honor_range = True
    yield "http://127.0.0.1:%d/file.seektest" % srv.server_port
    srv.shutdown()


def test_url_source_learns_size_and_coalesces_reads(range_server):
    s = UrlSource(range_server, budget=1 << 20)
    assert s.size == len(_BODY)
    assert s.read_at(0, 4) == b"SEEK"
    assert s.read_at(10, 6) == bytes(range(6, 12))     # same chunk: no new request
    assert s.read_tail(4) == b"TAIL"
    assert s.read_at(URL_CHUNK_BYTES - 2, 4) == _BODY[URL_CHUNK_BYTES - 2:URL_CHUNK_BYTES + 2]  # spans chunks
    # one size probe + chunk 0 + the tail chunk + chunk 1 for the spanning read
    assert len(_RangeHandler.requests) == 4


def test_url_source_refuses_a_server_that_ignores_range(range_server):
    _RangeHandler.honor_range = False
    s = UrlSource(range_server, budget=1 << 20)
    assert s.size is None
    assert s.read_at(0, 4) == b""   # does not download the body


def test_decode_url_seeking(range_server):
    from scigantic_headers import decode_url
    d = decode_url(range_server)
    assert d is not None and d.fields["tail"] == "TAIL"
