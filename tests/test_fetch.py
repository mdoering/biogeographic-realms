import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from realms.fetch import IncompleteDownload, fetch

BODY = bytes(range(256)) * 40


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        self.server.requests.append(self.headers.get("Range"))
        start = 0
        if self.headers.get("Range") and not self.server.ignore_range:
            start = int(self.headers["Range"].removeprefix("bytes=").removesuffix("-"))
        rest = BODY[start:]
        self.send_response(206 if start else 200)
        self.send_header("Content-Length", str(len(rest)))
        if start:
            self.send_header("Content-Range", f"bytes {start}-{len(BODY) - 1}/{len(BODY)}")
        self.end_headers()
        limit = self.server.send_limit
        self.wfile.write(rest if limit is None else rest[:limit])

    def log_message(self, *args):
        pass


@pytest.fixture
def server():
    srv = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    srv.requests, srv.send_limit, srv.ignore_range = [], None, False
    thread = threading.Thread(target=srv.serve_forever, daemon=True)
    thread.start()
    srv.url = f"http://127.0.0.1:{srv.server_address[1]}/file"
    yield srv
    srv.shutdown()


def test_complete_download_lands_at_dest(server, tmp_path):
    dest = tmp_path / "sub" / "file.bin"
    assert fetch(server.url, dest, expected_size=len(BODY)) == dest
    assert dest.read_bytes() == BODY
    assert not dest.with_name("file.bin.part").exists()


def test_short_download_raises_and_leaves_no_dest(server, tmp_path):
    server.send_limit = 1000
    dest = tmp_path / "file.bin"
    with pytest.raises(IncompleteDownload):
        fetch(server.url, dest, expected_size=len(BODY))
    assert not dest.exists()
    assert dest.with_name("file.bin.part").stat().st_size == 1000


def test_short_download_without_expected_size_is_caught_by_content_length(server, tmp_path):
    server.send_limit = 1000
    dest = tmp_path / "file.bin"
    with pytest.raises(IncompleteDownload):
        fetch(server.url, dest)
    assert not dest.exists()


def test_resumes_from_part_file(server, tmp_path):
    dest = tmp_path / "file.bin"
    dest.with_name("file.bin.part").write_bytes(BODY[:1000])
    fetch(server.url, dest, expected_size=len(BODY))
    assert server.requests == ["bytes=1000-"]
    assert dest.read_bytes() == BODY


def test_restarts_when_server_ignores_the_range(server, tmp_path):
    server.ignore_range = True
    dest = tmp_path / "file.bin"
    dest.with_name("file.bin.part").write_bytes(b"x" * 1000)
    fetch(server.url, dest, expected_size=len(BODY))
    assert dest.read_bytes() == BODY


def test_existing_complete_file_is_not_refetched(server, tmp_path):
    dest = tmp_path / "file.bin"
    dest.write_bytes(BODY)
    fetch(server.url, dest, expected_size=len(BODY))
    assert server.requests == []


def test_existing_file_of_wrong_size_is_an_error(server, tmp_path):
    dest = tmp_path / "file.bin"
    dest.write_bytes(b"short")
    with pytest.raises(IncompleteDownload):
        fetch(server.url, dest, expected_size=len(BODY))
    assert server.requests == []
