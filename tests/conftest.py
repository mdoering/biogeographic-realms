import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest


class _Handler(BaseHTTPRequestHandler):
    def _answer(self):
        self.server.requests.append((self.command, self.path, self.rfile.read(int(self.headers.get("Content-Length") or 0))))
        status, body = self.server.routes.get(self.path.split("?")[0], (404, {"message": "not found"}))
        payload = json.dumps(body).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    do_GET = _answer
    do_POST = _answer

    def log_message(self, *args):
        pass


@pytest.fixture
def api():
    """A local server answering canned JSON per path: `api.routes["/path"] = (status, body)`."""
    server = ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
    server.routes, server.requests = {}, []
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    server.base = f"http://127.0.0.1:{server.server_address[1]}"
    yield server
    server.shutdown()
