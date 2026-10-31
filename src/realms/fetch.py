"""Downloads that can't pass as complete unless they are."""

import http.client
import urllib.error
import urllib.request
from pathlib import Path

CHUNK = 1 << 20


class IncompleteDownload(Exception):
    pass


def _total_size(response, offset: int) -> int | None:
    content_range = response.headers.get("Content-Range")
    if content_range and "/" in content_range and not content_range.endswith("/*"):
        return int(content_range.rsplit("/", 1)[1])
    length = response.headers.get("Content-Length")
    return int(length) + offset if length is not None else None


def fetch(url: str, dest: Path, expected_size: int | None = None) -> Path:
    """Download `url` to `dest`, resuming from `<dest>.part`.

    `dest` only ever appears once the file has the expected size, or the size the server announced if none was given. A
    file already at `dest` is trusted unless it contradicts `expected_size`.
    """
    dest = Path(dest)
    if dest.exists():
        if expected_size is not None and dest.stat().st_size != expected_size:
            raise IncompleteDownload(f"{dest} has {dest.stat().st_size} bytes, expected {expected_size}")
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    offset = part.stat().st_size if part.exists() else 0

    headers = {"User-Agent": "biogeographic-realms"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    total = expected_size
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=120) as response:
            if response.status != 206:
                offset = 0
            total = expected_size if expected_size is not None else _total_size(response, offset)
            with part.open("ab" if offset else "wb") as out:
                while chunk := response.read(CHUNK):
                    out.write(chunk)
    except (http.client.HTTPException, urllib.error.URLError, TimeoutError, ConnectionError) as e:
        size = part.stat().st_size if part.exists() else 0
        raise IncompleteDownload(f"{url} stopped at {size} bytes: {e}") from e

    size = part.stat().st_size
    if total is not None and size != total:
        raise IncompleteDownload(f"{url} stopped at {size} of {total} bytes, kept {part}")
    part.replace(dest)
    return dest
