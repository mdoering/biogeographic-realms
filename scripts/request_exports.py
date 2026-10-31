"""Ask ChecklistBank for the exports that `fetch_reference.py` found missing: WRiMS and the GRIIS lists without one.

This writes to production, so it is for Markus to run with his own ChecklistBank (GBIF) account. Without `--submit` it only
lists what it would ask for. Each dataset is checked again first and skipped if a usable export now exists. Requests are
sequential. Once the jobs finish, run `fetch_reference.py` again to download them.

Usage: uv run python scripts/request_exports.py [--submit] [--user NAME]
"""

import argparse
import base64
import getpass
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

from realms.gbif import CLB_API, USER_AGENT
from realms.reference import latest_export

MANIFEST = Path("data/reference/manifest.json")
# extended carries distributions, environment and extinct, and fetch_reference.py takes it before a plain export
REQUEST = {"format": "COLDP", "extended": True, "synonyms": True}


def missing() -> list[tuple[str, int]]:
    manifest = json.loads(MANIFEST.read_text())
    return sorted((name, entry["dataset"]) for name, entry in manifest.items() if entry.get("status") == "missing")


def still_missing(name: str, dataset: int) -> bool:
    """A missing item that now has an export fetch_reference.py would use: extended for WRiMS, extended or plain for GRIIS."""
    if latest_export(dataset, extended=True):
        return False
    return name == "wrims" or latest_export(dataset) is None


def submit(dataset: int, auth: str) -> dict:
    request = urllib.request.Request(
        f"{CLB_API}/dataset/{dataset}/export",
        data=json.dumps(REQUEST).encode(),
        headers={**USER_AGENT, "Content-Type": "application/json", "Accept": "application/json", "Authorization": f"Basic {auth}"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read() or b"{}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--submit", action="store_true", help="really POST; without it nothing is written")
    parser.add_argument("--user", help="ChecklistBank (GBIF) username; asked for when not given")
    args = parser.parse_args()

    todo = [(name, key) for name, key in missing() if still_missing(name, key)]
    print(f"{len(todo)} datasets without a usable export: {', '.join(name for name, _ in todo)}")
    if not todo or not args.submit:
        if todo:
            print(f"dry run: would POST {json.dumps(REQUEST)} to {CLB_API}/dataset/<key>/export for each; add --submit")
        return

    user = args.user or input("ChecklistBank username: ")
    auth = base64.b64encode(f"{user}:{getpass.getpass('Password: ')}".encode()).decode()
    failed = 0
    for name, key in todo:
        try:
            job = submit(key, auth)
            print(f"ok       {name} (dataset {key}): job {job.get('key')} {job.get('status', '')}", flush=True)
        except urllib.error.HTTPError as e:
            failed += 1
            print(f"FAILED   {name} (dataset {key}): HTTP {e.code} {e.read()[:200].decode(errors='replace')}", flush=True)
            if e.code in (401, 403):
                sys.exit("not authorised: stopping")
        time.sleep(1)
    print(f"\n{len(todo) - failed} requested, {failed} failed. Follow them at {CLB_API}/export?datasetKey=<key>, then run fetch_reference.py")


if __name__ == "__main__":
    main()
