"""Which COL release does GBIF production use, and does it behave like the metadata says?

Takes species ids from the newest public extended release that don't exist in the release GBIF reports. If GBIF returns any of
them as its own match, the metadata understates the release and this exits non-zero. Requests are sequential.
"""

import json
import sys
from pathlib import Path

import numpy as np

from realms.gbif import CLB_API, ApiError, gbif_match, get_json, in_release, production_release

OUT = Path("results/release.json")
PAGES = 8
PAGE_SIZE = 200
MAX_OFFSET = 99_000


def main() -> None:
    production = production_release()
    _, info = get_json(f"{CLB_API}/dataset/{production}")
    print(f"GBIF production: {production} ({info.get('alias')}, {info.get('version')})")

    _, listing = get_json(f"{CLB_API}/dataset?releasedFrom=3&sortBy=created&limit=8")
    extended = [d for d in listing["result"] if d.get("origin") == "xrelease" and not d.get("private")]
    newest = extended[0]
    print(f"newest public extended release: {newest['key']} ({newest.get('alias')}, {newest.get('version')})")

    probe: list[dict] = []
    if newest["key"] != production:
        rng = np.random.default_rng(7)
        sample: dict[str, str] = {}
        for offset in rng.choice(MAX_OFFSET // PAGE_SIZE, size=PAGES, replace=False) * PAGE_SIZE:
            _, page = get_json(
                f"{CLB_API}/dataset/{newest['key']}/nameusage/search?rank=species&status=accepted&limit={PAGE_SIZE}&offset={offset}"
            )
            sample |= {r["usage"]["id"]: r["usage"]["name"]["scientificName"] for r in page["result"]}
        print(f"sampled {len(sample)} species ids", flush=True)
        for usage_id, name in sample.items():
            if not in_release(production, usage_id):
                match = gbif_match(name)
                probe.append({"id": usage_id, "name": name, "gbif_id": match.usage_id, "gbif_rank": match.rank, "gbif_match": match.match_type})

    known = [p for p in probe if p["gbif_id"] == p["id"]]
    for row in probe:
        print(f"  {row['id']:8} {row['name']:44} -> {row['gbif_id']} {row['gbif_rank']} {row['gbif_match']}")
    print(f"{len(probe)} ids absent from {production}, of which GBIF knows {len(known)}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        json.dumps(
            {"gbif_production": production, "alias": info.get("alias"), "version": info.get("version"), "newest_extended": newest["key"], "probe": probe},
            indent=2,
        )
    )
    if known:
        sys.exit(f"GBIF knows {len(known)} ids that release {production} doesn't have: the metadata understates the release")


if __name__ == "__main__":
    try:
        main()
    except ApiError as e:
        sys.exit(f"API failure, no conclusion: {e}")
