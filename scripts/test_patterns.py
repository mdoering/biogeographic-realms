"""Sample real records for each text term, to judge by reading how often it means a garden, zoo or captive organism.

The sample comes from GBIF's free text search, so it shows how precise a term is, not how common it is. Hits follow the same
kingdom and exclusion rules as the query.
"""

import random
import re
import urllib.parse
from pathlib import Path

from realms.gbif import GBIF_API, get_json
from realms.patterns import HUMAN_PATTERNS, is_human_held, like_match

OUT = Path("results/patterns.txt")
FIELDS = {"locality": "locality", "habitat": "habitat", "occurrenceremarks": "occurrenceRemarks"}
FETCH = 300
SHOW = 15
WIDTH = 110


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rng = random.Random(11)
    lines = []
    for term in sorted({p.term for p in HUMAN_PATTERNS}):
        query = urllib.parse.urlencode({"q": term.strip(), "limit": FETCH, "hasCoordinate": "true"})
        _, body = get_json(f"{GBIF_API}/v1/occurrence/search?{query}")
        hits = []
        for record in body["results"]:
            text = {column: record.get(field) for column, field in FIELDS.items()}
            if not is_human_held(text, record.get("kingdom")):
                continue
            for column, value in text.items():
                if like_match(value, term):
                    start = max(0, value.lower().index(term) - 45)
                    hits.append(f"{column[:4]} {str(record.get('kingdom'))[:4]} | {re.sub(r'\s+', ' ', value[start : start + WIDTH])}")
                    break
        lines.append(f"\n=== {term!r}: {len(hits)} of {len(body['results'])} fetched records count as held by people")
        lines += rng.sample(hits, min(SHOW, len(hits)))
    OUT.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
