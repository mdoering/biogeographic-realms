"""Generate the occurrence queries and have GBIF validate them. This submits no download.

Also the resolution 7 query with visit counts (`build_sql`'s `visit_resolution`, DECISIONS "Species-days", 2026-11-06): a visit is a species on one
date in one ISEA3H cell of `VISIT_RESOLUTION`, about 0.4 km2 (resolution 7 is about 23,000 km2, and each step divides by 3), near the 600 m and one day
below which van Merriënboer et al. (arXiv:2611.00284) found two GBIF records no more likely to be different species.
"""

from pathlib import Path

from realms.gbif import build_sql, validate_sql
from realms.patterns import HUMAN_PATTERNS

VISIT_RESOLUTION = 17

queries = {f"queries/occurrence_res{r}.sql": build_sql(r, HUMAN_PATTERNS) for r in (6, 7)}
queries["queries/occurrence_res7_visits.sql"] = build_sql(7, HUMAN_PATTERNS, visit_resolution=VISIT_RESOLUTION)
for name, sql in queries.items():
    validate_sql(sql)
    path = Path(name)
    path.parent.mkdir(exist_ok=True)
    path.write_text(sql + "\n")
    print(f"valid: {path}, {len(sql)} characters, {len(HUMAN_PATTERNS)} text patterns")
