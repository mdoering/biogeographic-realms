"""Generate the occurrence queries and have GBIF validate them. This submits no download."""

from pathlib import Path

from realms.gbif import build_sql, validate_sql
from realms.patterns import HUMAN_PATTERNS

for resolution in (6, 7):
    sql = build_sql(resolution, HUMAN_PATTERNS)
    validate_sql(sql)
    path = Path(f"queries/occurrence_res{resolution}.sql")
    path.parent.mkdir(exist_ok=True)
    path.write_text(sql + "\n")
    print(f"valid: {path}, {len(sql)} characters, {len(HUMAN_PATTERNS)} text patterns")
