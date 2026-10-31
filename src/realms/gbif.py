"""GBIF and ChecklistBank: which COL release is in use, and the occurrence query.

An API failure raises rather than reading as an answer, so a timeout can never pass as "this id isn't in the release".
"""

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Sequence
from dataclasses import dataclass

from realms.patterns import COLUMNS, HumanPattern

CHECKLIST_KEY = "7ddf754f-d193-4cc9-b351-99906754a03b"
GBIF_API = "https://api.gbif.org"
CLB_API = "https://api.checklistbank.org"
EBIRD = "4fa7b334-ce0d-4e88-aaae-2e0c138d049e"
INATURALIST = "50c9509d-22c7-4a22-a47d-8c48425ef4a7"
USER_AGENT = {"User-Agent": "biogeographic-realms"}
ATTEMPTS = 3


class ApiError(Exception):
    pass


class QueryInvalid(Exception):
    pass


@dataclass(frozen=True)
class Match:
    usage_id: str | None
    rank: str | None
    match_type: str | None


def get_json(url: str) -> tuple[int, dict | None]:
    """Status and JSON body for a 200 or 404. Anything else raises, after retrying server errors and network failures."""
    last = None
    for attempt in range(ATTEMPTS):
        if attempt:
            time.sleep(0.5 * attempt)
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=USER_AGENT), timeout=60) as response:
                return response.status, json.load(response)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 404, None
            last = f"HTTP {e.code}"
            if e.code < 500:
                break
        except (urllib.error.URLError, TimeoutError, ConnectionError, json.JSONDecodeError) as e:
            last = repr(e)
    raise ApiError(f"{url}: {last}")


def production_release(api: str = GBIF_API) -> int:
    """The ChecklistBank key of the COL release GBIF production matches occurrences against."""
    _, body = get_json(f"{api}/v2/species/match/metadata?checklistKey={CHECKLIST_KEY}")
    key = ((body or {}).get("mainIndex") or {}).get("clbDatasetKey")
    if not key:
        raise ApiError(f"no clbDatasetKey in the match metadata: {body}")
    return int(key)


def in_release(release_key: int, usage_id: str, api: str = CLB_API) -> bool:
    status, _ = get_json(f"{api}/dataset/{release_key}/nameusage/{urllib.parse.quote(usage_id)}")
    return status == 200


def gbif_match(name: str, api: str = GBIF_API) -> Match:
    query = urllib.parse.urlencode({"checklistKey": CHECKLIST_KEY, "scientificName": name})
    status, body = get_json(f"{api}/v2/species/match?{query}")
    if status != 200 or body is None:
        # no match is a 200 with matchType NONE, so a 404 means the endpoint itself has gone
        raise ApiError(f"species match answered {status} for {name!r}")
    usage = body.get("usage") or {}
    return Match(usage.get("key"), usage.get("rank"), (body.get("diagnostics") or {}).get("matchType"))


_NON_NATIVE_MEANS = "'introduced','introducedAssistedColonisation','vagrant'"
_NON_NATIVE_DEGREE = (
    "'managed','cultivated','captive','released','failing','casual','unestablished','reproducing','established',"
    "'colonising','spreading','naturalized','invasive','widespreadInvasive'"
)
_CLASH_KINGDOMS = "'animalia','plantae','fungi'"
# A pro parte name could be any of several taxa and a misapplied one was used for another taxon, so the accepted species GBIF gives is a guess.
_UNPLACEABLE = "'ambiguous_synonym','misapplied'"

# Layer 1 of the cleaning: the flags GBIF can check. Everything else is counted around this.
_CLEAN = (
    "basisofrecord <> 'LIVING_SPECIMEN' "
    f"AND (establishmentmeans.concept IS NULL OR establishmentmeans.concept NOT IN ({_NON_NATIVE_MEANS})) "
    f"AND (degreeofestablishment.concept IS NULL OR degreeofestablishment.concept NOT IN ({_NON_NATIVE_DEGREE})) "
    "AND NOT COALESCE(ARRAY_CONTAINS(issue,'TAXON_MATCH_FUZZY'),FALSE) "
    f"AND (v_kingdom IS NULL OR kingdom IS NULL OR LOWER(v_kingdom) NOT IN ({_CLASH_KINGDOMS}) OR LOWER(kingdom) NOT IN ({_CLASH_KINGDOMS}) "
    "OR LOWER(v_kingdom) = LOWER(kingdom)) "
    f"AND (taxonomicstatus IS NULL OR LOWER(taxonomicstatus) NOT IN ({_UNPLACEABLE}))"
)
_NC = "license = 'CC_BY_NC_4_0'"
_OPEN = "license IN ('CC0_1_0','CC_BY_4_0')"
_PART = "CASE WHEN level0gid IS NULL THEN 'sea' ELSE 'land' END"


def _sql_literal(term: str) -> str:
    return term.replace("'", "''")


def _pattern_sql(pattern: HumanPattern) -> str:
    if pattern.column not in COLUMNS:
        raise ValueError(f"unknown text column {pattern.column}")
    for term in (pattern.term, *pattern.unless):
        if "%" in term or "_" in term:
            raise ValueError(f"wildcard in text term {term!r}")
    text = f"LOWER(COALESCE({pattern.column},''))"
    if pattern.term != pattern.term.strip():
        # a term with outer spaces must also match at the very start and end of the text
        text = f"CONCAT(' ', {text}, ' ')"
    parts = [f"LOWER(kingdom) = '{pattern.kingdom.lower()}'"] if pattern.kingdom else []
    parts.append(f"{text} LIKE '%{_sql_literal(pattern.term)}%'")
    parts += [f"{text} NOT LIKE '%{_sql_literal(phrase)}%'" for phrase in pattern.unless]
    return f"({' AND '.join(parts)})"


def build_sql(resolution: int, patterns: Sequence[HumanPattern]) -> str:
    """The occurrence query: species x cell x land/sea x country, counting each filter rather than applying it."""
    if not 1 <= resolution <= 22:
        raise ValueError(f"ISEA3H resolution must be 1 to 22, not {resolution}")
    human = " OR ".join(_pattern_sql(p) for p in patterns)
    cell = f"GBIF_ISEA3HCode({resolution}, decimallatitude, decimallongitude, 0)"
    return f"""SELECT specieskey, countrycode, {cell} AS cell, {_PART} AS part,
 COUNT(*) AS n_raw,
 SUM(CASE WHEN {_CLEAN} THEN 1 ELSE 0 END) AS n,
 COUNT(DISTINCT CASE WHEN {_CLEAN} THEN datasetkey END) AS n_datasets,
 COUNT(DISTINCT CASE WHEN {_CLEAN} THEN "year" END) AS n_years,
 SUM(CASE WHEN ({_CLEAN}) AND {_NC} THEN 1 ELSE 0 END) AS n_nc,
 SUM(CASE WHEN ({_CLEAN}) AND {_OPEN} THEN 1 ELSE 0 END) AS n_open,
 COUNT(DISTINCT CASE WHEN ({_CLEAN}) AND {_OPEN} THEN datasetkey END) AS n_datasets_open,
 COUNT(DISTINCT CASE WHEN ({_CLEAN}) AND {_OPEN} THEN "year" END) AS n_years_open,
 SUM(CASE WHEN ({_CLEAN}) AND datasetkey = '{EBIRD}' THEN 1 ELSE 0 END) AS n_ebird,
 SUM(CASE WHEN ({_CLEAN}) AND datasetkey = '{INATURALIST}' THEN 1 ELSE 0 END) AS n_inat,
 SUM(CASE WHEN ({_CLEAN}) AND LOWER(taxonomicstatus) = 'synonym' THEN 1 ELSE 0 END) AS n_synonym,
 SUM(CASE WHEN ({_CLEAN}) AND ({human}) THEN 1 ELSE 0 END) AS n_human_text
FROM occurrence
WHERE hascoordinate = TRUE AND hasgeospatialissues = FALSE AND occurrencestatus = 'PRESENT' AND specieskey IS NOT NULL
 AND basisofrecord <> 'FOSSIL_SPECIMEN' AND (coordinateuncertaintyinmeters IS NULL OR coordinateuncertaintyinmeters <= 100000)
GROUP BY specieskey, countrycode, {cell}, {_PART}"""


def validate_sql(sql: str, api: str = GBIF_API) -> None:
    """Ask GBIF to check the query. This submits nothing and needs no login."""
    request = urllib.request.Request(
        f"{api}/v1/occurrence/download/request/validate",
        data=json.dumps({"format": "SQL_TSV_ZIP", "sql": sql}).encode(),
        headers={**USER_AGENT, "Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            if response.status == 201:
                return
            raise QueryInvalid(f"HTTP {response.status}: {response.read().decode(errors='replace')[:2000]}")
    except urllib.error.HTTPError as e:
        raise QueryInvalid(f"HTTP {e.code}: {e.read().decode(errors='replace')[:2000]}") from e
