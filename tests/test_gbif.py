import json
import urllib.parse

import pytest

from realms.gbif import ApiError, CHECKLIST_KEY, QueryInvalid, build_sql, gbif_match, in_release, production_release, validate_sql
from realms.patterns import HumanPattern

PATTERNS = (HumanPattern("locality", "botanic"), HumanPattern("habitat", "cultivated"))


def test_production_release_reads_main_index(api):
    api.routes["/v2/species/match/metadata"] = (200, {"created": "2026-07-18", "mainIndex": {"clbDatasetKey": "315557", "datasetKey": CHECKLIST_KEY}})
    assert production_release(api.base) == 315557
    assert f"checklistKey={CHECKLIST_KEY}" in api.requests[0][1]


def test_production_release_without_a_key_is_an_error(api):
    api.routes["/v2/species/match/metadata"] = (200, {"mainIndex": {}})
    with pytest.raises(ApiError):
        production_release(api.base)


def test_in_release_true_on_200_false_on_404(api):
    api.routes["/dataset/315557/nameusage/C56N"] = (200, {"id": "C56N"})
    assert in_release(315557, "C56N", api.base) is True
    assert in_release(315557, "GONE", api.base) is False


def test_in_release_raises_on_500(api):
    api.routes["/dataset/315557/nameusage/C56N"] = (500, {})
    with pytest.raises(ApiError):
        in_release(315557, "C56N", api.base)


def test_in_release_raises_when_the_server_is_unreachable():
    with pytest.raises(ApiError):
        in_release(315557, "C56N", "http://127.0.0.1:9")


def test_gbif_match_reads_usage_and_diagnostics(api):
    api.routes["/v2/species/match"] = (
        200,
        {"synonym": False, "usage": {"key": "8VR75", "name": "Agave L.", "rank": "GENUS"}, "diagnostics": {"matchType": "HIGHERRANK", "confidence": 98}},
    )
    match = gbif_match("Agave rafaeliorum", api.base)
    assert (match.usage_id, match.rank, match.match_type) == ("8VR75", "GENUS", "HIGHERRANK")
    query = urllib.parse.parse_qs(urllib.parse.urlparse(api.requests[0][1]).query)
    assert query == {"checklistKey": [CHECKLIST_KEY], "scientificName": ["Agave rafaeliorum"]}


def test_gbif_match_with_no_usage(api):
    api.routes["/v2/species/match"] = (200, {"diagnostics": {"matchType": "NONE"}})
    match = gbif_match("Nothing", api.base)
    assert (match.usage_id, match.rank, match.match_type) == (None, None, "NONE")


def test_sql_contains_resolution_and_zero_uncertainty():
    sql = build_sql(6, PATTERNS)
    assert "GBIF_ISEA3HCode(6, decimallatitude, decimallongitude, 0)" in sql
    assert "GBIF_ISEA3HCode(7" not in build_sql(6, PATTERNS)


def test_sql_escapes_apostrophes():
    sql = build_sql(7, (HumanPattern("locality", "jardin d'essai"),))
    assert "LIKE '%jardin d''essai%'" in sql


def test_sql_rejects_wildcards_in_a_term():
    with pytest.raises(ValueError):
        build_sql(7, (HumanPattern("locality", "100%"),))


def test_sql_rejects_unknown_pattern_column():
    with pytest.raises(ValueError):
        build_sql(7, (HumanPattern("recordedby", "zoo"),))


def test_sql_rejects_a_resolution_outside_the_grid():
    with pytest.raises(ValueError):
        build_sql(0, PATTERNS)


def test_sql_hard_filters_and_columns_present():
    sql = build_sql(7, PATTERNS)
    where = sql[sql.index("WHERE") :]
    for hard in (
        "hascoordinate = TRUE",
        "hasgeospatialissues = FALSE",
        "occurrencestatus = 'PRESENT'",
        "specieskey IS NOT NULL",
        "basisofrecord <> 'FOSSIL_SPECIMEN'",
        "coordinateuncertaintyinmeters <= 100000",
    ):
        assert hard in where
    # flags are counted, not filtered: raw counts must still see them
    assert "LIVING_SPECIMEN" not in where
    for column in ("n_raw", "n_datasets", "n_years", "n_nc", "n_datasets_open", "n_years_open", "n_ebird", "n_inat", "n_human_text"):
        assert f" AS {column}" in sql
    assert " AS n," in sql
    assert sql.rstrip().endswith("END")
    assert "GROUP BY specieskey, countrycode," in sql


def test_sql_patterns_only_feed_the_human_text_count():
    sql = build_sql(7, PATTERNS)
    assert sql.count("LIKE '%botanic%'") == 1
    assert "LOWER(COALESCE(locality,'')) LIKE '%botanic%'" in sql
    assert sql.index("LIKE '%botanic%'") > sql.index(" AS n_inat,")


def test_sql_limits_a_pattern_to_its_kingdom():
    sql = build_sql(7, (HumanPattern("locality", "botanic", kingdom="Plantae"),))
    assert "(LOWER(kingdom) = 'plantae' AND LOWER(COALESCE(locality,'')) LIKE '%botanic%')" in sql


def test_sql_excludes_negated_phrases():
    sql = build_sql(7, (HumanPattern("habitat", "kultiviert", unless=("nicht kultiviert",)),))
    assert "(LOWER(COALESCE(habitat,'')) LIKE '%kultiviert%' AND LOWER(COALESCE(habitat,'')) NOT LIKE '%nicht kultiviert%')" in sql


def test_sql_escapes_apostrophes_in_exclusions():
    sql = build_sql(7, (HumanPattern("habitat", "cultivé", unless=("n'est pas cultivé",)),))
    assert "NOT LIKE '%n''est pas cultivé%'" in sql


def test_sql_rejects_wildcards_in_an_exclusion():
    with pytest.raises(ValueError):
        build_sql(7, (HumanPattern("habitat", "x", unless=("a_b",)),))


def test_kingdom_clash_covers_only_animal_plant_fungi():
    sql = build_sql(7, PATTERNS)
    assert "LOWER(v_kingdom) NOT IN ('animalia','plantae','fungi')" in sql
    assert "LOWER(kingdom) NOT IN ('animalia','plantae','fungi')" in sql


def test_validate_sql_accepts_201(api):
    api.routes["/v1/occurrence/download/request/validate"] = (201, {"sql": "SELECT 1"})
    validate_sql("SELECT 1", api.base)
    method, _, body = api.requests[0]
    assert method == "POST"
    assert json.loads(body) == {"format": "SQL_TSV_ZIP", "sql": "SELECT 1"}


def test_validate_sql_raises_with_the_server_message(api):
    api.routes["/v1/occurrence/download/request/validate"] = (400, {"error": "bad column"})
    with pytest.raises(QueryInvalid, match="bad column"):
        validate_sql("SELECT nonsense", api.base)
