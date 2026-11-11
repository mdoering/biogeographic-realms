from pathlib import Path
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


def test_gbif_match_reads_the_accepted_usage_of_a_synonym(api):
    api.routes["/v2/species/match"] = (
        200,
        {
            "synonym": True,
            "usage": {"key": "5G6ZJ", "rank": "SUBSPECIES"},
            "acceptedUsage": {"key": "QLXF", "rank": "SPECIES"},
            "diagnostics": {"matchType": "EXACT"},
        },
    )
    match = gbif_match("Canis lupus familiaris", api.base)
    assert (match.usage_id, match.accepted_id, match.accepted_rank, match.match_type) == ("5G6ZJ", "QLXF", "SPECIES", "EXACT")


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
        "occurrencestatus = 'DETECTED'",
        "specieskey IS NOT NULL",
        "basisofrecord <> 'FOSSIL_SPECIMEN'",
        "coordinateuncertaintyinmeters <= 100000",
    ):
        assert hard in where
    # flags are counted, not filtered: raw counts must still see them
    assert "LIVING_SPECIMEN" not in where
    for column in ("n_raw", "n_datasets", "n_years", "n_nc", "n_open", "n_datasets_open", "n_years_open", "n_ebird", "n_inat", "n_synonym", "n_human_text"):
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


def test_gbif_match_raises_on_404(api):
    # the match API says "no match" with a 200: a 404 means the endpoint moved, not that GBIF doesn't know the name
    with pytest.raises(ApiError):
        gbif_match("Agave rafaeliorum", api.base)


def test_open_counts_name_the_open_licences():
    sql = build_sql(7, PATTERNS)
    assert " AS n_open" in sql
    assert "license IN ('CC0_1_0','CC_BY_4_0')" in sql
    # an unspecified licence must not count as open
    assert "NOT (license" not in sql


def test_kingdom_clash_needs_a_known_kingdom():
    # without this, the record of a species with no kingdom makes the clause NULL and drops out of every clean count
    assert "(v_kingdom IS NULL OR kingdom IS NULL OR LOWER(v_kingdom) NOT IN" in build_sql(7, PATTERNS)


def test_matches_that_cant_be_placed_leave_the_clean_count():
    # a null status must not make the clause NULL, which would drop the record from every clean count
    assert "(taxonomicstatus IS NULL OR LOWER(taxonomicstatus) NOT IN ('ambiguous_synonym','misapplied'))" in build_sql(7, PATTERNS)


def test_records_matched_through_a_synonym_stay_and_are_counted():
    sql = build_sql(7, PATTERNS)
    assert "taxonomicstatus" not in sql[sql.index("WHERE") :]
    assert "LOWER(taxonomicstatus) = 'synonym' THEN 1 ELSE 0 END) AS n_synonym" in sql


def test_sql_pads_terms_with_outer_spaces():
    sql = build_sql(7, (HumanPattern("locality", " zoo ", kingdom="Animalia", unless=("not a zoo",)),))
    assert "CONCAT(' ', LOWER(COALESCE(locality,'')), ' ') LIKE '% zoo %'" in sql
    assert "CONCAT(' ', LOWER(COALESCE(locality,'')), ' ') NOT LIKE '%not a zoo%'" in sql


def test_sql_without_visits_is_the_query_of_the_downloads_so_far():
    sql = build_sql(7, PATTERNS)
    assert "n_visits" not in sql and sql == build_sql(7, PATTERNS, visit_resolution=None)


def test_sql_counts_visits_as_distinct_dates_at_the_visit_resolution():
    sql = build_sql(7, PATTERNS, visit_resolution=17)
    rest = 'END, COALESCE("month",0), COALESCE("day",0), GBIF_ISEA3HCode(17, decimallatitude, decimallongitude, 0)) AS n_visits'
    assert sql.count('THEN COALESCE("year",0) ' + rest) == 2  # n_visits and n_visits_open
    assert "CAST(" not in sql and f"{rest}_open\nFROM occurrence" in sql
    end = sql.index("AS n_human_text") + len("AS n_human_text")
    assert sql[:end] + sql[sql.index("\nFROM occurrence") :] == build_sql(7, PATTERNS)  # nothing else changes


@pytest.mark.parametrize("visit_resolution", [7, 6, 23])
def test_sql_refuses_a_visit_resolution_not_finer_than_the_cells(visit_resolution):
    with pytest.raises(ValueError, match="visit resolution"):
        build_sql(7, PATTERNS, visit_resolution=visit_resolution)


def test_sql_of_the_downloads_before_the_occurrence_status_vocabulary_changed():
    # GBIF replaced Present by Detected on 2026-11-05; with status PRESENT the query is the one of downloads 0020744 and 0003457
    from realms.patterns import HUMAN_PATTERNS

    root = Path(__file__).resolve().parent.parent / "queries"
    assert build_sql(7, HUMAN_PATTERNS, status="PRESENT") + "\n" == (root / "occurrence_res7.sql").read_text()
    assert build_sql(7, HUMAN_PATTERNS, visit_resolution=17, status="PRESENT") + "\n" == (root / "occurrence_res7_visits.sql").read_text()


@pytest.mark.parametrize("status", ["Detected", "DETECTED' OR 1=1 --", ""])
def test_sql_refuses_an_occurrence_status_that_is_not_a_vocabulary_term(status):
    with pytest.raises(ValueError, match="vocabulary term"):
        build_sql(7, PATTERNS, status=status)
