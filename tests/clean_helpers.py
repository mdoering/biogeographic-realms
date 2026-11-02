"""Hand-made tables for the tests of `realms.clean` and `realms.clean_report`.

The grid has cells named by numbers, so that "lowest cell id" can be told from the lowest id as text ("9" before "100"). A cell is
`(land_km2, sea_km2)`, and a part with less than `ZERO_KM2` of area is none.
"""

import polars as pl

from realms.clean import Cleaned, Settings
from realms.grid import Grid
from realms.natives import References
from realms.occurrences import COUNTS, DISTINCT

# the design's layer 5 thresholds, which were the defaults of `Settings` until setting 132 of the calibration became them (Markus, 2026-11-01). A test of
# what the design's values do passes this, and does not lean on the defaults
DESIGN_5 = Settings(satellite_record_share=0.01, satellite_cell_share=0.05, under_rate_factor=0.01, group_min_records=1000)


def make_grid(cells: dict[str, tuple[float, float]], edges: list[tuple[str, str]] = (), centres: dict[str, tuple[float, float]] | None = None) -> Grid:
    """A grid of `cells` (name: (land_km2, sea_km2)) with `edges` for neighbours. Cells sit a degree of longitude apart on the equator unless
    `centres` gives a (lat, lon) for a name."""
    names = list(cells)
    centres = centres or {}
    table = pl.DataFrame(
        {
            "cell": names,
            "lat": [centres[c][0] if c in centres else 0.0 for c in names],
            "lon": [centres[c][1] if c in centres else float(i) for i, c in enumerate(names)],
            "area_km2": [land + sea for land, sea in cells.values()],
            "land_km2": [land for land, _ in cells.values()],
            "sea_km2": [sea for _, sea in cells.values()],
        }
    )
    adjacency = {c: set() for c in names}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    return Grid(6, table, adjacency)


# 10 and 13 are land only, 11 and 9 are coastal, 12 is open sea, 100 is land only. In a chain 10 - 11 - 12 - 13.
GRID = make_grid({"10": (1000.0, 0.0), "11": (500.0, 500.0), "12": (0.0, 1000.0), "13": (1000.0, 0.0)}, [("10", "11"), ("11", "12"), ("12", "13")])


def row(species, part, country="AA", n=10, **given):
    """A row of country counts. Every count is `n` but the distinct ones, which are 1, and the human-text one, which is 0."""
    values = {c: (1 if c in DISTINCT else n) for c in COUNTS} | {"n_human_text": 0} | given
    return {"species_id": species, "countrycode": country, "part": part, **values}


def counts(*rows) -> pl.DataFrame:
    schema = {"species_id": pl.String, "countrycode": pl.String, "part": pl.String, **{c: pl.Int64 for c in COUNTS}}
    return pl.DataFrame(list(rows), schema=schema)


def taxon(species, kingdom="Animalia", cls="Insecta", environment=None, extinct=False):
    return {"species_id": species, "kingdom": kingdom, "class": cls, "environment": environment, "extinct": extinct}


def refs(
    taxa=(), groups=None, wcvp=(), griis=(), wrims=(), tdwg3=(), admin1=(), griis_no_cells=(), wrims_no_country=(), no_wild_range=(), griis_name_only=(), griis_name_only_no_cells=()
) -> References:
    taxa = pl.DataFrame(list(taxa), schema={"species_id": pl.String, "kingdom": pl.String, "class": pl.String, "environment": pl.String, "extinct": pl.Boolean})
    griis_schema = {"species_id": pl.String, "country": pl.String, "subdivision": pl.String}
    wrims_schema = {"species_id": pl.String, "area_id": pl.String, "area": pl.String, "country": pl.String}
    if groups is None:
        groups = [(s, "Insects") for s in taxa["species_id"]]
    return References(
        taxa=taxa,
        groups=pl.DataFrame(list(groups), schema={"species_id": pl.String, "tax_group": pl.String}, orient="row"),
        wcvp=pl.DataFrame(list(wcvp), schema={"species_id": pl.String, "tdwg3": pl.String, "native": pl.Boolean}, orient="row"),
        griis=pl.DataFrame(list(griis), schema=griis_schema, orient="row"),
        wrims=pl.DataFrame(list(wrims), schema=wrims_schema, orient="row"),
        tdwg3_parts=pl.DataFrame(list(tdwg3), schema={"part": pl.String, "tdwg3": pl.String}, orient="row"),
        admin1_cells=pl.DataFrame(list(admin1), schema={"cell": pl.String, "country": pl.String, "subdivision": pl.String, "via": pl.String}, orient="row"),
        griis_no_cells=pl.DataFrame(list(griis_no_cells), schema=griis_schema, orient="row"),
        griis_name_only=pl.DataFrame(list(griis_name_only), schema=griis_schema, orient="row"),
        griis_name_only_no_cells=pl.DataFrame(list(griis_name_only_no_cells), schema=griis_schema, orient="row"),
        wrims_no_country=pl.DataFrame(list(wrims_no_country), schema=wrims_schema, orient="row"),
        no_wild_range=frozenset(no_wild_range),
    )


def kept(cleaned: Cleaned, column="records") -> dict:
    return {(s, p): v for s, p, v in cleaned.presence.select("species_id", "part", column).iter_rows()}


def gone(cleaned: Cleaned, layer=None, rule=None) -> set:
    removed = cleaned.removed
    if layer is not None:
        removed = removed.filter(pl.col("layer") == layer)
    if rule is not None:
        removed = removed.filter(pl.col("rule") == rule)
    return {(s, p, c) for s, p, c in removed.select("species_id", "part", "countrycode").iter_rows()}
