"""Check GBIF's ISEA3H grid files at resolutions 5, 6 and 7 against what the grid must be, and print the cell areas.

An ISEA3H grid has 10 * 3**r + 2 cells, of which 12 are pentagons and the rest hexagons, and the cells tile the whole globe. The
land is measured against Natural Earth, whose own area it must come close to.
"""

import collections
from pathlib import Path

import geopandas as gpd
import pyproj
import shapely

from realms.grid import load_grid

ROOT = Path("data/reference")
LAND = ROOT / "natural_earth" / "ne_10m_land.zip"
RESOLUTIONS = (5, 6, 7)
EARTH_KM2 = 510.06e6
AREA_TOLERANCE = 0.0005
LAND_TOLERANCE = 0.02


def natural_earth_km2() -> float:
    geod = pyproj.Geod(ellps="WGS84")
    return sum(abs(geod.geometry_area_perimeter(part)[0]) for shape in gpd.read_file(LAND).geometry for part in shapely.get_parts(shape)) / 1e6


def main() -> None:
    land_km2 = natural_earth_km2()
    print(f"Natural Earth land: {land_km2:,.0f} km2")
    for resolution in RESOLUTIONS:
        grid = load_grid(ROOT / "grids" / f"ISEA3H-{resolution}.gpkg", LAND)
        cells = grid.cells
        expected = 10 * 3**resolution + 2
        assert grid.resolution == resolution, f"{grid.resolution=}"
        assert cells.height == expected, f"resolution {resolution}: {cells.height} cells, expected {expected}"
        assert cells["cell"].n_unique() == expected, f"resolution {resolution}: ids are not unique"

        neighbours = collections.Counter(len(near) for near in grid.adjacency.values())
        assert neighbours == {5: 12, 6: expected - 12}, f"resolution {resolution}: cells by number of neighbours {dict(neighbours)}"

        area = cells["area_km2"]
        assert abs(area.sum() / EARTH_KM2 - 1) < AREA_TOLERANCE, f"resolution {resolution}: cells add up to {area.sum():,.0f} km2"
        assert abs(cells["land_km2"].sum() / land_km2 - 1) < LAND_TOLERANCE, f"resolution {resolution}: {cells['land_km2'].sum():,.0f} km2 of land"
        assert (cells["land_km2"] >= 0).all() and (cells["sea_km2"] >= 0).all()

        coastal = cells.filter((cells["land_km2"] > 0) & (cells["sea_km2"] > 0)).height
        print(
            f"resolution {resolution}: {expected:,} cells, 12 with 5 neighbours and {neighbours[6]:,} with 6; "
            f"median cell {area.median():,.0f} km2 (range {area.min():,.0f} to {area.max():,.0f}); "
            f"{coastal:,} coastal cells; total area {area.sum():,.0f} km2; "
            f"land {cells['land_km2'].sum():,.0f} km2 ({cells['land_km2'].sum() / land_km2 - 1:+.4%} against Natural Earth)"
        )
    print("ok")


if __name__ == "__main__":
    main()
