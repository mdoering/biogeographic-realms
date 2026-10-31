"""The species-by-cell matrix published with Costello et al. 2017, on 5 degree c-squares."""

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
from scipy.sparse import csr_matrix

EARTH_RADIUS_KM = 6371.0088
SIDE = 5.0


def csquare_bounds(code: str) -> tuple[float, float, float, float]:
    """(lat_min, lat_max, lon_min, lon_max) of a 5 degree c-square such as "3603:3".

    The first digit is the global quadrant (1 NE, 3 SE, 5 SW, 7 NW), the next is tens of latitude and the next two are tens
    of longitude. The digit after the colon picks the quarter: 1 low lat low lon, 2 low lat high lon, 3 high lat low lon, 4
    high lat high lon, where low and high are in absolute degrees.
    """
    head, quarter = code.split(":")
    quadrant, lat_tens, lon_tens, q = int(head[0]), int(head[1]), int(head[2:4]), int(quarter)
    lat_low = lat_tens * 10 + (SIDE if q in (3, 4) else 0.0)
    lon_low = lon_tens * 10 + (SIDE if q in (2, 4) else 0.0)
    lat = (-(lat_low + SIDE), -lat_low) if quadrant in (3, 5) else (lat_low, lat_low + SIDE)
    lon = (-(lon_low + SIDE), -lon_low) if quadrant in (5, 7) else (lon_low, lon_low + SIDE)
    return (float(lat[0]), float(lat[1]), float(lon[0]), float(lon[1]))


def csquare_area_km2(code: str) -> float:
    lat_min, lat_max, lon_min, lon_max = csquare_bounds(code)
    sine_span = math.sin(math.radians(lat_max)) - math.sin(math.radians(lat_min))
    return EARTH_RADIUS_KM**2 * math.radians(lon_max - lon_min) * sine_span


def csquare_adjacency(cells: list[str]) -> dict[str, set[str]]:
    """Cells sharing an edge, across the equator and the dateline too."""
    by_corner = {(b[0], b[2]): c for c, b in ((c, csquare_bounds(c)) for c in cells)}
    adjacency: dict[str, set[str]] = {c: set() for c in cells}
    for (lat, lon), cell in by_corner.items():
        east = lon + SIDE if lon + SIDE < 180.0 else -180.0
        for corner in ((lat, east), (lat + SIDE, lon)):
            neighbour = by_corner.get(corner)
            if neighbour is not None and neighbour != cell:
                adjacency[cell].add(neighbour)
                adjacency[neighbour].add(cell)
    return adjacency


@dataclass(frozen=True)
class Matrix:
    species: list[str]
    cells: list[str]
    X: csr_matrix  # species x cell

    def presence(self) -> pl.DataFrame:
        coo = self.X.tocoo()
        cells = np.array(self.cells, dtype=object)[coo.col].tolist()
        species = np.array(self.species, dtype=object)[coo.row].tolist()
        return pl.DataFrame({"cell": cells, "species_id": species}, schema={"cell": pl.String, "species_id": pl.String})

    def genera(self) -> pl.DataFrame:
        """The matrix holds names and no taxonomy, so the first word has to stand in for a genus id, and family stays null."""
        return pl.DataFrame(
            {
                "species_id": self.species,
                "genus_id": [name.split(" ", 1)[0] for name in self.species],
                "family_id": [None] * len(self.species),
            },
            schema={"species_id": pl.String, "genus_id": pl.String, "family_id": pl.String},
        )


def load_matrix(path: Path) -> Matrix:
    df = pl.read_csv(path, infer_schema_length=0)
    cells = df.columns[1:]
    values = df.select(pl.col(cells).cast(pl.Int8)).to_numpy()
    return Matrix(species=df[df.columns[0]].to_list(), cells=cells, X=csr_matrix(values, dtype=np.int8))
