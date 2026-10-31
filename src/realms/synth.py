"""Synthetic networks with a known answer."""

import numpy as np
import polars as pl


def planted_network(
    n_realms: int = 4,
    cells_per_realm: int = 30,
    endemic_species: int = 200,
    widespread_species: int = 40,
    occupancy: float = 0.3,
    seed: int = 1,
) -> tuple[pl.DataFrame, pl.DataFrame, dict[str, int]]:
    """Realms of cells with their own endemic species and family, plus widespread species across all of them.

    `endemic_species` is per realm. Each species occupies `occupancy` of its range, at least two cells. Returns the presence
    table (`cell, species_id`), the taxa table (`species_id, genus, family`) and the true realm of each cell.
    """
    rng = np.random.default_rng(seed)
    realms = {r: [f"r{r}-c{i}" for i in range(cells_per_realm)] for r in range(n_realms)}
    all_cells = [c for cells in realms.values() for c in cells]

    def occupy(species: str, cells: list[str]) -> list[tuple[str, str]]:
        n = min(len(cells), max(2, round(occupancy * len(cells))))
        return [(str(c), species) for c in rng.choice(cells, size=n, replace=False)]

    presence: list[tuple[str, str]] = []
    taxa: list[tuple[str, str, str]] = []
    for r, cells in realms.items():
        for s in range(endemic_species):
            species = f"e{r}-{s}"
            presence += occupy(species, cells)
            taxa.append((species, f"G{r}-{s // 10}", f"F{r}"))
    for s in range(widespread_species):
        species = f"w{s}"
        presence += occupy(species, all_cells)
        taxa.append((species, f"Gw-{s // 10}", "W"))

    truth = {c: r for r, cells in realms.items() for c in cells}
    return (
        pl.DataFrame(presence, schema={"cell": pl.String, "species_id": pl.String}, orient="row"),
        pl.DataFrame(taxa, schema={"species_id": pl.String, "genus_id": pl.String, "family_id": pl.String}, orient="row"),
        truth,
    )
