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


def grid_presence(n_species: int, side: int = 90, seed: int = 1) -> tuple[pl.DataFrame, pl.DataFrame, dict[str, set[str]]]:
    """Species with round ranges on a `side` x `side` grid of cells, for measuring how the pipeline scales.

    A species gets a random centre and a log-normal radius (capped at 60 cells), and takes each cell within that radius with
    probability 0.6, or its centre alone if none is taken. Returns the presence table (`cell, species_id`), the taxa table
    (`species_id, genus_id, family_id`, five species to a genus and a hundred to a family) and each cell's four neighbours.
    The draws come in a fixed order, so a given seed always gives the same ranges.
    """
    rng = np.random.default_rng(seed)
    n_cells = side * side
    centres = rng.integers(n_cells, size=n_species)
    radius = np.minimum(rng.lognormal(0.5, 1.0, size=n_species), 60)

    ranges = []
    for s in range(n_species):
        r = radius[s]
        row, col = divmod(centres[s], side)
        reach = int(np.ceil(r))
        grid_row, grid_col = np.meshgrid(
            np.arange(max(0, row - reach), min(side, row + reach + 1)),
            np.arange(max(0, col - reach), min(side, col + reach + 1)),
            indexing="ij",
        )
        within = (grid_row * side + grid_col)[(grid_row - row) ** 2 + (grid_col - col) ** 2 <= r * r]
        taken = within[rng.random(len(within)) < 0.6]
        ranges.append((taken if len(taken) else np.array([centres[s]])).astype(np.int32))

    cell_index = np.concatenate(ranges)
    species_index = np.repeat(np.arange(n_species, dtype=np.int32), [len(r) for r in ranges])

    cell_names = pl.select(("c" + pl.int_range(n_cells).cast(pl.String).str.zfill(len(str(n_cells - 1)))).alias("cell")).to_series()
    species = pl.int_range(n_species)
    taxa = pl.select(
        ("s" + species.cast(pl.String)).alias("species_id"),
        ("g" + (species // 5).cast(pl.String)).alias("genus_id"),
        ("f" + (species // 100).cast(pl.String)).alias("family_id"),
    )
    presence = pl.DataFrame({"cell": cell_names.gather(cell_index), "species_id": taxa["species_id"].gather(species_index)})

    names = cell_names.to_list()
    adjacency: dict[str, set[str]] = {name: set() for name in names}
    for i in range(n_cells):
        for j in (i + 1 if (i + 1) % side else None, i + side if i + side < n_cells else None):  # right and below
            if j is not None:
                adjacency[names[i]].add(names[j])
                adjacency[names[j]].add(names[i])
    return presence, taxa, adjacency
