"""Nested levels of realms instead of one realm count (DECISIONS, "Nested levels", 2026-11-05).

The merge heights ("Is there a natural k?") single out no k, but a floor (five blocks about equally distinct) and an order in which the smaller
realms join. So the realms are reported at nested levels, and a unit is a realm by how distinct it is, not by the 1M km2 size floor:

- `leaf_realms`: Infomap's units, merged as `candidates.merge_sequence` merges them (the pair with the lowest beta-sim, on the union of their taxa),
  for as long as some pair is below `cut`. What is left are the finest realms. Small islands join the realm they share most with.
- `nested_merges`: the leaf realms merged on in the same way, but land only with land and sea only with sea (`sea_realms`, by area), the two joined last.
  Unconstrained, the vast ocean realm, whose coasts hold taxa from everywhere, takes in Madagascar, Australia and South and East Asia before any
  land realm does. Marine and terrestrial realms are two systems in most schemes, and the data put the sea's join with land at the top anyway.
- `at_cuts`: the map of each level, cut from the nested merges at a beta-sim. A cut never joins land and sea.

`candidates.py`, the selection and the Costello gate are unchanged: this is a separate path.
"""

import numpy as np
import polars as pl

from realms.candidates import _Merger
from realms.endemism import Labels
from realms.matrix import CellTaxa

SEA_SUFFIX = ":sea"
MERGES = {"k": pl.Int64, "height": pl.Float64, "keep": pl.Int64, "drop": pl.Int64, "kind": pl.String}


def _relabel(current: Labels, keep: int, drop: int) -> Labels:
    return {c: (keep if r == drop else r) for c, r in current.items()}


def leaf_realms(units: Labels, ct: CellTaxa, cut: float) -> tuple[Labels, pl.DataFrame]:
    """`units` merged as `merge_sequence` merges them while the closest pair's beta-sim is below `cut`, and the merges made (`MERGES`, `kind` null).

    A row's `k` is the count before its merge. The merged unit keeps the label of `keep`, as in `merge_sequence`.
    """
    merger = _Merger(units, ct)
    current = dict(units)
    rows = []
    while merger.active.sum() > 1:
        keep, drop = merger.closest_pair()
        height = float(merger.beta_sim()[keep, drop])
        if height >= cut:
            break
        rows.append({"k": int(merger.active.sum()), "height": height, "keep": merger.units[keep], "drop": merger.units[drop], "kind": None})
        current = _relabel(current, merger.units[keep], merger.units[drop])
        merger.merge(keep, drop)
    return current, pl.DataFrame(rows, schema=MERGES)


def sea_realms(labels: Labels, areas: dict[str, float], share: float = 0.5) -> dict[int, bool]:
    """Whether each realm is a sea realm: at least `share` of its area in sea parts (units whose id ends in `SEA_SUFFIX`)."""
    total: dict[int, float] = {}
    sea: dict[int, float] = {}
    for unit, realm in labels.items():
        total[realm] = total.get(realm, 0.0) + areas[unit]
        if unit.endswith(SEA_SUFFIX):
            sea[realm] = sea.get(realm, 0.0) + areas[unit]
    return {r: total[r] > 0 and sea.get(r, 0.0) / total[r] >= share for r in total}


def nested_merges(labels: Labels, ct: CellTaxa, sea: dict[int, bool]) -> tuple[pl.DataFrame, dict[int, Labels]]:
    """`labels` merged down to one realm, land only with land and sea only with sea until one of each is left, and the map at every count.

    Each step merges the pair of the same kind with the lowest beta-sim on the union of their taxa (ties to the lowest positions, as `_Merger`);
    the last merge joins land and sea. A row (`MERGES`) has the count before the merge, its height, the two labels and `kind`: "land", "sea" or "both".
    """
    merger = _Merger(labels, ct)
    kind = np.array([sea[u] for u in merger.units])
    current = dict(labels)
    maps = {int(merger.active.sum()): dict(current)}
    rows = []
    upper = np.triu(np.ones((len(merger.units),) * 2, dtype=bool), k=1)
    while merger.active.sum() > 1:
        beta = merger.beta_sim()
        allowed = np.outer(merger.active, merger.active) & upper
        same = allowed & np.equal.outer(kind, kind)
        pool = same if same.any() else allowed
        keep, drop = np.unravel_index(np.argmin(np.where(pool, beta, np.inf)), beta.shape)
        keep, drop = int(keep), int(drop)
        joined = "both" if kind[keep] != kind[drop] else ("sea" if kind[keep] else "land")
        k = int(merger.active.sum())
        rows.append({"k": k, "height": float(beta[keep, drop]), "keep": merger.units[keep], "drop": merger.units[drop], "kind": joined})
        current = _relabel(current, merger.units[keep], merger.units[drop])
        merger.merge(keep, drop)
        maps[k - 1] = dict(current)
    return pl.DataFrame(rows, schema=MERGES), maps


def at_cuts(merges: pl.DataFrame, maps: dict[int, Labels], cuts: list[float]) -> dict[float, Labels]:
    """The map of each cut: the realms before the first merge of `merges` at a height of at least the cut, or joining land and sea, which a cut never
    does (that join is the construction's, and its height says nothing). The merges are in the order made."""
    found = {}
    for cut in cuts:
        first = merges.filter((pl.col("height") >= cut) | (pl.col("kind") == "both")).head(1)
        k = int(first["k"].item()) if first.height else min(maps)
        found[cut] = maps[k]
    return found


def levels_of(units: Labels, ct: CellTaxa, areas: dict[str, float], leaf_cut: float, cuts: list[float]) -> tuple[dict[float, Labels], pl.DataFrame]:
    """The map of every level: the leaf realms at `leaf_cut`, and `at_cuts` of their nested merges at each of `cuts` (each above `leaf_cut`).
    Also the merges: the leaf merges, then the nested ones (`MERGES`, with a `stage` column)."""
    if any(c <= leaf_cut for c in cuts):
        raise ValueError(f"the cuts of the levels above the leaves must be above the leaf cut {leaf_cut}: {cuts}")
    leaves, leaf_merges = leaf_realms(units, ct, leaf_cut)
    merges, maps = nested_merges(leaves, ct, sea_realms(leaves, areas))
    levels = {leaf_cut: leaves} | at_cuts(merges, maps, cuts)
    table = pl.concat([leaf_merges.with_columns(pl.lit("leaf").alias("stage")), merges.with_columns(pl.lit("nested").alias("stage"))])
    return levels, table


def realm_recovery(reference: Labels, replicate: Labels, areas: dict[str, float]) -> dict[int, float]:
    """For each realm of `reference`, the highest Jaccard index by area with a realm of `replicate`, on the cells both maps have."""
    common = [c for c in reference if c in replicate]
    inside: dict[tuple[int, int], float] = {}
    ref_area: dict[int, float] = {}
    rep_area: dict[int, float] = {}
    for c in common:
        a, b, w = reference[c], replicate[c], areas[c]
        inside[(a, b)] = inside.get((a, b), 0.0) + w
        ref_area[a] = ref_area.get(a, 0.0) + w
        rep_area[b] = rep_area.get(b, 0.0) + w
    best = {r: 0.0 for r in set(reference.values())}
    for (a, b), w in inside.items():
        best[a] = max(best[a], w / (ref_area[a] + rep_area[b] - w))
    return best


def _boot(shared: tuple, task: tuple[int, np.ndarray]) -> list[dict]:
    from realms.candidates import resolve_pieces
    from realms.cluster import run_infomap
    from realms.network import network_from

    ct, areas, adjacency, reference, leaf_cut, cuts, options, seed, trials = shared
    b, draws = task
    resampled = ct.resample(draws)
    partition = run_infomap(network_from(resampled), seed=seed + b, trials=trials, options=options)
    units = resolve_pieces(partition.leaf_labels(), adjacency, resampled)
    levels, _ = levels_of(units, resampled, areas, leaf_cut, cuts)
    return [{"boot": b, "cut": cut, "realm": realm, "jaccard": j}
            for cut, labels in levels.items() for realm, j in realm_recovery(reference[cut], labels, areas).items()]


def bootstrap_recovery(
    ct: CellTaxa,
    areas: dict[str, float],
    adjacency: dict[str, set[str]],
    reference: dict[float, Labels],
    leaf_cut: float,
    cuts: list[float],
    options,
    n_boot: int = 20,
    seed: int = 1,
    trials: int = 3,
    workers: int = 1,
) -> pl.DataFrame:
    """How well each realm of each level of `reference` comes back when the whole path (Infomap, `resolve_pieces`, `levels_of`) runs on resampled species.

    The draws are those of `select.bootstrap_stability` (species with replacement, all drawn up front from one generator at `seed`, Infomap at `seed`
    plus the bootstrap's number), so the result does not depend on `workers`. A row: `boot`, `cut`, `realm` and its `jaccard` (`realm_recovery`).
    """
    from realms.parallel import map_tasks

    rng = np.random.default_rng(seed)
    n_species = len(ct.species_ids())
    if not n_species:
        raise ValueError("bootstrap_recovery resamples species, so the matrix needs the species rank")
    draws = [rng.integers(n_species, size=n_species) for _ in range(n_boot)]
    shared = (ct, areas, adjacency, reference, leaf_cut, cuts, options, seed, trials)
    rows = [row for rows in map_tasks(_boot, list(enumerate(draws)), shared, workers) for row in rows]
    return pl.DataFrame(rows, schema={"boot": pl.Int64, "cut": pl.Float64, "realm": pl.Int64, "jaccard": pl.Float64})


# ---------------------------------------------------------------------------------------------------------------------------- the cut rules
# (DECISIONS, "Nested levels: the cut rules", 2026-11-06): which levels are robust to the sampling of species (rule 2), and whether each realm of a
# level holds more endemics than a random region of its area and pieces (rule 6).


def replay_heights(start: Labels, merges: pl.DataFrame, ct: CellTaxa) -> list[float]:
    """The height of each of `merges` (in the order made, from the map `start`) when the same groups merge in the same order on `ct`.

    A group with no cell in `ct` has no taxa, and its beta-sim with anything is 1."""
    merger = _Merger(start, ct)
    index = {u: i for i, u in enumerate(merger.units)}
    heights = []
    for keep, drop in merges.select("keep", "drop").iter_rows():
        i, j = index[keep], index[drop]
        heights.append(float(merger.beta_sim()[i, j]))
        merger.merge(i, j)
    return heights


def _replay(shared: tuple, task: tuple[int, np.ndarray]) -> list[float]:
    ct, start, merges = shared
    _, draws = task
    return replay_heights(start, merges, ct.resample(draws))


def height_replicates(ct: CellTaxa, start: Labels, merges: pl.DataFrame, n_boot: int = 199, seed: int = 1, workers: int = 1) -> np.ndarray:
    """`replay_heights` on `n_boot` resamples of the species (drawn up front from one generator at `seed`, as `bootstrap_recovery` draws them):
    resamples x merges."""
    from realms.parallel import map_tasks

    rng = np.random.default_rng(seed)
    n_species = len(ct.species_ids())
    if not n_species:
        raise ValueError("height_replicates resamples species, so the matrix needs the species rank")
    draws = [rng.integers(n_species, size=n_species) for _ in range(n_boot)]
    rows = map_tasks(_replay, list(enumerate(draws)), (ct, start, merges), workers)
    return np.array(rows, dtype=np.float64).reshape(n_boot, merges.height)


def robust_levels(merges: pl.DataFrame, replicates: np.ndarray, share: float = 0.95) -> pl.DataFrame:
    """Every place in the merge order where a cut could fall, and whether it is robust (rule 2).

    The level after the first `done` merges (`merges` in the order made, the land-sea join excluded) holds where every merge made has a lower height
    than every merge to come. It is robust where that holds in at least `share` of `replicates` (resamples x merges, `height_replicates`); its cut
    is the middle of the reference's gap: half way between the highest merge made and the lowest to come. A row: `done`, `realms` (the count after
    them), `below` and `above` (the gap's two ends), `cut`, `separated` (the share of replicates where a cut fits) and `robust`."""
    within = merges.filter(pl.col("kind") != "both")
    h = within["height"].to_numpy()
    reps = replicates[:, : len(h)]
    start = within["k"][0] if within.height else 1
    rows = []
    for done in range(1, len(h)):
        below, above = float(h[:done].max()), float(h[done:].min())
        separated = float(np.mean(reps[:, :done].max(axis=1) < reps[:, done:].min(axis=1))) if len(reps) else float("nan")
        rows.append({"done": done, "realms": int(start) - done, "below": below, "above": above,
                     "cut": (below + above) / 2 if above > below else None, "separated": separated, "robust": above > below and separated >= share})
    schema = {"done": pl.Int64, "realms": pl.Int64, "below": pl.Float64, "above": pl.Float64, "cut": pl.Float64, "separated": pl.Float64, "robust": pl.Boolean}
    return pl.DataFrame(rows, schema=schema)


def _rank_incidences(ct: CellTaxa, cells: list[str]) -> list:
    everything = ct.rows(cells)
    from realms.matrix import Incidence

    return [Incidence(everything.x[:, np.flatnonzero(ct.rank == r)], everything.totals[np.flatnonzero(ct.rank == r)]) for r in np.unique(ct.rank)]


def realm_excess(
    labels: Labels, ct: CellTaxa, adjacency: dict[str, set[str]], areas: dict[str, float], n_null: int = 99, seed: int = 1, tolerance: float = 0.95
) -> pl.DataFrame:
    """Each realm's endemic share against the same realm in `n_null` random maps of the same realm areas and pieces (`nullmodel.null_partitions`), the
    shares averaged equally over the ranks of `ct`, as `nullmodel.endemism_excess` averages them (rule 6).

    A taxon is endemic to a realm when at least `tolerance` of its cells lie inside. A row: `realm`, `observed`, `null_mean`, `null_q95`, `excess`
    (observed less the null mean) and `p`, the share of the null maps (counting the observed one) whose realm has at least the observed share."""
    from realms.endemism import realm_share
    from realms.nullmodel import null_partitions

    cells = sorted(labels)
    realms = sorted(set(labels.values()))
    position = {r: i for i, r in enumerate(realms)}
    tables = _rank_incidences(ct, cells)

    def shares(map_: Labels) -> np.ndarray:
        cell_realm = np.array([position[map_[c]] for c in cells])
        return np.mean([realm_share(t, cell_realm, len(realms), tolerance) for t in tables], axis=0)

    observed = shares(labels)
    nulls = np.array([shares(m) for m in null_partitions(labels, adjacency, areas, n=n_null, seed=seed)])
    p = (1 + (nulls >= observed).sum(axis=0)) / (n_null + 1)
    return pl.DataFrame({"realm": realms, "observed": observed, "null_mean": nulls.mean(axis=0), "null_q95": np.quantile(nulls, 0.95, axis=0),
                         "excess": observed - nulls.mean(axis=0), "p": p},
                        schema={"realm": pl.Int64, "observed": pl.Float64, "null_mean": pl.Float64, "null_q95": pl.Float64, "excess": pl.Float64, "p": pl.Float64})


# ---------------------------------------------------------------------------------------------------------------------------- which groups carry the signal
# (DECISIONS, "Which groups carry the realm signal", 2026-11-07): a class whose species are no more endemic to the realms than to random regions of the
# same areas only adds noise. Species in one cell are endemic to any map, so they are not scored, and a class too small to test is never judged.


def _endemic_species(x, totals: np.ndarray, cell_realm: np.ndarray, n_realms: int, tolerance: float) -> np.ndarray:
    """Whether each column of `x` (cells x species, binary) has at least `tolerance` of its `totals` cells in one realm."""
    from scipy.sparse import coo_matrix

    membership = coo_matrix((np.ones(len(cell_realm)), (np.arange(len(cell_realm)), cell_realm)), shape=(x.shape[0], n_realms)).tocsr()
    counts = (x.T.tocsr() @ membership).toarray()
    return counts.max(axis=1) / np.maximum(totals, 1) >= tolerance


def class_signal(
    labels: Labels, ct: CellTaxa, adjacency: dict[str, set[str]], areas: dict[str, float], species_class: dict[str, str],
    n_null: int = 99, seed: int = 1, tolerance: float = 0.95, min_species: int = 30, alpha: float = 0.05,
) -> pl.DataFrame:
    """For each class, the share of its species endemic to a realm of `labels` (at least `tolerance` of its cells inside), against the same share in
    `n_null` random maps of the same realm areas and pieces (`nullmodel.null_partitions`).

    Only species in at least two cells of `ct` are scored (`scored`); `species` counts them all. A class is `testable` with at least `min_species` scored
    species, and is to be dropped (`drop`) only if testable and its share is not above the null maps' (p > `alpha`), p counting the observed map."""
    from realms.matrix import SPECIES
    from realms.nullmodel import null_partitions

    cells = sorted(labels)
    realms = sorted(set(labels.values()))
    position = {r: i for i, r in enumerate(realms)}
    table = ct.rows(cells)
    columns = np.flatnonzero(ct.rank == SPECIES)
    ids = [ct.taxa[j].removeprefix("species:") for j in columns]
    x = table.x[:, columns].tocsc()
    totals = table.totals[columns]
    classes = np.array([species_class.get(i) or "" for i in ids])
    scored = totals >= 2

    def endemic(map_: Labels) -> np.ndarray:
        return _endemic_species(x, totals, np.array([position[map_[c]] for c in cells]), len(realms), tolerance)

    observed = endemic(labels)
    nulls = [endemic(m) for m in null_partitions(labels, adjacency, areas, n=n_null, seed=seed)]
    rows = []
    for name in sorted(set(classes) - {""}):
        mine = classes == name
        use = mine & scored
        n = int(use.sum())
        if n:
            obs = float(observed[use].mean())
            null = np.array([float(e[use].mean()) for e in nulls])
            p = float((1 + (null >= obs).sum()) / (n_null + 1))
            null_mean, null_q95 = float(null.mean()), float(np.quantile(null, 0.95))
        else:
            obs = null_mean = null_q95 = p = float("nan")
        testable = n >= min_species
        rows.append({"class": name, "species": int(mine.sum()), "scored": n, "observed": obs, "null_mean": null_mean, "null_q95": null_q95,
                     "excess": obs - null_mean, "p": p, "testable": testable, "drop": bool(testable and p > alpha)})
    schema = {"class": pl.String, "species": pl.Int64, "scored": pl.Int64, "observed": pl.Float64, "null_mean": pl.Float64, "null_q95": pl.Float64,
              "excess": pl.Float64, "p": pl.Float64, "testable": pl.Boolean, "drop": pl.Boolean}
    return pl.DataFrame(rows, schema=schema)
