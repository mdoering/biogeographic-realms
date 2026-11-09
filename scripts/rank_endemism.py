"""Rule 6 split by rank: each realm's endemism at species, genus and family rank on their own, at every level of `scripts/nested_cuts.py`.

`uv run python -I scripts/rank_endemism.py --counts data/gbif/<key> [--run native_only-eukaryotes] [--n-null 99]`

Rule 6 (`nested.realm_excess`) averages the three ranks equally, and so does the endemism excess of the configuration choice. This asks whether each
rank on its own is above the null, as the proposed half B pass mark C1 would require (PREREGISTRATION.md, section 9). It changes no method.

The levels are rebuilt as `scripts/nested_cuts.py` builds them (its `cuts.json`), and the null maps are drawn as `nested.realm_excess` draws them
(`nullmodel.null_partitions`, seed 1), so the rank average here reproduces `nested-cuts/excess.csv`, which is checked before anything is written.
A taxon is endemic to a realm when at least 95% of its units lie inside. Writes to `res<r>/<run>/rank-endemism/` (or `--out`), with `run.log`:

- `realms.csv`: per level, realm and rank, the taxa present and endemic, the endemic share against the same realm in the null maps (mean, q95,
  excess, p).
- `levels.csv`: per level and rank, the area-weighted share (the measure of `nullmodel.endemism_excess` at one rank) against the null maps, and the
  realms with p above `--alpha` at that rank.
- `level<n>.png`, one per reported level: a map per rank, each realm shaded by its endemic share at that rank, on a scale per rank; a realm not
  above the null at that rank is hatched. Drawn from the realm outlines of `nested/maps/level<i>.gpkg` (`scripts/nested_maps.py`).
- `report.md`: the tables per reported level.

Names come from `<run>/realm-names.json` if it is there (realm count -> realm -> name), else realms are numbered only. Half A only.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import shapely  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.patheffects import withStroke  # noqa: E402
from matplotlib.ticker import PercentFormatter  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from shapely.ops import polylabel  # noqa: E402

from realms.candidates import resolve_pieces  # noqa: E402
from realms.endemism import area_weighted, realm_share  # noqa: E402
from realms.explore import assigned_only, network_matrix, peak_gb, read_selected, read_units  # noqa: E402
from realms.hexmaps import DPI, GLOBE, OUTLINE, _polygons, _projected  # noqa: E402
from realms.maps import SURFACE, TEXT  # noqa: E402
from realms.matrix import RANKS, Incidence  # noqa: E402
from realms.nested import leaf_realms, nested_merges, sea_realms  # noqa: E402
from realms.nullmodel import null_partitions  # noqa: E402
from realms.occurrences import IncompleteIngest, resolution_of  # noqa: E402
from realms.review import markdown  # noqa: E402

log = logging.getLogger("rank_endemism")
TOLERANCE = 0.95
# the default sequential ramp, blue, light to dark (dataviz palette, steps 100 to 700)
BLUES = LinearSegmentedColormap.from_list("blues", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])


def counts(table: Incidence, cell_realm: np.ndarray, n_realms: int) -> tuple[np.ndarray, np.ndarray]:
    """Taxa present in and endemic to each realm, as `endemism.realm_share` counts them."""
    membership = coo_matrix((np.ones(len(cell_realm), dtype=np.int32), (np.arange(len(cell_realm)), cell_realm)),
                            shape=(table.x.shape[0], n_realms)).tocsr()
    inside = (table.xt @ membership).tocoo()
    present = np.bincount(inside.col, minlength=n_realms)
    endemic = np.bincount(inside.col[inside.data / table.totals[inside.row] >= TOLERANCE], minlength=n_realms)
    return present, endemic


def _percent(share: float) -> str:
    return f"{share:.0%}" if share >= 0.1 else f"{share:.1%}"


def _place(ax, x: float, y: float, text: str, placed: list) -> None:
    """The text at (x, y), moved below or above by one or two of its heights if it would overlap a label placed so far; left where it is if all do."""
    label = ax.text(x, y, text, ha="center", va="center", fontsize=7.5, color=TEXT, zorder=5, linespacing=1.1,
                    path_effects=[withStroke(linewidth=2.5, foreground=SURFACE)])
    renderer = ax.figure.canvas.get_renderer()
    box = label.get_window_extent(renderer=renderer).expanded(1.05, 1.1)
    for shift in (0.0, -box.height, box.height, -2 * box.height, 2 * box.height):
        moved = box.translated(0, shift)
        if not any(moved.overlaps(other) for other in placed):
            break
    else:
        moved, shift = box, 0.0
    if shift:
        label.set_position(ax.transData.inverted().transform((moved.x0 + moved.width / 2, moved.y0 + moved.height / 2)))
    placed.append(moved)


def plot_level(gpkg: Path, rows: pl.DataFrame, names: dict[str, str], sea: set[int], alpha: float, out: Path, title: str) -> None:
    """One map per rank, each realm shaded by its endemic share on a scale of that rank's own (0 to the largest share of a land realm, to the
    next 10%), since family shares are a fraction of species ones. A sea realm (`sea`) above that is drawn at the top of the scale and says so. The label gives the share, the null mean and endemic / present taxa; a realm not above the
    null at that rank is hatched and gives its p."""
    realms_geo = gpd.read_file(gpkg, layer="realms")
    if set(realms_geo["realm"]) != set(rows["realm"]):
        raise SystemExit(f"FAILED {gpkg} has realms {sorted(realms_geo['realm'])}, the level has {sorted(set(rows['realm']))}")
    shapes = {int(r): _projected(g) for r, g in zip(realms_geo["realm"], realms_geo.geometry, strict=True)}
    globe = _polygons(_projected(shapely.box(-180, -90, 180, 90)))[0]
    fig, axes = plt.subplots(len(RANKS), 1, figsize=(13, 7.2 * len(RANKS)), facecolor=SURFACE, layout="constrained")
    for ax, rank in zip(axes, RANKS, strict=True):
        ax.fill(*globe.exterior.xy, facecolor=GLOBE, edgecolor=OUTLINE, linewidth=0.5, zorder=0)
        x0, y0, x1, y1 = globe.bounds
        ax.set(xlim=(x0, x1), ylim=(y0, y1), aspect="equal")
        ax.axis("off")
        mine = rows.filter(pl.col("rank") == rank).sort("observed", descending=True)
        land = mine.filter(~pl.col("realm").is_in(list(sea)))
        norm = Normalize(0, float(np.ceil((land if land.height else mine)["observed"].max() * 10) / 10), clip=True)
        for realm, share in mine.select("realm", "observed").iter_rows():
            gpd.GeoSeries([shapes[realm]]).plot(ax=ax, color=BLUES(norm(share)), edgecolor=SURFACE, linewidth=0.5, aspect=None, zorder=1)
        failing = mine.filter(pl.col("p") > alpha)
        if failing.height:
            with plt.rc_context({"hatch.linewidth": 0.6}):
                gpd.GeoSeries([shapes[r] for r in failing["realm"]]).plot(ax=ax, facecolor="none", edgecolor="#0b0b0bb0", hatch="////", linewidth=0,
                                                                          aspect=None, zorder=2)
        ax.apply_aspect()
        placed: list = []
        for realm, share, null_mean, p, endemic, present in mine.select("realm", "observed", "null_mean", "p", "endemic", "present").iter_rows():
            biggest = max(shapely.get_parts(shapes[realm]), key=lambda s: s.area)
            spot = polylabel(biggest, tolerance=5e4)
            text = f"{names.get(str(realm), realm)}\n{_percent(share)} (null {_percent(null_mean)}), {endemic:,}/{present:,}"
            if share > norm.vmax:
                text += "\nabove the colour scale"
            _place(ax, spot.x, spot.y, text + ("" if p <= alpha else f"\nnot above the null, p = {p:.2f}"), placed)
        ax.set_title(f"{rank}: {mine.height - failing.height} of {mine.height} realms above the null (p ≤ {alpha})", loc="left", color=TEXT, fontsize=12)
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=BLUES), ax=ax, orientation="vertical", fraction=0.02, pad=0.01, shrink=0.7)
        bar.set_label(f"endemic share of {rank} taxa", color=TEXT)
        bar.ax.yaxis.set_major_formatter(PercentFormatter(1.0))
        bar.outline.set_visible(False)
    fig.suptitle(title + "\nA taxon is endemic with at least 95% of its units in the realm. Label: endemic share (null mean), endemic / present taxa."
                 "\nEach rank has its own colour scale, set by the land realms.", x=0.01, ha="left", color=TEXT, fontsize=12)
    fig.savefig(out, dpi=DPI // 2, facecolor=SURFACE)
    plt.close(fig)


def score(args, run_dir: Path, out: Path, cuts: dict, rule6: pl.DataFrame, names: dict, started: float) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The tables, written to `out`: each level rebuilt, its rank shares against the null maps, checked against rule 6."""
    found = read_selected(run_dir, cuts["config"], args.configs)
    units = read_units(run_dir)
    ct = network_matrix(run_dir, units, found.merged)
    scored = assigned_only(found.partition, ct, units)
    del ct
    leaves = resolve_pieces(found.partition.leaf_labels(), scored.adjacency, scored.ct)
    base, _ = leaf_realms(leaves, scored.ct, cuts["base_cut"])
    _, maps = nested_merges(base, scored.ct, sea_realms(base, scored.areas))
    ranks = [RANKS[r] for r in np.unique(scored.ct.rank)]
    if ranks != list(RANKS):
        raise SystemExit(f"FAILED the matrix has ranks {ranks}, not {list(RANKS)}")
    log.info(f"matrix: {len(scored.ct.cells)} units; " + ", ".join(f"{r} {int((scored.ct.rank == i).sum())}" for i, r in enumerate(RANKS)) + f"; {time.perf_counter() - started:.0f} s")

    realm_rows, level_rows = [], []
    for level in cuts["levels"]:
        n = level["realms"]
        labels = maps[n]
        cells = sorted(labels)
        realms = sorted(set(labels.values()))
        position = {r: i for i, r in enumerate(realms)}
        everything = scored.ct.rows(cells)
        tables = [Incidence(everything.x[:, cols], everything.totals[cols]) for cols in (np.flatnonzero(scored.ct.rank == i) for i in range(len(RANKS)))]
        cell_area = np.array([scored.areas[c] for c in cells])

        def placed(map_) -> np.ndarray:
            return np.array([position[map_[c]] for c in cells])

        observed_realm = placed(labels)
        observed = np.array([realm_share(t, observed_realm, len(realms), TOLERANCE) for t in tables])  # ranks x realms
        null_realms = [placed(m) for m in null_partitions(labels, scored.adjacency, scored.areas, n=args.n_null, seed=1)]
        nulls = np.array([[realm_share(t, cr, len(realms), TOLERANCE) for t in tables] for cr in null_realms])  # maps x ranks x realms

        # the rank average is rule 6: it has to reproduce nested-cuts/excess.csv
        mine = rule6.filter(pl.col("realms") == n).sort("realm")
        averaged = observed.mean(axis=0)
        p_avg = (1 + (nulls.mean(axis=1) >= averaged).sum(axis=0)) / (args.n_null + 1)
        if mine["realm"].to_list() != realms or not np.allclose(mine["observed"].to_numpy(), averaged, rtol=0, atol=1e-12) \
                or not np.allclose(mine["p"].to_numpy(), p_avg, rtol=0, atol=1e-12):
            raise SystemExit(f"FAILED at {n} realms the rank average does not reproduce nested-cuts/excess.csv")

        level_names = names.get(str(n), {})
        for i, rank in enumerate(RANKS):
            present, endemic = counts(tables[i], observed_realm, len(realms))
            p = (1 + (nulls[:, i, :] >= observed[i]).sum(axis=0)) / (args.n_null + 1)
            for j, realm in enumerate(realms):
                realm_rows.append({"realms": n, "cut": level["cut"], "realm": realm, "name": level_names.get(str(realm), ""), "rank": rank,
                                   "present": int(present[j]), "endemic": int(endemic[j]), "observed": float(observed[i, j]),
                                   "null_mean": float(nulls[:, i, j].mean()), "null_q95": float(np.quantile(nulls[:, i, j], 0.95)),
                                   "excess": float(observed[i, j] - nulls[:, i, j].mean()), "p": float(p[j])})
            whole = area_weighted(observed[i], observed_realm, cell_area)
            whole_null = np.array([area_weighted(nulls[b, i], cr, cell_area) for b, cr in enumerate(null_realms)])
            failing = [level_names.get(str(realms[j]), str(realms[j])) for j in np.flatnonzero(p > args.alpha)]
            level_rows.append({"realms": n, "cut": level["cut"], "rank": rank, "observed": whole, "null_mean": float(whole_null.mean()),
                               "null_q95": float(np.quantile(whole_null, 0.95)), "excess": float(whole - whole_null.mean()),
                               "p": float((1 + (whole_null >= whole).sum()) / (args.n_null + 1)), "realms_failing": "; ".join(failing)})
        log.info(f"{n} realms: rank average reproduces rule 6; failing p > {args.alpha}: "
                 + "; ".join(f"{r['rank']} [{r['realms_failing']}]" for r in level_rows[-len(RANKS):]) + f"; {time.perf_counter() - started:.0f} s")

    realm_table = pl.DataFrame(realm_rows)
    level_table = pl.DataFrame(level_rows)
    realm_table.write_csv(out / "realms.csv")
    level_table.write_csv(out / "levels.csv")
    with pl.Config(tbl_rows=200, float_precision=3, tbl_width_chars=220):
        log.info(f"per level and rank:\n{level_table}")
    return realm_table, level_table


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--run", default="native_only")
    parser.add_argument("--n-null", type=int, default=99)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--configs", type=Path, default=Path("explore/configs.toml"))
    parser.add_argument("--results", type=Path, default=Path("results/explore"))
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--maps-only", action="store_true", help="redraw the maps and the report from realms.csv and levels.csv already written")
    args = parser.parse_args()

    try:
        resolution = resolution_of(args.counts)
    except (IncompleteIngest, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    run_dir = args.results / f"res{resolution}" / args.run
    out = args.out or run_dir / "rank-endemism"
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
                        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out / "run.log", mode="a" if args.maps_only else "w")])
    cuts = json.loads((run_dir / "nested-cuts/cuts.json").read_text())
    rule6 = pl.read_csv(run_dir / "nested-cuts/excess.csv")
    names_file = run_dir / "realm-names.json"
    names = json.loads(names_file.read_text()) if names_file.exists() else {}
    log.info(f"---- rank endemism: resolution {resolution}, run {args.run}, configuration {cuts['config']}, base cut {cuts['base_cut']}, "
             f"levels {[lv['realms'] for lv in cuts['levels']]}, reported {cuts['reported']}, {args.n_null} null maps")
    started = time.perf_counter()

    if args.maps_only:
        realm_table, level_table = pl.read_csv(out / "realms.csv"), pl.read_csv(out / "levels.csv").with_columns(pl.col("realms_failing").fill_null(""))
    else:
        realm_table, level_table = score(args, run_dir, out, cuts, rule6, names, started)
    index_of = {}
    nested_realms = pl.read_csv(run_dir / "nested/realms.csv")
    for idx, count in nested_realms.group_by("level").len().iter_rows():
        index_of[int(count)] = int(idx)
    lines = [f"## Endemism by rank, {args.run} {cuts['config']} (resolution {resolution}, COL 315557)\n",
             "Rule 6 on each rank on its own. A taxon is endemic to a realm when at least 95% of its units lie inside; the share is the endemic taxa of "
             f"those present. The null is {args.n_null} random maps of the same realm areas and pieces (seed 1, as rule 6). The rank average of these "
             "shares is rule 6 and reproduces `nested-cuts/excess.csv`.\n"]
    for n in cuts["reported"]:
        rows = realm_table.filter(pl.col("realms") == n)
        cut = rows["cut"][0]
        lines.append(f"\n### {n} realms (cut {cut:.3f})\n")
        lines.append(markdown(level_table.filter(pl.col("realms") == n).select("rank", "observed", "null_mean", "null_q95", "excess", "p", "realms_failing"),
                              {"observed": ".3f", "null_mean": ".3f", "null_q95": ".3f", "excess": ".3f", "p": ".2f"}))
        wide = (rows.with_columns((pl.col("observed") * 100).round(1).cast(pl.String) + " (" + (pl.col("null_mean") * 100).round(1).cast(pl.String) + ") "
                                  + pl.col("endemic").cast(pl.String) + "/" + pl.col("present").cast(pl.String)
                                  + pl.when(pl.col("p") > args.alpha).then(pl.lit(" p=") + pl.col("p").round(2).cast(pl.String)).otherwise(pl.lit("")))
                .pivot(on="rank", index=["realm", "name"], values="observed").sort("realm"))
        lines.append("\nPer realm: endemic share % (null mean %), endemic / present taxa; p shown where above " + f"{args.alpha}.\n")
        lines.append(markdown(wide))
        if n in index_of:
            gpkg = run_dir / f"nested/maps/level{index_of[n]}.gpkg"
            png = out / f"level{index_of[n]}.png"
            sea = set(nested_realms.filter((pl.col("level") == index_of[n]) & (pl.col("sea_share") >= 0.5))["realm"].to_list())
            plot_level(gpkg, rows, names.get(str(n), {}), sea, args.alpha, png,
                       f"{args.run}, {cuts['config']}: {n} realms, cut {cut:.3f}, endemic share by rank, resolution {resolution}, COL 315557")
            lines.append(f"\n![{n} realms]({png.name})\n")
            log.info(f"wrote {png}")
    (out / "report.md").write_text("\n".join(lines) + "\n")
    log.info(f"rank endemism: {time.perf_counter() - started:.0f} s, peak {peak_gb():.1f} GB")


if __name__ == "__main__":
    main()
