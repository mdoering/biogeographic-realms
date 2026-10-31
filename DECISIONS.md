# Decision log

Every change to a method or threshold, dated, with its reason. The design is in [docs/design.md](docs/design.md).

## 2026-10-31

- Kingdom clash narrowed to Animalia, Plantae and Fungi. Verbatim names like `Metazoa` and real disagreements like Chromista
  against Plantae aren't wrong matches, whereas a plant record matched to an animal homonym is.
- Family split seed fixed at `20261031`, before anything was clustered.
- Weight of a genus or family link is `1 - log k / log K` (k cells occupied, K cells in total), checked against the formula in
  Infomap Bioregions 2. Species links keep weight 1, as in that paper.

### Reproducing Castro-Insua et al. 2018 on Costello's matrix

- Reproduced with their method: 1,759 cells with at least 5 species, Simpson dissimilarity, R `hclust` with `ward.D2`. ANOSIM R
  matches their Figure 1 within 0.015 at all 12 points read from it, including the peak at k=7 (0.442), and the k=7 map matches
  their Figure 2a realm for realm. One of their statements doesn't hold exactly: Ward doesn't beat average for every k above 6,
  as they tie at k=16, which their own Figure 1 also shows.
- **The published map isn't unique.** 55.5% of cell pairs have a dissimilarity of exactly 1, so Ward's tree depends on how ties
  break. On the same distances, SciPy agrees with R at a V-measure of only 0.755. Shuffling the order of the cells gives 0.583
  to 0.768 against R's map, and ANOSIM R from 0.318 to 0.433 against R's 0.442. Any single Ward tree on these data is one draw
  from many.
- Consequence: Ward on beta-sim can't serve as a single-run cross-check or fallback. During exploration, try a consensus over
  shuffled cell orders, and report the spread between orders either way.
- `hierarchical_labels` uses `cut_tree`, as `fcluster` returns fewer than k groups when merge heights tie.

### Trying Infomap and the selection rules on Costello's matrix

- **Infomap degenerates at its default scale.** Species and genus nodes give 16 modules, the largest holding 75% of the 1,759 cells.
  After the size floor, 5 units remain, so the merge sequence has a single candidate. It selected k=5 and fell short of the
  stability threshold (0.786 against 0.8). The design's fallback doesn't trigger, because 5 units aren't "fewer than 5".
- **The usable scale lies somewhere else.** A scan of Infomap settings on this matrix, as (leaf modules, share of the largest):
  species only, default (26, 0.57); species and genus, default (27, 0.73); two-level (16, 0.92 with genus); Markov time 0.75
  (372, 0.33 species only; 354, 0.36 with genus); Markov time 0.5 (over 1,600, 0.01). Nothing here was tried between 0.75 and 1.
- **Adding genus makes modules coarser.** Infomap Bioregions 2 warns that ancestral nodes collapse modular structure, and that
  shows here. Evolutionary depth may need a lower weight than the paper's whole-tree formula, or integration at one rank.
- To settle while exploring half A, using only endemism excess and stability: a Markov time sweep, or recursing Infomap inside
  large modules. Runtime and memory aren't a concern at this size, with about 100 s and 5 to 7 GB including 20 bootstraps.
- The area-weighted endemism score of a map dominated by one huge realm is high by construction (0.678 here), since nearly every
  species occurs in it and most occur nowhere else. Only the excess over the null, which has the same areas, is meaningful.
- Costello's own 30 realms, assigned from their published polygons by cell centre (1,584 of 1,759 cells), give a strict mean endemism
  of 19% on their own matrix (range 4% to 43%), against the published 42% (17% to 84%). This counts a species' presence in cells no polygon
  covers against its endemism. Their figure may come from cell groups before the boundaries were smoothed by hand, but that's unverified.
  So 42% can't serve as a reference value without its method.

### Text patterns

Terms were judged by reading samples of real records from GBIF's free text search (`scripts/test_patterns.py`), dropping or narrowing
any where over a third of the sample wasn't a garden, zoo or held organism.

- **Limited by kingdom.** Garden terms and `cultivated`-type terms count for plants only: animals in a botanic garden were overwhelmingly
  wild ("Santa Cruz Botanic Garden", "Jardin botanique de Montréal"), and "On cultivated Quercus" describes a wild insect on a host plant.
  Zoo terms count for animals only: plants at a zoo were mostly wild ("Boulevard du Jardin Zoologique", "along zoo fence").
- **Dropped.** `aquarium` ("East London, Aquarium Rocks" is a place, and the rest were rescued wild animals), `garden` ("Garden compost",
  "Rich garden soil"), `hortus` (`botanic` already catches "Hortus Botanicus", and the rest were places).
- **Negations excluded.** "nicht kultiviert", "no cultivada", "not a zoo animal", and "cultivated land" or "champ cultivé", which describe a
  weed's habitat.
- **Captivity in free text needs a statement of fact.** Remarks discuss it without meaning it ("do not mark this observation as captive",
  "not in captivity!"), so `captive` and `in captivity` only count in `locality` and `habitat`, and remarks need "captive born" or similar.
- `greenhouse`, `glasshouse` and `gewächshaus` count for everything, since an organism under glass says nothing about a native range.
- Accepted imprecision: hits remain imperfect ("Arboretum de Sibang" is natural forest, "Bairro Jardim Botânico" a neighbourhood). A
  species only loses its presence in a cell if every record of it there is hit, so a false hit costs little.

### COL release

- 2026-10-31: GBIF production is on 315557 (COL26.6 XR), per metadata and a probe: 19 sampled ids absent from it, GBIF knows none of them
  (`scripts/check_release.py`). Some of those names match a different id in 315557, so ids can't be joined across releases.

### Reference data

- Export flavours differ more than the design assumed. A plain COLDP export has only 6 columns and no distributions. A classification-only
  export is one flat `dataset-<key>.tsv` with ranks and `clb:taxGroup`. Only extended exports carry `col:environment`, `col:extinct`, flat ranks
  and `Distribution.tsv`. So COL uses an extended export plus the classification-only one for `taxGroup`, and WCVP an extended one.
- **Rank coverage in 315557 is good enough for the three-rank network:** of 2,529,340 accepted species, 99.0% have a genus and 97.2% a family.
  The weakest large kingdom is Fungi at 95.1%, and animals are at 96.9%. `environment` covers only 33.0%, as the design expected, and 6.9% are extinct.
  `taxGroup` covers 99.8% of species and suits the `group` of the family split and the per-group runs.
- **In WCVP an empty `col:establishmentMeans` means native.** Unflagged regions are right for potato (Andes), maize (Mexico and Guatemala), tomato
  (Peru), *Eucalyptus globulus* (Tasmania and Victoria), *Carpobrotus edulis* (the Cape) and a random 50 species read by eye. 239,132 regions are flagged
  `introduced` and 2,324 `uncertain`. 13,289 of 367,237 species have both native and introduced regions.
- **GRIIS.** 249 of the 304 lists only have a plain export, 32 an extended one, 23 have none. A plain species list is enough, as every species listed is
  alien in that country. Joining them to COL has to go through names: `col:ID`s are GBIF URLs (`https://www.gbif.org/species/8825970`) in the lists
  that have them. Lists also hold higher ranks, which need filtering out.
- **Missing, needing a POST on production to create:** WRiMS (26853), and GRIIS for AQ, AS, BL, BQ-Bonaire, BQ-Saba, FK, GQ, GS, GW, KP, KR, KY, KY-Grand,
  KY-Little, MF, PG, PS, SH-Ascension, SH-StHelena, SH-Tristan, SS, TR, VG. Nothing requested yet.

### Null model

- Pieces grow together, always extending the one furthest below its area, followed by a pass that hands border cells to neighbours below
  their area. The design's order of growing one piece at a time missed target areas by up to 44 of 144 cells. Now 92% of maps are within one
  cell, and the worst of 200 was 11 off. Null maps come out blobbier than real realms, which could understate the null's endemism slightly.
- A cell nothing touches joins the piece furthest below its area, not the nearest by centre, because the null works from adjacency alone.

### After review of the branch

- **The family split keeps families whole.** A lineage takes its most common `taxGroup`, as 20 real families carry two groups (Asteraceae among
  them) and 4,625 species have none. On the real table: no species without a half, no lineage in both halves.
- **Genus and family names are qualified by kingdom.** 2,272 genus names recur across kingdoms in 315557 (*Drosophila* in three), covering
  92,973 species. Unqualified, they would tie fly cells to fungus cells.
- **Endemism counts every cell a taxon occurs in, even outside the map.** Before, a species' cells outside the labelled set were ignored. This
  changed the Costello figure above from 23% to 19%.
- **Pieces and merges use every rank in the network, not species alone.** The design says "sharing the most species", written before evolutionary
  depth was added. Comparing species only would discard the genus and family signal the network was built to carry. `rank_presence` provides the
  combined table.
- **A separate piece stays if it shares nothing with any other realm,** as moving it would be arbitrary. The design only covers the case where
  there is a realm it shares something with.
- **`bootstrap_stability` reruns the caller's ranks, size floor and `k_max`.** It used to rerun a fixed pipeline whatever produced the maps.

Query changes, which the first download freezes:

- `n_open`, `n_datasets_open` and `n_years_open` count `license IN ('CC0_1_0','CC_BY_4_0')`, as the design's robustness run says, rather than
  anything not CC BY-NC. `n_nc` stays.
- A record of a species with no kingdom no longer drops out of the clean counts (1,753 species in 315557).
- `" zoo "` now matches at the start and end of a text, and `" zoo,"` and `" zoo."` are added, so "Berlin Zoo" and "Zoo, Berlin" count.

### Taxa by id, accepted only

Markus's call, replacing the kingdom-qualified names above. All figures are for 315557, from `scripts/check_reference.py` and the GBIF API.

- **Genus and family nodes are COL ids, found by walking `parentID`.** The kingdom qualifier only separated names across kingdoms. Within a
  kingdom, several ids hold species under 8,080 genus names (379,141 species, 28,942 of them outside the name's largest id) and 1,004 family
  names (1,215,886 species, 14,416 outside). The walk gives the same coverage as the flat columns (99.0% and 97.2%) and the same names for all
  2,529,340 species, and no accepted species sits under a genus or family that isn't accepted. `load_taxa` refuses a release where one does.
- **Ids also split duplicates, and many repeated names look like duplicates.** Curculionidae is `8HBQ5` (78,396 species) and `GD8GY` (158
  species, provisionally accepted, no authorship). Of the 1,040 family ids that aren't the largest of their name, 872 are provisionally accepted
  in the same order. Of the 8,996 such genus ids, 4,879 are provisionally accepted in the same family, while 672 are accepted in another family,
  the pattern of true homonyms. Some hold many records: `N7XLH` Eunotiaceae has 50,312 against 14,005 on the larger `VHBM4`. A duplicate becomes a
  narrow extra genus or family, which can inflate endemism at those ranks. **Open:** left as ids for now. Measure the effect on half A by
  rerunning with same-named usages in one lineage merged.
- **The family split stays by name.** There a duplicate on both sides would leak relatives across the hold-out, and a true homonym in one half
  costs nothing.
- **Species names repeat as well,** with 34,116 names on two or more accepted ids (69,363 usages). They stay separate ids. In a random 25 of those
  names GBIF holds 18,545 records on one id and 39 on the others, so little is at stake.
- **Only accepted species become nodes.** `build_network` and `presence_by_rank` drop any species missing from `load_taxa`, which keeps accepted
  and provisionally accepted species. A provisionally accepted species (174,878 of the 2,529,340, 0.29% of records directly) is a taxon and not a
  synonym.
- **Records matched through a synonym count for the accepted species.** That is 113.8M records, 3.1%, too many to lose. Of 3,000 sampled records
  matched to synonyms, ambiguous synonyms and misapplied names, 2,970 have a species in their COL classification, always accepted or provisionally
  accepted, with synonyms of a subspecies rolling up to its species. The other 30 match at genus level and fall out with `specieskey IS NOT NULL`.
  `n_synonym` allows a run without them.
- **Matches to ambiguous synonyms and misapplied names leave the clean count** (3.33M and 19k records, 0.09%), staying in `n_raw`. A pro parte name
  could be any of several taxa and a misapplied one was used for another taxon, so the species GBIF gives them is a guess.
- **To check on the pilot:** the sample reads the API's COL classification, not the SQL columns. Every `specieskey` must be in the accepted table,
  and `n_synonym` should be near 3.1% of records. Zero would mean the SQL `taxonomicstatus` isn't COL's.

### Scale (2026-10-31)

Speed-ups must not change a result, and each one is measured on the same synthetic grid: 8,100 cells, 1M species, 31.6M presences, 89M
links with genus and family. All figures are from this Mac (48 GB, 14 cores). The grid is `synth.grid_presence`, which draws in the order
of the scratch generator the baseline was measured with: the same 31,608,488 presences, and 89,059,784 links on 1,210,000 taxa.

Baseline, the string-keyed code before any speed-up:

| Stage | Time | Note |
|---|---|---|
| `build_network` | 37 s | peak RSS 5.0 -> 12.7 GB |
| `rank_presence` | 5 s | |
| `incidence`, one call | 22 s | peak 15.5 GB |
| `add_link` loop into Infomap | 70 s | one Python call per link |
| Infomap, 1 trial | 23 s | the synthetic network has no modules; real data will be slower |
| `resolve_pieces` (529 units) | 24 s | nearly all of it the `incidence` call |
| `_Merger` init / per merge | 26 s / 0.05 s | init is the `incidence` call |
| one null map (`_one_partition`) | 0.09 s | |
| `realm_share`, one call at 89M nonzeros | 0.29 s (5 realms) to 1.24 s (60) | scoring one null map at one rank |

One selection run (56 candidate maps, 99 null maps each, 3 ranks, 20 bootstraps) then costs about 1 h of `endemism_excess` rebuilding
incidence, 2 h of scoring null maps, and 1 h of bootstraps rebuilding network and tables from strings. Peak memory near 17 GB came from
string-keyed Polars frames and grows with the links, so resolution 7 (about 4x these links) would not fit in 48 GB.

**`CellTaxa`, one integer-coded matrix built once** (`scripts/benchmark.py --species 1000000`). Peak RSS is the whole process so far,
generator included.

| Stage | Time | Target | Peak RSS | Target |
|---|---|---|---|---|
| `grid_presence` | 10.8 s | | 2.1 GB | |
| `cell_taxa` | 3.4 s | 40 s | 3.6 GB | 8 GB |
| `network_from` | 0.6 s | 10 s | 5.8 GB | |
| `rows`, all 8,100 cells | 0.2 s | 2 s | 5.8 GB | |
| `resample`, one draw | 2.2 s | 15 s | 6.3 GB | |

Building the network took 37 s and now takes 4.0 s (`cell_taxa` and `network_from`), and what `incidence` took 22 s for, `rows` does in
0.2 s. The results are the same: `network_from` equals the old `build_network` bit for bit on planted, hand-made and grid cases, and
`scripts/costello_selection.py` still gives a byte-identical `selection.csv` and a `selection.log` that differs only in its stamps.

- Genus and family weights come from a table of `information_weight(k, n_cells)` for `k = 1..n_cells`, so they are the same floats as
  before and a column held by every cell still drops out of the network. `rows` and endemism keep such a column, as before.
- **`resample` drops a cell that no drawn species occupies,** as the old `_resample` did when the cell's rows fell out of the joined
  presence table. This matters beyond tidiness, because the weights count cells: keeping the empty cell would change them. On a
  400-cell grid with 300 species a draw lost 7 cells. The old function drew among every species in the presence table, the new one among
  the accepted species that have a column, which is the same set when `taxa` lists every species present, as it does for Costello.
- **`bootstrap_stability` takes its ranks from the matrix and needs the species rank.** The `ranks` argument is gone, so the network, the
  pieces and the merges of a bootstrap always use the ranks of the matrix that made the maps. It resamples species, so a matrix without them
  is refused, where `ranks=("genus",)` used to resample species and cluster on genera alone. Nothing in the pipeline runs that way.

**Bulk Infomap input, parallel bootstraps and null maps** (`scripts/benchmark.py --species 1000000`, same grid). Infomap itself is not run,
because the grid has no modules to find. Where a stage needs a partition, 529 blocks of about 4 x 4 cells stand in for the leaf modules (the
unit count of the baseline's `resolve_pieces` row), on cells of 100,000 km2 so that most clear the size floor. The maps are the 56 that `candidate_maps` builds from
them, k = 5 to 60.

| Stage | Time | Target | Was |
|---|---|---|---|
| Infomap input, `load_infomap` | 13.9 s | 10 s | 70 s |
| one bootstrap, Infomap replaced by the blocks | 44.6 s | 20 s | about 3 min |
| `realm_share`, 60 contiguous realms | 0.3 s | 0.5 s | 1.24 s |
| `realm_share`, 60 scattered realms (worst case) | 0.5 s | 0.5 s | |
| `Incidence.xt`, once per table | 0.4 s | | |
| `excess_table`, 56 maps x 99 null maps x 3 ranks, 8 workers | 238.8 s | 20 min | about 3 h |

`excess_table` ran in full, not on a subset, with a worker peak of 1.7 GB each and 11.7 GB in the parent. A run of 8 of the 56 maps took
46.9 s, which scales to 5.5 min, so a subset is a fair estimate here.

**Two targets are missed, and the time is in code this change doesn't touch.**

- **Infomap input** is 13.9 s (13.4 to 16.9 s across runs). Building the array takes 0.3 s. The other 13.6 s is `addLinksFromNumpy2D`
  inside Infomap, about 150 ns a link and growing with the link count (10M links take 0.9 s). Sorted links, as ours are, load faster than
  shuffled ones (10M take 0.9 s against 2.2 s). It is five times faster than the loop.
- **A bootstrap without Infomap** is 44.6 s, and 38 s of it is `merge_sequence`: 527 merges from 529 units down to 2, at 68 ms each, as
  every merge rebuilds the units x taxa matrix. `apply_size_floor` takes 2.7 s, resampling 2.3 s, `resolve_pieces` 0.8 s and the rest 0.6 s,
  so without `merge_sequence` the bootstrap is under 7 s. `merge_sequence` and `apply_size_floor` are untouched here, and the pilot will
  time them at full scale. The merges grow with the number of leaf modules Infomap finds, so real data may be better or worse than 529.

- **`realm_share` was already fast.** Against the same table the old product `x.T @ membership` takes 0.23 s with contiguous realms and 0.56 s
  with scattered ones, so most of the 1.24 s went when `rows` started giving int32 counts. The cached transpose helps where realms are
  scattered (the product goes from 0.56 s to 0.45 s) and not where they are contiguous, which is what real and null maps are (0.23 s
  before, 0.28 s now). The transpose costs 0.4 s once per table.
- **Results are the same.** `run_infomap` gives the loop-built reference partition, `workers=2` gives the same stability and excess as one
  worker, and `scripts/costello_selection.py --workers 4` still gives a byte-identical `selection.csv` and a `selection.log` that differs
  only in its stamps. Its stability step took 30 s on 4 workers against 85 s on one.
- **Draws are made up front,** all `n_boot` of them in order from the one generator, and bootstrap `b` seeds Infomap with `seed + b`, so
  which worker runs which bootstrap changes nothing. A null map comes from `seed` and its own map, as in `endemism_excess`.
- **Workers use one OpenMP thread,** unless `OMP_NUM_THREADS` is set. Infomap starts a thread per core and gains nothing from them: a
  three-trial run on Costello takes 3.9 s on 14 threads and 3.6 s on one, but burns 7 times the CPU. With 4 workers that crowds the
  processes out (8 bootstraps took 23 s, against 34 s serial), and with one thread each they take 14.7 s. The partition is the same on one
  thread. The parent process is left alone.
- **`--markov-time` is passed only when it isn't 1.0,** so the default command line is the one Infomap ran before. Below 1 gives more,
  smaller modules. `bootstrap_stability` takes the same value as the run that made the maps, as it does for `trials` and the size floor.
  `costello_selection.py` has no flag for it, because no run uses anything but the default yet.
