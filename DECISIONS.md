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
  KY-Little, MF, PG, PS, SH-Ascension, SH-StHelena, SH-Tristan, SS, TR, VG. Nothing requested yet. *(Since requested and finished, see "Reference joins" below.)*

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
  Settled on a Denmark download before the pilot's query was trusted, see "Column check download" below.

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

### Family split frozen (2026-10-31)

Written before anything is clustered on GBIF data, by `scripts/write_split.py` to `data/split/halves-315557.parquet` (columns `species_id`,
`lineage`, `half`, sorted by `species_id`; not committed). The freeze (M2) cites the hash below.

- **Inputs:** COL release 315557, seed `20261031`, `load_taxa` left-joined to `load_tax_groups` on `species_id` (`tax_group` as `group`, null
  where a species has none), then `assign_halves`. 2,529,340 species in 46,615 lineages, and the file is the same on a second run.
- **`SPLIT_SHA256` = `88f237b47d4964108c161e08acb00c9cf03601418dd2e0f16b49a43704845fa9`.** It hashes the content and not the Parquet file,
  whose bytes carry the writer's version, so a Polars upgrade would fail an unchanged split. The content is the UTF-8 rows
  `species_id TAB lineage TAB half LF`, sorted by `species_id` (byte order). `split.split_sha256` computes it for both the script and
  `split.load_split`, which raises `ValueError` when a file's hash differs. An independent Python sort and loop on the real file gives the
  same hash.
- **Balance:** half A has 1,264,661 species in 23,335 lineages, half B 1,264,679 in 23,280. By lineage rank: 15,866 families, 391 orders, 104
  classes and 30,254 species with no family, order or class, each its own lineage (15,170 in A, 15,084 in B).
- **Per group,** counting each species under its own group, so a lineage that spans groups is counted in each and the lineages sum to more than
  46,615:

| Group | Species A | Species B | Lineages A | Lineages B |
|---|---|---|---|---|
| (none) | 2,312 | 2,313 | 40 | 39 |
| algae | 10,138 | 10,138 | 287 | 288 |
| angiosperms | 185,871 | 186,764 | 229 | 229 |
| animals | 12,152 | 12,226 | 12,050 | 12,049 |
| arachnids | 57,131 | 57,131 | 413 | 412 |
| archaea | 382 | 382 | 23 | 24 |
| arthropods | 1,823 | 1,425 | 44 | 49 |
| ascomycetes | 56,613 | 56,613 | 2,248 | 2,249 |
| bacteria | 10,739 | 10,739 | 321 | 321 |
| basidiomycetes | 28,754 | 28,753 | 176 | 176 |
| bivalves | 14,295 | 14,296 | 174 | 176 |
| bryophytes | 10,381 | 10,380 | 136 | 137 |
| chordates | 46,475 | 46,452 | 935 | 937 |
| coleoptera | 188,848 | 188,872 | 107 | 108 |
| crustacean | 43,668 | 44,050 | 691 | 692 |
| diptera | 90,414 | 90,458 | 118 | 119 |
| eukaryotes | 3,439 | 2,605 | 672 | 672 |
| fungi | 108 | 141 | 91 | 1 |
| gastropods | 56,335 | 56,276 | 404 | 405 |
| gymnosperms | 881 | 880 | 13 | 13 |
| hemiptera | 58,699 | 58,698 | 143 | 143 |
| hymenoptera | 81,535 | 81,405 | 79 | 80 |
| insects | 260 | 453 | 38 | 3 |
| lepidoptera | 91,218 | 91,218 | 75 | 75 |
| molluscs | 271 | 272 | 17 | 21 |
| orthoptera | 15,492 | 15,491 | 42 | 43 |
| otheranimals | 80,789 | 80,788 | 1,756 | 1,756 |
| otherarthropods | 17,962 | 17,963 | 231 | 233 |
| otherfungi | 1,739 | 1,740 | 95 | 95 |
| otherinsects | 33,875 | 33,874 | 322 | 321 |
| othermolluscs | 7,676 | 7,618 | 318 | 320 |
| plants | 290 | 291 | 238 | 247 |
| prokaryotes | 0 | 1 | 0 | 1 |
| protists | 30,604 | 30,605 | 578 | 578 |
| pseudofungi | 1,048 | 924 | 1 | 37 |
| pteridophytes | 7,318 | 7,318 | 30 | 32 |
| trichoptera | 6,879 | 6,879 | 24 | 24 |
| viruses | 8,247 | 8,247 | 187 | 187 |

- **Where a group's halves differ by more than a few species,** one large lineage sits on one side: `Periconiaceae` (141 species) in B for
  `fungi`, `Peronosporaceae` (1,048) in A for `pseudofungi`, class `Insecta` (321) in B for `insects`, and `Cyclopidae` (1,085) and class
  `Copepoda` (388) in A for `arthropods`. A lineage is balanced under its most common group and takes its other-group species with it, so the
  856 `Asteraceae` species grouped as `eukaryotes` sit in A with the rest of the family. 23 lineages hold species of several groups (70,012
  species).

### ISEA3H grid, land and sea (2026-10-31)

`realms.grid` reads GBIF's `ISEA3H-{5,6,7}.gpkg` (`scripts/fetch_reference.py` puts them in `data/reference/grids/`) and splits each cell against
Natural Earth `ne_10m_land`. `scripts/check_grid.py` asserts, at all three resolutions, `10 * 3**r + 2` cells, exactly 12 cells with 5 neighbours and
all the others with 6, a total area within 0.05% of 510.06M km2, and a land total within 2% of Natural Earth's.

- **Median cell area:** 209,785 km2 at resolution 5 (2,432 cells), 69,888 km2 at 6 (7,292) and 23,298 km2 at 7 (21,872). The ranges are 147,798 to
  215,731, 60,071 to 70,617 and 16,307 to 23,965 km2, the smallest being the 12 pentagons.
- **Area is geodesic (WGS84), and the cells add up to 510,065,622 km2 at every resolution,** which is the area of the ellipsoid. Summing the polygons as
  the files draw them gives 510,057,396 km2 at resolution 6, 8,225 km2 short, because of the poles (below).
- **Land:** 146,836,143, 146,840,908 and 146,841,179 km2 at resolutions 5, 6 and 7, against 146,841,798 km2 for Natural Earth (the geodesic area of its
  polygons): -0.0039%, -0.0006% and -0.0004%. A cell's land is the geodesic area of its overlap with the union of the land polygons, kept within
  the cell's area. 812, 1,608 and 3,192 cells are coastal, with some land and some sea.
- **Neighbours share a boundary segment of positive length.** Comparing coordinates exactly finds too few: with the antimeridian already dealt
  with, only 7,143 of the 7,292 cells at resolution 6 had 6 neighbours, because the vertices two neighbours share differ by up to 1.4e-8 degrees,
  and two edges that are the same line are then no longer collinear. A boundary now counts when more than 1e-4 degrees of it lie within 1e-6
  degrees (about 10 cm) of the other's. At all three resolutions no two vertices of different cells lie between 1.4e-8 and 1 degree apart, so the
  tolerance has a wide margin.
- **Antimeridian:** the files cut the 86 cells that cross it (at resolution 6) into two pieces at -180 and 180, so a cell's neighbours were on
  either side of the line. Each cell is redrawn in the 360 degrees of longitude around its own centre, which joins the pieces, and each cell is also
  tried 360 degrees to either side of the other, which is the same as treating 180 and -180 as one line.
- **Poles:** at every resolution the pole lies on an edge between two cells, running from longitude L to L + 180 at one latitude, which is a geodesic
  over the pole but is drawn along the parallel. One of the two cells is drawn the long way round (a self-intersecting polygon, at resolutions 5 and 7),
  and where the cell is also cut at the antimeridian the cap over the pole is missing. With only the antimeridian dealt with, 4 cells per resolution broke
  the 12 and 6 count. Such an edge is now replaced by a detour along the two meridians to the pole and back, after which the cells tile the globe.
- **The id doesn't give the resolution.** Its first digit, after any minus sign, is 2, 5, 6 or 7 at resolution 6, and 2, 4, 5 or 7 at resolution 5, so the
  resolution is read from `res`. Ids have up to 19 digits, more than a float holds, so they are text from the moment they are read.

### Column check download (2026-10-31)

The clean filter reads `LOWER(taxonomicstatus)` for `ambiguous_synonym`, `misapplied` and `synonym`, and the occurrence API shows COL's status only
inside the COL classification (its top-level `taxonomicStatus` is the backbone's). If SQL's column were the backbone's, `n` and `n_synonym` would be
quietly wrong. A Denmark-only download, `queries/check_columns_dk.sql`, settles that and tests the cell codes against the grid files.

- **Download:** `0019716-261023074918325`, DOI 10.15468/dl.67gagy, submitted 2026-10-31 14:28 UTC. 309,066 rows, 67,114,870 records, saved as
  `data/gbif/0019716-261023074918325.zip`. Its one member is `<key>.csv`, which is tab-separated with a header. The API's COL `taxonomicStatus`
  facet for the same filter was saved at 14:29 UTC as `data/gbif/0019716-261023074918325.api-facet.json`.
- **Order:** the res 6 pilot (`0019721-261023074918325`, DOI 10.15468/dl.erz6hs) was submitted at the same time, so this check confirms its query
  after the fact and did not gate it. It holds, so the pilot's filter stands.
- **`scripts/check_columns.py data/gbif/<key>.zip`** asserts checks 1 to 3 and reports 4 and 5. It compares with the saved facet, not a fresh API
  call: the download is a fixed snapshot and a later call drifts with new records. Run it as `uv run python -I`, since it reads a downloaded file.
- **1. SQL `taxonomicstatus` is COL's.** The statuses are the API's five plus null, and not the backbone's `DOUBTFUL`, `HETEROTYPIC_SYNONYM` or
  `PROPARTE_SYNONYM`. That equality is the deciding test, and each count must also be within 1% of the API's. AMBIGUOUS_SYNONYM (67,757), MISAPPLIED
  (587) and PROVISIONALLY_ACCEPTED (163,644) match exactly. SYNONYM (1,690,718) is 0.01% below the API and ACCEPTED (64,876,895) 0.20% below, put
  down to snapshot timing. Null status covers 315,269 records, which is the API's 67,247,339 without `checklistKey` less its 66,932,070 with it. The
  snapshot holds only the second figure, so the script reports the null count and doesn't assert it.
- **2. Species:** all 41,892 distinct `specieskey` values (308,869 rows) are accepted species of 315557. 4,501,865 records, 6.7%, have no `specieskey`.
- **3. Cells:** all 8 `cell6` codes are `geogrid_id`s of `ISEA3H-6.gpkg` and all 18 `cell7` codes of `ISEA3H-7.gpkg`, compared as text. The first digit of
  an id is not the resolution (see the grid entry above), so the resolution comes from `Grid.resolution`, and the script checks that each file holds
  the resolution its name says.
- **4. Land and sea parts:** 11,099,401 of the 67,114,870 records (16.54%) are in the sea part, 56,015,469 in the land part. No record lies on a part of
  zero area, at either resolution: none in the sea part of a cell without sea, none in the land part of a cell without land. How many of the
  sea-part records are terrestrial species on the coast is still to be looked at, in the QA of the downloaded records.
- **A part of zero area is below 1000 m2, not exactly 0.** `Grid.cells` has float noise: cell `2888795476011250000` at resolution 6 lies wholly on land
  and has a `sea_km2` of 3.2e-8. The smallest real parts on the grids are 0.006 km2 of land and 0.26 km2 of sea, so the threshold splits them cleanly.
- **5. Growth:** 117,960 distinct (species, cell) pairs at resolution 6 and 163,091 at 7, a ratio of 1.38. Denmark spans 8 cells at resolution 6, so
  this is a lower bound.

### Ingest and QA of the res 6 download (2026-10-31)

`realms.occurrences` streams a SQL download of counts into Parquet, checks it, and reads it per half. `scripts/ingest.py data/gbif/<key>.zip
--resolution <r>` writes `data/gbif/<key>/part-00000.parquet` and onwards and `results/pilot/res<r>/qa.json` (both git-ignored), and exits non-zero
unless the three pass marks hold. Run it as `uv run python -I`, since it reads a downloaded file.

- **Run:** `0019721-261023074918325` (DOI 10.15468/dl.erz6hs), 19,531,246 rows, exit 0.
  - **Time:** the ingest took 4 s, the whole script 22 s (load the grid and the taxa, QA).
  - **Parquet:** 262 MB in 4 parts of up to 5,000,000 lines, from a 679 MB zip holding 1.15 GB of text, which is never written to disk.
  - **Peak memory:** 3.6 GB RSS for the ingest, read after it and before the grid and taxa load, and 8.2 GB by the end, after loading the taxa and
    the grid and running the QA over all 19.5M rows (not broken down further). The ingest holds one chunk of text and its table at a time, never the 19.5M rows.
  - **Check:** the rows (19,531,246), `sum(n_raw)`, `sum(n)` and `sum(n_synonym)` of the Parquet equal `awk`'s sums over the TSV member exactly.
- **Pass marks:** unaccepted species, unknown cells, and a synonym share of 2% to 5% of `n`. The ranges are not to be widened to pass a run.
  All three hold.
  - **Unaccepted species:** 0 species, 0 rows.
  - **Unknown cells:** 0 codes. 7,234 of the 7,292 cells of the resolution 6 grid hold records.
  - **Synonyms:** 100,329,919 of 3,381,603,728 clean records, 2.97%. The API's COL facets gave 3.6% of all records the same day, a different
    base, and the 0.6 point difference was not looked into.
- **Records:** 3,405,053,062 in all (`n_raw`) and 3,381,603,728 clean (`n`), 99.3%. 3,148,622,379 clean records are in land parts (15,002,114 rows,
  3,091 cells) and 232,981,349 in sea parts (4,529,132 rows, 6,675 cells). 229,803 rows, 3,347,345 clean records, have no `countrycode`.
- **Records on a part of zero area** (area below `ZERO_KM2`, 1000 m2), global counts only:
  - land parts of cells with no land: 35 cells, 2,152 rows, 24,676 records, 24,638 clean (0.0008% of the clean land-part records);
  - sea parts of cells with no sea: 955 cells, 405,640 rows, 2,337,244 records, 2,312,519 clean (0.99% of the clean sea-part records).
    These are records the SQL calls sea (no `level0gid`) in a cell that Natural Earth's land covers entirely. The cells with most of them lie in
    France, Switzerland, Germany and around the Great Lakes, so it is most likely inland water that GADM gives no country. Not looked into further.
  The records are kept for now, as parts like any other. Dropping them is not decided here.
- **`ZERO_KM2` is in `realms.grid`,** and both `scripts/check_columns.py` and `qa` import it, so a part of zero area means one thing.
- **Max over countries:** `load_presence` sums the counts of a (species, part) over its countries, except the four distinct counts `n_datasets`,
  `n_years`, `n_datasets_open` and `n_years_open`, which take the maximum over the countries. The query counts distinct datasets and years inside a
  (species, country, cell, part) group, and a dataset or a year can count in two countries of one border cell, so adding them overcounts. The true
  number lies between the maximum and the sum. The maximum is a lower bound, exact when one country holds the part, so a threshold of the form
  "at least k datasets or years" is conservative. `load_country_counts` keeps the countries apart, for GRIIS.
- **Seal on half B:** `load_presence` and `load_country_counts` raise `SealedHalf` for half "B" or "all" unless `unsealed=True`, before anything
  is read, and a species that is not in the split is left out of every half. `qa` and the script print and write global numbers only: no species,
  lineage or per-taxon figure. The unaccepted-species check gives how many species and rows and names none, where `check_columns.py` lists
  examples, so no figure can be traced to a taxon. Unknown cells are listed (up to 10 codes), since a cell is not a taxon.
  - **`unsealed=True` call sites so far:** none. The calibration on half B's plants is planned as the only one before the freeze, and will be listed here.
  - **`scan_counts` reads, which are not sealed:** `scan_counts` is public and returns every species, half B included, so the seal doesn't cover it.
    It has two readers, and both touch all species: `qa`, called from `scripts/ingest.py`, which produces global numbers only, and
    `_select` in `realms/occurrences.py`, which checks the seal before it reads and returns half A, or half B only when `unsealed=True`.
    `tests/test_occurrences.py` fails if any other module under `src/` or `scripts/` mentions `scan_counts`, so a new reader needs a ruling
    listed here.
- **Not done here:** how many sea-part records are terrestrial species on the coast (open since the column check) needs the species'
  environment, and wasn't part of this check.

### Pilot at resolution 6, and the scale decision (2026-10-31)

`scripts/pilot.py --counts data/gbif/0019721-261023074918325 --workers 12` ran the pipeline once on half A at resolution 6, from 20:49 to 21:48
(59 min), and wrote `results/pilot/res6/measure.json` and `pilot.log` (git-ignored). `--project-only` makes the projection again from a saved
`measure.json`. Nothing of half B is loaded, printed or written: `load_presence(..., "A")` with no `unsealed`, so the call-site list above is still
empty, and the module summaries below are of half A's parts.

**Data.** 638,157 of half A's 1,264,661 species have a clean record (`n > 0`) in 8,945,552 (species, part) rows, over 9,556 parts (3,066 land and
6,490 sea) with 28,843 pairs of adjacent parts, two of them with no neighbour. The matrix has 727,680 taxa (638,157 species, 82,306 genera, 7,217
families) and 14,944,355 links, 737,236 network nodes with the parts. **923 parts have an area below `ZERO_KM2`** (33 land, 890 sea), holding 199,562
rows, 2.2% of the rows and 0.066% of the clean records. They are kept as they are.

**Stages** (one trial of Infomap each, a process's peak is its own: Infomap, `candidate_maps`, the taxa load and the workers' parent run in child
processes, and the parent's high-water mark of 7.0 GB is `load_presence`):

| Stage | Time | Peak RSS |
|---|---|---|
| load taxa (child) | 12.9 s | 15.1 GB (12.2 to 15.1 over the runs) |
| load split, grid, presence | 0.2 s, 3.7 s, 1.2 s | parent 7.0 GB |
| parts, `cell_taxa`, `network_from`, save links | 0.2 s, 1.1 s, 0.1 s, 0.1 s | the matrix is 0.075 GB, the `.npy` 0.36 GB |
| Infomap 1 trial, Markov time 1.0 | 32.3 s (2.6 s of it loading the links) | 5.15 GB (2.1 GB after loading) |
| Infomap 1 trial, Markov time 0.75 | 30.2 s | 4.16 GB |
| Infomap 10 trials, Markov time 0.75 | 238 s, 7.9 times one trial | 5.11 GB |
| `candidate_maps`, Markov time 0.75 | 1,555 s | 3.03 GB |
| one bootstrap, 3 trials, 2 of them in 2 workers | 1,328 s and 1,009 s | workers 4.93 and 4.81 GB, their parent 0.77 GB |
| `excess_table` | fails, see below | |
| `excess_table`, without the zero-area parts | 13.5 s for one map | 0.95 GB |

- **`candidate_maps` at 0.75:** `resolve_pieces` 0.8 s, `apply_size_floor` 1,554 s, `merge_sequence` 0.15 s. The floor makes 7,039 merges of 0.22 s each,
  as every merge goes over a dense units x units array (0.4 GB, 0.2 s at 7,044 units), and leaves 5 units, so `merge_sequence` has nothing to
  merge and there is one map, k = 5. The bootstrap scores only that map, with a median stability of 0.31 over its 2 draws.
- **Infomap trials run one after the other:** 10 trials cost 7.9 times one, and its peak is 1 GB above that of one trial.
- **Markov time 1.0 is degenerate on half A:** Infomap puts 9,487 of the 9,556 parts (99.3% of the parts and of the area) in one module, and
  `candidate_maps` raises `DegenerateClustering`, 2 units above the size floor where 5 are needed. Recorded, not raised, and the later stages
  used 0.75 as the brief says.

**Module structure, one trial, seed 1, half A only.** Both partitions have one level. The place of a module is the median of its parts' centres.

| Markov time | Modules | Largest | Modules of one part | The largest, in parts (median centre) |
|---|---|---|---|---|
| 1.0 | 35 | 99.3% of the parts | 31 | 9,487 (5 N, 7 E: all the rest of the globe), 33 (20 S, 48 E), then 3 and 2, and 31 single parts |
| 0.75 | 7,044 | 15.7% of the parts | 6,874 | 1,497 (40 N, 33 E), 632 (5 S, 64 W), 39 (30 S, 20 E), 32 (19 S, 47 E), 23 (1 N, 113 E), 21 (21 N, 159 W), then 11 or fewer |

At 0.75 the six largest modules hold 2,244 parts, 23% of them, and 7,039 of the 7,044 are under the 1,000,000 km2 floor. The structure is made of a
few real blocks and thousands of single parts that the floor then merges one by one.

**Between them** (extra, `SURVEY` in the script: one trial each, and `candidate_maps` on each, fewest units first):

| Markov time | Modules | Largest | One part | `apply_size_floor` | Units after the floor | Maps |
|---|---|---|---|---|---|---|
| 1.0 | 35 | 99.3% | 31 | 0.3 s | degenerate (2) | 0 |
| 0.97 | 43 | 99.2% | 38 | 0.4 s | degenerate (2) | 0 |
| 0.95 | 62 | 89.6% | 44 | 0.6 s | 5 | 1 |
| 0.9 | 176 | 77.0% | 154 | 1.7 s | 7 | 3 |
| 0.85 | 705 | 56.3% | 651 | 8.4 s | 10 | 6 |
| 0.8 | 2,567 | 42.3% | 2,475 | 87.9 s | 10 | 6 |
| 0.75 | 7,044 | 15.7% | 6,874 | 1,554 s | 5 | 1 |

Infomap costs 30 to 40 s and 4.1 to 5.2 GB at every Markov time. `candidate_maps` costs about the units to the power 1.8 from 705 to 2,567 units and 2.8 from 2,567 to 7,044, as the dense arrays take over, and
nothing below 700 units. After the floor there are only 5 to 10 units at every Markov time from 0.95 down, so at most 6 maps (k = 5 to 10), where the selection
assumes 56 (k = 5 to 60).

**`excess_table` fails on the zero-area parts.** `null_partitions` divides by the area of each piece of a realm (`nullmodel.py`,
`-remaining[p] / targets[p]`), and a piece made only of parts with no area has none: `ZeroDivisionError: float division by zero`, at Markov times
0.75 and 0.9 alike. It looks like a limit of earlier code with the zero-area parts kept, which the brief says to keep and a later task folds away,
so it is **not patched here**. To time the stage and weigh its workers, it was run again on the same maps without the 923 zero-area parts (8,633
parts, labels, adjacency and areas cut down, the matrix as it was). That run is the `excess_table_without_zero_area_parts` row, it is not the
pilot's map, and its excess (0.061 for k = 5) is not a result. With one map only one of the 12 workers had a task. Earlier, on 3 maps and 4 workers
(Markov time 0.9), each worker peaked at 0.54 to 0.62 GB and took 16 s a map.

**`--workers 12`** is the 14 cores less two for the parent and the rest of the machine. Memory allows many more: the bootstrap's workers peak at
5 GB and the excess's at 0.6 GB, so 12 excess workers need about 9 GB beside a parent of 1 GB. The numbers came from the calibration run (Markov times 1.0
and 0.9, `--results results/pilot-calibration`), which ran before this one.

**Projection to resolution 7** (`project` in the script; each stage is scaled from its measure and the larger option is taken wherever it is
unsure). The headline has 2 times the links of half A, so 5 times (R3, growth 2.5) or 2.76 times (Denmark's 1.38) the links measured: 74.7M and
41.2M. Cells grow 3.0 times (7,292 to 21,872). Infomap, the resampling and their memory scale with links, with no fixed share, and a 10-trial
Infomap is the larger of 10 one-trial runs and the 10-trial run measured. `candidate_maps` has units, so cells, times the merges, each over a dense
units x units array (cells squared) or the links, whichever is slower, and the dense arrays need cells squared the memory. `excess_table` scales
with the larger of cells and links, and with 56 maps. **1.38 is a lower bound** on the growth, as Denmark spans 8 cells, and is only shown beside
2.5, which decides.

| | Growth 2.5 (decides) | Growth 1.38 |
|---|---|---|
| Links at resolution 7, headline | 74.7M | 41.2M |
| (a) headline Infomap peak, alone, with the matrix; limit 40 GB | 25.8 GB, 26.2 GB: **ok** | 14.2 GB, 14.4 GB: ok |
| (b) selection run, limit 24 h | **245 h: fails** | 245 h: fails |
| Of it: 10-trial Infomap | 0.42 h | 0.23 h |
| Of it: `candidate_maps` of the headline | 11.65 h (27.3 GB) | 11.65 h (27.3 GB) |
| Of it: 20 bootstraps, 1 worker (27.3 GB each) | 233 h, 11.65 h each, nearly all `candidate_maps` | 233 h |
| Of it: `excess_table`, 7 workers (4.8 GB each) | 0.15 h | 0.05 h |
| (b) without `candidate_maps` | 0.57 h | 0.28 h |
| (c) disk left after the zip and Parquet (2.35 GB; free now 55.7 GB); minimum 10 GB | 53.3 GB: **ok** | 54.4 GB: ok |

**Decision: stay at resolution 6.** (b) fails by a factor of 10, and all of the failure is `apply_size_floor`: at resolution 7 the 7,044 units of
Markov time 0.75 become about 21,000, with a merge of about 2 s over a 3.6 GB dense array, 11.7 hours and 27 GB for the headline's `candidate_maps`
and again for each bootstrap. Taking `candidate_maps` out leaves 0.57 hours, and (a) and (c) hold, so the memory and the disk are not what rules
resolution 7 out. Two things would change the answer, and neither is decided here:

- **The Markov time, or whatever leaves fewer units.** The same pipeline with Markov times 1.0 and 0.9 (176 units, 7 after the floor, 3 maps, 2.3 s
  of `candidate_maps`, a bootstrap of 96 s and 5.3 GB) gives 26.8 GB for (a), 3.6 h for (b) and 53.8 GB left for (c): go, and 1.25 h at 1.38.
  Markov time 1.0 gives one module with the whole globe in it, so no setting yet gives both a usable partition and a cheap one.
- **A cheaper `apply_size_floor`,** which recomputes a dense units x units beta-sim for every merge where one row changes. It was not touched.

**Correction (2026-11-01, Task 8a):** the "(b) without `candidate_maps`: 0.57 h" in the table above left out the bootstraps' own work. `project`
subtracted the headline's `candidate_maps` time (1,555 s) from a bootstrap worker's whole time (1,328 s), which is another call, and the difference
was negative, so it clamped it to zero and costed the 3-trial Infomap, the resample and `network_from` of 20 bootstraps at nothing. They are 76 s of
an 80 s worker (measured in the re-run below), which at 5 times the links is 0.1 h a bootstrap with one worker fitting: **2.68 h at growth 2.5 and
0.86 h at 1.38**, where it said 0.57 h and 0.28 h. The decision above is not changed by it, as the 245 h were nearly all `candidate_maps`.

**Measurement notes.** The peaks of Infomap, `candidate_maps` and the workers come from `ru_maxrss` inside each child, so they are exact. The
parent's high-water mark is not a stage's own. The bootstrap ran 2 draws in 2 workers, so its time is that of 2 running together and not of 12.
`run_infomap`'s `linked` set was left as it is: at 14.9M links it takes 0.2 s, and its list of Python ints is about 0.5 GB by size (not measured), a tenth of Infomap's own 5.2 GB. Five
times that at resolution 7 is still a tenth of Infomap's 26 GB.

#### After the incremental size floor (Task 8a) (2026-11-01)

`apply_size_floor` (commit `81a05e7`) no longer builds a units x units array, dense or sparse. It keeps a taxa -> units index, the transpose of the
units x taxa matrix with the unit of each entry changed as units merge (one of the two entries of a taxon held by both is dropped), and the richness
of each unit as a vector. One merge then reads the overlap of the one small unit with every other in a single gather, picks the partner by the same
rule as before, takes the next small unit from a heap, and the labels are resolved once at the end through `absorbed_by`. `merge_sequence` still uses
`_Merger`, which lost its `among` argument, as only the old floor used it. **The results are identical, down to the order of merges and every tie:**
the partner is the lowest `1 - shared / min(richness)` in float64 (1.0 where the poorer unit has no taxa), ties to the lowest position, and the
smallest area goes first, ties to the lowest label. Checked against the frozen reference on 10 planted cases, 5 with units of equal areas and
identical taxon sets, 3 of 400 units over 3 ranks (90% under the floor) and a unit without taxa, and against the old code on half A (Markov times
0.75, 0.8 and 0.85, the old run took 1,527 s at 0.75): `maps`, `sub_realms` in their order, and `units` are equal. The Costello selection gives
the same `selection.csv` and `selection.log` byte for byte. Flipping the partner's tie to the highest position fails all 19 floor tests against the
reference, and flipping the tie in area to the highest label fails 4 of the 5 tie tests.

**Stages** (the pilot run again, `results/pilot/res6/measure.json`, the earlier one kept as `measure-before-8a.json`; same inputs, the Infomap runs
agree with the first run in modules and in seconds):

| Markov time | Units after pieces | `candidate_maps` | `apply_size_floor` | Units after the floor |
|---|---|---|---|---|
| 0.95 | 62 | 0.6 s (1.1 s) | 0.14 s (0.55 s) | 5 |
| 0.9 | 176 | 0.9 s (2.3 s) | 0.20 s (1.65 s) | 7 |
| 0.85 | 705 | 1.6 s (9.9 s) | 0.20 s (8.4 s) | 10 |
| 0.8 | 2,567 | 2.6 s (90 s) | 0.56 s (87.9 s) | 10 |
| **0.75** | 7,044 | **4.6 s (1,555 s)** | **3.63 s (1,554 s)** | 5 |

(The figures in brackets are Task 8's.) At 0.75, `resolve_pieces` is 0.83 s, `apply_size_floor` 3.63 s for 7,039 merges (0.5 ms each) and `merge_sequence`
0.19 s, which is 4% of `candidate_maps`, so it stays on `_Merger` (it is 7 to 23% at the survey's other times, which is 0.14 to 0.19 s there, and it is
the one stage that has dense arrays, so a floor that leaves hundreds of units would bring it back: 68 ms a merge, Task 3). The peak of the call is
0.99 GB (3.03 GB before), under the 4.15 GB of one Infomap run. **The target of 60 s and one Infomap's peak are met.** `apply_size_floor` costs the
units to the power 0.8 from 705 to 2,567 units and 1.85 from 2,567 to 7,044. A bootstrap worker takes 80.2 s and 79.5 s (1,328 s and 1,009 s),
of which `candidate_maps` is 4.2 s and 3.2 s, and the rest, 76 s, is the 3-trial Infomap, the resample and `network_from`; the peaks are 4.83 and 4.43 GB.
Its median stability is 0.311 as before, the maps being the same. The 10-trial Infomap takes 249 s (5.07 GB), `excess_table` is where it was (the
zero-area failure, then 14.7 s and 0.95 GB for one map without those parts, which is a timing and memory proxy only).

**`project`, changed to fit** (script and `tests/test_pilot.py`): `candidate_maps` time is the measure times cells (c = 3.0) times the larger of the
cell growth and the link growth f, which is c x f = 15 at growth 2.5 (the brief's cells times links growth) and c x c = 9 at 1.38 (where the links
grow less than the cells and a merge's vectors, with a slot per unit, are the slower part); the work of a merge is the taxa of one unit and the units
holding them, which is links, and the merges are c times as many. Its memory follows the links (f), with no array of units squared. Each bootstrap
worker times its own `candidate_maps` inside the worker, so that its seconds and the rest of the worker's (links, f) are scaled apart, taking the
slowest worker's of each, and the headline's `candidate_maps` no longer stands in for the worker's. If neither excess stage was measured, `project` goes on
without it, lists it under `unmeasured` and does not say go. `--project-only` needs only `--results`. **Workers are as many as fit in 40 GB for all of
them together** (`BUDGET_GB`, next to the parent, and at most one fewer than the cores), and a single process (the headline Infomap, alone or with
the matrix) is held to the same 40 GB.

| | Growth 2.5 (decides) | Growth 1.38 (a lower bound) |
|---|---|---|
| Links at resolution 7, headline | 74.7M | 41.2M |
| (a) headline Infomap peak, alone, with the matrix; limit 40 GB | 25.9 GB, 26.3 GB: **ok** | 14.3 GB, 14.5 GB: ok |
| (b) selection run, limit 24 h | **3.05 h: ok** | 0.97 h: ok |
| Of it: 10-trial Infomap | 0.40 h | 0.22 h |
| Of it: `candidate_maps` of the headline | 0.02 h (5.0 GB) | 0.01 h (2.7 GB) |
| Of it: 20 bootstraps | 2.47 h, 1 worker of 24.1 GB, 0.12 h each of which 0.02 h `candidate_maps` | 0.69 h, 2 workers of 13.3 GB |
| Of it: `excess_table` | 0.16 h, 7 workers of 4.8 GB | 0.05 h, 13 workers of 2.6 GB |
| (b) without `candidate_maps` | 2.68 h | 0.86 h |
| (c) disk left after the zip and Parquet (2.35 GB; free now 54.6 GB); minimum 10 GB | 52.2 GB: **ok** | 53.3 GB: ok |

**Decision: go to resolution 7** on the same three tests, with growth 2.5 deciding: (a), (b) and (c) all hold, and 1.38 agrees. The earlier "stay at
resolution 6" came from `apply_size_floor` alone, and it is gone. What is left of (b) is the bootstraps' own Infomap, 2.1 of the 3.05 hours, at
one worker because five times the links makes a worker 24 GB. The decision is not sensitive to the new `candidate_maps` model: its share of (b) is
0.42 h, so it would have to cost 51 times more to reach 24 h. What it does not rest on is a measurement at resolution 7: every figure is the res 6
measure scaled, and the first real resolution 7 stage will be the check. At Markov time 0.75 half A still leaves 5 units after the floor and one map
(k = 5), where a selection run assumes 56, and `excess_table` is still not run on the pilot's own map (zero-area parts), so neither is settled here. The
bootstrap measure ran 2 draws at once, not 12.

## 2026-11-01

### Reference joins

`src/realms/natives.py` and `lists.py` join WCVP, GRIIS and WRiMS to COL 315557 ids, and `load_references(grid)` returns them as `References`.
Figures are from `scripts/check_natives.py` on the real exports (`results/natives.json`).

**The 24 missing exports now exist.** Markus requested WRiMS and the 23 GRIIS lists. `fetch_reference.py` again: 316 items, 0 missing,
and WRiMS (26853) has an extended export of 2026-10-31T23:47 (4.5 MB). Layer 3 therefore runs with WRiMS and all 304 GRIIS lists.
`load_references` still skips and names (`References.missing`) a WRiMS or GRIIS export the manifest lists as missing, and raises only for
COL, WCVP, the two polygon files and the grid.

**Licences, as found in the sources** (also in `manifest.json`):

- **WGSRPD level 3** (`github.com/tdwg/wgsrpd`, `geojson/level3.geojson`, pinned to commit `52da782`, 2021-07-15). The repository has no
  licence file and the GitHub API reports none. The shapefile metadata (`level3/level3.shp.xml`) says Access_Constraints "None" and
  Use_Constraints "None, but dataset is in first release. Thanks to ESRI for use of the Orginal Administration Boundaries", with ESRI's
  notice "Administrative boundary files are the intellectual property of ESRI and its licensors and are using herein with permission.
  (c) 1992-1997 ESRI, GMi. All Rights Reserved." Credit: Justin Moat, RBG Kew, for Brummitt 2001 (TDWG, Plant Taxonomic Database Standards
  No. 2). **That is no open licence.** The polygons stay in `data/` and are not committed or republished, and the published result
  should cite Brummitt 2001 and say where the regions came from.
- **Natural Earth admin-1** (`ne_10m_admin_1_states_provinces`, 5.1.1): public domain, per naturalearthdata.com/about/terms-of-use/.

**The name index.** `name_index` has 4,715,525 rows: the 2,354,462 accepted and 174,878 provisionally accepted species, and 2,186,185
names of species-rank synonyms pointing to one of them. A synonym of a subspecies or variety rolls up to that taxon's species (about 160,800
of 2.3M species-rank synonyms), as GBIF counts those records for the species. Synonyms whose accepted taxon is a genus (1,150) or higher are
left out, as are ambiguous synonyms (129,976) and misapplied names (6,093), since a pro parte name could belong to several species. The index
has a `status` column besides the brief's columns, which the matching and `NO_WILD_RANGE` need. It takes 4 s and 3.2 GB to build.

**Matching** is by exact `col:scientificName` (no authorship) and kingdom, with the clash narrowed to Animalia, Plantae and Fungi as in the
occurrence query: a plant never matches an animal homonym, while `Bacteria` or `Viruses` against COL's `Bacillati` or new virus kingdoms is
no clash. WCVP has no kingdom column, so its names count as Plantae. A plain GRIIS list has none either, and takes the kingdom row above
each entry.

- **Ties are broken by standing, not guessed.** The brief reads: a name that fits several accepted species is ambiguous. Read strictly,
  10,037 of WCVP's 371,830 names (2.7%, which is 8.6% of its distribution rows) and 4,660 of GRIIS's 31,575 (14.8%, which is 26.5% of its
  entries) are ambiguous. 7,948 of WCVP's are one accepted and one provisionally accepted species of the same name (*Zea mays* is `5CXX5` and `VJSFD`), as
  are most of GRIIS's, which are the duplicates the "Taxa by id" section already found. Now an accepted species wins over a provisionally
  accepted one and that over a synonym's, as COL and GBIF rank them, and only a tie within a standing stays ambiguous, never guessed:
  102 WCVP names (0.1% of its rows) and 136 GRIIS names (0.4% of its entries). Accepted over synonyms alone would leave 8,071 names (7.3%
  of rows) and 3,665 (22.5% of entries). **Markus can revert this** by dropping the standing order in `match_names`; the cost is the
  8.6% of WCVP's rows and the 26.5% of GRIIS's entries that then join to nothing.
- **Match rates** (distinct names; rows in brackets). Ambiguous names count as neither matched nor failed.

| Source | Names | Matched | Ambiguous | Unmatched |
|---|---|---|---|---|
| WCVP | 371,830 | 98.3% | 102 | 6,217 (1.7%) |
| WRiMS | 2,511 | 98.6% | 0 | 35 (1.4%) |
| GRIIS, 304 lists pooled | 31,575 | 91.0% (95.7% of 120,547 entries) | 136 | 2,704 (8.6%) |

  A further 617,893 of WCVP's 1,984,179 distribution rows (31%) are of genera (338,608) and infraspecific taxa (279,285), which no species
  name joins, and 12,800 more belong to the unmatched names. Of the unmatched GRIIS names, 743 hold non-ASCII text, because some lists
  carry authorship in `scientificName` with its accents broken (`ForsskÃ¥l`), and 488 hold punctuation. Nine of the 303 GRIIS lists that have entries have more than
  10% unmatched (CY-Cyprus, DZ, IL, KW, MH-Namdrik-2, SY, TC, TK, WF; the largest is IL with 762 names and 10.5%). No whole source does.
- **GBIF fallback for those nine only.** `gbif_fallback` asked the match API v2 (`checklistKey`, at most 5 requests a second, GET only)
  about their 132 unmatched names and took 97, as `EXACT` or `VARIANT` matches whose accepted usage is a species of 315557 and whose
  kingdom doesn't clash; a synonym of a subspecies, a higher rank or a species not in this release is refused. They are saved to
  `data/reference/gbif_fallback.parquet`, which `load_references` applies to any name still unmatched, so a rerun doesn't depend on
  GBIF. **The file is tied to the COL release:** its Parquet metadata holds the release key (`write_fallback`), and `load_references`
  raises if the file has none or another than the manifest's `col` release, since a species id means nothing in another one. It drops
  with a warning that counts them any placed species the current index lacks. GRIIS rose to 91.3% of names (2,607 unmatched) and WRiMS to 98.7%.

**WCVP.** Species rank only; WCVP's genus and infraspecific distributions aren't added (of 371,856 species, 3 have a subspecies' or
variety's distribution and none of their own). An empty `establishmentMeans` is native (the "Reference data" section above), `introduced`
isn't, and the 2,324 `uncertain` rows are dropped. 330 rows hold a level 1 or 2 code (`84`, `7`) and are dropped, as they can't say which
region. Both native and introduced regions stay, each with its level 3 code, since Task 11 has three modes. Where one species gets both for a
region (two WCVP names that COL holds as one species) it is native. Result: 1,341,114 rows of 356,880 species, 1,198,410 native and 142,704
introduced.

**GRIIS.** Only accepted species-rank entries are alien taxa: the 66,316 higher ranks (families 32,434, orders 15,608, `other` 5,761, classes
4,453, genera 2,230, phyla 2,111, kingdoms 878) and 1,924 subspecies are dropped, and so are synonyms. **Subspecies are dropped, not rolled up to
their species,** because an alien subspecies doesn't make its species alien: rolling up would mark the whole species as alien in the country
where only one infraspecific form is, and clean away the native records of the rest. A GRIIS list names a subspecies as a taxon in its own right, and the species-level join has nothing to put it on. The country and subdivision come from
the file name. Result: 112,116 pairs of 27,389 species. **A finding the plain lists can't show:** 3,457 of the 52,462 species of the 55
extended lists carry `native` distribution rows (`GRIIS-EC` and `GRIIS-SR` flag every species so, though *Cajanus cajan* is no native of
Ecuador, and the US list's own remarks say "Native to: Eurasia" for species flagged native). The flag looks like an import artefact, and as
the plain exports lack it, **every listed species is kept as alien**, as "Reference data" above argued.

**WRiMS.** 8,573 area rows are `introduced` (of 11,093: 380 `native`, 8 `uncertain` and 727 with no means, which are left out) and 99 of
them belong to genera and infraspecific taxa. Result: 8,347 pairs of 2,139 species in 686 areas. `References.wrims` has `species_id`,
`area_id`, `area` and `country`. **The areas don't map to sea parts.** They are Marine Regions places (MRGIDs) with no geometry in the
export: 244 of the 686 are the part of a sea inside one country's waters ("Italian part of the Tyrrhenian Sea"), 78 are exclusive economic
zones, and the rest are whole seas, countries, states and harbours. `country` is read from the name (a table of 100 demonyms and country
names in `MARINE_COUNTRIES`) for 6,352 of 8,347 pairs (76%), and is null for an open sea, a bay or a harbour. Placing a pair on sea parts needs
the polygons from Marine Regions' WFS, which isn't fetched. Whether GBIF's `countrycode` is set for marine records at sea is **not
checked**, so `country` may not be usable to clean them.

**Land parts and TDWG regions.** `tdwg3_parts(grid, level3)` gives each land part (`land_km2` above 1000 m2) every level 3 region whose polygon
shares some area with the cell's polygon, which is the whole cell and not its land, as WGSRPD draws its own coast. Cell polygons are drawn as
in `grid.py`, and regions are tried a turn of longitude to either side for the antimeridian. At resolution 6, 3,145 land parts: 3,065 have a
region (1.57 on average, at most 7) and 80 have none (specks of coast and ice the scheme leaves out, the largest the Ross Ice Shelf at 1,150
km2). At resolution 7: 8,279 parts, 163 without. The polygons are in the grid file, not in `Grid`, so the function reads the file named in
the manifest (`isea3h_<resolution>`).

**Sub-national GRIIS lists on the grid.** `admin1_cells(grid, admin1, subdivisions)` places a list's subdivision on Natural Earth's admin-1
polygon of that name (compared in any of its name fields, folded to letters and digits; the list's own country where several polygons carry
it, any country where one polygon alone does, since GRIIS says `BQ` for Saba and Natural Earth `NL`). `GRIIS-US-United` is the conterminous
states, which is the US less Alaska and Hawaii. A cell belongs by its centre, as the brief says, and **islands that no centre falls on take the
cells whose polygon overlaps them** instead (`via` = `overlap`), as the cells are 70,000 km2 at resolution 6: by centre alone only Alaska,
Svalbard and the US list got any cell, and Hawaii, the Azores, Madeira, Réunion and 17 others got none. Result at resolution 6: 186 cells for 24
of the 85 subdivisions (138 by centre, 48 by overlap; 457 cells at resolution 7). **61 subdivisions have no polygon in Natural Earth** and
their 7,468 pairs (6.7%) can't be placed: the 33 Marshall Islands atolls, the 7 Canary Islands (`IC`), 5 Tuvalu islands, Tonga's Niua and
Protected, St Helena, Tristan, Juan Fernández, Rapa Nui, Galápagos, Pitcairn, Cyprus, Clipperton, Cayman, Soqotra, Huvalu and the two French
southern lists. The country of a sub-national list is the alias's, which is the sovereign's for overseas places (`FR-Reunion`) and `IC` for
the Canaries, so Task 11 should match those on the cells and not on GBIF's `countrycode`.

**`NO_WILD_RANGE`: 16 names, 16 COL ids.** Humans, 13 domestic animals and 2 crops.
Each was checked against 315557: the name must be an accepted species there (`resolve_no_wild_range` checks it again at load, taking every
id of the name within its kingdom, and drops with a warning a name COL has as a synonym, like `Gallus domesticus`), and COL must not file a
wild form of the taxon under it. Entries, each with its reason in `lists.py`: *Homo sapiens*; *Bos taurus*, *Bos indicus*, *Bubalus bubalis*,
*Equus caballus*, *Equus asinus*, *Canis familiaris*, *Felis catus*, *Sus domesticus*, *Camelus dromedarius*, *Camelus bactrianus*, *Lama
glama*, *Vicugna pacos*, *Cavia porcellus*; *Triticum aestivum*, *Brassica napus*.

- **Dropped for covering wild populations** (`NO_WILD_RANGE_DROPPED`): *Ovis aries* (COL holds *Ovis orientalis*, the mouflon, as
  its synonym), *Capra hircus* (COL files the Chiltan markhor, a wild population, under it, as the synonyms *Capra falconeri chiltanensis*
  and *Capra hircus chialtanensis*; the first version of the list kept it as a minor case, and the rule has no minor cases),
  *Gallus gallus domesticus* (COL has no such species, *Gallus domesticus* being a synonym of *Gallus gallus*, which holds the
  red junglefowl), *Zea mays* (accepted subspecies *parviglumis* and *huehuetenangensis*, the wild teosintes), *Hordeum vulgare* (wild barley as
  subspecies *spontaneum*), *Glycine max* (wild soybean as subspecies *soja* and the synonym *Glycine soja*), *Solanum tuberosum* (186
  species-rank synonyms, wild potatoes among them, and WCVP gives it an Andean native range), and from the first look *Oryza sativa*,
  *Sorghum bicolor* and *Capsicum annuum*. **The brief named maize, barley, soybean and potato as examples, and these four are the ones
  dropped.** Markus can put them back in `NO_WILD_RANGE`; the cost is the teosinte, wild barley, wild soybean and wild potato records going
  with them.
- **Kept with a caveat:** *Canis familiaris* also holds the dingo, a feral dog with no natural range. COL holds the dog apart from the wolf (*Canis lupus*, left out like *Sus scrofa*, which holds the wild boar),
  so the case the brief imagined, a dog as *Canis lupus*, doesn't occur here.

**What `References` says it could not place,** so that Task 11 needn't recompute it: `griis_no_cells` holds the `griis` pairs of a
sub-national list whose subdivision has no admin-1 cells (7,468 pairs of 61 subdivisions at resolution 6, the figure above), and
`wrims_no_country` the `wrims` pairs whose area names no country (1,995 of 8,347). `load_references` also warns, with the count and the
first examples, of WCVP level 3 codes that the polygon file's `LEVEL3_COD` lacks, since no land part would have them. There are none in
the real files (the rows with a code that isn't three letters are dropped earlier, see WCVP above).

**Cost.** The check of all sources takes 56 s and peaks at 13.1 GB (16.2 GB on the first run, so
it varies from run to run), almost all of it `load_taxa`, which `load_references` calls and which was
already known to need about 15 GB. Everything new reads the exports through Polars scans and peaks at 3.2 GB (`name_index`).

