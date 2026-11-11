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
  - **`unsealed=True` call sites so far:** one, implemented and not yet run. `score_held_out` in `scripts/calibrate.py`, **in a run of its own**: `--score-half-b --setting N` (a required `N`, the
    setting Markus chooses after reading the dry run at resolution 7) scores setting `N` of the 81 on half B's **calibration species**: the plants WCVP gives a native region **and every other species of their
    taxonomic groups**, which are not all plants (`under_reported` needs the group's records per cell, so the score is like for like with A's). Of the 197,567 species of half B that it asks for, 2,586 are not
    Plantae in COL (see the calibration section for the breakdown, counted from the split and the taxonomy alone). **That run loads the references and those species of half B, and nothing else**: no known cases, no
    grid of settings, no regions, no cross-check, and nothing of half A, so A's data is not in memory and the setting is not spent with the choice. **The dry run, without the flag, never loads half B.**
    It is the only call planned before the freeze and **runs once, at resolution 7**. Before anything is loaded it refuses: any resolution but 7 (the ingest's, from its `_COMPLETE`); **unless the dry run at
    resolution 7 is complete** (`results/calibrate/res7/report.json`, which the dry run removes when it starts and writes last, at resolution 7 on half A with its `chosen` setting, from the same ingest by the zip's
    sha256); and **if a scoring was started or done**: any `results/calibrate/res*/held_out.started` or `res*/held_out.json` under `--results`, **or the fixed marker `data/split/held_out-315557.started`** next to
    the frozen split, which holds whatever `--results` says and is not in `results/`. **Both markers are written by exclusive creation before half B is loaded** (the time, the commit, `N`, the dry run's choice,
    whether `N` is that choice, the counts and the zip's sha256), and `held_out.json` after it (`N`, the dry run's choice, whether they are the same, the thresholds, the score and how many species were loaded), so a
    run that fails after the load (out of memory at resolution 7, or `calibrate` raising) is not retried by running the script again, as that would be a second unsealed load. A failed run needs Markus's explicit
    decision, and there is no override flag. **The content of `data/split/held_out-315557.started` is copied into this file, which is tracked, as soon as it is written**, so that the record of the one load does not
    rest on a file git ignores. The split is cut to those species before it loads, so no other species of half B is read. **Only counts and rates are printed or saved** (the score, and how many species were loaded and
    how many are not Plantae): nothing per taxon. `tests/test_calibrate.py` scans the syntax tree of every module under `src/` and `scripts/` and fails unless this is the only call that passes `unsealed=` with anything
    but the constant `False`, or `**` keywords to a loader, or calls `_select` (which takes `unsealed` by position) from outside `occurrences.py`, or holds the string `"unsealed"`.
  - **`scan_counts` reads, which are not sealed:** `scan_counts` is public and returns every species, half B included, so the seal doesn't cover it.
    It has two readers, and both touch all species: `scripts/ingest.py`, for `qa` and for the cell check of `--mark-complete`, which produce global numbers only, and
    `_select` in `realms/occurrences.py`, which checks the seal before it reads and returns half A, or half B or `"all"` only when `unsealed=True`.
    `tests/test_occurrences.py` fails if any other module under `src/` or `scripts/` mentions `scan_counts`, so a new reader needs a ruling
    listed here, and (since the final review) **if any module but `occurrences.py` reads the parts with a Parquet reader of its own** (`scan_parquet`, `read_parquet` and the like, on a path with `part-` in it,
    a name made from one, or the counts' directory), found by the syntax tree.
- **Not done here:** how many sea-part records are terrestrial species on the coast (open since the column check) needs the species'
  environment, and wasn't part of this check.

### Pilot at resolution 6, and the scale decision (2026-10-31)

`scripts/pilot.py --counts data/gbif/0019721-261023074918325 --workers 12` ran the pipeline once on half A at resolution 6, from 20:49 to 21:48
(59 min), and wrote `results/pilot/res6/measure.json` and `pilot.log` (git-ignored). `--project-only` makes the projection again from a saved
`measure.json`. Nothing of half B is loaded, printed or written: `load_presence(..., "A")` with no `unsealed`, so the pilot adds no call site to the list above
(which was empty then and has one entry since, the scoring of half B, not yet run), and the module summaries below are of half A's parts.

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

**GRIIS** (the ruling that every listed species is alien is superseded on 2026-11-01 by "GRIIS by its own establishment means" below). Only accepted species-rank entries are alien taxa: the 66,316 higher ranks (families 32,434, orders 15,608, `other` 5,761, classes
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

**Sub-national GRIIS lists on the grid** (placement by centre superseded on 2026-11-01 by "Sub-national GRIIS lists on the cells they overlap" below). `admin1_cells(grid, admin1, subdivisions)` places a list's subdivision on Natural Earth's admin-1
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


### Cleaning, layers 1 to 4

`src/realms/clean.py` (`clean`, `Settings`), `src/realms/clean_report.py` (`removal_summary`, the provisional twins) and `scripts/clean_counts.py`. `clean` takes the country counts of half A and
returns `Cleaned`: `presence` (species and part, summed over countries, with a `records` column), `removed` (one row for every row of the counts
that a layer took out or moved: `layer`, `rule`, species, part, country, `records`, `to_part`, and `presence_lost`), `unplaced` and `sides` (layer 4's counts of
species). Nothing is removed without a layer and a rule, `RULES` says which layer each rule is of and a rule that is not in it raises, and each can be switched off: a layer with `layers`, one rule with `Settings.rules_off`
(`extinct`, `no_wild_range`, `griis`, `griis_subnational`, `wrims`), WCVP with `wcvp_mode`. Layer 5 was accepted in `layers` and did nothing until
Task 12, which added it (see "Layer 5" below). All thresholds are starting values until the freeze.

**Records.** `records` is `n - n_human_text` with layer 1 on and `n_raw` with it off. Layer 1 takes out a row with no such record: rule
`failed_flags` if `n` is 0 (every record failed the flag checks) and `human_text` otherwise. An emptied row is logged with its `n_raw`, as that is what it
took with it. A row that loses some records and keeps others is not in `removed`: `n_raw - records` in `presence` is what it lost.
**The unit of `removed.records` is therefore `n_raw` for layer 1 and the row's `records` for every other layer** (the clean records `n - n_human_text`, or `n_raw` with
layer 1 off), so it must not be summed across layer 1 and the rest. From layer 2 on records are conserved: `presence.records` plus the `records` of the removals
of layers 2 to 4 (a null `to_part`; layers 2 to 5 since layer 5 exists) are the input's `n - n_human_text` with layer 1 on, **and its `n_raw` with layer 1 off** (the first version of this
sentence gave `n - n_human_text` for both), which `test_records_are_conserved` states for every WCVP mode and `side_conflict`, with layer 1 on and off. The other counts
(datasets, years, licence and source counts) are kept as the download has them, with the distinct ones taking the maximum over rows, as `load_presence` does.

**Zero-area parts: layer 0, rule `zero_area_part`.** Records in a part with an area below `ZERO_KM2` fold into the other part of the same cell, first of
all and with every layer off. It is not switchable: a zero-area part cannot carry a realm (`excess_table` divides by its area), so this repairs the data
and is no cleaning. A cell whose two parts are both below `ZERO_KM2` keeps its records, and a row folded into a part that holds none of the
species' records creates the row there. Rows with no records at all are moved but not logged.

**Layer 2.** `extinct` (COL's flag) and `no_wild_range` (`NO_WILD_RANGE`), each its own rule.

**Layer 3: order of the rules.** `griis`, then `griis_subnational`, `wrims`, and WCVP last. A row goes with the first rule that catches it, so
that `removed` never counts a row twice. (Since the final review, a row can be logged twice by one layer 3 rule: once, and once more for the records a fold of
layer 4 brought back to it. See "Layer 3 judges the rows that layer 4 moved" below.) WCVP is last on purpose: what `wcvp_*` reports is what no alien list had already removed, which is
the number that matters for the worry below. To measure one source alone, switch the others off. The runs of M1 ("layer 2 alone", "GRIIS alone") are
`layers={1, 2}` and `layers={1, 3}` with `wcvp_mode="off"` and `rules_off={"wrims", "griis_subnational"}`.

**WCVP, three modes, and why.** WCVP's native regions come from regional floras, and a flora stops at a TDWG region's border. Cutting a plant back to
them may predefine realm edges along TDWG regions, and a map drawn from such data could then "find" those edges. Layer 3 is the only rule that does it, so
it has three modes. Which one the analysis uses is decided in M1 by the circularity gate (the share of realm boundary within one cell of a TDWG level 3
boundary against rotated boundaries, with WCVP off, `introduced` and `native_only`), not by endemism excess, which trimming raises by construction.

- `native_only`, the design's rule and the default (`wcvp_outside_native`): a plant with at least one native region loses the parts touching none of them.
- `introduced` (`wcvp_introduced`): a plant loses the parts that touch an introduced region and no native one. Parts WCVP says nothing about stay,
  so documented introductions go without ranges being cut back to the floras' polygons.
- `off`: WCVP removes nothing, and GRIIS and WRiMS still apply.

Rulings that go with them, which Markus can revert in `_wcvp`:

- **Plants only:** `kingdom` is Plantae in `refs.taxa`. A plant WCVP does not list is not touched, and `native_only` does not touch a plant with
  introduced regions only.
- **A part touches a region** when it is in `tdwg3_parts` for it. A cell on a border touches both regions, and a part touching a native region counts as
  native whatever else it touches.
- **Sea parts have no TDWG region, so WCVP never removes a sea part.**
- **Nor does it remove a land part that no region reaches** (80 land parts at resolution 6: specks of coast and ice the scheme leaves out). The
  brief reads "parts that touch none of its native regions", which would take them, but WCVP says nothing there rather than that the plant is not native. This is
  an extension of the ruling above.

**GRIIS and sub-national lists.** A country's list (subdivision null) applies to rows by `countrycode`, rule `griis`, **except in the cells of that
country's sub-national lists**, which apply there by cell (`admin1_cells`), whatever the row's country code, rule `griis_subnational`. So the
country's list applies to the rest of the country, as the brief says, and a species on the conterminous United States' list is not removed from
Hawaii's cells, where Hawaii's list rules. Overseas lists carry the sovereign's code (`FR-Reunion`) while GBIF's records say `RE`, which is why they match on
cells and never on `countrycode`. Both rules apply to land and sea parts. **Superseded by Task 11a** ("Sub-national GRIIS lists on the cells they overlap", below): a
sub-national list now applies to the rows in every cell that overlaps its subdivision, and **where GBIF uses the list's country code in those cells, only to the rows with
that code**; the cells alone decide only where the code does not occur there (the overseas lists of `FR`). "Whatever the row's country code" and "never on `countrycode`"
describe the first version. The separate name `griis_subnational` is so that the sub-national lists' effect
can be told from the country lists'.

**What Task 10 could not place is never applied** (`griis_no_cells`, `wrims_no_country`) and is listed in `Cleaned.unplaced` with the records of each
species, for the species that have records: pairs at the level of a list or area. `removal_summary` puts them in as rows `unplaced_griis` and
`unplaced_wrims`, with a boolean column `unplaced` True so that a consumer filters on it (`~unplaced` for the removals); the rule names are kept. This replaces the
`applied` column the first version had, which said the same the other way round. `rows` are the pairs, and `records` the species' records in all, which is only an
upper bound on the gap. No pair of the gap is removed, and nothing is counted for a rule that is off.

**WRiMS** removes a row when the species is introduced for the row's `countrycode` (the country the area name gives) and the part is a sea part, rule `wrims`.
It rests on the Denmark check: GBIF sets `countrycode` on sea records inside EEZs.

**Layer 4.** A land species is exactly `terrestrial` in COL's environment, or has at least `land_share` of its records on land parts, and a sea species is the
same the other way round, with the same share (`land_share` is symmetrical, which is why it must be above one half). The share is taken after layers 1 to 3.
COL's environment is a sorted, comma-separated list, so "exactly terrestrial" is the value `terrestrial` alone, and a species with freshwater or brackish
in it is exempt, as is all of Aves. A row on the wrong side folds into the other part of its cell if that has area (`land_sea_fold`), or else into the
right-side part of the neighbouring cell that holds most of the species' records, the lowest cell id on a tie (`land_sea_neighbour`), or else is dropped
(`land_sea_dropped`). Three points the brief leaves open:

- **A neighbour counts only if the species has records on the part.** A neighbour with the right kind of part and none of the species' records is
  not one: folding there would invent a presence where the species was never seen, so the row is dropped. **"Most records" is what the neighbour holds
  before this pass moves anything**, so the records that the same pass folds into a part within its own cell do not count towards it. It could change which
  neighbour takes a row, and it is a choice of simplicity (one pass, no order among species' rows). It touches only `land_sea_neighbour`, which is
  20,781 records with `native_only` (20,945 with `introduced`, 21,178 with `off`), 0.001% of the clean records.
- **The lowest id is the lowest as a number** (ids are long integers kept as text, and `9` comes before `100`).
- **Layer 3 judges again every row that layer 4 moved, where it went** (added after the final review, 2026-11-01, see "Layer 3 judges the rows that layer 4 moved" below):
  before that, a fold could put back records that layer 3 had taken out of the part it folds into.
- **Where COL's environment and the records disagree, the environment decides** (`Settings.side_conflict="environment"`; `"leave"` leaves the species alone).
  The brief's "or" lets a species be both: of half A's species, 1,536 are exactly terrestrial with at least 95% of their records in sea parts, and 4,729 are
  exactly marine with that on land parts. They are mostly rare (the median is 2 records) shore molluscs, crustaceans and insects, whose records fall in a
  coastal cell's land part because GADM draws country polygons into the shallows, and a land/sea label that rests on GADM is the weaker word. The effect on the
  map is small either way: 183 species lose every presence with `environment`, none with `leave`. Markus can switch it.

### Layers 1 to 4 on half A at resolution 6 (2026-11-01)

`scripts/clean_counts.py` on `data/gbif/0019721-261023074918325`, half A, every layer on, 3 modes in one process (55 s, peak 18.9 GB: the references
peak at 15.5 GB, mostly `load_taxa`, and the counts' load adds to what they leave). Outputs in `results/clean/res6/<mode>/` (git-ignored), with the log in `results/clean/res6/clean.log` and the whole
report in `report.json`. **The counts are of COL 315557 (COL26.6 XR).** The input is 10,126,245 rows of 647,326 species and 1,909,146,347 clean records
(`n - n_human_text`), in 8,919,665 species-and-part presences.

| Layer | Rule | Rows | Presences | Records | Share of clean records |
|---|---|---|---|---|---|
| 0 | `zero_area_part` (moved) | 216,241 | 199,419 | 1,255,754 | 0.07% |
| 1 | `failed_flags` | 166,251 | 153,218 | 1,188,700 (`n_raw`) | |
| 1 | `human_text` | 26,852 | 25,752 | 52,328 (`n_raw`) | |
| 2 | `extinct` | 58,303 | 57,416 | 622,394 | 0.03% |
| 2 | `no_wild_range` | 8,036 | 7,052 | 628,291 | 0.03% |
| 3 | `griis` | 325,520 | 262,907 | 84,405,810 | 4.42% |
| 3 | `griis_subnational` | 78,247 | 76,318 | 47,923,761 | 2.51% |
| 3 | `wrims` | 3,125 | 2,930 | 259,332 | 0.01% |
| 3 | `wcvp_outside_native` (`native_only`) | 355,511 | 338,112 | 7,176,959 | 0.38% |
| 4 | `land_sea_fold` (moved) | 794,698 | 744,151 | 9,954,205 | 0.52% |
| 4 | `land_sea_neighbour` (moved) | 5,627 | 5,574 | 20,781 | 0.00% |
| 4 | `land_sea_dropped` | 12,989 | 12,874 | 37,974 | 0.00% |

The rows of layers 0 to 2 and of the three alien-list rules are the same in all three modes, and layer 4's differ by under 2% (it follows what layer 3 leaves). Presences are the
(species, part) pairs left with nothing by the rule, so a country's presences do not add up (use `removal_summary(by=())`), and a fold's presence
reappears at the part it went to. Per layer, `native_only` lost 199,419 (0), 178,943 (1), 64,468 (2), 680,267 (3) and 762,599 (4) presences.
**The result is 7,391,207 presences of 611,163 species and 1,768,091,826 records** (of 8,919,665 presences, 83%).
Records are conserved: 1,768,091,826 left + 1,250,685 (layer 2) + 139,765,862 (layer 3) + 37,974 (layer 4, dropped) = 1,909,146,347, the clean records of the input,
as the folds move records and the rows of layers 0 and 1 carry none of that count.

**Layer 3 is mostly records, not presences, and mostly GRIIS.** The two GRIIS rules take 132.3M records (6.9% of the clean records) in 403,767 rows
and 339,225 presences (3.8%). Of the country lists' records, 51.0M are chordates and 28.9M angiosperms, and AU (15.2M), GB (13.1M), SE (10.4M) and CA (9.0M) are the
largest countries. The sub-national lists' 47.9M are almost all `US` (47.6M), and 42.2M of the 47.9M are chordates: *Sturnus vulgaris* (15.5M), *Passer
domesticus* (12.5M), *Columba livia* (6.3M) and *Streptopelia decaocto* (3.7M) lead. **The gap the lists leave:** 4,027 unplaced GRIIS pairs and 1,140 unplaced WRiMS pairs name species with records in half A, of 7,468 and 1,995
in all.

**The effect of each WCVP mode** (plants, `kingdom` Plantae, are 169,163 species and 2,568,687 presences of half A with WCVP off):

| `wcvp_mode` | Rule rows | Rule presences | Rule records | Presences left (all species) | Plant presences left | Species left |
|---|---|---|---|---|---|---|
| `off` | 0 | 0 | 0 | 7,704,561 | 2,568,687 | 612,687 |
| `introduced` | 110,910 | 99,962 | 5,271,123 | 7,620,169 (-1.1%) | 2,484,295 (-3.3%) | 612,648 (-39) |
| `native_only` | 355,511 | 338,112 | 7,176,959 | 7,391,207 (-4.1%) | 2,255,333 (-12.2%) | 611,163 (-1,524) |

`native_only` removes 313,354 plant presences net, 3.7 times as many as `introduced`, but only 36% more records (7.18M to 5.27M): a presence `introduced` takes
has 62 records on average, and each of the 228,962 more that only `native_only` takes has 8. Those are parts that WCVP lists neither as native nor as introduced for
the plant (it removes from 50,866 species, `introduced` from 5,063), which is where the incompleteness of WCVP and the edges of its regions would act, and it is
what the circularity gate in M1 tests. The three modes agree on every other rule.

**The zero-area fold.** 923 parts (890 sea, 33 land: the pilot's figures, and the same in every mode), 216,241 rows, 1,255,754 clean records and
199,419 presences. Sea: 890 parts, 215,089 rows, 1,243,297 records of 46,370 species. Land: 33 parts, 1,152 rows, 12,457 records of 827 species. (Task 7's 2.3M
and 24.6k are over all species, with half B.) None of the 8,624 parts of the result has an area below `ZERO_KM2`, and no cell of the grid has both parts below it,
so no record is left without a part to go to.

**The Caspian.** Natural Earth does give its cells sea area: 13 cells in a box around it (lon 44.5 to 56.5, lat 35 to 49.5, on the cells' centres) have 1,000 km2 or
more of sea, 10 of them 9,269 to 66,317 km2 and 3 of 2,299 to 3,616 km2, and 3 more with 34 to 393 km2 are left out. So the fold of zero-area parts never touches
them (0 records folded), and the records GBIF puts in their sea parts stay there as far as that fold goes. There are few. **On one measure throughout, the clean
records `n - n_human_text`:** 8,160 of the 622,268 in these cells (1.3%) are in the sea parts and 614,108 in the land parts, and the cells that are mostly water hold their
records in the land part of the shore (cell 644184800048851186: 3,900 records on its 3,803 km2 of land and 3 on its 66,317 km2 of sea). After all layers 1,755 records are in
these sea parts with `native_only` (1,656 with `introduced`, 1,617 with `off`), against the 8,160 that came in. (The first version of this paragraph set the input's `n`
against the output's `records`, and gave 8,161, 622,663 and 614,502.) So the Caspian is a sea part in the network, as intended, and a sparse one, which the sparse-part rules
of M1 will have to judge.

**Provisionally accepted twins** (reported, nothing changed). 33,811 provisionally accepted species ids share a name and kingdom with an accepted species
(33,789 with one accepted twin, which Task 10 joins to, and 22 with several, which it leaves ambiguous). Half A has records on 2,835 of them: 59,773
rows and 2,728,658 clean (`n`) records, 0.14% of half A's; 2,789 species and 58,565 rows pass layer 1. None of them has several accepted twins.
So a GRIIS, WCVP or WRiMS record on such a name reaches the accepted twin and not the provisional id, whose records are then not cleaned by it.

**Layer 4 in numbers.** After layer 3, `native_only`: of 611,346 species 454,559 are land, 80,916 sea, 27,543 exempt (freshwater, brackish or Aves) and the rest have
no side. 6,265 species were in the conflict above. The folds move 9.95M records (0.52%), of 744,151 presences, and drop 37,974 records of 12,874 presences,
which end 183 species' presences altogether.

**GRIIS lists whole countries, and some of its largest removals are in the species' native range (observed, nothing changed).** A country's list removes a
species from every row of that country. Among the biggest single removals are species that are native there: *Gymnorhina tibicen*, the Australian Magpie, in AU (2.8M records),
*Anas platyrhynchos*, the Mallard, in SE (1.5M) and DK (0.8M), and *Columba livia* in IN (1.1M), where the Rock Dove is wild, while in CA (1.7M) and GB it is feral. Presumably GRIIS lists a
species wherever some population is introduced (released Mallards, feral pigeons), which the rule turns into the whole country. It is the rule as the brief gives it, and the
design's softening is the sub-national lists, which help only where they exist. Nothing in `References` says where in a country a species is native, so `clean` cannot soften this.
GRIIS's effect on its own is one of the runs of M1, and these are the first rows to look at there.

### Layer 5, the rare-outlier heuristics (2026-11-01)

`src/realms/outliers.py` (`core_parts`, `detached_satellites`, `under_reported`, `thin_long_range`, `flag_outliers`, `suspect_cells`), called by `clean` as layer 5. The rules
`detached_satellite`, `under_reported` and `thin_long_range` are in `RULES` at layer 5 and in `SWITCHABLE`, so `rules_off` switches each off, and `layers` the whole layer. The thresholds are
the `Settings` of Task 11, at the design's values, and all starting values until the freeze. Records stay conserved from layer 2 to layer 5 (`test_records_are_conserved` now also
runs with layer 5 on, at the design's thresholds and at looser ones, as the random tables of 8 cells are too small to meet the design's). Rulings, with the reason for each:

- **Cells, not parts.** A species' records in a cell are summed over its parts and countries, the rules judge that sum, and a flagged (species, cell) takes every row of the species there under the
  rule's name. The design speaks of cells and records ("a gap of 500 km", "5% of its cells", "a cell with at least 1,000 records of the group"), and a coastal cell's two parts would otherwise be
  judged apart, with half the records each, which makes every shore species look half as common.
- **Distances are great-circle arcs between cell centres, in kilometres**, never degrees. **The gap is `satellite_gap_km` plus the grid's median neighbour spacing**, measured from the grid's
  adjacency (292.5 km at resolution 6, 160.7 km at 7, so a gap of 792.5 and 660.7 km). The design's 500 km is a gap between ranges and cells are joined by their centres, so the cell's own width is
  added, and nothing is tied to a resolution. Test: the same 667 km between a patch and the range is detached on a grid of 1-degree cells and joined on one of 3-degree cells.
- **Core cells.** A species' core is its top cells by records, up to and including the cell that takes the running total to `core_share` (90%) of its records, and **every cell under
  `heuristic_min_records` (20) records.** Cells with equal records are ordered by cell id (as text), so **among tied cells at the line, which are protected is reproducible but arbitrary.** (A first
  version extended the core to every cell tied with the smallest core cell. A long tail of equal small cells that crosses the line, such as singletons, then became core whole: 1,000 records in 17 cells of
  50 and 150 singletons made all 167 cells core where the rule gives 67, and neither rule could touch the tail. Reverted, with a test for that case.)
  `detached_satellite` and `under_reported` never remove a core cell (and never a patch that holds one). **`thin_long_range` is the exception and removes a whole species, core cells and
  all**: the brief says it takes the species out of the network, and the design says such a species "stays out of the network". A species with under 5 records is all core, so the guardrail
  would otherwise make the rule do nothing.
- **Every rule judges the rows layers 1 to 4 left, and not one another's removals.** The core is taken once, from what layer 4 leaves, and a cell that two rules flag goes with the first of
  `thin_long_range`, `detached_satellite`, `under_reported`. What a rule flags then does not depend on the order, and switching a rule off hands its cells to the next rule that flags them
  (not back). The thin rule's species have under 5 records, so the others never flag them unless the thresholds are changed.
- **Strictness follows the design's words.** "Under" is strict: a patch with under 1% of the records and under 5% of the cells, a species with under 5 records, cells over 1,000 km apart. "At most" and
  "at least" are not: at most 5 records in the cell, at least 1,000 records of the group, at least 20 records, at least 2 datasets and 3 years. A share against a count is compared to 1e-9
  (0.05 times 60 cells is 3.0000000000000004 in floating point, and 3 cells of 60 are not under 5%).
- **Satellite patches.** Only species with at least 20 records, on a graph of the grid's cells within the gap of each other (a k-d tree on unit vectors), and a species' patches are the connected
  components of its occupied cells, by `scipy.sparse.csgraph`, for all species at once. A patch goes when it holds under 1% of the records and under 5% of the cells, unless it has at least 2
  datasets and 3 years. **A patch's datasets and years are the largest over its cells (and so over its parts and countries), which is a lower bound**, as `load_presence`'s counts are: distinct
  counts do not add up across rows. It leans towards removing. **`n_datasets` and `n_years` count every record that passes layer 1's flag checks, so they include those that only the text patterns
  caught (`n_human_text`), and the fold of zero-area parts can bring in a record-less row's counts (maximum of the two). Both can only raise a patch's counts, which leans towards keeping.** The two
  biases are not corrected, and which is larger is not known.
- **Reporting rate.** `taxGroup` is the group (`refs.groups`, 37 groups). The group total in a cell is over all the species of the group there, small ones included, in what layers 1 to 4 left; the
  species' median rate is over its core cells, whatever their group totals. A species with no group is not judged by this rule (1,434 species and 25,827 records of half A). Its cell goes if the
  group has at least 1,000 records in the cell, the species at most 5 there and at least 20 in all, and a rate under 1% of its median core rate.
- **`suspect_cells`** lists and removes nothing. A cell's share is the species layer 5 took out of it over those it had before (what `presence` holds plus what layer 5 removed there), and it is
  suspect when the share is at least 3 times the median of its neighbours' shares (the cells in `grid.adjacency`) and at least 10 percentage points above it. **A neighbour with no species has no share and
  is left out of the median, and a cell with no neighbour that has one is not judged.** It has no floor on the number of species, as the design has none; see the result below for what that does to the
  top of the list.

Layer 5 is fast: all the patches of 252,498 species are found in one pass (no loop over species), in about 3.5 s at resolution 6, and the layer takes about 7.5 s and 2 to 3.5 GB on top of what layers 1 to 4 hold
(see below). At resolution 7 the grid has 21,872 cells and 59 links per cell within the gap, against 7,292 and 27.9, so with 2.5 times the rows the patch search grows about 5 times (an estimate, not
a measurement on resolution 7 counts).

### Layer 5 on half A at resolution 6 (2026-11-01)

`scripts/clean_counts.py` on `data/gbif/0019721-261023074918325`, half A, all five layers, 3 modes in one process: 76 s, peak 18.9 GB (the peak is decided by the references' load and has been 15.5 to 18.9 GB over runs
of the same code). `clean` takes 14.7, 14.1 and 13.3 s per mode where it took 8.9, 7.7 and 7.0, so **layer 5 costs about 5 to 7 s per mode**. Measured on its own (before the correction of the core
below, which changes only the core step): 8.3 s without layer 5 and 15.3 s with it on the same inputs, and with it on `clean`'s resident set did not rise above the peak of layers 1 to 4 (14.4 against 14.5 GB, from a
7 GB start). `flag_outliers` alone on the 7.4M rows of layers 1 to 4 takes 7.5 s and 1.8 to 3.5 GB above what is resident. Far under the 10 minutes and 30 GB set. Outputs as for layers 1 to 4, in
`results/clean/res6/<mode>/`, plus `suspect_cells.csv`.

**Layers 1 to 4 did not change.** Against a copy of the outputs of the run before layer 5 existed, in all three modes: the rows of `removed` with a layer below 5 are identical (2,051,400, 1,811,857 and
1,708,998 rows), the presence of the old run is the new one plus what layer 5 removed, row for row and record for record, the unplaced tables and the summaries of layers 0 to 4 (per rule, per rule and
group and country) are identical, and so are layer 4's species counts, the zero-area fold and the input counts in `report.json`.

**The core was corrected after review, and these are the figures after it.** The first version extended the core to every cell tied with the smallest core cell, which made a long tail of equal small cells
that crosses the line (singletons) core whole and out of reach of the satellite and under-reported rules. With the literal rule layer 5 takes more: `native_only` 492,302 rows against 462,523 (satellite 239,486
rows against 216,750, under-reported 197,606 against 190,563, thin unchanged at 55,210, as all of a thin species' cells are core either way), 858,568 records against 818,980, and 498 suspect cells against 550.

**What layer 5 took with `native_only`** (of 7,391,207 presences, 611,163 species and 1,768,091,826 records that layers 1 to 4 left):

| Rule | Rows | Presences | Records | Species | (species, cell) pairs |
|---|---|---|---|---|---|
| `detached_satellite` | 239,486 | 237,819 | 470,530 | 45,439 | 235,767 |
| `under_reported` | 197,606 | 192,396 | 320,939 | 55,677 | 188,922 |
| `thin_long_range` | 55,210 | 55,114 | 67,099 | 22,456 | 54,946 |
| all of layer 5 | 492,302 | 485,329 (6.6%) | 858,568 (0.049%) | 101,348 (16.6%) | 479,635 |

**The result is 6,905,878 presences of 588,707 species and 1,767,233,258 records.** Records are conserved: 1,767,233,258 left + 858,568 (layer 5) = 1,768,091,826, what layers 1 to 4 left. The records
of a removed (species, cell) are few: the median of a satellite is 1 record and the mean 2.0 (largest 1,341). So layer 5 removes presences, which are what the network is built from, and hardly any
records. **`thin_long_range` takes 22,456 species out whole** (3.7% of the species, 10% of the 218,380 that have under 5 records), and they are the only species layer 5 leaves with no presence at all: a species
with 20 records keeps its core. In the other modes (rows, presences, records, species, species taken out by the thin rule, suspect cells, and what is left):

| `wcvp_mode` | Rows | Presences | Records | Species | Thin species | Suspect cells | Left: presences | species | records |
|---|---|---|---|---|---|---|---|---|---|
| `native_only` | 492,302 | 485,329 | 858,568 | 101,348 | 22,456 | 498 | 6,905,878 | 588,707 | 1,767,233,258 |
| `introduced` | 556,109 | 548,582 | 955,875 | 109,105 | 23,311 | 506 | 7,071,587 | 589,337 | 1,769,041,809 |
| `off` | 565,667 | 557,749 | 975,460 | 109,267 | 23,300 | 512 | 7,146,812 | 589,387 | 1,774,293,399 |

`introduced` and `off` have 13% and 15% more rows than `native_only`: presumably the plant presences that WCVP's native regions would have cut are caught by layer 5 instead, in part (by rule, `off` has
279,596 satellite rows against 239,486, 228,794 under-reported rows against 197,606, and 57,277 thin ones against 55,210). The groups are in `summary_by_group.csv`.

**Suspect cells: 498 with `native_only`.** The 10 with the highest share (the order the design asks for) are all near-empty cells of open ocean, where one to seven species in all were removed (centre, removed of
species, the neighbours' median share):

| Cell | Centre (lat, lon) | Removed of species | Neighbours' median |
|---|---|---|---|
| 613444144062176855 | 13.4, 62.2 | 7 of 7 | 25.1% |
| 616349856157052712 | 16.3, 157.1 | 4 of 4 | 32.9% |
| 7227605003122292665 | -27.6, -122.3 | 3 of 3 | 8.3% |
| 2817244927072977949 | -17.2, 73.0 | 2 of 2 | 0% |
| 7220778330100356614 | -20.8, -100.4 | 2 of 2 | 25.5% |
| 2806606672077933218 | -6.6, 77.9 | 1 of 1 | 7.7% |
| 2811415353082192981 | -11.4, 82.2 | 1 of 1 | 0% |
| 2817019090077513526 | -17.0, 77.5 | 1 of 1 | 0% |
| 2817664893084291419 | -17.7, 84.3 | 1 of 1 | 14.2% |
| 2820296930084253319 | -20.3, 84.3 | 1 of 1 | 20.0% |

So the list by share is dominated by cells with so few species that a share means little (the median suspect cell has 17 species, and 465 of the 498 have at least 2, 311 at least 10, 129 at least 50). The design sets no
floor and none was added; a floor on the number of species is the natural fix, and Markus can set it. **Ranked by the species removed, the top 10 are real, well-sampled places** (centre, removed of species, share, neighbours' median):
24.5, 120.6 (2,132 of 12,722, 16.8%, 5.4%), 42.8, 132.6 (931 of 5,385, 17.3%, 4.4%), 32.0, -64.1 (930 of 3,100, 30.0%, 3.6%), 32.3, -17.4 (853 of 3,249, 26.3%, 6.8%), 79.2, 11.3 (629 of 2,968, 21.2%, 6.5%),
52.1, 105.7 (618 of 2,558, 24.2%, 6.6%), 59.1, -94.7 (591 of 3,365, 17.6%, 5.8%), 59.9, 68.3 (546 of 3,243, 16.8%, 4.3%), 64.7, -16.5 (542 of 2,353, 23.0%, 5.4%) and -25.2, 134.8 (448 of 2,916, 15.4%, 3.9%).
By their centres they are near Taiwan, Primorsky Krai, Bermuda, Madeira, Svalbard, Lake Baikal, Hudson Bay, western Siberia, Iceland and central Australia. **Layer 5 takes 15 to 30% of the species of these
cells, and 48 to 92% of what it takes there is the satellite rule** (Bermuda 854 of 930, Madeira 725 of 853, Svalbard 535 of 629, Iceland 478 of 542; the Taiwan cell is the exception at 1,024 of 2,132, with 861 under-reported and 247 thin).
That is what the rules are built for, a species with a single record or dataset far from its range or in a cell where its group is heavily sampled, and it is also where real disjunct occurrences on oceanic islands
are, which realm boundaries depend on. The design's guard is 2 datasets and 3 years, and a single-source island record is removed. **Whether that over-cleans the islands is not known** and it is a calibration
question for M1 (the thresholds are starting values), which the cells of this list are the first place to look at. Nine of the top 10 by share and eight of the top 10 by species removed are the same as `native_only`'s
in both other modes.

### Calibration of layer 5 on WCVP plants, and the known cases (2026-11-01)

`src/realms/calibrate.py` (`native_status`, `tally`, `rates`, `settings_grid`, `baseline`, `calibrate`, `readings`, `choose`, `score`, `kingdom_counts`, `region_loss`, `top_regions`, the known cases `CASES`, `evaluate_case`, `known_cases`) and `scripts/calibrate.py`. **No threshold
changed:** the `Settings` defaults are the design's, and the setting chosen below is a dry run at resolution 6. The binding calibration is repeated at resolution 7, where it is also the one scoring of half B. Rulings, with the reason for each:

- **What the calibration measures.** With layer 3 off, layer 5 has to catch the plant presences outside WCVP's native regions by its own evidence. `native_status` marks each plant presence (species, part) `native` where the land part touches
  one of the species' native TDWG level 3 regions and `outside` where it touches none, as layer 3's `native_only` does (a part touching a native and an introduced region is native). Only plants with a native region and only land parts with a TDWG region get a status.
  `rates` are the share of outside presences removed and the share of native presences lost, **counted in presences** (species, part), which is what the network is built from. A test checks that `outside` is exactly what layer 3 takes in `native_only`.
- **Three readings of the 5% limit, and the one that chooses.** The thresholds only change layer 5, so layers 1, 2 and 4 run once and layer 5 runs per setting. What counts as lost can be read three ways, and `readings` gives all three:
  `layer5` (the presences layer 5 takes of those layers 1, 2 and 4 leave: what the thresholds do and nothing else); `run_all` (the whole of layers 1, 2, 4 and 5, over the raw presences, `n_raw` above 0, counting every raw presence that is not in the result);
  and `run_drops` (the same, but **only drops count**: a presence is lost if it is a raw presence, absent after cleaning, and none of its rows was moved by a fold, that is layer 0's zero-area fold or layer 4, a removal-log row with a `to_part`.
  A moved presence is moved, not lost. The shares are over all the raw presences of the kind, outside removed counted the same way). **The plan words the limit over the whole run, so the controller ruled `run_drops` (provisional until Markus confirms it before the resolution 7 run).**
  `choose` reads it by default, with the same ties (lower native loss, then the grid order) and the same fallback, and `calibrate` has its columns as `run_*`. The first version chose on `layer5`, and its choice is reported below for comparison.
  A test builds a native presence that layer 4 folds and one that layer 1 drops, and checks that the first is not lost and the second is, and the columns of `calibrate` are checked against `clean`'s own removal log for all 81 settings.
- **The seam.** `clean.outlier_rows` is the public seam: `_layer_5` calls it, so the calibration runs layer 5 per setting with the same code and not a copy. **Checked on the real counts:** `calibrate`'s counts for settings 0, 40 and 80 are identical to `clean` with layers 1, 2, 4 and 5 at that
  setting, and the synthetic test does it for all 81, with the rows of two countries for one satellite part and one part of a range (`clean` judges rows by country, `calibrate` presences). `Baseline` records the settings it was made with, and `calibrate(..., start)` refuses a `start` whose `land_share`,
  `side_conflict` or layer 2 rules differ from the grid's. **Those rules are derived from `RULES`** (`EARLY_RULES`: the switchable rules of layers 1, 2 and 4, which is `extinct` and `no_wild_range`), and not from a hand-written set, so a new switchable rule of those layers is covered.
- **The species** are the plants WCVP gives a native region and **every other species of their taxonomic groups**: `under_reported` compares a species with its group's records in a cell, so the groups are kept whole and the group totals are those of the real run. **Not all of those are plants.**
  Counted from the split, `refs.taxa`, `refs.groups` and `refs.wcvp` alone (no occurrences were read), the calibration set has 197,509 species in half A (194,975 Plantae, 2,534 not) and **197,567 in half B (194,981 Plantae, 2,586 not Plantae)**, 395,076 in all (5,120 not).
  All the species that are not Plantae are in the group `eukaryotes`; half B's are 1,411 Protozoa, 504 Animalia, 599 with no kingdom in COL (two of them WCVP plants), 67 Chromista and 5 Fungi. The group `algae`, which holds Chromista and Protozoa too, is not in the set, as no WCVP plant is in it.
  These are species that the split and the taxonomy give, whether they have records or not. In the dry run, 161,145 species of half A have records. **Checked:** at the design's setting, `clean` on all of half A (all species, layers 1, 2, 4 and 5) gives the same 525,201 / 107,326 / 1,808,902 / 56,953 as the calibration set alone.
- **The grid of 81 settings** is as the brief has it: satellite record share 0.5, 1 and 2%, satellite cell share 2, 5 and 10%, under-reported rate factor 0.5, 1 and 2%, group minimum 300, 1,000 and 3,000, in that order (the first varies slowest), with the design's setting in the middle (index 40).
- **The choice** is the setting that removes the most outside presences while losing at most 5% of the native ones (`choose`, a tolerance of 1e-9 on the share against the count), by the reading above. With no setting within 5%, the lowest native loss is taken and the run says so.
- **Half B is scored behind `--score-half-b`, and was not scored.** (Since the final review a run of its own, `--score-half-b --setting N`, after the dry run at resolution 7: see the call-site list above and
  "The scoring of half B, a run of its own" below.) `score_held_out` is implemented and tested on a made-up download and hand-made tables only. It is listed in the call-site list of the seal above. It scores half B's **calibration species, which are the plants and the other species of their groups**
  (the species that are not Plantae are kept because the reporting rate needs the group totals, so that the score is like for like with A's; the numbers are above). The flag refuses to run at any resolution but 7, and refuses if any `results/calibrate/res*/held_out.started` or `res*/held_out.json` exists (tested, before anything is loaded). The marker is written before the load and the result after it, so a run that fails in between leaves the marker and a retry is refused (tested with a loader that fails); it says that a retry needs Markus's decision, and there is no override flag.
  It **prints and saves counts and rates only, and how many species were loaded and how many are not Plantae: nothing per taxon** (the function returns exactly that, and a test and an assertion in the script check the keys). **No unsealed call was run in this task, and nothing of half B was loaded, counted or looked at:** the script loads
  `load_country_counts(..., "A")` only, and the halves of the known cases and the counts of calibration species come from the split file and the taxonomy, which have no occurrences.
- **The animal cross-check against GRIIS is pending, until Task 11a** (done in "Task 11a on half A at resolution 6" below). The GRIIS lists are name-only for now and include species that are native where listed (5.8% of the entries that can be checked), so a check against them would count native presences as errors.
- **Known cases** (`CASES`) run on half A fully cleaned (all five layers, `wcvp_mode="native_only"`, the design's thresholds). **Areas are TDWG level 3 regions, and not country codes:** an area is the set of cells whose land touches a region (or the regions of a level 1 region), and
  both parts of such a cell are in it. The criteria were set before the run and are not tuned: tiger, *Panthera tigris*, no presence in Europe (level 1 region 1) and some left elsewhere; potato, *Solanum tuberosum*, no presence outside its WCVP native regions and some inside;
  starling, *Sturnus vulgaris*, no presence in North America (level 1 region 7) and some left elsewhere, with what removed it listed by rule; a Hawaiian endemic bird, *Drepanis coccinea* (COL's accepted name is *Vestiaria coccinea*, the iiwi, id 7G2WX, which holds the synonym's records),
  at least 95% (the native-loss limit) of its presences in the cells of HAW kept, of at least one; Lemuridae (24 species), no presence outside the cells touching MDG and some inside. A failure is reported as a failure.
  - **Waiting for step 9:** the tiger's family, Felidae, is in half B. The other four families are in half A (Solanaceae, Sturnidae, Fringillidae, Lemuridae).

### Calibration dry run and known cases on half A at resolution 6 (2026-11-01)

`scripts/calibrate.py --counts data/gbif/0019721-261023074918325 --resolution 6`, half A, outputs in `results/calibrate/res6/` (`settings.csv`, `regions.csv`, `known_cases.json`, `report.json`, `calibrate.log`). **This is a dry run: the binding calibration is
repeated at resolution 7 (`--resolution 7 --score-half-b`), and nothing below is a threshold decision.** 242 s in all and 19.0 GB at the peak (the references decide the peak, and it has varied from 15.5 to 19.0 GB between runs of this and the cleaning code): load 21 s, the known cases 19 s (a full clean of half A),
layers 1, 2 and 4 on the plants 3 s, and layer 5 at 81 settings 191 s (2.4 s a setting, with both readings), against the 30 minutes allowed. At resolution 7 the patch search grows about 5 times (see layer 5), so the 81 settings may take around 15 minutes (an estimate, not a measurement).

**The plants:** 2,557,367 presences after layers 1, 2 and 4, 2,334,103 with a status: 1,808,902 native and 525,201 outside (22.5%). The raw presences (`n_raw` above 0) with a status are 1,820,925 native and 565,525 outside.
(Since `native_status` counts kingdom Plantae only, after the final review, 2 native presences fewer in each: 2,334,101, 1,808,900 and 1,820,923. Nothing else of this section moves; see "Smaller fixes after the final review" below.)

**The 81 settings, by the whole run counting drops only.** Outside presences removed run from 21.8% to 31.8% and native presences lost from 3.1% to 6.2%. **54 of 81 settings lose at most 5%** (by layer 5 alone, 79 do). The 27 that do not are exactly those with a rate factor of 2%.
**The chosen setting is 75:** record share 2%, cell share 10%, rate factor 1%, group minimum 300. The top 10 within the limit (the last column is layer 5's reading of the same setting):

| Setting | Record share | Cell share | Rate factor | Group min | Outside removed | Rate | Native lost | Loss | Layer 5 alone: outside / native |
|---|---|---|---|---|---|---|---|---|---|
| **75 (chosen)** | 2% | 10% | 1% | 300 | 161,436 | 28.5% | 83,468 | 4.58% | 22.2% / 3.37% |
| 76 | 2% | 10% | 1% | 1,000 | 161,334 | 28.5% | 83,389 | 4.58% | 22.1% / 3.36% |
| 77 | 2% | 10% | 1% | 3,000 | 160,937 | 28.5% | 83,139 | 4.57% | 22.1% / 3.35% |
| 48 | 1% | 10% | 1% | 300 | 158,397 | 28.0% | 81,681 | 4.49% | 21.6% / 3.27% |
| 49 | 1% | 10% | 1% | 1,000 | 158,295 | 28.0% | 81,602 | 4.48% | 21.5% / 3.26% |
| 50 | 1% | 10% | 1% | 3,000 | 157,891 | 27.9% | 81,352 | 4.47% | 21.5% / 3.25% |
| 66 | 2% | 5% | 1% | 300 | 153,422 | 27.1% | 80,154 | 4.40% | 20.6% / 3.18% |
| 67 | 2% | 5% | 1% | 1,000 | 153,310 | 27.1% | 80,075 | 4.40% | 20.6% / 3.17% |
| 21 | 0.5% | 10% | 1% | 300 | 153,201 | 27.1% | 79,092 | 4.34% | 20.6% / 3.12% |
| 22 | 0.5% | 10% | 1% | 1,000 | 153,098 | 27.1% | 79,013 | 4.34% | 20.5% / 3.12% |
| *40 (the design's)* | 1% | 5% | 1% | 1,000 | 152,607 | 27.0% | 79,622 | 4.37% | 20.4% / 3.15% |
| *0 (the mildest)* | 0.5% | 2% | 0.5% | 300 | 123,558 | 21.8% | 57,332 | 3.15% | 14.8% / 1.89% |
| *80 (the first version's choice, over the limit now)* | 2% | 10% | 2% | 3,000 | 178,738 | 31.6% | 112,062 | 6.15% | 25.5% / 4.98% |

**The three readings, for the chosen setting, the design's, and the first version's choice** (outside presences removed, native presences lost; the whole-run shares are over the raw presences):

| Setting | Layer 5 alone | Whole run, counting everything | Whole run, counting drops only |
|---|---|---|---|
| **75 (chosen)** | 116,350 of 525,201 (22.15%); 60,890 of 1,808,902 (3.37%) | 161,450 of 565,525 (28.55%); 83,567 of 1,820,925 (4.59%) | 161,436 (28.55%); 83,468 (4.58%) |
| 40 (the design's) | 107,326 (20.44%); 56,953 (3.15%) | 152,621 (26.99%); 79,721 (4.38%) | 152,607 (26.99%); 79,622 (4.37%) |
| 80 (chosen by layer 5 alone before) | 133,954 (25.51%); 90,006 (4.98%) | 178,752 (31.61%); 112,161 (6.16%) | 178,738 (31.61%); 112,062 (6.15%) |

- **A moved presence hardly matters:** counting folds as lost adds 99 native presences to setting 75's whole-run loss (83,567 against 83,468) and 14 outside ones. The whole run loses about 1.2 percentage points of native presences before layer 5 (to layers 1, 2 and 4).
- **The choice moved from 80 to 75.** Under the whole-run reading the 5% limit binds on the rate factor: every setting with 2% is over it and every one with 1% or 0.5% is within it, so the rate factor is bracketed by the grid. **The record share and the cell share are not:** the chosen setting has the largest value of both, so a larger share would remove more, at a cost that the grid cannot show.
  Averaged over the other thresholds, the rate factor does most (0.5%, 1%, 2%: 23.8%, 26.5%, 29.9% of outside presences removed, at 3.34%, 4.31%, 5.93% of native lost), then the cell share (2%, 5%, 10%: 25.0%, 27.1%, 28.1%, at 4.34%, 4.56%, 4.68%), and the record share (0.5%, 1%, 2%: 26.3%, 26.8%, 27.0%, at 4.46%, 4.54%, 4.58%)
  and the group minimum (300, 1,000, 3,000: 26.7%, 26.7%, 26.6%) hardly matter.
- **Layer 5 alone removes at most about a quarter of the outside presences, and with layers 1, 2 and 4 the whole run about 28%, and costs 3.4% (layer 5) or 4.6% (whole run) of the native ones**, so it is no substitute for layer 3 where WCVP exists. `outside` is a lower bound on what is not native, as WCVP is coarse (TDWG level 3, cells of 265 km): a presence in a cell that touches a native region is
  native here, however far from the species' real range, and a garden record inside the range cannot be told apart.

**Native loss per TDWG level 3 region, at the chosen setting (75), by layer 5's reading** (a native presence counts in each region its part touches where the species is native, so the rows do not add up to the total; the whole run's loss per region was not made). 324 regions have at least 100 native presences, and the median loss is 2.2%. The 20 highest:

| # | Region | Level 1 | Land parts | Native presences | Lost | Loss |
|---|---|---|---|---|---|---|
| 1 | Marianas (MRN) | 6 | 4 | 231 | 44 | 19.0% |
| 2 | Caroline Is. (CRL) | 6 | 12 | 570 | 107 | 18.8% |
| 3 | France (FRA) | 1 | 18 | 32,069 | 4,622 | 14.4% |
| 4 | New South Wales (NSW) | 5 | 20 | 21,056 | 2,616 | 12.4% |
| 5 | California (CAL) | 7 | 15 | 23,839 | 2,907 | 12.2% |
| 6 | Victoria (VIC) | 5 | 10 | 8,767 | 964 | 11.0% |
| 7 | Marquesas (MRQ) | 6 | 2 | 248 | 27 | 10.9% |
| 8 | Cape Verde (CVI) | 2 | 4 | 584 | 63 | 10.8% |
| 9 | Spain (SPA) | 1 | 18 | 33,466 | 3,559 | 10.6% |
| 10 | Seychelles (SEY) | 2 | 3 | 326 | 32 | 9.8% |
| 11 | Germany (GER) | 1 | 11 | 15,351 | 1,468 | 9.6% |
| 12 | Azores (AZO) | 2 | 3 | 732 | 69 | 9.4% |
| 13 | Sri Lanka (SRL) | 4 | 5 | 1,937 | 174 | 9.0% |
| 14 | Maldives (MDV) | 4 | 6 | 129 | 11 | 8.5% |
| 15 | Finland (FIN) | 1 | 13 | 7,472 | 605 | 8.1% |
| 16 | Italy (ITA) | 1 | 15 | 24,986 | 1,992 | 8.0% |
| 17 | Ogasawara-shoto (OGA) | 3 | 3 | 103 | 8 | 7.8% |
| 18 | Sweden (SWE) | 1 | 15 | 14,092 | 1,084 | 7.7% |
| 19 | New York (NWY) | 7 | 6 | 5,320 | 392 | 7.4% |
| 20 | Portugal (POR) | 1 | 5 | 6,040 | 445 | 7.4% |

**Islands do lose more natives in some places, but not as a rule.** Nine of the 20 are island regions (Marianas, Caroline Is., Marquesas, Cape Verde, Seychelles, Azores, Sri Lanka, Maldives, Ogasawara-shoto, mostly with few natives, 103 to 732, but Sri Lanka has 1,937),
and the other 11 are large, heavily sampled regions of Australia, western and northern Europe and North America, which lose the most presences in number: France 4,622, Spain 3,559, California 2,907, New South Wales 2,616, Italy 1,992 and Germany 1,468. Summed over the regions with at least 100 native presences, those touching at most 5 land parts
(99 regions, a proxy for small regions that includes mainland states and not only islands) lose 2.7% (6,044 of 224,081), and the larger ones (225) 2.9% (77,462 of 2,662,995). The loss follows sampling density at least as much as isolation (presumably, as `under_reported` compares a species with the records of its group in the cell, and native species of a heavily sampled cell with few records look under-reported).
Of the cells Task 12 found layer 5 emptying, Taiwan (4.9%) and Madeira (5.6%) lose about twice the median, Iceland (2.4%) and Svalbard (1.0%) lose no more than it, and Bermuda has the highest loss of all, 26 of 79 native presences (32.9%), but is under the 100 that a region needs to be listed.

**Known cases, on half A fully cleaned (all five layers, `native_only`, the design's thresholds; they do not depend on the calibration): 1 of the 4 run passes, 3 fail, and 1 waits.**

| Case | Result | Numbers (presences, species and part) |
|---|---|---|
| Tiger, *Panthera tigris*, no presence in Europe | **waiting for step 9** | family Felidae is in half B |
| Potato, *Solanum tuberosum*, none outside its native regions | **pass** | 779 presences before, 92 after; outside its native regions 588 before, 0 after; native 81 before, 71 after. Removed by `wcvp_outside_native` 328, `griis` 204, `griis_subnational` 88, `under_reported` 15, `failed_flags` 14, `detached_satellite` 13, `zero_area_part` 22, `human_text` 2 |
| Starling, *Sturnus vulgaris*, no presence in North America | **fail** | 473 presences in North America before, 75 after (1,751 and 908 in all). Removed there by `griis` 175, `griis_subnational` 149, `zero_area_part` 64, `detached_satellite` 8, `under_reported` 2 |
| Hawaiian endemic, *Drepanis coccinea* (iiwi), 95% of the presences in the cells of HAW kept | **fail** | 9 presences in the area before, 8 after (10 and 9 in all): 89%, which is as near 95% as 9 presences allow. The one lost, by `under_reported`, has 1 record, in a cell whose centre is at sea (18.5 N, 153.9 W) and which touches HAW; the cells with the 26,856 and 1,391 records are kept |
| Lemuridae (24 species), only on Madagascar | **fail** | 172 presences in the cells of MDG before, 138 after; outside them 56 before, **9 after**: 4 on parts with no TDWG region (sea parts), 2 in the Comoros (726 records), and 3 in cells that touch several regions: Denmark and Sweden (19 records), Poland and Belarus (6) and North Carolina, South Carolina and Virginia (6) |

- **Starling:** the 75 presences left are 30 sea parts and 45 land parts, and none of the land parts is in a cell of a sub-national GRIIS list. GRIIS lists the starling for the US only through two sub-national lists, Alaska and `United` (the lower 48, which `natives.py` places on the cells whose **centre** is in a state, 115 cells), and the country rule has no US entry.
  So a coastal cell whose centre is at sea is not covered (6 of the 15 cells that touch California are not, and California alone has 6 presences left with 639,344 records), and the starling's US rows there stay. The rows of a left part also include other countries' (a border cell with US rows is kept by them). This is layer 3's placement of sub-national lists and not a calibration question, and **nothing was
  changed to make the case pass.** GRIIS is name-only until Task 11a.
- **Lemuridae:** the 9 left are in two groups. **4 of the 9 are parts with no TDWG region**, sea parts, which no explanation here covers (a one-off look: two are the sea part of the cell of the Comoros presences, 58 records, and two are in a cell of the Mozambique Channel at 14.5 S, 43.4 E, 9 records). **The Comoros and captive explanation covers 5 of the 9**: the 2 Comoros land presences,
  where lemurs live as an introduced population, and the 3 in Europe and North America, which look like captive animals (the Duke Lemur Center is in North Carolina); neither was checked in the data. Layer 5 spares a species with under 20 records (all its cells are core) and a patch with at least 2 datasets and 3 years; which of the two applies to each was not looked at.
  Removing them would need knowledge of the range, which layer 3 has for plants and GRIIS for aliens, and not a heuristic. Reported as a failure.
- **Hawaii** shows how little a presence-count criterion can say for a species of 9 presences: one single-record presence in the open sea is 11%. The range itself (the cells with 26,856, 1,391, 1,286 and 236 records) is kept. The criterion is as it was set, and the case fails by it.

### Changes after the review of Task 13 (2026-11-01)

- **The reading of the 5% limit** changed from layer 5 alone to the whole run of layers 1, 2, 4 and 5 counting drops only (see the rulings above), by the controller's ruling, provisional until Markus confirms it before the resolution 7 run. The chosen setting at resolution 6 moved from 80 to 75, and the island table is the one of setting 75.
  The three readings are reported for the chosen and the design setting. The first version's numbers (setting 80, 25.5% outside removed at 4.98% native lost by layer 5 alone, 6.15% over the whole run) are in the readings table.
- **The load of half B** was described as plants only, and it is the plants and every other species of their taxonomic groups; the code docstring, the call-site list, the ruling on half B, the README and the test name now say so, with the number of species that are not Plantae (2,586 of 197,567), counted from the taxonomy alone. The group mates stay in, so that the score is like for like with A's.
- **The once guard:** `--score-half-b` refuses any resolution but 7 and refuses if any `res*/held_out.json` exists, and not only the one of its own resolution. A second round made it key on a marker written **before** the load as well: `held_out.json` is written after the unsealed load, so a run that fails in between (out of memory, `calibrate` raising)
  left nothing and a retry would have been a second load. `res*/held_out.started`, with the time and `git rev-parse HEAD`, is written first (exclusive creation), and either file refuses the run with a message that a retry needs Markus's explicit decision. No override flag.
- **`Baseline` records its settings** and `calibrate(..., start)` checks them against the grid's; the rules layers 1, 2 and 4 read are derived from `RULES`; the test world has the rows of two countries for one part, so that the seam is checked on rows by country.

**Not done here:** the cross-check on A's animals against GRIIS (Task 11a); the scoring of half B's calibration species and the held-out rates as the remaining error (resolution 7, with `--score-half-b`); the cases of half B (step 9); Markus's confirmation of the reading of the 5% limit; a grid that brackets the record share and the cell share.

### GRIIS by its own establishment means (2026-11-01)

`natives.griis_names` now reads each export's distribution table (`Distribution.tsv`, `col:taxonID` and `col:establishmentMeans`) and gives each entry a `means`;
`griis_pairs` applies the introduced ones; `References` says what the lists hold. This **supersedes the ruling of "Reference data" above that every listed species is kept as
alien**, which rested on the plain exports having no flag to read. Task 11 found it removing native records (*Gymnorhina tibicen* in AU, 2.8M records). Rulings, with the reason for each:

- **An entry applies only if it is `introduced`; an empty value counts as introduced.** A GRIIS list is a list of aliens, so an entry with no value is not a statement that it is native. *(Replaced on 2026-11-11: an empty value is left out. See "The GRIIS lists after ChecklistBank's reimport".)* `native` and `uncertain` never apply, as `uncertain` is left out in WCVP.
  An entry with no distribution row at all in a list that has a table would count as empty too, and there is none (0 of 52,462 species entries). The means of an entry that several rows or several names of one COL species give is settled **native first, then introduced, then empty, and `uncertain` only where nothing else is said**,
  so a species that is native anywhere in a list is no alien of it, as in WCVP. 11 species of the 55 lists have both native and introduced rows (AU 7, IE, IT, PT-Azores and IC-Tenerife 1 each). A value other than `native`, `introduced` and `uncertain` raises, as in WCVP.
- **A list with no distribution table is not applied, and is reported.** Only extended exports carry one. `References.griis` holds the pairs of the lists that have a table, and `References.griis_name_only` (with `griis_name_only_no_cells`) those of the lists that have none, which `clean` applies only under
  **`Settings.griis_without_status`, `skip` (the default and the rule) or `apply`**: `apply` applies them name-only as every list was before, **and is an interim setting until Markus has requested extended exports for the 250 lists that only have a plain one** (`scripts/request_exports.py --griis-extended` still reports 250, so
  Step 0 was not run: nothing was fetched or submitted). `References.griis_lists` has, per list, `table`, `entries`, and the entries that are `introduced`, `unstated` (an empty value, counted as introduced), `native` and `uncertain`; a list with no table has its counts of the kinds at 0 and says only `entries`.
  The unplaced pairs of the lists with no table are reported only where they would be applied. `scripts/clean_counts.py` and `scripts/calibrate.py` take `--griis-without-status`.
- **What the lists say** (`scripts/check_natives.py`, species-rank entries; the brief's 61,226 distribution rows and 3,534 native ones are of all ranks). 303 of the 304 lists have entries. The 55 with a table hold 52,462 entries: 48,972 introduced, 28 with no value (all GRIIS-TR), **3,457 native** and 5 uncertain. The other 248 hold 68,085 entries that say nothing.
  Native entries by list: GRIIS-EC 731 of 731 and GRIIS-SR 176 of 176 (all of them, so **these two lists apply nothing**), US-Hawaii 440 of 5,542, US-United 429 of 8,349, IC-Tenerife 423 of 1,119, IC-Gran-Canaria 287 of 777, SG 211 of 695, AU 193 (and 1 uncertain) of 4,463, GB 155 of 2,783, IT 125 of 1,701, JP 84 of 2,331, GR 43 (4 uncertain) of 670.
  Result: **46,519 pairs of 20,556 species in 53 lists apply, and 62,375 pairs of 14,170 species in 248 lists say nothing**, against the 112,116 pairs of 27,389 species of before (the 3,222 pairs fewer are the native and uncertain entries).
- **The `native` flag is not reliable in the US list, and the rule is applied all the same (reported, nothing changed).** In `GRIIS-US-United` *Passer domesticus* and *Columba livia* are `native` with no remark, while *Sturnus vulgaris* and *Streptopelia decaocto* are `introduced`, and 40 of its 431 native rows carry the remark "Native to: Eurasia" or "Native to: Europe" (23 of Hawaii's 451).
  Neither bird is native to the US, so the House Sparrow and the Rock Dove are no longer removed from the lower 48 by GRIIS (12.4M and 6.8M records fewer than before even with the plain lists applied, 13.1M and 9.2M with them skipped: the first two rows of the table below). The brief's rule is followed, since the flag is the list's own statement; whether to override it for the US is Markus's call. The `degreeOfEstablishment` of every native row is `native` as well, and that of introduced rows is empty.
- **Reports added, no rule changes:** `realms.calibrate.griis_cross_check` (layer 5 against GRIIS on half A's animals, deferred from Task 13; `scripts/calibrate.py` writes `griis_cross_check.json`), `scripts/griis_review.py` with `realms.clean_report.griis_totals`, `griis_change` and `griis_feral_candidates` (what GRIIS took, what changed, and the likely feral forms of native species), and `griis_lists` in the reports of `check_natives.py` and `clean_counts.py`. Their figures are in the section below.

### Sub-national GRIIS lists on the cells they overlap (2026-11-01)

`natives.admin1_cells` places a sub-national list on **every cell whose polygon overlaps the subdivision's admin-1 polygons with some area** (touching along an edge is not overlapping), and `via` says whether the cell's centre is in the polygon (`centre`) or only the overlap puts it there (`overlap`).
It replaces placement by cell centre, with overlap only for the islands no centre fell on, which left coastal and border cells out: the starling's US presences by the sea and at the Canadian border were not covered (Task 13: 75 of 473 North American presences left). Rulings:

- **A sub-national list applies to the rows in those cells, and where GBIF uses the list's country code in them, only to the rows with that `countrycode`** (`griis_subnational`). A cell is bigger than the places, so a border cell holds a neighbour's rows too, and a CA row in a cell that overlaps a US state is Canada's: it stays, and Canada's own list may take it (`griis`).
  **"GBIF uses the code" is read per list, from the counts, in that list's cells:** the code occurs as a `countrycode` among the rows (all rows of the download of half A, before layer 1 takes any out) of the cells of the list. The brief reads "the code occurs as a `countrycode` in the counts", which I first took over the whole counts, and **changed to the cells of the list on the real data**: the sovereign's code is used elsewhere for an overseas place and not in it.
  `FR` is in the counts and in none of the cells of Reunion, Guadeloupe, Martinique and Saint-Pierre, whose rows say `RE`, `GP`, `MQ` and `PM`, so over the whole counts those four lists would apply to nothing (Reunion alone, 84,691 records of `RE` that the list takes). Read per list they apply by the cells alone, and the lists of the US, PT, TO and FM, where the code is used, apply by code. **On the real data only the FR overseas lists take the cells-alone way** (Reunion and Guadeloupe have a table, Martinique and Saint-Pierre are plain): **`IC`, the Canary Islands, has no Natural Earth polygon, so it has no cells at all** (IC-Tenerife and IC-Gran-Canaria are among the 6 subdivisions without cells, below). It is the case of the tests, where no row says `IC` and GBIF says `ES`, and not one of the data.
  **The cost of the cells-alone case is that a neighbour's rows in the same cells go too:** in the default run `MU` 1,056 records (Reunion's cell), `AG` 1,930, `KN` 1,906, `DM` 1,004 and `MQ` 6,893 (Guadeloupe's cells), which is what happened to them before this change (`AG` 4,267, `KN` 3,138, `MQ` 7,516). The code is not told from the polygon's `iso_3166_2` (`FR-RE`), which would spare them; that is a refinement not made.
  `NO-Svalbard` is the other way: `NO` occurs in its cells (960 records, 3 rows) and GBIF says `SJ` (347,303), so by code it would take almost nothing; its list is plain and not applied under `skip` in any case.
- **A country's own list still applies to the rest of the country**, meaning its rows outside the cells covered by that country's sub-national lists (`covered`, now larger). Tested with a cell by centre and a cell by overlap only: both are left to the sub-national list, and the rows in the rest of the country are the country list's.
  **`covered` is made of the lists that are applied only** (review of Task 11a, 2026-11-01). The first version took every placed list, so under `skip` a sub-national list with no distribution table, which is not applied, still exempted its cells from the country's own list, and the country list's rows there were taken by no rule. A list counts as applied when it has pairs in `References.griis` (and under `apply` in `griis_name_only` too). Tested with a country list with a table and a plain sub-national list. **No real-data figure changes:** `clean_counts.py` run again on half A in the three modes, with `skip` and with `apply`, gives the same rows, presences and records for every rule and the same presences left as the first run (no `FR` row is in the cells of the plain FR lists, so the country list had nothing there to take). The GRIIS pairs of `load_references` are made by `natives.griis_split`, which `griis_pairs` and `griis_name_only_pairs` call too, so the tested functions are the ones the references use.
  It is by definition that a border cell that overlaps a subdivision is the sub-national list's for that country's own rows, so a species on the country list only keeps its rows there.
- **Result at resolution 6:** 263 cells for the same 24 subdivisions (138 by centre, 125 by overlap only), against 186 (138 and 48); at resolution 7 it was not run. **Task 10's unplaced subdivisions are unchanged: 61 have no polygon in Natural Earth**, and their pairs fell from 7,468 to 6,820 (1,630 in 6 subdivisions of lists with a table: IC-Gran-Canaria, IC-Tenerife, KY-Grand, KY-Little, SH-StHelena, SH-Tristan, and 5,190 in 55 of lists with none), the 648 being native and uncertain entries left out.

### Task 11a on half A at resolution 6 (2026-11-01)

`scripts/check_natives.py`, `scripts/clean_counts.py` (three WCVP modes, once with `--griis-without-status skip` and once with `apply`, into `results/clean/res6` and `results/clean-apply/res6`), `scripts/griis_review.py` and `scripts/calibrate.py` (both settings, `results/calibrate/res6` and `results/calibrate-apply/res6`) on `data/gbif/0019721-261023074918325`.
The runs before this task are the baseline (`results/clean-before/res6`, which are those of "Layers 1 to 4 on half A", layer 5 included). Half A only. **The counts are of COL 315557 (COL26.6 XR).** The GRIIS rules take **the same rows in all three WCVP modes** (as before: they run first).

**GRIIS before and after** (the rules `griis`, the country lists, and `griis_subnational`, together; presences are the (species, part) pairs a rule left with nothing, counted once across both):

| Run | `griis` rows / records | `griis_subnational` rows / records | Both: rows | Presences | Records | Species |
|---|---|---|---|---|---|---|
| before | 325,520 / 84,405,810 | 78,247 / 47,923,761 | 403,767 | 339,225 | 132,329,571 | 10,011 |
| `skip` (the default) | 160,442 / 50,648,154 | 85,707 / 30,668,327 | 246,149 | 215,670 | 81,316,481 | 8,007 |
| `apply` | 317,611 / 72,807,032 | 88,591 / 30,775,593 | 406,202 | 341,369 | 103,582,625 | 9,771 |

- Presences by rule: `griis` 262,907 before, 132,928 `skip`, 255,649 `apply`; `griis_subnational` 76,318, 82,742 and 85,720. **`skip` takes 39% fewer records than before and `apply` 22% fewer**: `apply` is before less the native entries and plus the cells; `skip` is that less the 248 plain lists (56% of the entries that GRIIS has).
  The `griis` records of the biggest countries fall in `skip` from AU 15.2M to 7.5M, GB 13.1M to 9.5M and SE 10.4M to 0 (GRIIS-SE has no table), while CA stays 9.0M.
- What is left of the clean records and presences, with every layer (layer 3 takes WCVP last, so it takes more where GRIIS takes less): `native_only` 6,963,228 presences of 588,877 species and 1,814,972,688 records in `skip` (6,909,903, 588,738 and 1,796,325,353 in `apply`, against 6,905,878, 588,707 and 1,767,233,258 before),
  `introduced` 7,141,188, 589,497 and 1,817,072,357 (7,075,025, 589,368 and 1,798,119,294; before 7,071,587, 589,337 and 1,769,041,809) and `off` 7,254,566, 589,544 and 1,825,239,385 (7,147,634, 589,416 and 1,803,049,778; before 7,146,812, 589,387 and 1,774,293,399).
  `wcvp_outside_native` takes 425,026 rows and 10.4M records in `skip` (349,740 and 6.8M in `apply`, 355,511 and 7.2M before), and `wrims` 3,765 (2,972 and 3,125) rows, for the same reason. Gaps the lists leave: 860 unplaced GRIIS pairs of species with records in `skip` (3,677 in `apply`, 4,027 before) and 1,140 WRiMS pairs.
  **These figures are from before layer 3 judged the rows that layer 4 moved; the figures after it are in that section below** (`native_only` in `skip`: 6,952,838 presences of 588,865 species and 1,814,830,625 records).
- **The sub-national rule's records** fell from 47.9M to 30.7M (`skip`), the US's from 47.62M to 30.47M. That is the native entries of the US lists, not the placement: the CA records that the cells of the lower 48 took, **54,020 before, are 0 now**, and the plain Alaska list is added back in `apply` (US 30.55M).

**The 20 species whose GRIIS removals changed most** (records, the two rules together; the change in `skip`, with `apply` beside it where the species is in both lists):

| Species | Before | `skip` | Change | `apply` |
|---|---|---|---|---|
| *Passer domesticus* | 16,739,383 | 3,640,762 | -13,098,621 | 4,375,509 |
| *Columba livia* | 12,588,666 | 3,369,763 | -9,218,903 | 5,775,906 |
| *Gymnorhina tibicen* | 2,933,086 | 143,212 | -2,789,874 | 143,212 |
| *Anas platyrhynchos* | 2,509,043 | 233,481 | -2,275,562 | |
| *Branta canadensis* | 2,521,461 | 629,185 | -1,892,276 | 1,825,981 |
| *Cygnus olor* | 3,644,601 | 1,901,042 | -1,743,559 | 3,021,674 |
| *Phasianus colchicus* | 2,685,988 | 1,266,015 | -1,419,973 | 1,864,178 |
| *Cygnus cygnus* | 1,241,827 | 0 | -1,241,827 | |
| *Ocyphaps lophotes* | 1,228,523 | 11 | -1,228,512 | 11 |
| *Sturnus vulgaris* | 20,145,777 | 21,317,168 | +1,171,391 | 21,448,162 |
| *Branta leucopsis* | 608,470 | 0 | -608,470 | |
| *Wallabia bicolor* | 570,889 | 32 | -570,857 | 32 |
| *Streptopelia decaocto* | 5,182,525 | 5,704,102 | +521,577 | 5,796,888 |
| *Geopelia placida* | 520,643 | 0 | -520,643 | 0 |
| *Geopelia humeralis* | 463,837 | 0 | -463,837 | 0 |
| *Phascolarctos cinereus* | 451,337 | 0 | -451,337 | 0 |
| *Mareca strepera* | 392,912 | 0 | -392,912 | 0 |
| *Mareca penelope* | 354,612 | 0 | -354,612 | 0 |
| *Apis mellifera* | 414,627 | 70,777 | -343,850 | 115,139 |
| *Alectura lathami* | 277,990 | 0 | -277,990 | 0 |

The Australian natives (the magpie, the crested pigeon, the swamp wallaby, the koala, the brush-turkey) are no longer removed in AU: GRIIS-AU calls them native (193 entries), and what is left of the magpie (143,212) is GB, BE and NZ, where it is introduced. The ducks and swans fall under `skip` because GRIIS-SE and GRIIS-DK have no table, and the two species that rise,
the starling and the collared dove, rise because the overlapping cells add coastal and border cells to the US list, where both are `introduced`. With `apply` the top 20 is the same but for *Anas platyrhynchos* and *Cygnus cygnus* (SE and DK name-only), who drop out of it. Files: `results/clean/res6/griis_review/griis_change.csv` and `results/clean-apply/res6/griis_review/griis_change.csv`.

**The starling, *Sturnus vulgaris*, in North America (`scripts/calibrate.py`): it does not pass in either setting, and why.**

| Setting | Presences in North America before | After | Removed there by | Left in |
|---|---|---|---|---|
| Task 13 | 473 | 75 | `griis` 175, `griis_subnational` 149, `zero_area_part` 64, `detached_satellite` 8, `under_reported` 2 | 30 at no region (sea parts), British Columbia 7, Alaska 6, California 6, Florida 5, Texas 5, Greenland 4, Maine 4, ... |
| `skip` (the default) | 473 | **31** | `griis` 175, `griis_subnational` 188, `zero_area_part` 64, `detached_satellite` 4, `under_reported` 11 | no region 14, Alaska 11, British Columbia 5, Greenland 4, Yukon 2, Bahamas 1, Florida 1, Newfoundland 1 |
| `apply` | 473 | **6** | `griis` 175, `griis_subnational` 225, `zero_area_part` 64, `detached_satellite` 2, `under_reported` 1 | Greenland 4, no region 2 |

The placement fix took the starling from 75 presences to 6 with all lists and to 31 under `skip`: California, Texas, Maine and Mexico are gone, as the coastal and border cells are now the US-United list's, and the CA rows in them go by Canada's country list. **What is left is not a placement question:**
under `skip` **the starling is named by the plain lists of Alaska, Saint-Pierre, the Bahamas and Bermuda, which have no distribution table and so are not applied**: the 31 are in the cells of Alaska (their `US` rows stay, and the cells reach British Columbia and Yukon, whose `CA` rows go), in the Florida and Bahamas cell (the `US` rows go, the `BS` rows stay: 1,853 records on the land part), in the Saint-Pierre cell (the `PM` rows stay, the `CA` ones go), and four Greenland presences with their sea parts.
Under `apply` those lists apply and **the 6 left are Greenland's** (4 land presences and 2 sea parts, `GL` rows of 1 to 3 records each): GRIIS-GL does not name the starling, so no GRIIS rule can take them, and layer 5 spares them. The case will pass when the extended exports of the plain lists arrive, if they say `introduced` for the starling, except for Greenland, which GRIIS has no entry for.
The calibration itself is untouched (the same chosen setting, 75, and the same readings: layer 3 is off there). Of the other cases, the potato passes in both settings (93 and 92 presences left), and the Hawaiian bird and Lemuridae fail exactly as before.

**The cross-check of layer 5 against GRIIS on half A's animals** (`realms.calibrate.griis_cross_check`; layers 1, 2, 4 and 5 on everything half A has with layer 3 off, so layer 5 sees the group totals of the real run; an animal presence is *listed* where GRIIS lists its species as introduced in that country or sub-national place, by layer 3's own rules, from the lists with distributions only, so it is the same in both settings):

| Setting | Presences listed as introduced taken by layer 5 | Other animal presences taken | (species, cell) pairs, listed / other | Records, listed / other |
|---|---|---|---|---|
| 75 (chosen): 2%, 10%, 1%, 300 | 5,663 of 43,388 (**13.05%**) | 319,321 of 4,186,023 (**7.63%**) | 14.13% / 7.89% | 0.021% / 0.038% |
| 40 (the design's): 1%, 5%, 1%, 1,000 | 5,350 of 43,388 (**12.33%**) | 280,753 of 4,186,023 (**6.71%**) | 13.35% / 6.93% | 0.018% / 0.033% |

There are 4,229,411 animal presences after layers 1, 2 and 4; 43,388 of them (1.0%, 2,977 species, 64.5M records) are listed. **Layer 5 takes 1.7 to 1.8 times as many of the listed presences as of the others, and 13% of them at most, so the heuristics find some of what GRIIS says is introduced and far from all of it.**
By records the listed presences lose less than the others (0.02% against 0.04%): they are the large, long-established populations of the best sampled countries (the US, GB, AU), which layer 5 spares by design as core cells, and a presence it takes has a handful of records. So it is no substitute for GRIIS where GRIIS has a table, and it is not tuned on it. The listed set is also biased towards the 55 countries and places with an extended list, which is not the world.
The unit that is listed is the (species, cell) pair, so another country's rows in a listed species' border cell count as listed with it, though no GRIIS rule took them. Report only: nothing was tuned on it, and the base rates of the two sets differ (a listed species is one GRIIS knows because it is widespread). Files: `results/calibrate/res6/griis_cross_check.json`.

**Report only, for Markus and M1: the 20 largest GRIIS removals, by records, of species that also have clean, unlisted records in at least three neighbouring countries.** Two countries are **neighbours where GBIF has rows of both in the same cell or in adjacent cells of the grid** (from all of half A's rows, `realms.clean_report.neighbour_countries`; it takes in sea neighbours, as the cells are 265 km). "Clean" is `n - n_human_text` above 0 and "unlisted" is not taken by a GRIIS rule.
These are the likely feral or released forms of native species. GRIIS calls them introduced, so **no rule changes**. `scripts/griis_review.py`, files `griis_neighbours.csv` and, with 100 records at least in each neighbour, `griis_neighbours_min100.csv` in `results/clean/res6/griis_review` and `results/clean-apply/res6/griis_review`. All 20 rows for each of the two settings follow.

**With `skip`, the default** (`results/clean/res6/griis_review/griis_neighbours.csv`; records removed by the two GRIIS rules in that country, then how many neighbouring countries hold clean, unlisted records of the species, with the three largest as country:records):

| # | Species | Country | Records removed | Neighbours | The largest neighbours' unlisted records |
|---|---|---|---|---|---|
| 1 | *Sturnus vulgaris* | US | 16,740,718 | 7 | RU:95874 BM:22434 BS:2468 |
| 2 | *Streptopelia decaocto* | US | 4,257,018 | 7 | BS:21677 CU:16389 RU:10210 |
| 3 | *Sturnus vulgaris* | CA | 3,290,292 | 3 | US:22801 PM:627 GL:13 |
| 4 | *Passer domesticus* | CA | 2,077,204 | 3 | US:13225272 PM:578 GL:4 |
| 5 | *Columba livia* | CA | 1,675,621 | 3 | US:6978356 PM:209 GL:1 |
| 6 | *Cygnus olor* | US | 1,527,775 | 4 | RU:27868 CA:97 MX:17 |
| 7 | *Streptopelia decaocto* | FR | 1,007,611 | 24 | ES:571789 GB:558339 NL:430752 |
| 8 | *Ondatra zibethicus* | NL | 898,530 | 14 | CH:4064 PL:2104 AT:840 |
| 9 | *Sturnus vulgaris* | AU | 840,097 | 3 | NF:2290 ID:2 NC:2 |
| 10 | *Acridotheres tristis* | AU | 808,680 | 5 | NC:1986 ID:932 SB:471 |
| 11 | *Spilopelia chinensis* | AU | 704,001 | 4 | NZ:31290 ID:29827 NC:1822 |
| 12 | *Passer domesticus* | AU | 611,436 | 5 | NF:2014 NC:1089 PG:657 |
| 13 | *Cistus ladanifer* | ES | 559,068 | 6 | PT:10977 FR:5371 MA:103 |
| 14 | *Alopochen aegyptiaca* | NL | 555,805 | 12 | DK:19071 CZ:10821 SE:7877 |
| 15 | *Columba livia* | ES | 542,079 | 13 | GB:589200 PT:257033 FR:238891 |
| 16 | *Phasianus colchicus* | GB | 534,745 | 13 | NL:366845 DK:201015 DE:106883 |
| 17 | *Dactylis glomerata* | GB | 470,172 | 15 | FR:473886 NL:172114 DK:93628 |
| 18 | *Columba livia* | AU | 434,438 | 4 | ID:4791 NF:1100 NC:290 |
| 19 | *Quercus suber* | ES | 428,837 | 8 | PT:258918 FR:11965 IT:2051 |
| 20 | *Branta canadensis* | BE | 380,770 | 11 | GB:695480 NL:368057 DK:113480 |

**With `apply`** (`results/clean-apply/res6/griis_review/griis_neighbours.csv`; the lists of DK, SE and IN are plain, so they apply only here):

| # | Species | Country | Records removed | Neighbours | The largest neighbours' unlisted records |
|---|---|---|---|---|---|
| 1 | *Sturnus vulgaris* | US | 16,762,779 | 5 | RU:95874 CA:664 KY:6 |
| 2 | *Streptopelia decaocto* | US | 4,263,041 | 5 | RU:10210 KY:663 MX:215 |
| 3 | *Cygnus olor* | US | 1,527,775 | 4 | RU:27868 CA:97 MX:17 |
| 4 | *Anas platyrhynchos* | SE | 1,487,602 | 10 | NL:1053414 NO:787556 DE:489491 |
| 5 | *Cygnus cygnus* | SE | 1,184,996 | 11 | DK:334975 NO:261715 FI:91773 |
| 6 | *Columba livia* | IN | 1,103,587 | 10 | TH:97859 NP:31985 CN:14815 |
| 7 | *Streptopelia decaocto* | FR | 1,007,611 | 23 | ES:571789 GB:558339 NL:430752 |
| 8 | *Cygnus olor* | SE | 974,504 | 10 | NL:831258 DK:684569 NO:341163 |
| 9 | *Ondatra zibethicus* | NL | 898,530 | 8 | GB:129 FR:120 DE:23 |
| 10 | *Branta canadensis* | SE | 876,176 | 6 | NL:368057 AX:328 DE:151 |
| 11 | *Sturnus vulgaris* | AU | 840,097 | 3 | NF:2290 ID:2 NC:2 |
| 12 | *Anas platyrhynchos* | DK | 781,870 | 12 | GB:1209662 NL:1053414 NO:787556 |
| 13 | *Spilopelia chinensis* | AU | 704,001 | 3 | NZ:31290 NC:1822 TL:1552 |
| 14 | *Passer domesticus* | AU | 611,436 | 5 | NF:2014 NC:1089 PG:657 |
| 15 | *Cistus ladanifer* | ES | 559,068 | 6 | PT:10959 FR:5371 MA:103 |
| 16 | *Alopochen aegyptiaca* | NL | 555,805 | 8 | CZ:10821 CH:4382 DE:261 |
| 17 | *Branta leucopsis* | SE | 548,226 | 11 | NL:385046 DK:289910 NO:94174 |
| 18 | *Columba livia* | ES | 542,079 | 13 | GB:589200 PT:257033 FR:238891 |
| 19 | *Phasianus colchicus* | GB | 534,745 | 11 | NL:366845 DE:106883 IM:6184 |
| 20 | *Dactylis glomerata* | GB | 470,172 | 14 | FR:473886 NL:172114 DK:93628 |

These include the examples of the brief only under `apply`: **the Mallard in SE (1.49M, row 4) and DK (0.78M, row 12)** and **the Rock Dove in IN (1.10M, row 6)**, with the Whooper and Mute Swans in SE. Under `skip` (the default) the Swedish, Danish and Indian rows are not removed, and the list is the US, Canadian, French, Dutch, Australian, Spanish, British and Belgian cases above.
In `skip` the House Sparrow and the Rock Dove appear in CA (rows 4 and 5) with the US as their largest neighbour, 13.2M and 7.0M unlisted records, because `GRIIS-US-United` calls them native (see above), not because the species is native there. **The definition is generous**: any three neighbours count, however few records, so the list also has the starling in AU, whose neighbours hold NF 2,290, ID 2 and NC 2.
`griis_neighbours_min100.csv` asks 100 records of each neighbour: under `apply` it leaves out the starling in the US and AU, the Mute Swan in the US and the Muskrat in NL, and keeps the Mallard, the swans, the Rock Dove in IN (9 neighbours) and the Pheasant in GB. A neighbour that holds unlisted records may also hold them because it has no list, or no table, and not because the species is native there: a native range is not what this list shows.

### The ingest's completion marker, and the resolution read from it (2026-11-01)

From the final review of the branch. **An ingest is complete when its directory holds `_COMPLETE`**, and the loaders refuse one that does not. The reason: the resolution 7 ingest is larger, and an out-of-memory kill or a
`SIGKILL` in the middle of it runs no cleanup, so it could leave parts that `load_presence` and `load_country_counts`, and with them the one unsealed load of half B, would read without complaint.

- **What it holds:** `rows`, `resolution`, `zip` (the name), `zip_sha256`, `written` and `how` (`ingest`, or `mark-complete`). It is written to a temporary file and renamed.
- **Who writes it:** `scripts/ingest.py`, after the ingest **and** its QA, and only if the three pass marks hold (`realms.occurrences.mark_complete`), after checking that the parts hold the rows the
  ingest counted. **`ingest` removes the marker before anything else**, then the earlier parts, so a run that is killed half way leaves parts and no marker, and one that fails its QA leaves none either.
- **Who reads it:** `_select`, so both loaders, right after the seal and before any row is read (`read_complete`): no marker, a marker without one of its keys, or parts that hold other rows than it says
  (from the Parquet metadata, a global count that reads no row) raise `IncompleteIngest`. `scan_counts` does not check it, as `qa` reads the parts before the marker exists.
- **The res 6 ingest got its marker by `--mark-complete`** (`uv run python -I scripts/ingest.py data/gbif/0019721-261023074918325.zip --resolution 6 --mark-complete`), which ingests nothing: it checks that
  the parts hold the rows of `results/pilot/res6/qa.json` (19,531,246, its `rows` and its `ingest.rows`), that `qa.json` passes the marks, and that every cell of the parts is in the grid of the stated resolution
  (all 7,234 are in `ISEA3H-6`), and then writes the marker with the zip's sha256 (`9fbba9b1ecd81fdd47bd9febfef0fac18b01792ed7ea29d07459a352d154bc03`, the same as `shasum -a 256`), `how` `mark-complete`. 4 s, 0.7 GB.
- **`--resolution` of `scripts/clean_counts.py` and `scripts/calibrate.py` is read from the marker** (`realms.occurrences.resolution_of`), before the references are loaded. The flag stays, as a statement
  that must agree with the marker, and a mismatch stops the script at once; before, `--resolution` defaulted to 6 and a mismatch failed only after the 16 GB load of the references (if at all). The grid file
  is the marker's resolution's, and its own resolution is checked against it too. `--score-half-b` checks the marker's resolution, so stating 7 for an ingest of another resolution does not get past it.

### Layer 3 judges the rows that layer 4 moved (2026-11-01)

From the final review of the branch. **Layer 4 put back presences that layer 3 had emptied.** `clean` runs layer 4 after layer 3, and layer 4 folds a species' rows on the wrong side into the
other part of the cell, or into a neighbour's, including into a part that layer 3 had just emptied. Two hand-made cases: a marine crab that WRiMS lists for AA loses its 50 AA records in the sea part of a
coastal cell, and layer 4 then folds its 5 AA records from the cell's land part back into the sea part; a terrestrial plant loses a land part to `wcvp_outside_native`, and its records in that cell's
sea part are folded back into it. On the res 6 data (`skip`, read from the outputs by the reviewer): of the presences WRiMS emptied, 1,057 of 3,403 were back at the end (64,196 records folded in), of
those `wcvp_outside_native` emptied 6,994 of 394,350 (65,359), and of those `wcvp_introduced` emptied 5,485 of 141,317. So `presence_lost`, and the `presences` of `removal_summary`, said these rules
had removed presences that were there at the end.

**The ruling, and the fix:** after layer 4, the rules of layer 3 run again on every row that layer 4 moved, at the part it went to, and a row they take is logged under the layer 3 rule (layer 3), after
the fold (layer 4) that brought it. **Layer 3's judgement is about the place**, so a fold must not bring back what it took out of a place, and that holds for a fold into the neighbouring cell too. It runs
only when layer 3 is on, with the same `Settings` (`rules_off`, `wcvp_mode`, `griis_without_status`), in the same order of rules. Every rule decides a row by its key alone (species, country code and part,
with what the references and the counts say of them), so a row the first pass kept is kept again, and what the second pass takes is the records the folds brought.

- **A row can now be logged twice by layer 3**, once before the fold and once for the records the fold brought to it, by the same rule. The rows and the records of the log still add up (each row of the
  log is records of its own), and the presence is counted once by `removal_summary` and `griis_totals`.
- **`presence_lost` says what the species had left right after the rule took the row.** A fold of another country's rows, which layer 3 does not take where they went, can still bring the presence back:
  WRiMS lists the crab for AA, and the BB records folded into the sea part are the sea part's, as they would have been had GBIF put them there. That is the rule as it is, and it is tested. With
  `native_only` in `skip`, 51 presences are like that at the end (WRiMS 20, `griis` 25, `griis_subnational` 6; `apply` 65), against 8,091 before. None is of WCVP, which judges every country of a part alike.
- **Tests:** the reviewer's two cases (`test_a_fold_does_not_put_back_what_wrims_took_from_the_part_it_folds_into`, `..._what_wcvp_took_...` in both WCVP modes), a fold into the neighbouring cell, a
  fold of another country that stays, the same settings (`rules_off`, and a plain sub-national list under `skip` and `apply` reached by a neighbour fold), and layer 3 off. In the random world, which now has 12
  shore species per seed from a stream of their own and `Sub` on cell 5 as well: **layer 3 run again on what `clean` leaves takes nothing**, over 6 seeds, the three WCVP modes, `skip` and `apply`, and layers
  1 to 5 and 1 to 4 (the rows are made from the input and the log alone, which also checks that the log accounts for every species, country and part to the record); a row is logged twice by layer 3 only where a
  fold reached it, and by one rule; and the second pass is reached for `griis`, `wrims` and both WCVP rules (`griis_subnational` by the hand-made neighbour). 49 of the tests fail on the code before.

**Res 6, half A, before and after,** `scripts/clean_counts.py --counts data/gbif/0019721-261023074918325 --resolution 6`, once with `skip` (`results/clean/res6`) and once with `apply`
(`results/clean-apply/res6`), 79 and 78 s, peaks 21.4 and 22.7 GB. **The rows of layers 0, 1, 2 and 4 are identical in every mode and both settings**, every layer 3 row of before is still there, and every new one is
at a part a fold of layer 4 reached. Layer 3's rules, rows / presences / records (the presences a rule left with nothing):

| Rule | Mode | `skip` before | `skip` after | `apply` before | `apply` after |
|---|---|---|---|---|---|
| `wrims` | every mode | 3,765 / 3,403 / 318,895 | 5,604 / 3,919 / 392,634 | 2,972 / 2,784 / 254,830 | 4,400 / 3,149 / 318,580 |
| `wcvp_outside_native` | `native_only` | 425,026 / 394,350 / 10,384,403 | 446,813 / 399,264 / 10,473,999 | 349,740 / 332,589 / 6,839,526 | 367,657 / 337,141 / 6,912,036 |
| `wcvp_introduced` | `introduced` | 164,313 / 141,317 / 8,184,933 | 177,287 / 142,565 / 8,254,445 | 106,168 / 95,411 / 4,949,270 | 115,517 / 96,375 / 5,002,467 |
| `griis_subnational` | `native_only`, `introduced` | 85,707 / 82,742 / 30,668,327 | 85,708 / 82,742 / 30,668,329 | 88,591 / 85,720 / 30,775,593 | 88,592 / 85,720 / 30,775,595 |
| `griis` | every mode | 160,442 / 132,928 / 50,648,154 | the same | 317,611 / 255,649 / 72,807,032 | the same |

(`griis_subnational` with WCVP off takes 2 rows and 3 records more, in both settings.) The second pass itself, in `skip`: `wrims` 1,839 rows (1,281 of them on a row it had taken before the fold), 73,739 records,
1,802 presences of 261 species; `wcvp_outside_native` 21,787 rows (16,385), 89,596 records, 20,787 presences of 5,439 species; `wcvp_introduced` 12,974 rows (11,441), 69,512 records, 12,077 presences of 1,823
species; `griis_subnational` one row of 2 records. In `apply`: `wrims` 1,428 rows, 63,750 records; `wcvp_outside_native` 17,917 rows, 72,510 records; `wcvp_introduced` 9,349 rows, 53,197 records.

**What is left** (presences, species, records; `skip`, then `apply`):

| Mode | `skip` before | `skip` after | `apply` before | `apply` after |
|---|---|---|---|---|
| `native_only` | 6,963,228 / 588,877 / 1,814,972,688 | **6,952,838 / 588,865 / 1,814,830,625** | 6,909,903 / 588,738 / 1,796,325,353 | 6,901,117 / 588,728 / 1,796,206,828 |
| `introduced` | 7,141,188 / 589,497 / 1,817,072,357 | 7,133,382 / 589,485 / 1,816,939,983 | 7,075,025 / 589,368 / 1,798,119,294 | 7,069,297 / 589,358 / 1,798,011,141 |
| `off` | 7,254,566 / 589,544 / 1,825,239,385 | 7,253,184 / 589,532 / 1,825,166,230 | 7,147,634 / 589,416 / 1,803,049,778 | 7,146,545 / 589,406 / 1,802,986,488 |

With `native_only` in `skip`, 10,612 presences are gone and 222 are new: every new one is a presence that layer 5 took before and keeps now, as its input changed (a species' own cells, and its group's totals).
Records are conserved: 1,814,830,625 left + 1,250,685 (layer 2) + 92,183,116 (layer 3) + 38,302 (layer 4, dropped) + 843,619 (layer 5) = 1,909,146,347, the clean records of the input.

**Layer 5 takes less**, as some of what it took were the odd presences the folds made and layer 3 now takes first. With `native_only` in `skip`, rows / records / species, before and after: `detached_satellite`
241,808 / 471,711 / 45,515 and 233,516 / 457,426 / 43,986; `under_reported` 200,614 / 326,151 / 55,808 and 196,485 / 319,154 / 55,284; `thin_long_range` 55,158 / 67,031 / 22,436 and 55,164 / 67,039 / 22,438;
all of layer 5 497,580 / 864,893 / 101,444 and 485,165 / 843,619 / 100,495, and 501 suspect cells and 500. With `introduced` 562,739 and 556,590 rows, 508 and 507 suspect cells; with `off` 571,167 and 570,852, 514 and 514.

**The calibration does not change, as it runs with layer 3 off**, and the known cases do not either. `scripts/calibrate.py --resolution 6`, run again in both settings (278 and 277 s, peaks 21.9 and 20.0 GB):
`settings.csv`, `regions.csv`, `known_cases.json` and `griis_cross_check.json` are byte-identical to the runs before, in `skip` and in `apply`, and so are the chosen setting (75), the readings, the regions and
the top settings in `report.json`. The known cases come from a full clean with layer 3 on, and none of their figures moved (the starling and the Hawaiian bird are birds, which layer 4 never folds; for the potato and the
lemurs the second pass changed no presence and no count of `removed_by_rule`, which counts a presence once). `scripts/griis_review.py` (both settings, with and without `--min-neighbour-records 100`): `griis_change.csv` and both
neighbour tables are identical, and `griis_totals.json` has the one `griis_subnational` row of 2 records more.

### The scoring of half B, a run of its own (2026-11-01)

From the final review of the branch, which found no breach of the seal. The plan allows exactly one `unsealed=True` load of half B before the freeze, the scoring of the calibration at resolution 7, and it has not run.

- **A run of its own** (I2). `scripts/calibrate.py` chose the setting and scored half B in one process: half A's data (about 15 GB at resolution 6, more at 7) was still in memory when half B was loaded, and half B
  was spent on the dry run's choice before Markus had seen the dry run. Now the dry run (no flag) never loads half B and ends with `report.json`, and **`--score-half-b --setting N`**, with a required `N`, loads the
  references and half B's calibration species and nothing else, and refuses unless the dry run at resolution 7 is complete (`res7/report.json` there, with its `chosen`, from the same ingest). The dry run
  removes its `report.json` when it starts, so a dry run that fails leaves none. `held_out.json` records `N`, the dry run's choice and whether they are the same. The guards that were there stay: resolution 7
  only (now the ingest's, from its marker), the markers, exclusive creation before the load, no override flag. The binding sequence is the dry run at 7, Markus's choice, then the scoring.
- **A marker that does not depend on `--results`** (I3). The once guard looked for `held_out.started` and `held_out.json` under `--results`, which is git-ignored and can be pointed elsewhere or cleared. **A second,
  fixed marker, `data/split/held_out-315557.started`**, next to the frozen split (its path is fixed in the script, relative to the repository), is now created with the other by exclusive creation before the load,
  and either one refuses a run. `data/` is git-ignored as well, so **the marker's content is copied into this file, which is tracked, as soon as it is written**: that is the record that cannot be cleared by accident.
- **The seal's tests read the syntax tree** (M2). The test of the one call site was a regular expression on `unsealed\s*=\s*True`, which a line break, `unsealed=flag`, `**kwargs` or a positional call of
  `_select` would pass. It now scans the syntax tree of every module under `src/` and `scripts/` (with `occurrences.py`): a call that passes `unsealed=` with anything but the constant `False`, `**` keywords to
  `load_presence`, `load_country_counts` or `_select`, a call of `_select` outside `occurrences.py` (or one inside that does not hand on its caller's `unsealed`), and the string `"unsealed"`. Exactly one call
  is found, in `score_held_out`. A second scan fails if any module but `occurrences.py` reads the parts with a Parquet reader of its own. Both scans are tested on code that unseals or reads the parts in each way
  they look for.
- **Tests,** on the hand-made world twice over (the species of half A and twins of them in half B, in one made-up download at resolution 7, with the loaders of the grid, references, split and regions stubbed):
  the dry run loads half A only; the scoring loads half B's calibration species once, unsealed, and nothing of half A, runs no known case, grid, region or cross-check, leaves the dry run's files as they were,
  and writes both markers and `held_out.json` with `N` and whether it is the dry run's choice; it refuses without a finished dry run at 7 on the same ingest (missing, no `chosen`, another resolution or half,
  another zip), a setting outside the 81, a missing `--setting`, `--setting` without the flag, any resolution but 7, a marker of either kind (the fixed one alone too), and a retry after a load that failed.

### Smaller fixes after the final review (2026-11-01)

- **`native_status` counts kingdom Plantae only** (M5), as layer 3's WCVP rule does (`clean._wcvp`): a species is a plant of the calibration where COL puts it in Plantae and WCVP gives it a native region.
  Before, a species WCVP lists and COL puts in no kingdom (3 such WCVP natives in half A, 2 in half B) got a status that layer 3 would never act on. `evaluate_case` takes `taxa` for `native_only`.
  **At resolution 6 on half A, 2 native presences fewer** (after layers 1, 2 and 4, 1,808,900; raw, 1,820,923; with a status 2,334,101), and nothing else moves: the same counts of outside presences removed
  and native presences lost for all 81 settings, the same 54 within the limit, the same chosen setting (75), the same top 10, and the shares differ by less than 1e-7. `scripts/calibrate.py --resolution 6`
  run again (279 s, 16.7 GB). Tested: a WCVP plant with no kingdom or another one gets no status, and the outside presences are exactly what layer 3 takes. `calibration_species`, the set that is loaded, is
  unchanged (the two of half B are in it as members of a plant group).
- **`run_infomap` finds the cells with links by `np.bincount`** (T3) instead of a Python set of the link's cell indices (89M links at resolution 7: about 1.3 s and a 3 GB list a worker, against 0.05 s, as Task 3's review measured). Tested against the reference copy, with
  an unlinked cell and with no link at all. **The Costello gate prints nothing:** `uv run python scripts/costello_selection.py --workers 8 | tee results/costello/selection.log`, then `diff` of `selection.csv`
  and of `selection.log` without its time and memory stamps against `results/costello/before/`.
- **Stale text** (M3): the pilot's sentence on the call-site list (it had none then, and has one now); the first ruling on sub-national lists (whatever the row's country code, never on `countrycode`), marked as
  superseded by Task 11a; the conserved total with layer 1 off, which is `n_raw` (here and in `Cleaned`); the comment on `held` in layer 4, which excludes the same-cell folds; `"all"` among the halves that
  `_select` returns only unsealed; and the `scan_counts` docstring, which said to look at it only in global numbers while `_select` reads it per species.
- **`docs/design.md` describes the one load of half B exactly** (T13): not "B's plants" but B's calibration species, the WCVP plants and the other species of their taxonomic groups, 197,567 species of
  which 2,586 are not Plantae, read only through `--score-half-b` at resolution 7. The spec is edited to say what the code does, which the review of Task 13 ruled (the group mates stay in, so that the score
  is like for like with A's); no method changes.

### Open points for Markus, before M2 (2026-11-01)

From the final review, recorded and not changed in code.

- **I4: how `group_min_records` scales from half A to the headline.** `under_reported` needs a group to hold at least `group_min_records` records in a cell, and the calibration chooses it (300, 1,000 or 3,000) on
  half A's group totals, which are about half of what the headline run on all species has in the same cell. A threshold chosen on A means about twice as much on all species: kept as a number, it lets
  under-sampled cells of the headline count as sampled where A's would not (more cells judged), and doubled, it keeps the meaning it had on A. At resolution 6 the calibration hardly depends on it (300,
  1,000 and 3,000 take 26.7%, 26.7% and 26.6% of the outside presences), but the headline and half B's confirmation should say which reading they take. To decide: a number of records, or a share of
  the group's records (which scales by itself).
- **M1: freeze the per-list "GBIF uses the code" decision.** `griis_subnational` applies a sub-national list by country code where the code occurs among the rows of the list's cells, and by the cells
  alone where it does not, and that is read from whatever counts `clean` is given (half A now; all species in the headline; half B's calibration species in the scoring). A list can change sides between
  runs (a code that occurs in a cell only through a species of the other half), and so can what it takes. Before M2 the decision should be made once, from all of the download's rows (global numbers, which
  QA already reads) or from A, written down per list, and passed to `clean` rather than read from its input.

### GBIF's near-miss cell codes at resolution 7 (2026-11-01)

The resolution 7 download `0020744-261023074918325` (DOI 10.15468/dl.baekdm, 28,052,905 rows) failed the QA's second pass mark on **7 cell codes in 35 rows (40 records, all of them clean)** that are not in
`ISEA3H-7.gpkg`. The res 6 download had none. This changes the method of the ingest, and the reason is that the QA's rule, any code the grid doesn't have fails the run, is right for a wrong grid and wrong
for a few cells whose centre GBIF and the grid file give a few hundred metres apart.

- **What the codes are.** A `geogrid_id` without its sign, zero-padded to 19 digits, is a 2-digit code, then |lat| in millionths of a degree in 8 digits and |lon| in 9. **Decoded against `lat_c` and `lon_c`, every
  id of the three grids agrees to 1e-6 degrees** (2,432, 7,292 and 21,872 cells at resolutions 5, 6 and 7; `tests/test_occurrences.py` runs it on the real files and skips if they are not there). The 2-digit code
  is the resolution plus 22 times the hemisphere (05, 27, 49, 71 at 5; 06, 28, 50, 72 at 6; 07, 29, 51, 73 at 7: north and east, south and east, north and west, south and west), and is the only thing
  that says the hemisphere. A negative id (12 cells to a grid, at |latitude| 0, 31.72 or 58.28) reads like the others: the sign adds nothing to the digits.
- **The seven codes** each have a grid id with the same 2-digit code whose centre is 0.11 to 2.29 km away, at a resolution where a cell is about 150 km across: GBIF's `GBIF_ISEA3HCode` and the file disagree on the
  centres of a few cells. Nothing else is off: the same grid at resolution 6 is exact.
- **The rule** (the plan's, written down before the run): an unknown code is **remapped** when all of these hold, and otherwise stays and fails the QA with its count as before: (1) it decodes to a place (|lat| at most 90,
  |lon| at most 180, an ASCII digit string of at most 19 digits), (2) a grid id with the same 2-digit code has its centre within `REMAP_KM = 5` km, (3) the next such id is more than `REMAP_RATIO = 3` times as far.
  **The code decides, not the nearest centre of any code**: a cell of another code, even 0.4 km away across the equator, is another cell and not a rounding of this one (tested). The hemisphere of a code is read from
  the grid itself (a code whose cells lie on both sides of the equator or the prime meridian gives no place, so it is never remapped) rather than a table, so it holds for any resolution. Distances are great-circle on a
  sphere of 6,371.0088 km, within 0.5% of the geodesic, which is all that a limit of 5 km and a ratio of 3 need.
- **Bounded, and listed.** At most `MAX_REMAPS = 100` distinct codes may be remapped in a download; one more stops the ingest with the count and removes its parts (more than that is another grid, not a
  rounding difference, and nothing is remapped). Every remap is printed and written to `qa.json` under `remaps` (code, target, rows, `n_raw`, `n`, km, and the totals and the rows merged into others), and nothing
  else is changed or dropped: a row with no cell, and a code that is no near miss, stay as they are and fail the QA.
- **Where.** At the ingest (`ingest_remapped`), so the parts never hold such a code and the loaders never see one; `ingest` as it was stays, and remaps nothing. `scripts/ingest.py` reads the cells and their centres
  (`realms.grid.read_centres`, a few milliseconds, no geometry) before it ingests and passes them in: it still loads the grid (a minute or more at resolution 7) after, so the ingest's peak memory is not the grid's.
  The rows of the remapped codes are held back while the parts are written (35 rows here), given their target cell, and written merged as the last part.
- **Rows that collide are merged as `load_presence` merges them.** A row of the download is one (species, country, cell, part), so a remapped row can have the key of a row already written, or of another remapped
  row: the counts add up (`n_raw`, `n`, `n_nc`, `n_open`, `n_ebird`, `n_inat`, `n_synonym`, `n_human_text`) and the four distinct counts (`n_datasets`, `n_years` and their `_open` forms) take the larger, which is a lower
  bound, as there. The row of an earlier part that is merged is taken out of that part (the part is streamed to a file of its own and renamed), so **no key is in two rows**. Tested against the loaders: a download with
  the near-miss codes replaced by their targets gives the same `load_presence` and `load_country_counts`, and the merged row is checked by hand, for a collision across parts, within a part, with a null country, and of
  two codes onto one cell. `merged_rows` is the rows of the download less the rows of the parts.
- **Not to be loosened:** the pass marks are as they were. The remap is the one exception to "a code the grid doesn't have stops the run", and it is bounded as above.

**The run** (`uv run python -I scripts/ingest.py data/gbif/0020744-261023074918325.zip --resolution 7`, exit 0, QA passes, `_COMPLETE` written; 7f0bbb6e4d0d012387a26436d2962043a38f41c54bdcaede620438a420ab7003 is the
zip's sha256, the same as `shasum -a 256`):

- **Rows:** 28,052,905 in the download, **28,052,882 in the Parquet**: 23 rows merged into others. The last part (`part-00006.parquet`) holds the 35 rows of the seven cells, 12 of which met no row and 23 of which
  merged with a row of the target. Every key is in one row (checked over all parts, 0 more than once).
- **Sizes:** the zip is 1,027,468,365 bytes (1.03 GB); the Parquet 391,561,410 bytes (392 MB) in 7 parts (the first run, before the change, wrote 391,387,694 bytes).
- **Time and memory:** the ingest took 6 s (5.2 s before the change) and the script 23 s; **peak 4,558 MB for the ingest** (read before the grid and the taxa are loaded; 3,879 MB before the change, and 4,353 MB
  in an earlier run of the merge as it is now: the figure moves by a few hundred MB between runs), 14,951 MB with the QA (`/usr/bin/time -l`: 15.7 GB maximum resident, 10.2 GB peak footprint).
  A first version of the merge read a whole part back to rewrite it and peaked at 5,423 MB; it now streams the part through the anti-join to a file and renames it.
- **QA:** unaccepted species 0, **unknown cells 0**, synonyms 100,324,056 of 3,381,828,966 clean records, 2.97%. `rows` fell by the 23, and **`n_raw` (3,405,278,365), `n` (3,381,828,966) and `n_synonym` are the same as in
  the failed run to the record**: nothing was lost. Zero-area parts as before (land 98 cells, 7,065 rows, 170,707 clean records; sea 2,534 cells, 633,274 rows, 2,902,211). 21,129 of the 21,872 cells hold records.
- **The seven remaps** (code in the download, the cell of the grid it went to, km between the code's own centre and that cell's, and the rows and the records of the code, all and clean):

| Code in the download | Remapped to | km | Rows | Records | Clean |
|---|---|---|---|---|---|
| `2902405731102034479` | `2902411730102033531` | 0.675 | 1 | 1 | 1 |
| `5100000000079533471` | `5100000000079532484` | 0.110 | 15 | 16 | 16 |
| `5121434201079619262` | `5121436312079621336` | 0.318 | 4 | 5 | 5 |
| `704809617102037511` | `704821168102036680` | 1.288 | 1 | 1 | 1 |
| `714370309102070552` | `714390873102070952` | 2.287 | 2 | 2 | 2 |
| `7300815160144987711` | `7300817749144998464` | 1.230 | 11 | 13 | 13 |
| `7302405731079534479` | `7302411730079533531` | 0.675 | 1 | 2 | 2 |
| | | | 35 | 40 | 40 |

  These are the targets and distances that were found by hand, to within 0.01 km, here by the rule: for each code the nearest cell of its code is the target and **the next is 155 km or more away**, 68 to 1,455 times
  as far (the rule asks for 3), so the 3-times test is not close for any of the seven.
- **The `_COMPLETE` of the res 6 ingest still verifies** (`read_complete` and `resolution_of(..., 6)`: 19,531,246 rows, zip sha256 `9fbba9b1...54bc03`, `how` `mark-complete`); res 6 was not ingested again.

**Growth from resolution 6 to 7, from global counts only** (distinct (species, cell) pairs, which is what the pilot's growth was counted in, and rows; the parts of both downloads, no species named):

| | Res 6 (`0019721-261023074918325`) | Res 7 (`0020744-261023074918325`) | Ratio |
|---|---|---|---|
| Rows | 19,531,246 | 28,052,882 (download: 28,052,905) | 1.436 |
| Distinct (species, cell) pairs | 15,340,244 | 23,504,905 | **1.532** |
| Distinct (species, cell, part) pairs | 17,572,638 | 26,369,126 | 1.501 |
| Species with a row | 1,290,444 | 1,290,499 | 1.000 |
| Cells with a row | 7,234 | 21,129 | 2.921 |

The pairs grow 1.53 times for 3 times the cells. The pilot's two figures for the growth were Denmark's 1.38 (a lower bound, as Denmark spans 8 cells at resolution 6) and 2.5, which decided the go to
resolution 7: the measure is between them and nearer the lower, so the decision stands on the more demanding number. Nothing is changed here; the pair growth is a global count and says nothing of the links, which
the pilot scaled from half A.

**Tests** (`tests/test_occurrences.py`, `tests/test_grid.py`): the plan's four (a near miss goes to the cell with its code and nearest centre; no cell within tolerance, an ambiguous near miss, still fail the QA; the remaps are
counted in `qa.json`), the limits at 5 km (4.9 and 5.1) and at 3 times (2.5, 2.9, 3.2 and 4.0), a decoy cell of another code with the same absolute coordinates and one across the equator, a negative id, a code that
is not a place (not digits, too long, past the pole or the antimeridian, a code the grid doesn't have), a code in the grid, two cells of one centre, a row with no cell, the cap, a failed run leaving no parts, the
merges above, the printed remaps and `read_centres`. Mutations were caught: a ratio of 2, a limit of 6 km, ignoring the hemisphere, adding the distinct counts, not removing the merged rows from their part, and
joining with nulls unequal.

### Calibration dry run and known cases on half A at resolution 7 (2026-11-01)

*History: the run of the 81 settings, with its rule and its rulings still open. Markus's rulings, the knee rule, the 144-setting grid and the dry run that replaces this one are in "The knee rule, the 144-setting grid and the resolution 7 dry run" at the end of the log; this run's files are kept in `results/calibrate/res7-grid81/`.*

`uv run python scripts/calibrate.py --counts data/gbif/0020744-261023074918325`, half A only, with no `--score-half-b`, so nothing of half B was loaded. It ran
in 666 s (layer 5 at 81 settings 499 s) and peaked at 21.7 GB. The outputs are in `results/calibrate/res7-grid81/` (git-ignored; moved there when the 144-setting run replaced them in `res7/`, which is what `--score-half-b --setting N` checks
for), and the run waited on Markus's choices below.

- **Within the 5% limit** (whole run of layers 1, 2, 4 and 5, drops only, the provisional reading): 58 of the 81 settings. By layer 5 alone, all 81 are within the limit.
- **Chosen by the rule** (the most outside presences removed within the limit): setting 7, with record share 0.5%, cell share 2%, rate factor 2% and group minimum 1,000. It removes 24.28% of the
  outside presences (211,947 of 872,784) and loses 4.99% of the native ones (151,232 of 3,029,872). It sits at an edge of the grid again (the smallest record and cell shares, the largest rate
  factor), so the grid doesn't bracket the optimum.
- **The design setting** (40: 1%, 5%, 1%, 1,000): 22.73% and 3.85%. **Resolution 6's choice** (75: 2%, 10%, 1%, 300): 23.75% and 4.00%. Setting 7 buys 0.5 points more outside removal than 75 for
  1.0 point more native loss. By layer 5's own reading the rule would take 78 (2%, 10%, 2%, 300), which loses 5.33% over the whole run.
- **What drives it:** the rate factor of `under_reported`. Averaged over the other thresholds, 0.5%, 1% and 2% give 20.1%, 22.3% and 25.3% removed outside, against 3.0%, 3.8% and 5.1% lost
  native. Record share, cell share and group minimum move both by under 2 points, and the group minimum by under 0.2.
- **Native loss by region** (TDWG level 3, at least 100 native presences, chosen setting, layer 5's reading): the top is heavily sampled mainland, France 17.1%, California 14.0%, New South Wales 13.5%,
  Spain 13.1%, Victoria 12.9%, Germany 12.0%, besides a few island groups (Marianas 16.7%, Caroline Is. 13.7%, Cape Verde 12.9%). Regions of at most 5 land parts lose 3.1% in all, larger ones 3.5%,
  so islands don't lose more as a rule. The sampled-mainland pattern points at `under_reported` taking rare natives in well-sampled cells, which the rate factor controls.
- **Known cases** (all five layers, `native_only`, `griis_without_status="skip"`):

  | Case | Result |
  |---|---|
  | Potato, outside its WCVP native regions | Pass: 1,024 presences outside before, 0 after, 110 of 127 inside kept |
  | Starling, North America | Fail: 1,045 before, 43 after |
  | Hawaiian endemic, *Drepanis coccinea* | Pass: 11 of 11 kept (at resolution 6 one of 9 was lost) |
  | Lemuridae, only on Madagascar | Fail: 56 outside before, 10 after |
  | Tiger, Europe | Waits for step 9 (Felidae is in half B) |

- **GRIIS cross-check, animals, chosen setting:** layer 5 takes 12.77% of the animal presences GRIIS lists as introduced (9,067 of 70,987), against 6.34% of the others (405,587 of 6,399,166). At
  the design setting the figures are 11.29% and 5.83%.

**For Markus, before `--score-half-b`:** confirm the 5% reading; choose N (the rule's 7, resolution 6's 75, or the design's 40); and decide whether to widen the grid towards smaller rate factors first.
The scoring refuses unless this run's `report.json` exists, and it records whether N is the dry run's choice.

### Two GRIIS datasets with one alias, and every GRIIS list with a table (2026-11-01)

**The defect.** ChecklistBank gives the alias `GRIIS-KY` to two datasets: 28941, "Global Register of Introduced and Invasive Species - Cayman Islands" (the country list, 85 species entries), and 28940,
"... - Cayman Brac, Cayman Islands" (a sub-national list of 8). **28940's alias is wrong at the source** (its siblings 28939 `GRIIS-KY-Grand` and 28938 `GRIIS-KY-Little` are right); Markus can fix it in ChecklistBank. Both wrote
`exports/extended/griis/GRIIS-KY.zip`, so `fetch_reference.py` stopped with "incomplete download ... GRIIS-KY.zip has 5108 bytes, expected 10615" and never rewrote the manifest, whose `griis/GRIIS-KY`
entry pointed at 28940. **Until now Cayman Brac's 8 species stood for the Cayman Islands' country list**, so they were taken as alien in the whole of KY and the country list's own 85 were not applied at all. The
runs that read that manifest (check_natives, clean_counts, calibrate, the dry run on half A) have the wrong list for KY.

**How it is handled.** `reference.griis_datasets` never returns two datasets with one alias. Where aliases collide it reads the titles (`<register> - <Country name>` for a country list, `<register> - <Place>, <Country name>` for a
sub-national one, with whichever dash the title uses): exactly one of the colliding datasets must be the country list, whose part after the dash is what every other title ends with, and it keeps the alias.
The others become `GRIIS-<country code>-<Place without spaces>` with the place as subdivision: **28940 is now `GRIIS-KY-CaymanBrac`**, subdivision "Cayman Brac". `natives` matches a subdivision folded to letters and digits, so `CaymanBrac` and
`Cayman Brac` are one name (tested). **Anything else raises a `ValueError` that names every key and title**, and nothing is picked: two country lists, two places and no country list, a place of another country,
a title with no part after the dash or none at all, two places that fold to one alias, a place of two parts (`North, East, Zedland`), a collision of an alias that already has a subdivision, and an alias made up that another dataset has. Natural Earth has one polygon for
KY, "Cayman Islands", and none for Cayman Brac, so **`KY-CaymanBrac` is unplaced** on the grid, as `KY-Grand` and `KY-Little` were. That is expected. Tests: `tests/test_reference.py` (the KY case, eight unresolvable collisions, the made-up alias, no collision unchanged,
a protected-areas list that is not a collision) and two in `tests/test_natives.py`; the suite has 1,030 tests, all passing.

**The fetch.** The stale `exports/extended/griis/GRIIS-KY.zip` held 28940's export, so it would have been read as the country list under the new manifest and it blocked the fetch (a file at its path with another size raises), and it was removed. It was the only file removed:
the exports are read only through the manifest's `path` (`natives._manifest_path`, `check_natives.py`), so the 250 plain exports that are no longer referenced are never read, among them `plain/griis/GRIIS-KY.zip`, which is 28941's own. The first rerun of `scripts/fetch_reference.py` completed:
**317 items, 0 missing**, and no incomplete download. `manifest.json` lists `griis/GRIIS-KY` (28941, 10,615 bytes) and `griis/GRIIS-KY-CaymanBrac` (28940, 5,108 bytes) apart, with `GRIIS-KY-Grand` and `GRIIS-KY-Little`. **All 305 GRIIS lists are extended exports**, with 305 different datasets
and paths. The `key` in each zip's `metadata.yaml` is the manifest's dataset and each has a `Distribution.tsv`.

**GRIIS after this, and before** (`scripts/check_natives.py -I`, species-rank entries, resolution 6, COL 315557, 16.3 GB peak, 54 s):

| | Before | Now |
|---|---|---|
| Lists | 304 in the manifest, 303 with entries | 305, all with entries |
| Lists with a distribution table | 55 | **305** |
| Lists without one | 248 | 0 |
| Species entries in lists with a table | 52,462 | **132,479** |
| Introduced | 48,972 | 127,338 |
| Native (left out) | 3,457 | 5,077 |
| Uncertain (left out) | 5 | 13 |
| No value, counted as introduced | 28 (all GRIIS-TR) | 51 (GRIIS-TR, NU, TO-Haapai) |
| Pairs applied (`References.griis`) | 46,519 of 20,556 species in 53 lists | **115,371 of 26,734 species in 303 lists** |
| Name-only pairs (`griis_name_only`) | 62,375 of 14,170 species in 248 lists | 0 |
| Unplaced subdivisions (no admin-1 cells) | 61: 6 with a table (1,630 pairs) and 55 with none (5,190 pairs) | **62 with 6,569 pairs**, all with a table: the 61 and `KY-CaymanBrac` |
| Placed subdivisions | 24 with 263 cells | 24 with 263 cells |

The 2 lists that apply nothing are `GRIIS-EC` (731 of 731 entries native) and `GRIIS-SR` (176 of 176), as before. Native entries remain the largest in GRIIS-EC, US-Hawaii (440 of 5,542), US-United (429 of 8,349), IC-Tenerife (423 of 1,119), IC-Gran-Canaria (287 of 777) and IC-Palma (252 of 661).
**`griis_without_status="apply"` no longer has anything to apply**: no list is without a table. GRIIS names matched to COL 92.0% after GBIF's fallback (7.5% unmatched, 140 ambiguous; 91.3% matched and 2,607 unmatched before); 3 lists are over 10% unmatched (IL, MH-Namdrik-2, TO-Haapai),
GBIF's match API was asked about 103 names and placed 78 (132 and 97 before), and `gbif_fallback.parquet` was rewritten by the script.

**The extended exports are not only the old lists with a table added.** The 250 plain exports were of older versions of the datasets: GRIIS-CY's plain one is version 4 and the extended one version 9, and 559 of the plain one's entries are ranked `other`, which the species rule drops (23 species entries
against 860 now). The 250 lists hold 80,017 species entries now against 68,170 in their plain exports (+11,847), and 85 more come from the Cayman Islands list taking the slot Cayman Brac's 8 held: the 132,479 of now against the 120,547 of before.
Nothing is cleaned or calibrated here. The numbers of the section "GRIIS by its own establishment means" above, of Task 11a and of the dry run before this one, were made with the old manifest, and the GRIIS rules of the cleaning and the calibration are to be run again on this one.

### The knee rule, the 144-setting grid and the resolution 7 dry run (2026-11-01)

*History in its rule and its choice: the walk along the Pareto frontier itself, with `MAX_MARGINAL = 4`, which chose setting 36 and stopped at a dent of the grid. Markus replaced it the same day by a walk along the convex hull at break-even, in "The knee walks the convex hull, at break-even" at the end of the log. His other decisions (reading c, B's calibration species, the grid), the fixture, the dry-run table and the known cases stand; the choice, the frontier table, the GRIIS cross-check at the chosen setting and the region table below are those of setting 36 and are redone for the new choice there.*

**Markus's decisions on 2026-11-01**, made after the res 7 dry run of the 81 settings (above, kept as history):

- The 5% limit is read over the whole run of layers 1, 2, 4 and 5, counting removals only (reading c). A moved presence is not lost. This is no longer provisional.
- Half B's scoring loads B's calibration species, meaning the WCVP plants plus the species of their tax groups (2,586 of them not Plantae). That is confirmed.
- **The choice rule changes to a knee rule.** Build the Pareto frontier of the settings on (native presences lost, outside presences removed), both under reading c. Sort by native lost ascending, then outside removed descending, then
  setting index ascending. Keep each row that removes more outside presences than every row before it. Walk the frontier from the lowest native loss. Accept a step while it costs at most `MAX_MARGINAL = 4.0` native presences lost per extra
  outside presence removed, and the setting's native loss stays within the 5% limit. Stop at the first step that fails either test. The chosen setting is the last one accepted. On the 81-setting res 7 table this rule gives setting 75 (2%, 10%,
  1%, 300); the step after it, to setting 8, costs 9.1 natives per non-native.
- **The grid widens to 144 settings:** satellite record share {0.5, 1, 2}%, satellite cell share {2, 5, 10, 20}%, under-reported rate factor {0.5, 1, 1.5, 2}%, group minimum {300, 1,000, 3,000}.

**The code** (`src/realms/calibrate.py`). `GRID` is the 144 settings, in the order of the grid (the first threshold varies slowest, the group minimum fastest), with the design's setting at **index 64** (the second value of each).
`MAX_MARGINAL = 4.0`. `frontier(table, reading)` is the Pareto frontier as above, with each row's step (`extra_outside`, `extra_native`) and `marginal`, the native presences it costs per extra outside presence removed. `knee(...)` walks it and
returns the chosen `setting`, whether the limit was `met`, and the frontier with three flags per row: `within_limit`, `affordable` (the step into it costs at most `MAX_MARGINAL`, compared as counts) and `accepted` (it and every row before it are both). `choose` returns the knee's
choice by default, and `choose(rule="most_removed")` the rule of before, for comparison. If not even the lowest-loss setting is within 5%, that one is taken and the caller is told. `chosen.setting` in `report.json` stays an integer, an index into the new grid, so the scoring's guard
reads it as before. `report.json` has, besides `chosen`, `design_setting` and `chosen_by_layer5_reading` (the knee by layer 5's reading): **`knee_rule`** (the whole frontier with its flags, `step_into` the chosen setting, `step_after` it and `stopped_by`, `cost` or `limit`) and
**`chosen_by_old_rule`**. Tests: the five of the brief (`test_frontier_keeps_only_settings_that_remove_more_for_their_loss`, `test_knee_rule_stops_at_the_first_step_dearer_than_max_marginal`, `test_knee_rule_never_passes_the_native_loss_limit`,
`test_knee_rule_on_the_resolution_7_table_of_81_chooses_75`, `test_the_grid_has_144_settings_in_a_stable_order`) and three more, among them one for the report of the dry run; the suite has 1,037 tests, all passing. The fixture
`tests/fixtures/res7_grid81_settings.csv` is the aggregate counts of half A's 81 settings at resolution 7 (thresholds and the counts of both readings: no taxon, no cell), and the rule picks 75 on it, as computed by hand; the old rule picks 7.

**The dry run at resolution 7** (`uv run python scripts/calibrate.py --counts data/gbif/0020744-261023074918325`, half A only, no `--score-half-b`; the references with all 305 GRIIS lists as extended exports and Cayman Brac apart from the Cayman Islands). 1,080 s
(layer 5 at 144 settings 905 s), peak 21.6 GB. The 81-setting outputs were moved to `results/calibrate/res7-grid81/` first; the new ones are in `results/calibrate/res7/` (git-ignored), which the scoring of half B reads.
**No `held_out*` file exists under `results/calibrate/` or `data/split/`, before or after the run, and half B was not loaded.**

- **Within the 5% limit** (reading c): 112 of the 144 settings; by layer 5 alone, all 144.
- **Chosen by the knee rule: setting 36**, with satellite record share 0.5%, satellite cell share 20%, under-reported rate factor 0.5% and group minimum 300. It sits in a corner of the grid: the lowest record share, rate factor and group minimum and the highest cell share. Its three readings:

  | Reading | Outside presences removed | Native presences lost |
  |---|---|---|
  | whole run, removals only (`run_drops`, the one that chooses) | 20.49% (178,804 of 872,784) | 3.03% (91,721 of 3,029,872) |
  | whole run, everything counted (`run_all`) | 20.49% (178,818 of 872,784) | 3.03% (91,836 of 3,029,872) |
  | layer 5 alone | 14.07% (115,086 of 817,911) | 1.71% (51,510 of 3,012,793) |

  The knee rule by layer 5's reading also takes 36.
- **The frontier** has 81 of the 144 settings, 70 of them within the limit. The walk accepted 14, and stopped at the 15th, whose step costs 17.4 natives for each outside presence. Around the chosen setting (whole-run counts; a step's marginal cost is `extra native / extra outside`):

  | Setting (record share, cell share, rate factor, group min) | Outside removed | Native lost | Step: extra outside, extra native | Marginal | |
  |---|---|---|---|---|---|
  | 26 (0.5%, 10%, 0.5%, 3,000) | 177,462 | 91,279 | +2,531, +644 | 0.25 | accepted |
  | 38 (0.5%, 20%, 0.5%, 3,000) | 178,363 | 91,429 | +901, +150 | 0.17 | accepted |
  | 37 (0.5%, 20%, 0.5%, 1,000) | 178,749 | 91,684 | +386, +255 | 0.66 | accepted |
  | **36 (0.5%, 20%, 0.5%, 300)** | **178,804** | **91,721** | **+55, +37** | **0.67** | **accepted, chosen** |
  | 62 (1%, 5%, 0.5%, 3,000) | 178,869 | 92,850 | +65, +1,129 | **17.37** | refused (cost) |
  | 61 (1%, 5%, 0.5%, 1,000) | 179,266 | 93,104 | +397, +254 | 0.64 | not reached |
  | 60 (1%, 5%, 0.5%, 300) | 179,323 | 93,142 | +57, +38 | 0.67 | not reached |

  Two more steps within the limit are dearer than 4, **17** (10.0) and **44** (25.4), and two beyond it, 70 (4.1) and 119 (62.7); the walk meets 62 first.
- **For Markus: what stops the walk.** The step to 62 is a single dear step, and the frontier goes on cheaply after it (61, 60, 110, 109, 108, 74, ... at 0.7 or less each). Counted from 36, **every later frontier setting within the limit costs at most 3.0 natives per extra
  outside presence removed over the whole way (61: 2.99, 60: 2.74, the others 1.6 or less), except 62 itself (17.4)**. The last frontier setting within the limit is **138** (2%, 20%, 1.5%, 300): 224,288 outside presences removed (25.70%), 143,410 native lost (4.73%), which is 1.14 natives
  per extra outside presence from 36 and is also what the old rule chose. The frontier is made of rungs of three settings (group minimum 3,000, 1,000, 300) of one family of the other three thresholds, and the step from one family to the next lands on the new family's first rung, here 62, which removes 65 more outside presences than 36 and loses 1,129
  more native ones: a dent in the frontier, and the first-failure rule stops at the dent. The rule is implemented as Markus gave it, and the choice stands at 36 until he says otherwise. The chosen setting removes less than the design's (64, below) while losing less: from 36 to 64 it is 1.28 natives per extra outside presence.
- **The old rule** (the most outside presences removed within 5%, for comparison): **setting 138** (2%, 20%, 1.5%, 300), with 25.70% (224,288) outside removed and 4.73% (143,410) native lost over the whole run (run_all 25.70% and 4.74%, layer 5 alone 19.75% and 3.46%).
  Its step from 36 is the one above. **The design's setting** (64: 1%, 5%, 1%, 1,000, on the frontier): 22.73% (198,417) and 3.85% (116,786) (run_all 22.74% and 3.86%, layer 5 alone 16.52% and 2.56%).
  The knee rule's choice on the 81-setting table, 75 (2%, 10%, 1%, 300), is **123** of the 144, on the frontier: 207,289 removed (23.75%), 121,246 lost (4.00%). The old rule's choice there, 7 (now 10), is not on the frontier.
- **What drives it:** the rate factor of `under_reported` again. Averaged over the other thresholds, rate factors 0.5, 1, 1.5 and 2% give 20.4%, 22.6%, 24.2% and 25.6% outside removed against 3.05%, 3.83%, 4.52% and 5.15% native lost, and **only 4 of the 36 settings at 2% are within the limit**.
  Cell shares 2, 5, 10 and 20% give 21.8% to 24.0% and 3.97% to 4.23%, record shares 22.7% to 23.6% and 4.06% to 4.20%, and the group minimum 23.1% to 23.2% and 4.12% to 4.15%.
- **Known cases** (all five layers, `native_only`, at the design's thresholds, so they do not depend on the choice; what changed since the 81-setting run is the refreshed GRIIS lists):

  | Case | Result |
  |---|---|
  | Potato, outside its WCVP native regions | Pass: 1,024 presences outside before, 0 after, 110 of 127 inside kept (as before) |
  | Starling, North America | Fail: 1,045 before, **7 after (43 before the refresh)**: Greenland 5 and no region 2; GRIIS took 404 and 527 sub-national, the zero-area fold 104 |
  | Hawaiian endemic, *Drepanis coccinea* | Pass: 11 of 11 kept (as before) |
  | Lemuridae, only on Madagascar | Fail: 56 outside before, **8 after (10 before)**, 178 of 217 inside kept; left in North Carolina and Virginia 2 each, Comoros, Denmark, Poland 1 each, 3 with no region |
  | Tiger, Europe | Waits for step 9 (Felidae is in half B) |

- **GRIIS cross-check, animals** (layers 1, 2, 4 and 5, layer 3 off; the lists with a distribution table are now all 305): at the chosen setting layer 5 takes **8.79%** of the animal presences GRIIS lists as introduced (9,103 of 103,595) against **4.70%** of the others (299,079 of 6,366,558), 1.87 times as much;
  records 0.03% against 0.03%. At the design's setting 10.50% (10,876) against 5.81% (369,940), 1.81 times. Before the refresh, at the 81-setting run's choice: 12.77% of 70,987 against 6.34% (2.01 times), and 11.29% against 5.83% at the design's. The listed presences are 103,595 now, as the lists with a table grew from 55 to 305.
- **Native loss by region** (TDWG level 3, at least 100 native presences, 325 of 366 regions, chosen setting, layer 5's reading, whose total loss is 1.71%; the median region loses 0.99%). Regions of at most 5 land parts (56) lose 1.43% (1,152 of 80,372), larger ones (269) 1.50% (61,964 of 4,123,626), so islands do not lose more as a rule,
  though island groups head the list besides the sampled mainland. The top 20 (a presence counts in each region its part touches, so the rows do not add up):

  | Region | Land parts | Native presences | Lost | Loss |
  |---|---|---|---|---|
  | Marianas | 6 | 251 | 43 | 17.1% |
  | Caroline Is. | 15 | 642 | 94 | 14.6% |
  | France | 38 | 57,815 | 7,230 | 12.5% |
  | Cape Verde | 3 | 434 | 53 | 12.2% |
  | Andaman Is. | 7 | 271 | 31 | 11.4% |
  | Ogasawara-shoto | 2 | 109 | 12 | 11.0% |
  | Marquesas | 5 | 267 | 28 | 10.5% |
  | Libya | 94 | 1,144 | 100 | 8.7% |
  | Madeira | 2 | 629 | 50 | 7.9% |
  | Germany | 27 | 32,534 | 2,560 | 7.9% |
  | Maldives | 7 | 163 | 12 | 7.4% |
  | Tubuai Is. | 3 | 111 | 8 | 7.2% |
  | Norfolk Is. | 2 | 224 | 16 | 7.1% |
  | Seychelles | 2 | 211 | 15 | 7.1% |
  | Chatham Is. | 2 | 322 | 22 | 6.8% |
  | Oman | 28 | 2,375 | 160 | 6.7% |
  | Sri Lanka | 8 | 2,431 | 163 | 6.7% |
  | Spain | 35 | 54,592 | 3,206 | 5.9% |
  | Austria | 11 | 14,680 | 859 | 5.9% |
  | California | 33 | 38,840 | 2,270 | 5.8% |

**Not done here:** the scoring of half B (Step 6, after review). The setting the rule chose is 36, and `res7/report.json` has `chosen.setting` 36 from the same ingest, so `--score-half-b --setting N` would find its guard satisfied.

### The knee walks the convex hull, at break-even (2026-11-01)

**Markus's decision, 2026-11-01**, after he had seen the walk along the Pareto frontier stop at setting 36 on a dent of the grid: **a walk along the upper convex hull, with `MAX_MARGINAL = 1.0`**, that is break-even, one native presence lost for each non-native presence removed.
Reading c, B's calibration species and the 144-setting grid stand as decided above.

**Why the hull.** The frontier of the 144 settings is smooth as a whole (81 settings, 70 of them within the limit, and from 36 every later within-limit setting but one costs at most 3.0 natives per extra outside presence), but it has dents of the grid: single dear steps, 62 (17.4), 17 (10.0)
and 44 (25.4) within the limit, where the lowest-loss rung of a new family of the other three thresholds joins it, after which the cost falls again. A walk from one frontier setting to the next ends at the first of them, which was 62, 65 more outside presences for 1,129 more natives, and with
it 5 points of outside removal short of settings that cost about one native for one non-native. The hull steps over a dent (its point is below the line between its neighbours), and along a hull the cost of a step only rises, so that stopping at the first step that costs more than the price is exact:
there is no cheaper step after it.

**The rule.** (1) Keep the settings that lose at most 5% of the native presences (reading c). (2) Build their Pareto frontier as before: sort by native lost ascending, outside removed descending, setting ascending, and keep each row that removes more outside presences than every row before it.
(3) Take the upper convex hull of that frontier in (native lost, outside removed): a point on or below the segment between its neighbours is dropped, collinear points with it. (4) Walk the hull from the lowest native loss and accept each step while (extra native lost) / (extra outside removed) is at most `MAX_MARGINAL`,
compared as counts with a tolerance of 1e-9. The chosen setting is the last vertex accepted. The settings past the limit are dropped **before** the hull is built, so a point beyond the limit cannot pull a point within it off the hull. If no setting is within the limit, the one that loses the fewest native
presences is taken (the first row of the frontier: of those that tie, the one that removes most), and the caller is told, as before.

**The code** (`src/realms/calibrate.py`). `frontier` is split into `pareto` (the rows), `with_marginals` (the step into each row from the one before, and its marginal cost) and `frontier` (the two in turn), and `hull(points)` is new. `knee` returns a `Knee` with the chosen `setting` and `met`, the whole Pareto frontier
as `steps` (with `within_limit` and `on_hull` per row, so a dent such as 62 stays in sight), the `hull` (with each vertex's step, `affordable` and `accepted`), and **the `limit`, `reading` and `max_marginal` it was made with**, which `report.json` copies and no longer hard-codes. `step_into()` and
`step_after()` are of the hull. `MAX_MARGINAL = 1.0`; the rule is `hull_knee` (the default of `choose`), and `choose(rule="most_removed")` stays for comparison. The first form of the knee, which walked the frontier and which this replaces, has no rule name any more. `report.json`'s `knee_rule` has `rule`
`"hull_knee"`, `reading`, `max_marginal`, `native_loss_limit`, `limit_met`, `frontier` (all 81 rows with `on_hull`), `hull`, `step_into`, `step_after` and `stopped_by` (`"cost"`, or `None` if the walk reached the end of the hull, where the limit ended it).
Tests (the suite has 1,051, all passing): the hull drops a dent point and collinear points; the walk is along the hull and a dent does not stop it; the price, the limit and the fallback; `step_after` and `stopped_by` are `None` at the end of the hull; both tables of counts, **the 81-setting one gives 72** (0.02, 0.10, 0.005, 300; hull
[2, 20, 47, 74, 73, 72, 75, 7], the step after, to 75, costs 1.27) **and the 144-setting one 132** (hull [2, 38, 86, 134, 133, 132, 135, 138], the step into 132 from 133 costs 0.71, the step after, to 135, 1.29); at a price of 4 the walk gives 75 on the 81-setting table and 138 on the 144-setting one. The fixtures are `tests/fixtures/res7_grid81_settings.csv`
and `tests/fixtures/res7_grid144_settings.csv`, aggregate counts of half A with the thresholds, no taxon or cell.

**The one-shot guard.** `dry_run_choice` (the check that `--score-half-b` makes before it loads anything) now also refuses a dry run of another grid: its chosen setting must be an index of the 144 and have the four thresholds that the grid's setting of that number has (a number means something only in the grid it was made on),
and `settings.csv` next to the report must exist and have one row for each of the 144 settings. Tested for a table of 81 rows, no table, other thresholds, no thresholds at all, and a number outside the grid (above and below). Checked on the real `res7/`: it accepts 132; it refuses a copy of the report with 132's thresholds changed, one with the number
changed to 75 and one next to the 81-row table of the first run. **It does not refuse another `--setting N`:** `score_half_b` takes any N of the 144, and records `dry_run_chosen` and `setting_is_dry_run_choice` (false for any N but 132) in both markers and in `held_out.json`, which also has the thresholds, which is the design of the run
(N is Markus's choice after he has read the dry run). Nothing was loosened, and the guard was not changed on that point.

**The dry run again** (`uv run python scripts/calibrate.py --counts data/gbif/0020744-261023074918325`, half A only): 1,101 s (layer 5 at 144 settings 922 s), peak 22.3 GB. **The table of the 144 settings is the same as the first run's** (checked against the fixture), so only what depends on the choice differs below; the known cases are at the design's thresholds and are as in
the section above. No `held_out*` file exists under `results/calibrate/` or `data/split/`, before or after, and half B was not loaded. `results/calibrate/res7/report.json` has `chosen.setting` 132.

- **Chosen by the hull walk at break-even: setting 132**, satellite record share 2%, satellite cell share 20%, under-reported rate factor 0.5% and group minimum 300. Its three readings:

  | Reading | Outside presences removed | Native presences lost |
  |---|---|---|
  | whole run, removals only (`run_drops`, the one that chooses) | 21.99% (191,887 of 872,784) | 3.26% (98,636 of 3,029,872) |
  | whole run, everything counted (`run_all`) | 21.99% (191,901 of 872,784) | 3.26% (98,751 of 3,029,872) |
  | layer 5 alone | 15.71% (128,453 of 817,911) | 1.94% (58,555 of 3,012,793) |

  The hull walk by layer 5's reading also takes 132. **The hull** (whole run, removals only; 8 vertices of the 81 frontier settings; a step's cost is `extra native / extra outside`):

  | Setting (record share, cell share, rate factor, group min) | Outside removed | Native lost | Step: extra outside, extra native | Cost | |
  |---|---|---|---|---|---|
  | 2 (0.5%, 2%, 0.5%, 3,000) | 163,804 | 86,936 | | | start |
  | 38 (0.5%, 20%, 0.5%, 3,000) | 178,363 | 91,429 | +14,559, +4,493 | 0.31 | accepted |
  | 86 (1%, 20%, 0.5%, 3,000) | 185,912 | 95,360 | +7,549, +3,931 | 0.52 | accepted |
  | 134 (2%, 20%, 0.5%, 3,000) | 191,478 | 98,345 | +5,566, +2,985 | 0.54 | accepted |
  | 133 (2%, 20%, 0.5%, 1,000) | 191,835 | 98,599 | +357, +254 | 0.71 | accepted |
  | **132 (2%, 20%, 0.5%, 300)** | **191,887** | **98,636** | **+52, +37** | **0.71** | **accepted, chosen** |
  | 135 (2%, 20%, 1%, 300) | 210,308 | 122,384 | +18,421, +23,748 | **1.29** | refused (cost) |
  | 138 (2%, 20%, 1.5%, 300) | 224,288 | 143,410 | +13,980, +21,026 | 1.50 | not reached |

  The walk ends at 132 on the price: the step to 135 loses 1.29 natives for each outside presence it removes, and the end of the hull, 138, is the last setting within the limit (and the old rule's choice). **The dent rows are in `report.json`'s frontier and not on the hull:** 62 (1%, 5%, 0.5%, 3,000) is 65 outside for 1,129 natives from 36, and 36 itself, the choice of the walk
  along the frontier (20.49% and 3.03%), is under the hull too. From 36 to 132 the cost is 0.53 natives per extra outside presence (13,083 outside for 6,915 natives). At a price of 4 the walk would reach 138. **The design's setting (64: 1%, 5%, 1%, 1,000)**, 22.73% (198,417) and 3.85% (116,786) (layer 5 alone 16.52% and 2.56%), is on the frontier and off the hull. **The old rule's** 138
  (2%, 20%, 1.5%, 300): 25.70% (224,288) and 4.73% (143,410).
- **GRIIS cross-check, animals** (layers 1, 2, 4 and 5, layer 3 off): at the chosen setting layer 5 takes **9.12%** of the animal presences GRIIS lists as introduced (9,451 of 103,595) against **5.71%** of the others (363,252 of 6,366,558), 1.60 times; records 0.03% against 0.04%. At the design's setting 10.50% against 5.81%, 1.81 times, as before.
- **Native loss by region** (TDWG level 3, at least 100 native presences, 325 of 366 regions, setting 132, layer 5's reading, whose total loss is 1.94%; the median region loses 1.16%). Regions of at most 5 land parts (56) lose 1.59% (1,279 of 80,372), larger ones (269) 1.70% (69,937 of 4,123,626), so islands still do not lose more as a rule, though island groups head the list
  with the sampled mainland. The top 20 (a presence counts in each region its part touches, so the rows do not add up):

  | Region | Land parts | Native presences | Lost | Loss |
  |---|---|---|---|---|
  | Marianas | 6 | 251 | 45 | 17.9% |
  | Caroline Is. | 15 | 642 | 108 | 16.8% |
  | Andaman Is. | 7 | 271 | 40 | 14.8% |
  | Marquesas | 5 | 267 | 36 | 13.5% |
  | Cape Verde | 3 | 434 | 55 | 12.7% |
  | France | 38 | 57,815 | 7,232 | 12.5% |
  | Ogasawara-shoto | 2 | 109 | 12 | 11.0% |
  | Libya | 94 | 1,144 | 118 | 10.3% |
  | Chatham Is. | 2 | 322 | 32 | 9.9% |
  | Oman | 28 | 2,375 | 197 | 8.3% |
  | Tubuai Is. | 3 | 111 | 9 | 8.1% |
  | Seychelles | 2 | 211 | 17 | 8.1% |
  | Maldives | 7 | 163 | 13 | 8.0% |
  | Madeira | 2 | 629 | 50 | 7.9% |
  | Sudan | 138 | 4,331 | 344 | 7.9% |
  | Germany | 27 | 32,534 | 2,560 | 7.9% |
  | Sri Lanka | 8 | 2,431 | 187 | 7.7% |
  | Norfolk Is. | 2 | 224 | 17 | 7.6% |
  | Fiji | 14 | 1,261 | 93 | 7.4% |
  | Society Is. | 4 | 424 | 30 | 7.1% |

**Not done here:** the scoring of half B. Setting 132 is the setting the rule chooses, and `res7/report.json` has it, from the same ingest, so `--score-half-b --setting 132` would find its guards satisfied.

### Half B's held-out rates: the one unsealed call (2026-11-01)

The one `unsealed=True` call before the freeze ran once:
`uv run python scripts/calibrate.py --counts data/gbif/0020744-261023074918325 --score-half-b --setting 132`, at resolution 7, after the dry run on half A above chose
setting 132 by the hull walk at break-even (Markus, 2026-11-01). It loaded the references and half B's calibration species only, meaning the WCVP plants plus the species of their tax groups:
197,567 species, 194,981 Plantae and 2,586 not. It printed and saved counts and rates only, nothing per taxon. It took 37 s and peaked at 17.9 GB.

The fixed marker, `data/split/held_out-315557.started`, copied here as the durable record (`results/calibrate/res7/held_out.started` is identical):

```
started 2026-11-01T16:19:27+00:00
commit 792815c605c4481bba6f2fd1c9cbd512cdc5074d
setting 132
dry_run_chosen 132
setting_is_dry_run_choice True
counts data/gbif/0020744-261023074918325
zip_sha256 7f0bbb6e4d0d012387a26436d2962043a38f41c54bdcaede620438a420ab7003
```

**Held-out rates of setting 132** (record share 2%, cell share 20%, rate factor 0.5%, group minimum 300). Half A's are from the dry run that chose it:

| Reading | Half B: outside removed | Half B: native lost | Half A: outside removed | Half A: native lost |
|---|---|---|---|---|
| Whole run, layers 1, 2, 4 and 5, removals only (the binding reading) | 24.53% (153,880 of 627,285) | **3.19%** (82,129 of 2,572,615) | 21.99% | 3.26% |
| Layer 5 alone | 17.56% (102,135 of 581,492) | 1.97% (50,456 of 2,560,544) | 15.71% | 1.94% |

The setting holds on the held-out half. Native loss is 3.19% against A's 3.26%, within the 5% limit, and it removes slightly more of B's outside presences than A's. Nothing points to the thresholds being fitted to A's plants. These are the remaining error rates of the cleaning that the design asks to report.

From now on half B stays sealed until M3. Both markers stop any second scoring. A retry would need Markus's explicit decision, and none is planned.

## 2026-11-02

### Setting 132 as the defaults, frozen GRIIS country codes, and the first cleaning at resolution 7 (2026-11-02)

**The thresholds.** Markus chose setting 132 of the 144 on 2026-11-01 (the hull walk at break-even on half A at resolution 7, scored once on half B: 3.19% of B's native presences lost against A's 3.26%, see "Half B's held-out rates"). **The defaults of `Settings` are now setting 132**, so that every run of M1 uses the
calibrated thresholds without having to say so, and a run that forgets to pass them does not fall back on values the calibration has since replaced:

| Threshold | Was (the design's) | Is (setting 132) |
|---|---|---|
| `satellite_record_share` | 0.01 | **0.02** |
| `satellite_cell_share` | 0.05 | **0.20** |
| `under_rate_factor` | 0.01 | **0.005** |
| `group_min_records` | 1000 | **300** |

Nothing else changed: `land_share` 0.95, `satellite_gap_km` 500, `satellite_min_datasets` 2, `satellite_min_years` 3, `under_max_records` 5, `thin_max_records` 5, `thin_km` 1000, `core_share` 0.9 and `heuristic_min_records` 20 stay, as do `wcvp_mode` `native_only`, `griis_without_status` `skip` and `side_conflict` `environment`.
The design's four are kept where they are wanted: `realms.calibrate.DESIGN` (setting 64 of the grid) and `tests/clean_helpers.DESIGN_5`, and every test that exercises what the design's values do (the layer 5 rules in `test_outliers`, the toy world of `test_calibrate`) passes it explicitly.
A cleaning at resolution 6 with the new defaults would take more than the logged res 6 runs, which were at the design's thresholds (record share 1%, cell share 5%, rate factor 1%, group minimum 1,000); none was made here.

**The calibration does not move with the defaults.** Checked, and tested where it can be:

- The grid sets the four thresholds itself (`replace(base, **GRID)`), so `settings_grid()` is the same 144 settings for any base: the test gives it bases with the design's, the new defaults and arbitrary thresholds, and compares it setting by setting with the thresholds of the res 7 dry run's `settings.csv` (`tests/fixtures/res7_grid144_settings.csv`).
  Setting 132 is the default and is also the grid's number 132; the design's is still number 64.
- Everything else that took the defaults by the way now names its thresholds: `dry_run` took the design's setting by `grid.index(Settings())` (which would now have found 132 and mislabelled the design as the chosen setting), ran the known cases on `Settings(wcvp_mode=...)` (at the design's thresholds, as its report says), and `baseline` and `griis_cross_check` had `Settings()` as their default argument.
  They use `DESIGN` now. A test reads the syntax trees of `src/realms/calibrate.py` and `scripts/calibrate.py`: every `Settings(...)` there names all four thresholds, every `clean(...)` passes its settings, and no default argument is a bare `Settings()`.
- `--score-half-b` is untouched (it loads no layer 3 and scores a setting of the grid), and **`results/calibrate/` and `data/split/held_out-315557.started` were not touched, nor was the dry run run again.** The dry run's known cases and GRIIS cross-check (both with layer 3 on) were made before the frozen table; a rerun would apply the codes below to the sub-national lists
  (the lists named under "Before and after", and only those), and the 144-setting table, which runs with layer 3 off, would not change.

**The frozen country codes of the sub-national lists.** A sub-national GRIIS list (the US states, Svalbard, the French overseas places, the Dutch Caribbean) is for the rows in the cells that overlap its subdivision, and the cells are bigger than the places, so a border cell holds a neighbour country's rows too. The question is which rows of those cells the list is for, by their country code,
and it is a decision, not a measurement: it should not depend on which rows a run happens to hold, or the later run on all species, or a subset such as the calibration's plants, could apply a list to other codes than half A does, a change of method hidden in the data. So it is made once and frozen. It was made three times in this task, each time because the one before applied some list to the wrong rows:

1. **First form (replaced).** Per place, whether GBIF uses the list's own country code in its cells: **by code** if that code occurs at all among the rows in the cells, **by cells alone** if not. "At all" is too weak. Svalbard is a list of `NO`, and GBIF codes its records `SJ`: in the 28 cells north of 76 degrees N and between 5 and 35 degrees E of half A at resolution 7 there are
   SJ 15,170 rows (337,765 records), no code 167, RU 40, GL 22 and NO 13, and the 13 NO rows put the list on "by code", so that it applied to the NO rows of its cells and to no SJ row (in the list's own 13 cells, which are what the cleaning uses: SJ 13,440 rows, NO 124, RU 15; the cleaning of half A took 0 rows by Svalbard's list).
2. **The dominant code (replaced; ruled in the review of Task 1, 2026-11-02).** Each place gets an effective country code, **the code with most rows** among half A's rows in the list's cells (rows of the counts, one for each species, country code and part; a row with no code is no candidate), a tie going first to the list's own country code and then to the lower code, and the list applies to the rows in its cells that carry it. That put Svalbard on `SJ`
   and the French overseas places on `RE`, `GP`, `MQ` and `PM`, but it also mis-assigned the small islands that share a cell with a bigger neighbour: **BQ/Bonaire went to `CW`** (Curacao: 1,946 of 3,854 rows, against BQ's 1,762) **and BQ/Saba and BQ/StEustatius to `MF`** (Saint Martin: 3,165 of 9,793, against BQ's 1,719), so that these three lists applied to Curacao's and Saint Martin's rows and not to BQ's,
   the opposite of what they are for (BQ rows 135 before, 0 then; 64 CW and 102 MF rows taken instead).
3. **The own code where GBIF uses it materially, else the dominant code (the rule; ruled in the second review round, 2026-11-02).** A list's effective code is **its own country code where the rows with that code are at least `OWN_CODE_SHARE` = 10% of all the rows in the list's cells** (the rows with no code included in the denominator; a starting value, see the gap below), **otherwise the code with most rows, with ties as in 2**
   (the list's own country code first, then the lower code). A place with no row with a code in its cells has its own country's code. The list applies to the rows in its cells that carry its effective code, and a row with no country code is judged by no sub-national list.

- **`scripts/griis_codes.py --counts DIR`** derives it once (`realms.clean.derive_list_codes`) from **half A's counts** at the ingest's resolution and writes `reference/griis_uses_code-res<r>.csv` (`country`, `subdivision`, `code`, sorted), tracked in git, for res 6 and res 7. The file holds the decision only: the shares are of the counts and are in the script's log and `clean_counts.py`'s report. It comes from half A because that is the half the exploration works on and the all-species run must not decide it:
  half B stays sealed until M3, and the decision has to exist before M1 starts. It prints every place with its code, the share of the rows in its cells that carry it, the list's own code's share and the next code, warns of each place whose code has under 80% (`LOW_SHARE`), leaves a table that is there and the same alone (compared sorted), and refuses one that differs or cannot be read unless `--overwrite` says so:
  a different table changes what every later run removes. Both files were replaced by each of the two fixes, which is what `--overwrite` is for.
- **`clean(..., list_codes=table)`** (it was `uses_code`) takes each place's code from the table instead of from `counts`, in both passes of layer 3 (the second judges the rows layer 4 moved). It raises if the table has a place twice, a missing or empty value, wrong columns or types, or no code for a list that is applied. `list_codes=None` derives the codes from `counts` in the same way, for tests.
  A country's list (`griis`) is left out of the rows a sub-national list applies to, that is of the pairs (cell, effective code) of the applied lists, and no longer of every row of the list's own country in its cells: a `NO` row in Svalbard's cells is Norway's and goes to Norway's list, while the BQ rows in the Dutch Caribbean cells are still the BQ lists'.
- **Who passes it:** `scripts/clean_counts.py` and the dry run of `scripts/calibrate.py` (its known cases and its GRIIS cross-check) load the table of their resolution and **refuse to start, before they load or write anything, if it is missing or unreadable**, saying to run `scripts/griis_codes.py`. `explore.py` (Task 4) is to do the same. `FROZEN_DIR`, the repository's `reference/`, is resolved from the scripts' own location, as `SEAL_MARKER` is,
  and not from `src/`. `clean_counts.py`'s `report.json` has every place with its code, rows, rows with the code, share, the own code's rows and share, and the next code (`list_codes.places_detail`), and the places under 80% (`low_share`).
- `natives.load_admin1_cells(grid, reference_dir)` is the admin-1 cells of `load_references` on their own (16 GB less), which `load_references` now calls; the script takes 15 s and 8 GB.

**Where the places sit around the 10%** (the own code's share of the rows in the list's cells, computed for all 24 places before the rule was written). **At resolution 7, the analysis resolution, the places split cleanly:** the sovereign-code places have next to none of their own code in their cells (FR/Guadeloupe, Martinique, Reunion and Saint-Pierre 0, NO/Svalbard 0.91%), and every other place has 17.55% or more
(BQ/Saba and BQ/StEustatius 17.55%, BQ/Bonaire 45.72%, then US/Alaska 87.23% and the rest above 90%). **The gap is 0.91% to 17.55%: the highest own-code share below 10% is Svalbard's 0.91%, and the lowest at or above it is Saba's and St Eustatius's 17.55%.** 10% is in the gap (its middle would be 9.2%), so it is kept as ruled.
**Two places sit between 5% and 20%: BQ/Saba and BQ/StEustatius at 17.55%**, well above 10% but not far above, and they are the places the rule exists for. The margin below the threshold is wide (Svalbard 0.91%, an eleventh of it) and the margin above it is 1.76 times.
**At resolution 6 the split is not clean**, because the cells are bigger and take in more of the neighbours: FR 0 (but Guadeloupe and Martinique, both `GP` at res 6), Svalbard 1.96%, then **BQ/Saba 9.93% and BQ/StEustatius 10.99%, either side of the 10%**, then KM/Mwali, Ndzuwani and Ngazidja 29.80%, BQ/Bonaire 31.65% and the rest above 80%. 10% is inside the gap between 9.93% and 10.99% but only just,
so at res 6 Saba goes to its dominant code `VI` and St Eustatius keeps `BQ`, and the res 6 table is not the res 7 one (it also puts Martinique on `GP`, not `MQ`). Resolution 7 is the analysis resolution and res 6 a rehearsal, and the threshold is not tuned to the res 6 knife edge: that Saba at res 6 is a coin toss is a reason to look at the places again if the analysis resolution ever changes.

**How the frozen tables were verified.** (1) The tests re-derive them: `test_the_frozen_res7_table_equals_the_derivation_from_half_a` and its res 6 sibling run the script's own `derive()` on the real half A (skipped where `data/` is missing) and compare it with the committed file. (2) An independent check with a throw-away script, not committed and sharing no code with `realms.clean`: plain polars on the parquet parts of half A and the admin-1 cells, counting rows by cell and code,
and the rule (the own code with a tenth of all the rows, else most rows, then the list's own code, then the lower code) applied in python to every place: **0 of the 24 places differ at res 7, and 0 of 24 at res 6**; and it reproduces the review's Svalbard counts for the 28-cell box exactly (SJ 15,170, no code 167, RU 40, GL 22, NO 13). (3) For the first form of the table, since replaced, the inline rule of `_layer_3` from before the change was recomputed on half A and agreed at both resolutions.
The own-code shares above are from a separate read-only script on the same counts, and the log of `scripts/griis_codes.py` and the report of `clean_counts.py` give the same figures.

**The codes at resolution 7 (the table that cleaning applies), with the own code's share and the share of the rows in each list's cells that carry the applied code** (all the rows in the cells are the denominator, those with no code included). **Five places have an applied code with under 80%:** the three Dutch Caribbean places, whose own code is a minority of a cell shared with a bigger neighbour by design of the rule, and Martinique (MQ 64.9%, Dominica's 2,982 rows next)
and Saint-Pierre (PM 67.2%, Canada's 1,172), where the code is the territory's own and the neighbours are what is spared:

| List (country / subdivision) | Own code's share of the rows in its cells | Applies to | Rows in its cells | Rows with that code | Its share | Next code | |
|---|---|---|---|---|---|---|---|
| BQ/Bonaire | 45.72% | BQ | 3,854 | 1,762 | 45.7% | CW (1,946) | under 80% |
| BQ/Saba | 17.55% | BQ | 9,793 | 1,719 | 17.6% | MF (3,165) | under 80% |
| BQ/StEustatius | 17.55% | BQ | 9,793 | 1,719 | 17.6% | MF (3,165) | under 80% |
| FM/Chuuk | 100.00% | FM | 967 | 967 | 100.0% |  |  |
| FM/Kosrae | 100.00% | FM | 884 | 884 | 100.0% |  |  |
| FM/Pohnpei | 100.00% | FM | 2,401 | 2,401 | 100.0% |  |  |
| FM/Yap | 100.00% | FM | 1,078 | 1,078 | 100.0% |  |  |
| FR/Guadeloupe | 0.00% | GP (not its own) | 10,332 | 10,274 | 99.4% | DM (44) |  |
| FR/Martinique | 0.00% | MQ (not its own) | 17,021 | 11,043 | 64.9% | DM (2,982) | under 80% |
| FR/Reunion | 0.00% | RE (not its own) | 12,118 | 12,112 | 100.0% | MU (6) |  |
| FR/Saint-Pierre | 0.00% | PM (not its own) | 3,570 | 2,398 | 67.2% | CA (1,172) | under 80% |
| KM/Mwali | 96.25% | KM | 1,466 | 1,411 | 96.2% | YT (55) |  |
| KM/Ndzuwani | 96.25% | KM | 1,466 | 1,411 | 96.2% | YT (55) |  |
| KM/Ngazidja | 97.94% | KM | 2,967 | 2,906 | 97.9% | YT (55) |  |
| NO/Svalbard | 0.91% | SJ (not its own) | 13,579 | 13,440 | 99.0% | NO (124) |  |
| PT/Azores | 100.00% | PT | 18,754 | 18,754 | 100.0% |  |  |
| PT/Madeira | 98.49% | PT | 6,740 | 6,638 | 98.5% | ES (102) |  |
| SH/Ascension | 100.00% | SH | 506 | 506 | 100.0% |  |  |
| TO/Haapai | 100.00% | TO | 1,246 | 1,246 | 100.0% |  |  |
| TO/Tongatapu | 100.00% | TO | 2,555 | 2,555 | 100.0% |  |  |
| TO/Vavau | 100.00% | TO | 1,526 | 1,526 | 100.0% |  |  |
| US/Alaska | 87.23% | US | 136,326 | 118,915 | 87.2% | CA (17,290) |  |
| US/Hawaii | 98.05% | US | 35,883 | 35,182 | 98.0% | UM (701) |  |
| US/United | 90.41% | US | 2,073,489 | 1,874,606 | 90.4% | CA (165,357) |  |

**The first cleaning at resolution 7** (`uv run python scripts/clean_counts.py --counts data/gbif/0020744-261023074918325`, half A, layers 1 to 5 at setting 132, `skip`, the three WCVP modes, the frozen res 7 table), run three times: with the first form of the table, with the dominant codes, and with the own code where it has a tenth of the rows, which is the one below.
**131 s, peak 17.8 GB** (30, 28 and 27 s for the modes' `clean`). Outputs in `results/clean/res7/<mode>/` (git-ignored) and `res7/report.json`. Half A at resolution 7 is 14,645,575 rows of 647,350 species, **1,909,274,465 clean records** (`n - n_human_text`) in 13,506,700 species-and-part presences of 637,903 species.
Records are conserved in every mode: what is left plus layers 2 to 5's removals is the 1,909,274,465.

`native_only`, with the same rule in the last two columns for res 6 (`results/clean/res6/report.json`, **which is not like for like**: it was made on 2026-11-01 before the GRIIS manifest had every list with a table, 55 of 303 lists then and **305 now**, and at the design's layer 5 thresholds. Its layer 3 differs by the manifest and its layer 5 by the thresholds as well as by the resolution, so the two are for orientation only):

| Layer | Rule | Rows | Presences | Records | Share of clean records | Res 6: rows | Res 6: records |
|---|---|---|---|---|---|---|---|
| 0 | `zero_area_part` (moved) | 339,515 | 326,012 | 1,648,535 | 0.09% | 216,241 | 1,255,754 |
| 1 | `failed_flags` (`n_raw`) | 229,154 | 218,561 | 1,476,800 |  | 166,251 | 1,188,700 |
| 1 | `human_text` (`n_raw`) | 32,584 | 31,961 | 62,260 |  | 26,852 | 52,328 |
| 2 | `extinct` | 63,424 | 62,852 | 622,411 | 0.03% | 58,303 | 622,394 |
| 2 | `no_wild_range` | 12,497 | 11,662 | 628,291 | 0.03% | 8,036 | 628,291 |
| 3 | `griis` | 520,496 | 459,801 | 73,046,493 | 3.83% | 160,442 | 50,648,154 |
| 3 | `griis_subnational` | 152,805 | 150,847 | 30,775,910 | 1.61% | 85,708 | 30,668,329 |
| 3 | `wrims` | 5,789 | 4,329 | 318,067 | 0.02% | 5,604 | 392,634 |
| 3 | `wcvp_outside_native` | 501,187 | 474,423 | 7,457,767 | 0.39% | 446,813 | 10,473,999 |
| 4 | `land_sea_fold` (moved) | 947,986 | 916,383 | 9,690,029 | 0.51% | 803,805 | 10,064,176 |
| 4 | `land_sea_neighbour` (moved) | 7,811 | 7,773 | 96,559 | 0.01% | 5,664 | 21,059 |
| 4 | `land_sea_dropped` | 20,411 | 20,306 | 55,203 | 0.00% | 13,108 | 38,302 |
| 5 | `detached_satellite` | 392,355 | 391,233 | 695,999 | 0.04% | 233,516 | 457,426 |
| 5 | `under_reported` | 155,322 | 152,558 | 252,536 | 0.01% | 196,485 | 319,154 |
| 5 | `thin_long_range` | 55,724 | 55,683 | 66,731 | 0.00% | 55,164 | 67,039 |

The rows of layers 0 to 2, `griis`, `griis_subnational` and `wrims` are the same in all three modes, and what differs by mode is the WCVP rule, which layer 4 follows, and layer 5, which judges what layers 1 to 4 leave:

| Rule | `native_only` | `introduced` | `off` |
|---|---|---|---|
| `wcvp_outside_native` / `wcvp_introduced` | 501,187 rows, 474,423 presences, 7,457,767 records | 162,935 rows, 145,838 presences, 5,081,307 records | none |
| `land_sea_fold` | 947,986 rows, 9,690,029 records | 954,352 rows, 9,713,735 records | 962,170 rows, 9,758,989 records |
| `land_sea_neighbour` | 7,811 rows, 96,559 records | 8,032 rows, 97,128 records | 8,256 rows, 97,947 records |
| `land_sea_dropped` | 20,411 rows, 55,203 records | 20,370 rows, 55,033 records | 20,352 rows, 54,828 records |
| `detached_satellite` | 392,355 rows, 391,233 presences, 695,999 records | 461,718 rows, 460,470 presences, 797,150 records | 467,623 rows, 466,324 presences, 808,927 records |
| `under_reported` | 155,322 rows, 152,558 presences, 252,536 records | 170,052 rows, 167,096 presences, 274,484 records | 178,850 rows, 175,654 presences, 289,640 records |
| `thin_long_range` | 55,724 rows, 55,683 presences, 66,731 records | 57,896 rows, 57,853 presences, 69,342 records | 57,878 rows, 57,835 presences, 69,322 records |
| **Presences left** | **10,711,332** of 588,594 species | **10,948,641** of 589,412 species | **11,073,392** of 589,460 species |
| **Records left** | **1,795,355,057** (94.0%) | **1,797,605,977** (94.2%) | **1,802,660,576** (94.4%) |

**Before and after the codes**, against the two earlier runs of this task (the first form of the table, then the dominant code, then the rule above; the counts of layers 0 to 2, `wrims` and `thin_long_range` do not change):

| | First form | Dominant code | Own code at a tenth (now) |
|---|---|---|---|
| `griis_subnational`, every mode: rows, presences, records | 152,792, 150,833, 30,776,040 | 152,836, 150,802, 30,775,301 | **152,805, 150,847, 30,775,910** |
| `griis`, every mode: rows, presences, records | 520,495, 459,801, 73,046,492 | 520,263, 459,717, 73,032,445 | **520,496, 459,801, 73,046,493** |
| `wcvp_outside_native` (`native_only`) / `wcvp_introduced` (`introduced`), rows | 501,209 / 162,935 | 501,269 / 163,020 | **501,187 / 162,935** |
| Presences left: `native_only`, `introduced`, `off` | 10,711,315, 10,948,639, 11,073,388 | 10,711,411, 10,948,691, 11,073,487 | **10,711,332, 10,948,641, 11,073,392** |
| Records left: `native_only`, `introduced`, `off` | 1,795,354,599, 1,797,605,821, 1,802,660,415 | 1,795,368,518, 1,797,619,453, 1,802,675,180 | **1,795,355,057, 1,797,605,977, 1,802,660,576** |

Place by place (`griis_subnational` rows and records, `native_only`; the three modes agree):

| Place | First form | Dominant code | Own code at a tenth (now) |
|---|---|---|---|
| **BQ/Bonaire** | 61 BQ rows, 6,044 records | 64 CW rows, 5,742 records | **61 BQ rows, 6,044 records** |
| **BQ/Saba** and **BQ/StEustatius** (one set of cells) | 74 BQ rows, 3,410 records | 102 MF rows, 3,103 records | **74 BQ rows, 3,410 records** |
| NO/Svalbard | 0 rows | 53 SJ rows, 485 records | 53 SJ rows, 485 records |
| FR/Martinique | 177 MQ rows, 7,460 records, and 37 rows of neighbours (13 LC, 13 DM, 11 VC: 556 records) | 177 MQ rows | 177 MQ rows |
| FR/Saint-Pierre | 61 PM rows, 1,501 records, and 3 CA rows (59 records) | 61 PM rows | 61 PM rows |
| FR/Reunion, FR/Guadeloupe, and every other place | as in the first form | as in the first form | as in the first form |
| the country lists of CW and MF (`griis`) | 213 and 127 rows (15,637 and 2,911 records) | 107 and 0 rows (4,500 and 0), the rest now the BQ lists' | 213 and 127 rows, as in the first form |

**The BQ places are restored** (135 BQ rows, 9,454 records, taken by the three lists, and none of CW's or MF's), and Svalbard and the French overseas places keep what the dominant code gave them. Against the first form the whole change of `griis_subnational` is Svalbard's 53 SJ rows (485 records) gained and Martinique's and Saint-Pierre's 40 neighbour rows (615 records) spared:
+13 rows and -130 records, and the one row more of `griis` (a NO row in Svalbard's cells that Norway's list now takes). The two GRIIS rules together take 103.8M records (5.44% of the clean records) in 673,301 rows. The few dozen rows that other rules gain or lose against the first form (`wcvp_outside_native` -22, `detached_satellite` -8 rows) are the cascade of rows that a list now spares or takes.

- **Per layer, `native_only` lost** (the presences a rule left with nothing, at the moment it took them) 326,012 (0), 250,499 (1), 74,514 (2), 1,089,300 (3), 944,462 (4) and 599,474 (5) presences, of 13,506,700, and **left 10,711,332 (79.3%)**. **These do not add up to the 2,795,368 presences that are gone:** a fold's presence reappears at the part it went to, so that a presence that layer 4 (or the zero-area fold)
  moved is counted among the lost at the part it left and is back at the part it went to, 488,893 presences in all. Records removed by layer: 1,250,702 (2), 111,598,237 (3), 55,203 (4, dropped) and 1,015,266 (5); the folds of layers 0 and 4 move 1,648,535 and 9,786,588 and lose none.
- **Layer 3 is again mostly records, and mostly GRIIS**: the two GRIIS rules take 103.8M records (5.44% of the clean records), `wcvp_outside_native` 7.5M (0.39%), WRiMS 0.3M. The pairs that could not be placed on the grid (62 sub-national lists with no admin-1 cells, 6,569 pairs) name species with records in half A in 3,558 rows (`unplaced_griis`), and 1,137 for WRiMS (`unplaced_wrims`).
- **Layer 4**, `native_only`: of 611,171 species 455,776 are land, 80,869 sea and 27,549 exempt (freshwater, brackish or Aves) and the rest (46,977) have no side; 6,232 are in the conflict (1,515 terrestrial with their records at sea, 4,717 marine with theirs on land). The zero-area fold takes 2,426 parts (2,334 sea, 92 land: 339,515 rows, 1,648,535 records), against 923 at resolution 6.
- **Layer 5 at setting 132**, per mode: `native_only` 603,401 rows, **599,474 presences**, 1,015,266 records (0.05% of the clean records) of 113,885 species, 22,334 of which it leaves with no presence at all (`thin_long_range`); `introduced` 689,666 rows, 685,419 presences, 1,140,976 records, 128,272 species (23,215 lost); `off` 704,351 rows, 699,813 presences, 1,167,889 records, 128,434 species (23,207 lost).
  Beside res 6, which was at the design's thresholds, setting 132 takes many more detached satellites (392,355 rows against 233,516) and fewer under-reported rows (155,322 against 196,485). That is the direction of the thresholds (a larger record share and cell share make more patches satellites, a lower rate factor flags fewer cells, a lower group minimum lets more cells be judged), but
  the resolution differs too, so it is not a measurement of them. The **suspect cells** (cells where layer 5 took far more of the species than next door, listed in `suspect_cells.csv`, nothing removed) are 2,074 in `native_only` and 2,094 in the others, against 500 at res 6.
- **The Caspian** at resolution 7: 27 cells with 1,000 km2 of sea or more; 1,625 of the 216,439 clean records in them are on the sea parts, and 390 after all layers (`native_only`); 68 records are folded by the zero-area rule.
- **Provisionally accepted twins** in half A: 2,834 species, 84,499 rows, 2,728,896 `n`; none has several accepted twins. Nothing is changed for them.

**Not done here:** no cleaning at resolution 6 (its logged results are of the design's thresholds and of the manifest before the 305 lists), no calibration, and no look at half B. A rerun of the sealed dry run would apply the codes above in its known cases and GRIIS cross-check (layer 3 on), and the table of its 144 settings (layer 3 off) would not change. `scripts/explore.py` (Task 4) is to pass the frozen table too.
The suite has 1,079 tests, all passing (1,051 before Task 1: the new ones are the frozen table and the effective codes, the script, the grid's independence of the defaults and `load_admin1_cells`).

### The exploration pipeline, and its first two stages on half A at resolution 7 (2026-11-02)

`scripts/explore.py` (library `src/realms/explore.py`, tests `tests/test_explore.py`) is the cached chain that every M1 run uses: clean, network units, `CellTaxa` and the network, Infomap, `candidate_maps`, screen, select. It adds no method: every threshold is that of the stage it calls (setting 132, `UnitSettings` as they stand, the design's 60 to 5 realms with the 1M km2 floor, 10 trials, 20 bootstraps of 3 trials and 99 null maps). Rulings made here, because the brief left them open:

- **R2: `explore/configs.toml` is tracked and holds one configuration, `default`, which is `InfomapOptions()`**, with `trials_screen = 3` and `n_null_screen = 19`. A `[configs.<name>]` table holds the fields of `InfomapOptions` that differ from Infomap's defaults, and a key that is not a field is an error. Task 9 adds the sweep configurations.
- **What a stage is cached on.** A stage has a stamp with the sha256 of its inputs and is skipped when its files are there and the stamp agrees. Cleaning hashes `Settings` (the WCVP mode, `rules_off` and every threshold), the layers, the content of the frozen GRIIS code table, the counts' `_COMPLETE` marker and the references' `manifest.json`. The units hash the cleaning's hash and `UnitSettings`; screen and select hash the units' hash, the configuration's `InfomapOptions`, the parameters (seed, `k_max`, `k_min`, floor, trials, bootstraps, nulls) and the stage. A change anywhere reruns what follows it, and only that. `--stage X` runs every stage up to X; `--force` reruns X only. The manifest is in the hash although the brief did not list it: a reference export fetched since changes the cleaning, and a hash that missed it would reuse a stale one.
- **The units are the network's cells.** `candidate_maps`, `bootstrap_stability` and `excess_table` get the units' areas and adjacency, and the matrix is built from the units' presence. A part left out of the network (`Units.excluded`) is in none of them, and gets its realm from `units.impute` in the outputs only, flagged `imputed` in `parts.parquet`. A unit that Infomap gave no module (it has no link) is treated the same way for the outputs, as `unlinked`; whether one occurs at res 7 is for the screen to show, and the planted tests cover the case.
- **Outputs.** `candidates.parquet` has `unit, k, realm, sub_realm, sub_realm_of` for every k: `sub_realm` is the unit's label before the size floor where the floor merged it into another (null otherwise), and `sub_realm_of` is the realm it sits in at that k. `parts.parquet` has `part, k, realm, unit, imputed, selected` for every k and every part, with the selected k marked (`selected`), the unit null for an imputed part. Tasks 5 and 12 draw maps from `parts.parquet`. Screen and select write the same map files (a screen at `trials_screen` trials, a select at 10), so they are those of the stage that ran last, `structure.json` says which (`stage`, `trials`), and a screen that runs again removes the select result it overwrote (`selected.json`, `selection.csv` and the select stamp). (Until the final review of M1: select now writes `<config>/select/`, see "Changes after the final review of M1".)
- **A degenerate clustering is a result.** `DegenerateClustering`, and a network with no link at all, are written to `structure.json` (`degenerate`, `degenerate_reason`) and to `selected.json`, with empty maps; nothing is raised.
- **Memory and workers.** Cleaning and the units each run in a child process (`log.json` of the run has the child's peak and the stage's seconds and `wall_seconds`, with the loading), so that the parent holds none of the 14 to 18 GB of the references. With `--workers N` the first two tasks of the bootstraps and of the null maps run in two workers, to measure a worker's peak, and the rest in as many as fit in 40 GB next to the parent at that peak times 1.25. Each worker's seconds and peak are in the configuration's `log.json`. The results do not depend on the number of workers (the tasks are pure functions, as `parallel.map_tasks` has it). With one worker, nothing is recorded per worker, and the stage's peak is the process's so far.
- **Half A only.** The only loader of counts is `load_country_counts(..., "A")`. Nothing here passes `unsealed` or touches half B, `scripts/calibrate.py`, `results/calibrate/` or `data/split/held_out-315557.started`. A test asserts the half that is loaded, and that the module has no word for another.

**The real-data stages, `--wcvp-mode native_only` at resolution 7** (the ingest `0020744-261023074918325`, half A, setting 132, `griis_without_status` skip; `results/explore/res7/native_only/`). Screen and select were not run; those are Task 9.

| Stage | Stage time | With loading | Peak memory (the child's) |
|---|---|---|---|
| clean | 31 s | 61 s | 16.0 GB, then 18.0 GB on the forced rerun |
| units | 4 s | 15 s | 6.2 GB |

The rerun with `--stage units` skipped both stages, and `--stage clean --force` reran cleaning and left the units current (same hash). That was a gap, not the intent: units built on the old cleaning stayed "current". Since the final review of M1, a stage that writes its outputs again deletes the stamps of the stages that depend on it, so the units would run again (see "Changes after the final review of M1").

- **Cleaning** left 10,711,332 presences of 588,594 species and 1,795,355,057 records on 22,147 parts (7,683 land), the same as `scripts/clean_counts.py` logged for `native_only` above.
- **Units.** 22,147 parts with records; 10,782 are not sparse and stay as their own unit, 11,365 parts are sparse and 5,005 of them are in 810 coarse units that entered the network. **11,292 units are in the network: 10,482 fine (6,365 land and 4,117 sea) and 810 coarse (114 land, 696 sea)**, holding 15,487 parts and 1,794,416,861 clean records (99.95%). 13,060 units were judged (10,782 fine and 2,278 coarse): of the coarse ones 1,416 stay sparse, 52 fail the composition gate and 810 enter the network.
- **Parts left out: 6,660, 938,196 clean records (0.052%), by reason**:

  | Reason | Parts (land / sea) | Units tried | Records | Species (distinct) |
  |---|---|---|---|---|
  | `sparse` (a coarse unit still under 0.5 coverage or 50 species) | 6,042 (692 / 5,350) | 1,416 coarse | 635,268 | 14,203 |
  | `composition` (fewer than 2 major groups sampled) | 618 (125 / 493) | 300 fine, 52 coarse | 302,928 | 6,478 |

  17,585 species have a record in some part that is left out. 290 of the 352 composition units have one major group sampled and 62 none.
- **By region** (WGSRPD level 1 of the cell centre; parts / records): sparse: sea 5,350 / 587,929, Antarctic 171 / 12,618, Asia-Temperate 125 / 6,784, Africa 111 / 5,943, Pacific 110 / 12,926, Asia-Tropical 70 / 3,457, Northern America 63 / 2,081, Southern America 27 / 2,286, Australasia 12 / 941, Europe 3 / 303. Composition: sea 493 / 262,432, Asia-Temperate 98 / 31,649, Africa 20 / 7,795, Northern America 3 / 378, Asia-Tropical 2 / 228, Southern America 2 / 446. `units/removed_by_region.csv` has it. The sea is 5,843 of the 6,660 parts (88%) and 91% of the records left out (850,361 of 938,196); on land the sparse parts are mostly in the Antarctic, Asia-Temperate, Africa and the Pacific, and the composition gate takes 98 Asia-Temperate parts of its 125.
- **The major groups of the gate** (`units.major_groups`, from COL's `taxGroup`, and for chordates the `class`), computed once in the cleaning stage and cached as `cleaned/major.parquet`: 569,770 of the 588,594 species with a presence are counted, 18,824 are not (algae, protists, bacteria and the other groups of `NOT_COUNTED`, and chordates with no class or one that is in no list; no chordate class warning was raised). Species per group: insects 210,675, plants 162,112, other arthropods 51,544, fungi 42,279, other invertebrates 36,218, molluscs 31,542, land vertebrates 20,223, fish 15,177. The mapping of `taxGroup`:

  | Major group | COL `taxGroup` |
  |---|---|
  | plants | angiosperms, gymnosperms, pteridophytes, bryophytes, plants |
  | insects | coleoptera, lepidoptera, diptera, hymenoptera, hemiptera, orthoptera, trichoptera, otherinsects, insects |
  | other arthropods | arachnids, crustacean, otherarthropods, arthropods |
  | molluscs | gastropods, bivalves, othermolluscs, molluscs |
  | fungi | ascomycetes, basidiomycetes, otherfungi, fungi |
  | other invertebrates | otheranimals, animals |

  The chordates (`taxGroup` chordates, or phylum Chordata whatever the group) are split by class: **land vertebrates** Mammalia, Aves, Reptilia, Amphibia; **fish** Teleostei, Elasmobranchii, Chondrostei, Holocephali, Myxini, Holostei, Petromyzonti, Coelacanthi, Cladistii, Dipneusti, Cephalaspidomorphi; **other invertebrates** (the invertebrate chordates) Ascidiacea, Thaliacea, Appendicularia, Leptocardii. A group is sampled in a unit where it has at least 5 species, and a unit needs 2 sampled groups (`UnitSettings`).
- **The res 5 edge tie-break.** A sparse part joins the res 5 cell that holds its res 7 cell's centre. 12.9% of the res 7 centres (2,826 of 21,872) lie on an edge or vertex of the res 5 grid, since res 6 centres sit on the vertices of res 5. They go to the **lowest cell id as text**, because the res 5 cells are sorted by id before `locate`, which takes the lowest index: so `"10"` is before `"9"` and a negative id before a positive one. That settles exact ties only; for a centre only almost on an edge the noise in the file's vertices (about 1e-8 degrees) decides, so the assignment there is stable but not meaningful, and a coarse unit's border is not to be read closely. This is unchanged since Task 2; it is logged here because the first real units were made now.

### Maps of the parts: how a cell is drawn (2026-11-02)

`hexmaps.part_geometries` and `plot_parts` draw the land and sea parts on an Equal Earth map (EPSG:8857), for the candidate realm maps.

- **Same cell, same land.** A cell is `grid._planar` (so a polar cell is closed to the pole and draws as a cap), the land is the overlap that `grid._land_overlaps` finds, which `_land_km2` now shares, and the sea part is the cell minus that land. The two tile the cell. A part with under `ZERO_KM2` is empty. For the map the drawing is then cut at the antimeridian, each piece moved back into -180..180, so no piece spans the globe.
- **Edges are straight in longitude and latitude**, as the grid files give vertices, cut into edges of at most 0.25 degrees before projecting. That is not the geodesic edge `grid.py` measures areas with, so on the map a polar part can be up to about 15% larger than its `area_km2` (38 of 3,246 parts at resolution 5 differ by more than 3% from the geodesic area, nearly all within 10 degrees of a pole or thin slivers). It only changes how the part looks, never the data; neighbours share the same edge, so nothing opens between them.
- **Two tools of shapely/GEOS 3.14 are avoided**: `segmentize` returns an empty polygon for some parts (own densifying instead), and `clip_by_rect` raises on some valid coastal parts (the cut at the antimeridian uses an intersection with a box). Where a projected part comes out self-touching, because a coast that follows a cell's edge bends across it, it is repaired with `make_valid`.
- **Colours** come from `maps._separating_colours`, with realm neighbours from part adjacency (`grid.part_adjacency`, passed as `adjacency`; by default the parts that share boundary on the map). `plot_cells` is untouched.

### Realm boundaries, and the TDWG-edge gate by rotated regions (2026-11-02)

`src/realms/boundaries.py` (tests `tests/test_boundaries.py`) measures realm boundaries on the grid and tests whether a WCVP mode draws its land boundaries along TDWG level 3 edges more than chance. Task 10 runs it on the three WCVP modes.

- **Edges.** `edge_lengths` keys every pair of neighbouring cells that `grid._neighbours` finds (so the adjacency is that of `load_grid`), and gives the geodesic km (WGS84) of the edge they share. An edge is the geodesic between two vertices of the grid file, as `grid.py` takes edges for areas, and pyproj measures a geodesic exactly, so no densifying is needed and the 0.1% bound holds without it. The shared edge is found as the boundary of one cell inside a band of `SNAP` around the other's, at each turn of longitude; a cell over a pole shares one meridian at one turn and the other at another, and the detour `_planar` draws over the pole measures as the geodesic over it. At resolutions 5 and 7 the edges of every cell add up to its geodesic perimeter within 2.4e-6, and there are 3N - 6 pairs. The band adds about 2e-6 of length (centimetres at the ends). Cached in `results/boundaries/edges-<stem>.parquet` with the gpkg's sha256 beside it.
- **A boundary** is between neighbouring cells whose parts of one kind carry different realms; a cell's land and sea parts never are one. `edge_share` is the share of boundary km with either cell in the near set, NaN for a map with no boundary (which then fails the gate).
- **Near.** `near_cells` is the cells whose polygon intersects two or more regions, touching included, plus their neighbours: the design's "within one cell width". Regions are also tried a turn of longitude to either side, and the pieces of a region cut at the antimeridian count as one region.
- **Rotations.** A draw is a longitude shift uniform in [0, 360), then a tilt about the horizontal axis through the equator at a uniform longitude, by an angle that puts the rotated pole uniformly on the spherical cap within `max_tilt_deg` (the cosine of the angle is uniform), so small tilts are not favoured. Regions are rotated as vertices on the unit sphere after cutting their edges to 0.25 degrees (they are straight in longitude and latitude in the file); a ring that runs out and back along one line, as Antarctica's does down the antimeridian to the pole, is cleaned of that spike; longitudes are unwrapped, and a ring that then winds round the globe is closed along the pole on the side of a point inside it; the result is cut at the antimeridian and made valid. Holes are rotated as rings of their own and taken out. Rotated regions keep their area on the sphere within 0.11% (the cut edges), Antarctica's included.
- **Boundary over land** is the union of the regions' boundaries (shared edges once, snapped at `SNAP`; the file's cut at the antimeridian left out), cut into pieces of at most 0.05 degrees, a piece over land where its midpoint is in Natural Earth land. Rotating the midpoints is then all a draw costs (0.2 s). A draw is kept when `|rotated - unrotated| <= land_tolerance * unrotated`; the seed fixes the draws; after 50 n draws too few kept is an error.
- **Starting values: `max_tilt_deg = 30`, `land_tolerance = 0.5`** (the brief's 0.2 cannot be met, below). `land_tolerance=None` keeps every draw and reads no land.
- **The gate is paired.** For each rotation r, `s_m(r) - s_base(r)` is the null raise of mode m, with `near_cells` computed once per rotation for all modes; the observed raise is `s_m(0) - s_base(0)`; a mode passes when its raise is at most the 95th percentile of its null raises (numpy's linear interpolation). Rotation moves region edges onto the sea, where no land boundary runs, which lowers every map's share alike; taking the difference to the `off` map on the same rotated regions cancels that, and the question becomes only whether the mode is closer to the TDWG edges than `off` is. Only land boundaries enter, as the regions are of land.

**Measured at resolution 7 with the real level 3 regions** (peak RSS 0.46 GB): `edge_lengths` 3.9 s cold (65,610 pairs), 0.01 s from the cache; `near_cells` 0.1 s (4,904 cells unrotated); rotating all 369 regions 0.6 s; the union of the boundaries is 1,256,588 km, of which 938,112 km lie over land unrotated.

**The land tolerance: 0.5, not 0.2.** The regions are land, so their boundaries lie almost wholly over land, coasts included (the WGSRPD coast falls mostly inside Natural Earth's). Rotated, they fall on land only in proportion: over 500 draws the rotated length was 0.31 to 0.83 of the unrotated one (median 0.49, 95th percentile 0.66). The brief's 20% keeps about 1 draw in 500, and `rotations(regions, 5, seed=1)` raised after its 250 draws (55 s); no tilt limit changes that, as the longitude shift alone moves the continents onto the ocean. The starting value is therefore 0.5, which keeps the more land-like half: `rotations(regions, 99, seed=1)` kept 99 of 224 draws (44%), with kept ratios 0.50 to 0.77 (median 0.56), in 50 s and 0.37 GB. Loosening the tolerance does not make the gate easier on WCVP: a rotation that puts the regions on little land leaves few land boundaries near them, which shrinks the null raises and makes the gate stricter, pushing towards `off`; letting in more of the less land-like draws moves the test towards a fair one, not towards passing.
- **A pairing-only reading.** `land_tolerance=None` keeps every draw (99 in 0.001 s) and leaves the land bias to the pairing alone. Task 10 reports the gate both ways, at 0.5 and with None.
- **Open for Markus:** the boundary over land is mostly coastline (about 1.13M of the 1.26M km of the union), not the internal borders that make cells near, so the filter tracks where the coasts land more than where the borders do. Measuring only the shared borders would be the alternative; the measure stays as built until Markus decides.

### Beta-sim over three ranks, and its consensus over shuffled cell orders (2026-11-02)

`betasim.rank_beta_sim` and `betasim.consensus` (tests in `tests/test_betasim.py`) are the cross-check and the fallback to Infomap. `beta_sim`, `hierarchical_labels` and `anosim_r` are unchanged, so the Costello gate is too.

- **`rank_beta_sim(ct, units)`** is the mean of three beta-sim matrices, on the species, the genus and the family columns of the cell-taxon matrix, for `units` in the order given (a subset or reordering of `ct.cells`; an unknown id raises). It is float32, symmetric, with a zero diagonal. A unit with no taxon at a rank (0/0) is at distance 1 from every other unit at that rank, as in `beta_sim`, even from another empty one, so such a unit is not read as similar to another. A matrix with no column at all at one of the three ranks raises, as the average would otherwise count a missing rank as full dissimilarity. It is built in row blocks of about 64 MB, so only the result is held at full size.
- **`consensus(dist, method, ks, shuffles=20, seed=1, ids=None)`.** Beta-sim has a great many tied distances (on Costello's matrix 55.5% of the pairs are at 1), so one Ward or UPGMA tree depends on the order of the cells. Each of `shuffles` seeded random orders of the units is clustered (`method` is `"ward"`, SciPy's ward.D2, or `"average"`, UPGMA) and the one tree is cut at every k at once (`cut_tree`, exactly k groups, back in the original order). For each k, **co-assignment** `co[i, j]` is the share of shuffles that put units i and j in one group, and the **groups** are an average linkage on `1 - co`, cut at k. **`spread`** is the mean V-measure (`sklearn.metrics.v_measure_score`) over all pairs of the shuffles' own groupings at that k: 1.0 when the order never mattered. Labels are numbered 0..k-1 by first appearance in the input order. `shuffles` must be between 2 and 255 (a spread needs two groupings, and the counts are one byte).
- **The canonical order is by unit id.** With `ids` (the unique id of each row of `dist`, which Task 11 passes as the network unit ids) the units are first put in the order of their sorted ids, and the shuffles, every linkage and the tie-breaks of the last linkage all run in that order. Ids never tie, so the groups are then **exactly the same for any order of the input** (`dist` and `ids` permuted together): the same partition, the same `spread`, tested on a matrix with 94% of the pairs at 1, where a single tree changes with the order and a consensus without ids does too. Only the label numbers follow the input order. This is done by indexing, without a second copy of the square. **Without `ids` the input order is the canonical one**, and the groups depend on it through the shuffles drawn for it and the tie-breaks, so Task 11 always passes `ids`. An earlier version of this entry said the order did not matter when the structure is clear; that held only for matrices with few ties (a test on one with 5% of pairs at 1 could not tell a consensus from a single tree) and is withdrawn.
- **What ids do not remove: the shuffle noise.** The seed still moves the groups. Measured on 300-unit presence matrices (6, 12 and 30 planted realms, 95.3%, 98.1% and 98.9% of pairs at 1; six seeds each, cut at the number of realms, V-measure between the groupings of two seeds):

  | | 20 shuffles, mean (min) | 100 shuffles, mean (min) |
  |---|---|---|
  | Ward, 6 realms | 0.870 (0.843) | 0.931 (0.911) |
  | UPGMA, 6 realms | 0.845 (0.779) | 0.913 (0.869) |
  | Ward, 12 realms | 0.825 (0.754) | 0.900 (0.852) |
  | UPGMA, 12 realms | 0.904 (0.842) | 0.979 (0.968) |
  | Ward, 30 realms | 0.859 (0.838) | 0.931 (0.883) |
  | UPGMA, 30 realms | 0.875 (0.818) | 0.915 (0.862) |

  So on a matrix that tied, two seeds agree at V 0.75 to 0.95 with 20 shuffles and 0.85 to 0.99 with 100, and the mean spread of the shuffles themselves was 0.50 to 0.72, so a single tree there is much less stable than the consensus. **The default stays 20.** Task 11 reports the seed sensitivity (the consensus at a few seeds, V-measure between them) beside the groups, and raises `shuffles` if it is low; the cost is linear in `shuffles` for the linkages, and the count matrices do not grow.
- **Cost, measured at 11,292 units** (res 7; synthetic grid presence, 150,000 species, 181,500 taxa columns, 13.8M links, 107 x 107 grid cells cut to the first 11,292; ties are rare in it, which does not change the cost of a linkage): `rank_beta_sim` 18 s, 1.8 GB peak; one linkage 0.8 s for Ward and for average (about 0.07 s to condense, a 510 MB float64 vector); `consensus(..., ks=range(5, 61), shuffles=20)` 180 s for Ward and 172 s for average (186 s for Ward with `ids` given, which adds only an index), **8.5 GB peak RSS** for either, of which 7.1 GB are the 56 uint8 count matrices (127 MB each). Task 11 calling both methods takes about 6 minutes and 9 GB at a time, well inside the limits, with no shortcut needed beyond what is built in: one linkage per shuffle cut at all ks, and one n x n uint8 matrix per k built incrementally. The real species rank has more links than the synthetic one (588k species), which only `rank_beta_sim` feels, so that may take a minute or two.

### Merging the duplicate genus and family ids of the XR (2026-11-02)

`cell_taxa(..., merge_duplicate_names=True)` merges the ids that share a name in one lineage, for the open item of "Taxa by id, accepted only". The default is off and leaves every matrix as it was, so the Costello gate is unchanged (`selection.csv` and `selection.log` byte-identical to `results/costello/before/`). Task 11 reruns the chosen configuration with it on and compares.

- **The rule (Markus's ruling, narrower than the brief's "share a name in one kingdom").** Genus ids merge when they share **kingdom, genus name and family name**; family ids when they share **kingdom, family name and order name**. A row with no id, name, kingdom or parent name never merges, because its lineage can't be established. Merging on kingdom and name alone would fuse true homonyms, a genus name under two families, which are the 672 genus ids "accepted in another family" above. Every id of a group is mapped to **the group's lowest id, compared as text**, so the survivor is the lowest and not the largest; the species columns' `genus_col` and `family_col` point at the survivor. `duplicate_names(taxa)` lists the groups as `rank, kingdom, name, parent, ids`. Tests: two ids with one name merge, one name in two kingdoms does not, one genus name in two families does not, and the default is unchanged.
- **Real data, 315557 (2,529,340 accepted species).**

  | | genus | family |
  |---|---|---|
  | ids that hold species | 240,578 | 16,906 |
  | duplicate groups (kingdom, name, parent) | 4,861 | 811 |
  | ids in those groups | 10,002 | 1,645 |
  | ids merged away | 5,141 | 834 |
  | species moved to another id | 102,554 | 334,851 |
  | of those, species outside their group's largest id (the ones that change a column's size) | 17,265 | 13,217 |
  | same-name groups by kingdom alone | 8,080 | 1,004 |
  | ids outside the largest of those | 8,996 | 1,040 |
  | kingdom-alone groups with ids under two or more parent names | 1,249 | 54 |
  | kingdom-alone groups where no two ids share a parent name (left entirely apart) | 3,221 | 193 |

  The kingdom-alone figures reproduce the 8,080 genus and 1,004 family names and the 8,996 and 1,040 non-largest ids of "Taxa by id" exactly, so the groups are the same population. The parent rule keeps 4,861 of 8,080 genus names and 811 of 1,004 family names as groups, and merges away 5,141 of the 8,996 genus ids and 834 of the 1,040 family ids; the rest are homonyms or have no parent name, and stay apart. The 5,141 are comparable with the 4,879 provisionally accepted in the same family. The "species moved" row counts relabelled species: when the lowest id is not the largest (1,706 genus groups, 218 family groups, for example a family whose small copy sorts first) the big id's species move to a smaller id's name, which changes no column's size. The row below it is the real change, 17,265 genus and 13,217 family species that join a larger taxon, against 28,942 and 14,416 outside the largest by kingdom alone.

### Units that Infomap puts in no module leave the scores (2026-11-02)

A unit with no link (every species it holds is not accepted) is a row of the cell-taxon matrix without a link, and Infomap leaves it in no module (`Partition.unassigned`). It was dropped from the maps but stayed in the matrix, `units.areas` and `units.adjacency` that `candidate_maps`, `bootstrap_stability` and `excess_table` were given, so it counted in the null maps and the totals of the scores for a map that did not contain it. The scored cells now equal the mapped ones. `explore.assigned_only` restricts the matrix, areas and adjacency to the assigned units after Infomap and before `candidate_maps`, in both the screen and the select stage. The row is gone from the matrix, the unit is in no neighbour set, and a taxon column that no remaining unit holds goes with its genus and family columns renumbered (none does in practice: a unit with no link holds no accepted species). `cluster.py`, `candidates.py`, `select.py`, `nullmodel.py`, `network.py`, `matrix.py` and `maps.py` are unchanged, so the Costello gate is too.

- **Where it is written.** `structure.json` has `unassigned` (the count, as before) and `unassigned_units` (the sorted ids), and a `warning` if they are more than 1% of `units`. `partition.parquet` still lists them with a null leaf and top, and `parts.parquet` imputes their parts from the neighbours (`labelled_units` makes them `unlinked` excluded parts), flagged as imputed, with no unit. They are in no row of `candidates.parquet`.
- **The bootstraps run on the restricted matrix.** `bootstrap_stability` runs Infomap again on each resample of the species. A resample drops the units that hold none of the drawn species from the matrix (`CellTaxa.resample`), and every unit that remains holds a species, which links at weight 1, so no unit is unassigned inside a bootstrap and `Partition.unassigned` is empty there. Units absent from a bootstrap partition do not break it: `v_measure` compares the cells that both maps contain. One consequence: the information weight of a genus or family depends on the number of cells, so the network that the bootstraps (and the null maps) see has one weight scale fewer cells than the one that produced the maps by `len(unassigned)`. That is the weight of the scored cells, which is what the design wants, and it only matters if there are unassigned units at all.
- **Tests.** `test_units_without_a_module_leave_the_scores_but_not_the_map` runs the two stages on the planted world with one and with three units that hold only an unaccepted species, spies on `bootstrap_stability` and `excess_table`, and checks the files.

### Infomap's scale on half A at resolution 7 (2026-11-02)

**Chosen: `mt080` (Markov time 0.80, everything else Infomap's defaults) at k=7**, stability 0.806, excess 0.1568, on `native_only` at resolution 7 (release 315557). Seven distinct non-degenerate configurations were bootstrapped at 10 trials, and `mt080` is the only one with a map at or above stability 0.8 (k=7, 9 and 10). Half A only, on the 11,292 network units; nothing from B. Tasks 10 and 11 use `mt080`.

**The 0.8 bar decides this, not excess.** Excess is flat: the maps of the six configurations other than `regularized` run from 0.125 to 0.188 in excess, `regularized` from 0.157 to 0.210, and the three stable `mt080` maps sit at 0.1535 to 0.1580. Many unstable maps beat them (for example `mt075` k=5 at 0.1875 and `regularized` k=11 at 0.2102). Stability is what separates `mt080` from the rest, and it clears the bar by 0.006 at k=7 (0.806), with k=8 at 0.793 just under it. See the caveats.

**Screen** (3 trials, 19 null maps; `results/explore/res7/native_only/<config>/screen.csv`). Ten configurations: Markov time 1.0 (`default`) down to 0.75, and at 1.0 `variable_markov_time`, `regularized`, `entropy_corrected` and `recurse_share=0.5`. Units above the floor = the top of the k range.

| config | top modules | single-unit | largest module, units / area | units above floor | maps (k) | best excess (k) |
|---|---|---|---|---|---|---|
| default | 15 | 3 | 34.6% / 52.0% | 6 | 2 (5-6) | 0.1442 (5) |
| mt095 | 21 | 3 | 39.5% / 56.6% | 8 | 4 (5-8) | 0.1436 (5) |
| mt090 | 109 | 74 | 39.7% / 57.0% | 8 | 4 (5-8) | 0.1617 (6) |
| mt085 | 454 | 418 | 38.2% / 56.2% | 8 | 4 (5-8) | 0.1660 (8) |
| mt080 | 1,118 | 1,059 | 30.6% / 29.3% | 14 | 10 (5-14) | 0.1886 (5) |
| mt075 | 2,526 | 2,450 | 24.5% / 20.9% | 12 | 8 (5-12) | 0.1889 (5) |
| variable | 2 | 0 | 99.6% / 99.7% | 1 | 0 | degenerate: "1 units above the size floor, need at least 5" |
| regularized | 7,208 (8,261 leaf, 2 levels) | 8,151 leaf | 8.0% top / 6.6% | 18 | 14 (5-18) | 0.2141 (10) |
| entropy | 15 | 3 | 34.6% / 52.0% | 6 | 2 (5-6) | 0.1442 (5) |
| recurse50 | 15 | 3 | 34.6% / 52.0% | 6 | 2 (5-6) | 0.1442 (5) |

No unit was unassigned in any configuration. The structure: from 1.0 to 0.85 the map is one module of 52-57% of the area, a few large ones and dust, and on the screen only 0.80 and 0.75 break the giant module up (at 10 trials `mt080` does not; see below); `regularized` has no giant module but 68.7% of its area in single-unit top modules; `variable` is one module. `default`, `entropy` and `recurse50` give the **same partition** (V-measure 1.0000, labels identical, top and leaf): `recurse50` never fires because the largest module holds 34.6% of the cells, under 50%, and `entropy_corrected` changes nothing on this network (not investigated further). The other pairs on the screen partitions are all under 0.91 (mt095-mt090 0.904, mt090-mt085 0.905, mt080-mt075 0.803). So ten configurations give eight distinct partitions, and `variable` is degenerate: seven could be bootstrapped.

**Shortlist (ruling R3), and why selection was widened.** R3 ranks non-degenerate configurations by their best screening excess over all k, orders ties by the smaller largest module, takes at most three and skips a copy of a listed partition: `regularized` 0.2141 (k=10), `mt075` 0.1889 (k=5), `mt080` 0.1886 (k=5). `mt075` and `mt080` are 0.0003 apart, so R3's smaller-largest-module clause orders them (24.5% against 30.6% of units, 20.9% against 29.3% of the area). The next was `mt085` at 0.1660. `entropy` and `recurse50` were skipped as copies of `default`, which ranked sixth. That shortlist rested on 3-trial partitions, and they did not all hold at 10 trials: `mt080` went from 1,118 modules with the largest at 29.3% of the area to 853 modules with the largest at 55.8% (816 single-unit modules, 37.5% of units), `mt095` from 21 to 47 modules and `mt090` from 109 to 122, while `default` (15), `mt085` (454), `mt075` (2,526) and `regularized` (8,261) did not change. The 10-trial `mt080` partition is close to `mt085`'s (V 0.881) and `mt090`'s (V 0.806), and further from `mt095` (0.797), `default` (0.684), `mt075` (0.658) and `regularized` (0.495), so the screen's order, and "mt080 breaks the giant module up", could not be trusted. A select costs 14 to 19 minutes, not the ~3 hours the plan assumed for "at most 3", so selection was run for `default`, `mt095`, `mt090` and `mt085` as well, which covers every distinct non-degenerate configuration. R3's shortlist stays on record as the screen ranking.

**Select** (10 trials, 20 bootstraps of 3 trials, 99 null maps; `selection.csv`). Stability is the median bootstrap V-measure; "scr" is the screening excess at the same k. Eligible means stability at least 0.8. The seven configurations, with their most stable map:

| config | maps (k) | most stable k | stability | excess | `fell_short` |
|---|---|---|---|---|---|
| default | 2 (5-6) | 6 | 0.674 | 0.1253 | yes |
| mt095 | 4 (5-8) | 8 | 0.759 | 0.1611 | yes |
| mt090 | 3 (5-7) | 6 | 0.770 | 0.1451 | yes |
| mt085 | 4 (5-8) | 8 | 0.730 | 0.1596 | yes |
| **mt080** | 6 (5-10) | 9 (7, 9, 10 are at or above 0.8) | 0.814 | 0.1535 | **no, k=7 chosen** |
| mt075 | 8 (5-12) | 12 | 0.723 | 0.1749 | yes |
| regularized | 14 (5-18) | 10 | 0.682 | 0.2068 | yes |

The three of the shortlist at every k:

| config | k | stability | excess | scr |
|---|---|---|---|---|
| regularized | 5 / 6 / 7 / 8 | 0.634 / 0.638 / 0.619 / 0.599 | 0.1586 / 0.1644 / 0.1575 / 0.1661 | 0.1447 / 0.1491 / 0.1484 / 0.1556 |
| | 9 / **10** / 11 / 12 | 0.641 / **0.682** / 0.670 / 0.673 | 0.2031 / 0.2068 / 0.2102 / 0.2057 | 0.2074 / 0.2141 / 0.2136 / 0.2000 |
| | 13 to 18 | 0.662 to 0.671 | 0.1920 to 0.2023 | 0.1951 to 0.2025 |
| mt075 | 5 / 6 / 7 / 8 | 0.600 / 0.665 / 0.704 / 0.717 | 0.1875 / 0.1780 / 0.1459 / 0.1506 | 0.1889 / 0.1846 / 0.1382 / 0.1479 |
| | 9 / 10 / 11 / **12** | 0.709 / 0.708 / 0.722 / **0.723** | 0.1462 / 0.1693 / 0.1781 / 0.1749 | 0.1406 / 0.1642 / 0.1757 / 0.1815 |
| mt080 | 5 / 6 | 0.763 / 0.738 | 0.1722 / 0.1739 | 0.1886 / 0.1724 |
| | **7** / 8 / 9 / 10 | **0.806** / 0.793 / 0.814 / 0.802 | **0.1568** / 0.1564 / 0.1535 / 0.1580 | 0.1374 / 0.1310 / 0.1324 / 0.1576 |

**The choice, by the design's rule across the seven.**
1. Eligible (stability at least 0.8): only `mt080` k=7, 9 and 10. The other six configurations have no eligible map. Their higher excess (`regularized` up to 0.2102, `mt075` 0.1875 at k=5) therefore does not count, and neither does `mt080`'s own k=5 and k=6, which have a higher excess (0.1722 and 0.1739) than k=7 and are out only by stability (0.763 and 0.738).
2. Highest excess among the eligible: `mt080` k=10, 0.1580.
3. The rule takes the coarsest k within 0.01 of the top excess: k=7 (0.1568, 0.0012 below) and k=9 (0.1535, 0.0045 below) are both within 0.01, and **k=7 is the coarser**.
4. Not reached: something was eligible, so no map falls short, and the beta-sim fallback does not apply (no bootstrapped configuration is degenerate; `variable` is, and was not bootstrapped).
Result: `fell_short` is false for `mt080`.

**The chosen map** (`mt080/candidates.parquet` and `parts.parquet`, k=7; 11,292 units, 303.0M km2: 125.2M km2 land, 177.9M km2 sea).

| realm | units (land / sea) | land km2 | sea km2 | share of area | imputed parts | separate pieces (largest, others) |
|---|---|---|---|---|---|---|
| 1 | 4,805 (197 / 4,608) | 0.99M | 175.62M | 58.3% | 5,940 of 14,516 (41%) | 53 (4,568 units; 52 pieces of 1-49) |
| 0 | 2,264 (2,182 / 82) | 44.45M | 0.97M | 15.0% | 489 of 3,027 (16%) | 33 (2,117; 32 of 1-33) |
| 5 | 1,019 (1,004 / 15) | 22.70M | 0.24M | 7.6% | 103 of 1,215 (8%) | 11 (1,006) |
| 3 | 1,061 (1,036 / 25) | 19.46M | 0.12M | 6.5% | 35 of 1,097 (3%) | 9 (1,046) |
| 4 | 999 (928 / 71) | 15.24M | 0.92M | 5.3% | 81 of 1,133 (7%) | 19 (971) |
| 2 | 722 (717 / 5) | 14.61M | 0.00M | 4.8% | 2 of 724 (0.3%) | 4 (719) |
| 6 | 422 (415 / 7) | 7.72M | 0.00M | 2.5% | 10 of 435 (2%) | 4 (419) |

The largest realm (1) holds 58.3% of the area and 42.6% of the units, and is almost all sea (4,608 of its 4,805 units). Every realm has one main piece holding 93.5% (realm 0) to 99.6% (realm 2) of its units (realm 1: 95.1%), plus small separate pieces (by adjacency of units, so an island is its own piece). 6,660 of the 22,147 parts are imputed (6,042 sparse, 618 failing the composition gate), 30% of the parts overall but **41% of the ocean realm's parts** (5,940 of its 14,516), so the ocean realm is the one most shaped by the imputation. 1,032 units were merged into a realm by the size floor from 843 sub-realms. The shares of area here are of the units (303.04M km2), the basis the select stage scores on; the parts, imputed ones included, cover 446.46M km2, of which 143.41M km2 (32%) is imputed and 132.68M km2 of that in the ocean realm, so on the parts' basis the ocean realm is 69.3% of the map (`review/realms-k7.csv` gives both).

**Caveats.**
- **A single seed.** `Params.seed` is 1 in code, for Infomap, the bootstraps (seed plus the bootstrap number) and the null maps, and a second seed was not run. Nothing says how far another seed would move a stability, or the partition (Infomap at 10 trials changed three of the screen's partitions).
- **The bootstraps use 3 trials against a 10-trial reference.** A bootstrap partition is the weaker clustering of a resample, which can lower V-measure by itself, in every configuration alike, so it may also have kept others under the bar. Both bear on `mt080`'s stability, which is 0.806 at k=7, 0.006 over the bar, with 0.814 and 0.802 at k=9 and 10 and 0.793 at k=8. With 20 bootstraps a different seed could put k=7 under 0.8 and `mt080` among the configurations that fell short. The choice rests on a margin that is narrow, and on the bar rather than on excess (above).
- **The screen's reading of `mt080` did not hold** (above): the giant module did not break up at 10 trials, and the chosen map has a realm of 58% of the area (the sea). The k range also changed, from 5-14 to 5-10.
- **Screening excess is a rough guide.** Where the partition is the same in both stages (`regularized`, `mt075`), select's 99 nulls differ from the screen's 19 by up to 0.015 (`regularized`, k=6; `mt075` 0.008): that is the noise of the screen, as `screen.csv` holds only k and excess and no null spread. `mt080`'s difference is up to 0.025, from the partition. The screen's order of `mt075` and `mt080` (0.0003) was noise; `regularized`'s lead of 0.025 was above that noise.
- **The k=5 peaks of the screen** (`mt075` 0.1889, `mt080` 0.1886, at the smallest candidate) did not carry into stability: k=5 is 0.600 for `mt075` and 0.763 for `mt080`.
- **`regularized`** has the highest excess (about 0.20 to 0.21 for k=9 to 18) but 68.7% of its area is in single-unit top modules, so most of its units reach the realms through the size floor, and its stability is 0.60 to 0.68. Not a map to select.
- **`mt080` k=7 has the lowest excess of any configuration's selected map but two** (0.157; `default` 0.125 and `mt090` 0.145 are lower): a stable map was chosen over higher-excess unstable ones, as the rule says.
- `default`, `entropy` and `recurse50` are one partition (above).
- Task 10 (the WCVP circularity gate) and Task 11 (the beta-sim cross-check) use `mt080`, the k=7 map; there were no unassigned units, so no `warning`.

**Cost.** The screen took 109 to 143 s per configuration (30 s for `variable`, which was degenerate), 4 workers, mostly the Infomap trials. The `log.json` peak of 6.5 to 14.9 GB is a high-water mark of the one process across the run (it never falls), not the memory of a configuration. Select took 850 s (`regularized`), 949 s (`mt075`) and 1,036 s (`mt080`) at a peak of 5.3 to 5.8 GB, and 1,056 to 1,149 s (`mt085`, `mt090`, `default`, `mt095`) at 7.1 GB: 14 to 19 minutes, 4 workers, of which the bootstraps are 495 to 724 s.

### The circularity gate on half A (2026-11-02)

**Headline mode: WCVP `native_only`.** Every WCVP mode passes the TDWG-edge gate, under both readings, so the rule takes the strongest, `native_only`. Half A only, at resolution 7 (release 315557), `mt080` (Markov time 0.8) on all five runs; nothing from B. `scripts/circularity.py` (helpers in `src/realms/circularity.py`, tests `tests/test_circularity.py`) reads each run's `parts.parquet` and writes `results/explore/res7/circularity/`: `runs.csv`, `pairs.csv`, `moves.csv`, `near_by_region.csv`, `gate.csv` and `run.log`. It took 191 s, of which 54 s went on the rotations at 0.5 and 68 s on each gate.

**The runs** (`--stage select`, `mt080`). Stability and falling short are context here, not a criterion:

| run | selected k | stability | fell short | parts | land / sea boundary km |
|---|---|---|---|---|---|
| `off` | 7 | 0.730 | yes | 22,159 | 34,567 / 69,608 |
| `introduced` | 7 | 0.762 | yes | 22,157 | 40,728 / 70,231 |
| `native_only` | 7 | 0.806 | no | 22,147 | 41,473 / 61,495 |
| no layer 2 (`native_only-no-extinct-no-no_wild_range`) | 12 | 0.801 | no | 22,167 | 49,374 / 193,953 |
| no GRIIS (`native_only-no-griis-no-griis_subnational`) | 9 | 0.699 | yes | 22,149 | 36,390 / 69,332 |

Under `mt080`, `off` and `introduced` fall short of 0.8, as does the run without GRIIS. The candidate ks run from 5-8 (`off`) to 5-12 (no layer 2). k = 5 to 8 is common to all five, and **k = 7** is the common k: `circularity.common_k` takes the shared k that the most runs selected (three of the five), and the files confirm it.

**Comparison by part, not by unit.** Cleaning changes the network units (11,292 to 11,340), so maps are compared over the parts both runs map (`parts.parquet`, imputed parts included, `realms.score.v_measure`). Every pair shares 22,144 to 22,156 parts, and 0 to 22 parts are in only one run (table). A moved boundary is the symmetric difference of the two runs' `boundary(labels, kind, edges)` edges, in geodesic km. An edge that both runs draw does not move, even where the realm on either side changed.

| pair | k | V, selected k | V, k = 7 | shared parts | only in a / b | land moved km, selected / k = 7 | sea moved km, selected / k = 7 |
|---|---|---|---|---|---|---|---|
| off / introduced | 7 / 7 | 0.890 | 0.890 | 22,156 | 3 / 1 | 37,549 / 37,549 | 30,045 / 30,045 |
| off / native_only | 7 / 7 | 0.873 | 0.873 | 22,145 | 14 / 2 | 39,383 / 39,383 | 32,841 / 32,841 |
| off / no layer 2 | 7 / 12 | 0.644 | 0.749 | 22,145 | 14 / 22 | 39,479 / 41,962 | 201,843 / 90,150 |
| off / no GRIIS | 7 / 9 | 0.949 | 0.951 | 22,144 | 15 / 5 | 16,428 / 16,543 | 23,348 / 23,348 |
| introduced / native_only | 7 / 7 | 0.953 | 0.953 | 22,146 | 11 / 1 | 17,034 / 17,034 | 20,441 / 20,441 |
| introduced / no layer 2 | 7 / 12 | 0.695 | 0.819 | 22,146 | 11 / 21 | 30,654 / 27,406 | 199,261 / 87,741 |
| introduced / no GRIIS | 7 / 9 | 0.890 | 0.887 | 22,145 | 12 / 4 | 36,479 / 34,082 | 29,247 / 29,247 |
| native_only / no layer 2 | 7 / 12 | 0.712 | 0.841 | 22,147 | 0 / 20 | 27,194 / 23,946 | 189,605 / 78,085 |
| native_only / no GRIIS | 7 / 9 | 0.885 | 0.883 | 22,146 | 1 / 3 | 34,560 / 31,218 | 34,493 / 34,493 |
| no layer 2 / no GRIIS | 12 / 9 | 0.659 | 0.761 | 22,147 | 20 / 2 | 34,231 / 35,256 | 197,221 / 85,904 |

**Moves by TDWG level 1 region** (land, at the selected k; `moves.csv` has every pair at both ks). An edge goes to the level 1 region whose level 3 polygons hold the centre of its first cell (`cell_a`, the lower id as text). If that centre is in no region, the second cell's is used, and otherwise the edge goes to `none`. A centre on an edge goes to the lowest level 3 code. Unlike `units.part_regions`, no nearest region is taken, so `none` holds coastal edges whose two cell centres both lie on the sea: 6,268 of the 21,872 cell centres are in a region.

| pair | moved land km | top regions |
|---|---|---|
| off / introduced | 37,549 | Asia-Temperate 12,619 (34%), Northern America 12,182 (32%), none 5,717 (15%) |
| off / native_only | 39,383 | Asia-Temperate 13,077 (33%), Northern America 11,981 (30%), none 6,275 (16%) |
| introduced / native_only | 17,034 | Asia-Temperate 10,292 (60%), none 3,966 (23%), Asia-Tropical 1,626 (10%) |
| native_only / no layer 2 | 27,194 | Asia-Temperate 7,933 (29%), none 6,657 (24%), Northern America 5,882 (22%), Southern America 3,254 (12%) |
| native_only / no GRIIS | 34,560 | Northern America 11,981 (35%), Asia-Temperate 11,192 (32%), none 4,200 (12%) |

**The gate** (`boundaries.gate`). `off` is the base, the regions are TDWG level 3, there are 99 rotations with seed 1, and each map is taken at its selected k, land boundaries only. All five runs went in, so the two layer variants have rows too, for information. `off`'s share near TDWG edges is 0.8328.

| mode | share | raise | null q95, `land_tolerance=0.5` | passes | null q95, `None` | passes |
|---|---|---|---|---|---|---|
| introduced | 0.8287 | -0.0041 | +0.1674 | yes | +0.1633 | yes |
| native_only | 0.7952 | -0.0376 | +0.1617 | yes | +0.1614 | yes |
| no layer 2 | 0.7646 | -0.0682 | +0.1313 | yes | +0.1264 | yes |
| no GRIIS | 0.8381 | +0.0053 | +0.0507 | yes | +0.0507 | yes |

**The two readings agree for every mode.** The no-GRIIS q95 is the same under both (0.0506671). That was checked, and it is a coincidence of order statistics: the 94th and 95th of the 99 sorted null raises come from the same two draws in both sets (seed 1 gives the same first draws, and the 0.5 filter keeps both). Neither WCVP mode raises its share near TDWG edges at all: both are below `off`, `native_only` by 0.038.

**The headline, by the rule (brief Step 3):** `native_only`, then `introduced`, then `off`. `native_only` passes at the default 0.5 reading, so it is the headline mode (`circularity.headline`). The layer variants are not candidates. Endemism excess plays no part.

**What layer 2 alone changes** (`native_only` against no layer 2). This is the largest change of any switch. Without layer 2 (extinct species and species with no wild range kept), the selected map is k = 12, not 7 (stability 0.801). Its V-measure with `native_only` is 0.712, or 0.841 at k = 7. The sea boundary triples (61,495 to 193,953 km at the selected k; 78,085 km of sea boundary moves even at k = 7). The land boundary that moves is 27,194 km, spread over Asia-Temperate, Northern America and Southern America, with 1,052 km in Antarctica against 90 km in every other run. Layer 2 matters most on the sea.

**What GRIIS alone changes** (`native_only` against no GRIIS). Without GRIIS (country and sub-national lists), the selected map is k = 9 at stability 0.699, short of the bar. Its V-measure with `native_only` is 0.885 (0.883 at k = 7), and 34,560 km of land boundary moves, 35% of it in Northern America. The no-GRIIS map is the closer one to `off` (V 0.949, 16,428 km moved). The Northern American boundary shows why. `off` and no GRIIS draw 2,369 km of land boundary there, against 11,315 km for `introduced` and 11,114 km for `native_only`, so the extra Northern American boundary, about 8,700 to 8,900 km net (8,745 km for `native_only`, 8,946 km for `introduced`), appears only when WCVP and GRIIS clean together. Neither does it alone.

**Caveats.**
- **The land measure is mostly coastline.** The rotations' land filter tracks where the coasts land more than the internal TDWG borders (Task 6). It stays as built and is on Markus's list of open points.
- **The gate has little power here.** At resolution 7, 60% of the cells with a land part (4,630 of 7,697) are within one cell of a TDWG edge, so any map's land boundary is mostly "near": the share is 0.80 to 0.84 for every run (`gate.csv`). The spread of the null comes from an ad hoc check, not from an output file: a throwaway script drew the same 99 rotations (`rotations(regions, 99, seed=1)`, tolerance 0.5) and took `edge_share` of each map's land boundary against `near_cells` of each rotated copy, as `gate` does internally (its 95th percentiles match `gate.csv`). On that check, `off`'s rotated share runs from 0.20 to 0.65 (5-95%, median 0.44), and the null raises are wide: for `native_only` from -0.142 to +0.162 (5-95%), median +0.028. A mode would need to raise its share by about 0.16 to fail. The WCVP modes lower it, so the verdict is not close, but the gate could miss a regional alignment that a global share dilutes.
- **A regional signal, outside the rule** (`near_by_region.csv`; an edge is near when either of its cells is, as in the gate). Against `off`, `native_only` adds 10,363 km of land boundary in Northern America and removes 1,618 km (net 8,745 km). The added km are 92.6% near TDWG edges (`introduced`: 10,564 km added, 90.9% near). The base rate is edge-level: of every edge between two land parts attributed to Northern America, 78.9% of the km are near, and 83.7% of those with both centres in Northern America. So the added boundary is 9 to 14 points above the base rate, and less near than `off`'s own Northern American boundary (2,369 km, 96.5% near; the same for no GRIIS). The 1,618 km removed were 100% near. Here the TDWG level 3 regions are US states and Canadian provinces. This does not change the verdict, which is defined on the whole map, but the WCVP-with-GRIIS boundary in Northern America is the place to look for a TDWG artefact. Its shape was not examined further (it runs from 23 to 77 degrees N).
- **A single seed.** The rotations (seed 1) and every run's Infomap, bootstraps and null maps (`Params.seed` 1) were run once.
- **`off` and `introduced` fall short of 0.8** under `mt080` (0.730 and 0.762), and so does no GRIIS (0.699). Their maps are less stable than the bar asks. The V-measures between runs (0.87 to 0.95 among the k = 7 maps) are about as high as the runs' own bootstrap stability (0.70 to 0.81), so part of every "move" may be Infomap's own variation rather than the cleaning. The pairs were not compared against a same-cleaning, other-seed baseline.
- **Parts present in only one run.** 0 to 22 parts per pair. V-measure leaves them out, but they can put edges in one run's boundary only, so they count among the moved km.
- **Different ks.** Pairs with no layer 2 (k = 12) or no GRIIS (k = 9) compare maps of different k at the selected k, which inflates the moves. The k = 7 columns compare like with like.

### Beta-sim and duplicate-name cross-checks (2026-11-02)

**The headline map is kept: by ruling, the cross-checks report agreement only** (`native_only`, `mt080`, k = 7). Neither check is degenerate, so the fallback does not apply, and the cross-check reports agreement only. Half A only, resolution 7 (release 315557); nothing from B. `scripts/crosscheck.py` (helpers in `src/realms/crosscheck.py`, tests `tests/test_crosscheck.py`) writes `results/explore/res7/crosscheck/`: `betasim_v.csv`, `betasim_seeds.csv`, `betasim_area.csv`, `betasim_labels.parquet`, `names_v.csv`, `names_selection.csv`, `names_duplicates.csv` and `run.log`.

**Step 1: beta-sim consensus on the headline's units.** The matrix is the one the select stage scored, built by `explore.network_matrix` and restricted by `explore.assigned_only` with the unassigned units of `partition.parquet` (none). It has 11,292 units and 666,842 taxon columns. `rank_beta_sim` took 14 s and 2.6 GB. `consensus` ran with `ids` (the unit ids), k = 5 to 12, 20 shuffles, and seeds 1, 2 and 3, for Ward and for UPGMA (`average`). Each run took about 73 s for Ward and about 153 s for UPGMA, at a 4.5 GB peak, 692 s in all. **Infomap's unit labels** come from the parts of `parts.parquet` at each k that are not imputed. Each such part carries its unit, every part of a unit has the same realm (`crosscheck.unit_labels` raises otherwise), and the result equals `candidates.parquet` at every k, which the script checks. V-measure is taken at unit level, over the 11,292 units that both hold.

| k | 5 | 6 | **7** | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|
| Ward, V against Infomap | 0.550 | 0.572 | **0.586** | 0.633 | 0.608 | 0.629 | | |
| UPGMA, V against Infomap | 0.437 | 0.409 | **0.463** | 0.475 | 0.469 | 0.422 | | |
| Ward, spread | 0.891 | 0.916 | 0.924 | 0.929 | 0.947 | 0.949 | 0.943 | 0.938 |
| UPGMA, spread | 0.739 | 0.995 | 0.995 | 0.995 | 0.950 | 0.996 | 0.983 | 0.992 |

These are seed 1. Seeds 2 and 3 give Ward 0.580 and 0.586 at k = 7, and UPGMA 0.463 for both. Infomap has no map beyond k = 10.

- **Seed sensitivity** (V between the consensus groupings of two seeds, k = 7): Ward 0.975 (seeds 1/2), 0.999 (1/3) and 0.974 (2/3); UPGMA 1.000 for all three pairs. Over every k from 5 to 12 the lowest is Ward's 0.950 at k = 11, and UPGMA is 1.000 throughout. No pair fell below 0.9, so there was no rerun with 100 shuffles. On the real matrix the consensus is much steadier than on the tied synthetic matrices of the beta-sim entry (0.75 to 0.95).
- **Area agreement at k = 7** (realms matched one to one by the Hungarian method on the table of shared unit area; `crosscheck.area_agreement`): **Ward 0.593** (seed 2: 0.583), **UPGMA 0.723** (all seeds).
- **Where they differ.** Infomap's realm 1 at k = 7 is the ocean: 97.5% of its parts are sea parts, and it holds 4,805 units and 58.3% of the area. Both methods respect its boundary. UPGMA is almost a two-way split. One group is the ocean (58.3% of the area, 99% of it in Infomap's ocean). The other lumps all six land realms (41.4%). The five remaining groups hold 23 units between them, which is the chaining of average linkage. UPGMA's high area agreement is therefore the ocean, plus the largest land realm matched to the lumped group. Ward splits the ocean into three groups (35.1%, 16.3% and 6.5% of the area, each 97-100% within Infomap's ocean) and has four land groups. One of those matches an Infomap realm (realm 3, 0.064 of the group's 0.071), and the other three each mix two or three Infomap land realms. Homogeneity (beta-sim groups within one Infomap realm) and completeness are 0.630 and 0.548 for Ward, and 0.330 and 0.776 for UPGMA.

**Step 2: duplicate names merged** (`--merge-duplicate-names`, the run `native_only-merged-names`, `mt080`, `--stage select`, `--workers 4`). The flag reaches `cell_taxa` through `explore.network_matrix`, which builds the one matrix that the screen and the select stage share. It enters the hash of the clean stage (whose `taxa.parquet` then keeps `kingdom`, `order`, `family` and `genus`) and of each configuration, only when it is on. `test_merging_duplicate_names_enters_the_hashes_only_when_it_is_on` pins the default hashes to the values they had before the option existed. The real run on ids reports every stage skipped (`--stage select --only mt080 --wcvp-mode native_only`: clean `e06a3837…`, units `25bcc4b6…`, mt080 screen and select skipped). The merged run is a directory of its own, so clean and units ran there too. The units stage does not use the matrix, but its hash is that of cleaning. Clean took 76 s at a 21.2 GB peak, against 18.0 GB on ids. The cleaned presence and the units (`members`, `areas`, `excluded`) are identical to the run on ids, so only the matrix differs. Among half A's 588,594 cleaned species, 973 genus groups (2,003 ids) merge **1,030 of 73,398 genus ids** away, and 238 family groups (487 ids) merge **249 of 5,964 family ids** away. The select stage took 1,070 s at an 8.5 GB peak.

| k | stability, ids | stability, merged | diff | excess, ids | excess, merged | diff |
|---|---|---|---|---|---|---|
| 5 | 0.763 | 0.667 | -0.096 | 0.1722 | 0.1574 | -0.0148 |
| 6 | 0.738 | 0.716 | -0.021 | 0.1739 | 0.1569 | -0.0170 |
| **7** | **0.806** (selected) | 0.705 | -0.101 | 0.1568 | 0.1613 | +0.0045 |
| 8 | 0.793 | **0.718** (selected) | -0.076 | 0.1564 | 0.1613 | +0.0050 |
| 9 | 0.814 | | | 0.1535 | | |
| 10 | 0.802 | | | 0.1580 | | |

- **Selected k:** 7 on ids (stability 0.806, not short). With names merged it is **k = 8, and it falls short**: no map reaches 0.8, and k = 8 is the most stable at 0.718. Infomap found 952 modules (901 of them single-unit), against 853 (816). The largest module holds 36.3% of the units and 54.6% of the area, against 37.5% and 55.8%. The candidate maps run from k = 5 to 8, against 5 to 10.
- **V-measure by part** (`circularity.compare`, imputed parts included): **0.883** between the selected maps (k = 7 on ids, k = 8 merged), and **0.888** at k = 7. Both runs hold all 22,147 parts, none of them in only one run.
- **Excess** moves by -0.017 to +0.005, inside the flat range Infomap's scale found. Stability falls by 0.02 to 0.10 at every k the two runs share.

**What this says about the headline's robustness.**
- **Its first split is robust; its land realms are Infomap's.** Two unconstrained similarity methods on the same units agree with Infomap on the ocean-land boundary and on almost nothing finer. UPGMA does not resolve the land at all at k = 7, and Ward resolves it differently. So the seven-realm map is not something any reasonable clustering of these units would draw. The land divisions come from the network method, which is what the design chose, but they have no independent confirmation from beta-sim at this k.
- **Its stability is marginal.** Changing 1,030 genus and 249 family columns, about 1.4% and 4.2% of them, moves the map about as much as switching the WCVP mode does (V 0.888 at k = 7, against 0.87 to 0.95 between the cleaning runs in the circularity gate). It also takes every k below the 0.8 bar. The headline's 0.806 is 0.006 above it. The map that merged names give at k = 7 is close to the headline (V 0.888), but whether k = 7 clears the bar depends on details this small.
- **Excess does not separate the two runs.** The merged run's selected map has an excess of 0.1613 (k = 8), against 0.1568 for the headline and 0.1535 to 0.1580 for its three stable maps. Merging makes the realms neither clearly more endemic nor clearly less.

**Caveats.**
- **Unit level against area level.** The beta-sim V-measures count units. A coarse unit (810 of them, from the res-5 grid) and a sparse ocean unit weigh the same as a fine land unit, and Infomap's ocean realm holds 43% of the units. The area agreement is the area-weighted reading, and the two disagree in rank: UPGMA scores lower than Ward in V, and higher in area. The merged-name comparison is by part, imputed parts included, and is not area-weighted either. The area agreement and the shares of this entry are of the units (303.04M km2, the ocean realm 58.3%); the parts with the imputed ones cover 446.46M km2 (143.41M km2 imputed, 132.68M km2 of it in the ocean realm), on which the ocean realm is 69.3%.
- **A single Infomap seed.** The beta-sim consensus was measured at three seeds, but Infomap ran once per run (`Params.seed` 1, 10 trials, with bootstraps and null maps from seed 1). The difference between the merged run and the run on ids includes Infomap's own run-to-run variation, which was not measured, so the drop in stability cannot be put on the merging alone. A same-matrix, other-seed baseline would separate the two.
- **The beta-sim k runs only to 12**, and the comparison with Infomap only to 10, Infomap's last candidate. Ward's V rises from 0.550 at k = 5 to 0.61-0.63 at k = 8-10, so the methods may agree better at finer k than the headline's.

### Maps at 19 and 30 realms without the size floor (2026-11-02)

**The gap.** The design (Choosing the headline level) promises candidate maps at 19 and 30 realms for M4's comparison, because merging Infomap's finest modules one pair at a time gives one nested map per realm count from 60 down to 5. With the size floor applied first, as `candidates.candidate_maps` does, that promise does not hold for the headline. Under `mt080` on `native_only` (half A, resolution 7, release 315557), only 10 of Infomap's 853 modules clear the 1M km2 floor. The other 843 become sub-realms, so the merge sequence starts at k = 10, and the select stage has no map at 19 or 30 realms.

**Ruling (controller, Task 12).** The maps at 19 and 30 realms are drawn with `candidates.merge_sequence` straight on Infomap's modules before the size floor. That is the same beta-sim pair-merging, on the matrix the select stage scored (`explore.network_matrix`, restricted by `explore.assigned_only`). `review.unfloored_maps` takes the finest modules (`Partition.leaf_labels` from `partition.parquet`), resolves their separate pieces with `candidates.resolve_pieces` as `candidate_maps` does first, and merges with no floor in between. The parts are imputed as the select stage imputes them (`explore._part_rows`). `candidates.py` is unchanged. `scripts/review_maps.py` checks the reconstruction: the same modules, floored and merged, give `candidates.parquet` at every k from 5 to 10, and its 843 sub-realms. On the figures and in the tables these maps are labelled "without the size floor". They were not bootstrapped or scored for excess, and **they are never used to choose anything**. M4 decides how to compare them. That includes whether a map without the floor is the right counterpart to schemes of 19 and 30 realms, or whether the floor, or the merge order, should change for that comparison. The merge took 11 s for both maps from 853 modules.

**What they look like** (`results/explore/res7/review/`: `k19-no-floor.png`, `k30-no-floor.png`, `no_floor_realms.csv`, `no_floor_parts.parquet`, `no_floor_units.parquet`).
- At **k = 19**, the ten largest realms carry the labels of the ten units above the floor (`merge_sequence` keeps the label of the unit merged into), from the ocean (308M km2 of parts, 69% of the map) down to New Guinea (1.0M km2). The areas of this entry are of parts, imputed ones included (446.46M km2 in all at k = 7, 143.41M km2 of it imputed and 132.68M km2 of that in the ocean realm); the scale and cross-check entries give the units' area (303.04M km2), on which the ocean realm at k = 7 is 58.3%, not 69.3%. The other nine are modules the floor would have absorbed: Madagascar (0.82M km2), New Zealand (0.27M km2; in realm 0 at k = 7), and seven islands or island groups of 0.3 to 42 thousand km2 (by their mean cell centre: Fiji, Socotra, New Caledonia, Hawaii, Kerguelen, the Mascarenes and the Seychelles).
- At **k = 30**, the eleven more realms include an Antarctic coastal realm (3.1M km2, 206 land parts, 164 of them imputed; at k = 7 these parts are in realm 0), a sea realm in the tropical eastern North Pacific (1.3M km2, 57 sea parts, 41 imputed; the largest sub-realm of the ocean realm at k = 7), and more single islands. The smallest realm has 153 km2 of parts at k = 30 and 257 km2 at k = 19.
- Without the floor, beta-sim merging keeps small, distinct islands apart longest. So between 10 and 30 realms the count grows mostly by islands, not by splitting the large realms. This is a property of the candidates, recorded for M4. It is not a judgement on the maps.

### Changes after the final review of M1 (2026-11-02)

The final review of the branch (99bfa17..5453f6f) found the pipeline's caching weaker than DECISIONS said, and three readers of a select result that had drifted apart. The fixes keep every hash of M1, so no run was redone.

- **Select has a directory of its own.** Screen and select wrote the same map files into `<config>/`, and after a select the screen's stamp still said current over select's 10-trial files (`native_only/mt080` had both stamps and select's `structure.json`). Select now writes `<config>/select/` (map files, `selection.csv`, `selected.json`, `stage-select.json`, `log.json`); the screen keeps `<config>/`. `scripts/migrate_select.py` (`explore.migrate_select`) moved the 12 select results of M1 at resolution 7 (`native_only` default, mt095, mt090, mt085, mt080, mt075, regularized, and mt080 of off, introduced, `native_only-merged-names`, `native_only-no-extinct-no-no_wild_range`, `native_only-no-griis-no-griis_subnational`), split their `log.json`, and deleted their 12 screen stamps, since the screen's map files were overwritten and are gone. `screen.csv` stays as the screen's record, with a `MIGRATED.txt` that says so; a screen asked for again reruns. The screen of `native_only/mt080` was run again (`--stage screen --only mt080 --workers 4`, 109 s): same hash, `screen.csv` byte-identical to the old one, and its own 3-trial `structure.json` (1,118 modules, against select's 853). The other eleven screens were not run again. The screens of entropy, recurse50 and variable never had a select and were left alone.
- **A rewritten stage takes its dependants' stamps.** The hashes are of a stage's inputs, never of what the stage before wrote, so `--stage clean --force` after a fix to `clean.py` left the units, screen and select current on the old cleaning. Chaining the hashes on the outputs would have changed every hash of M1 and forced a full rerun. Instead, a stage that writes its outputs deletes the stamps of the stages that depend on it in the run (`explore.invalidate_after`: cleaning those of the units and of every configuration's screen and select, the units those of every configuration; screen and select none). And each stage has a `METHOD_VERSION` (`explore.METHOD_VERSION`, all 1), which enters its hash only when it is not 1: raise it when the stage's behaviour changes with its inputs the same. A test pins the four at 1 and the M1 hashes. What this still does not catch: a parent rerun outside the pipeline (files replaced by hand), which only content hashes would.
- **One reader of a select result.** `explore.read_selected(run_dir, config, configs_path, seed=, merge_duplicate_names=)` checks that the units stamp is of the cleaning stamp's hash, recomputes the select hash from the units stamp, the configuration's options, the parameters of `configs.toml` at the seed and the merge flag (by default, from the run's name), and refuses a result that is not current, degenerate, or whose `parts.parquet` marks another k than `selected.json`. `scripts/circularity.py`, `scripts/crosscheck.py` and `scripts/review_maps.py` use it instead of their own readers, which checked less (circularity's not even the units). `explore.manifest_path` and `explore.part_rows` are public now; `natives._manifest_path`, which takes the manifest's content and returns None, is unchanged. Rerun on the moved results, the three scripts wrote byte-identical outputs (every CSV, `betasim_labels.parquet`, the `no_floor` parquets and `tables.md`), except `review/summary.csv`, whose `largest_realm_share` moved in the 16th significant digit: it summed the realms' areas in the order of a polars `group_by`, which was not fixed, so it varied from run to run in the last bit (30 repeats of that sum on the same `parts.parquet` gave four values). The group order is now fixed (`maintain_order=True`), so the next `review_maps.py` run settles the last digit and later runs repeat it.
- **A second seed.** `scripts/explore.py --seed N` sets `Params.seed`, which was already in the hash of screen and select, so seed 1's hashes are unchanged, and writes to `<config>-s<N>/` for N other than 1, so it never overwrites seed 1's result. No second seed was run: whether to run one on mt080 is Markus's call (the single seed is the main caveat to the 0.806 stability).
- **The half and the memory.** `Plan.half` defaults to A and refuses any other ("sealed until M3"); `run_name` and `clean_inputs` take the half, and the split file's sha256, only when it is not A, so a later run on another half can never share a stamp or a directory with half A's, and the M1 hashes stay. Nothing loads another half. `--budget-gb` (default 40, `BUDGET_GB`) is the memory the parent and the workers may hold together. **Measured clean peaks at resolution 7** (the child's, `log.json` of each run): native_only 18.0 GB, `native_only-no-griis-no-griis_subnational` 17.9 GB, `native_only-merged-names` 21.2 GB, off 23.3 GB, introduced 23.7 GB, `native_only-no-extinct-no-no_wild_range` 24.1 GB; the units 6.0 to 6.4 GB; a select with four workers 5.3 to 8.8 GB. The 14 to 18 GB and the 16 GB said earlier were too low.
- **mt080 was chosen on `native_only`** (Infomap's scale), before the circularity gate settled the headline mode; under mt080 the off and introduced runs fall short of the stability bar. If the headline mode changes, the configuration has to be chosen again on that mode.
- **Smaller.** `selected.json`'s `units` counts the units that were scored (those with a module), the same 11,292 as before since none is unassigned. `boundaries.LAND` and `CACHE` are resolved from the repository, not the working directory. `hexmaps._realm_neighbours` and `review.realm_neighbours` stay two functions (they differ on an asymmetric adjacency and in their keys, and merging them could move `plot_parts`' colours) and say so of each other. A test shows that the two pools of `measured_workers` (two workers to measure, then as many as fit) give the same `screen.csv` as one worker.

### One ocean realm is accepted (2026-11-02)

Markus reviewed the M1 maps and accepted the single ocean realm of the k = 7 map (`native_only`, mt080: 58.3% of the units' area, 69.3% of the parts' with the imputed ones) as what the data says. No method change is made to split it: a split would have to come from the endemism excess, the stability or the calibration, not from the expectation of several marine realms. The ocean realm's sub-realms (k = 30 shows a tropical eastern North Pacific one) stay available for M4's comparison.

### A second seed on mt080 (2026-11-02)

At Markus's request, `native_only/mt080` was run again at seed 2 (`--seed 2`, `mt080-s2/`; screen 2 min, select 15 min), which changes Infomap's trials, the bootstrap resamples and the null maps.

- **The maps barely move.** Infomap finds 889 modules against seed 1's 853 (unit V 0.964). At every k both seeds have (5 to 9), the maps agree at V 0.966 to 0.973 over all 22,147 parts. At k = 7 the realms are the same seven: 122 parts (0.55%, 20 of them imputed) change realm, most between realms 4 and 1, 5 and 0, and 0 and 4.
- **The stability does not hold the bar.** At seed 2 no map reaches 0.8: k = 5 to 9 score 0.683, 0.640, 0.752, 0.753 and 0.765 (seed 1: 0.763, 0.738, 0.806, 0.793, 0.814). Select picks k = 9 and marks it fell short. The excess is lower too, 0.143 to 0.160 against 0.153 to 0.174. Seed 2 has no k = 10 map, because only nine units clear the size floor.
- **What it means.** Between two seeds the stability of k = 7 moves by 0.054, while the reference maps agree at V 0.97. So most of the difference comes from the 20 bootstraps, not from the map. A gap of that size is larger than mt080's margin over the bar (0.006), and larger than the gaps that separated mt080 from mt090 (0.770), mt095 (0.759) and the default (0.674) at seed 1. The choice of mt080 in "Infomap's scale" therefore rests on one draw of the bootstraps. How to read the stability across seeds is left to Markus; no rule is changed here.

## 2026-11-03

### Stability over five seeds (2026-11-03)

Markus's ruling after the second seed (above): stability is judged over seeds, and every configuration is judged the same way. The rule is fixed here, before any further seed is run.

- **Seeds.** Every one of the seven distinct non-degenerate configurations of "Infomap's scale" (`default`, `mt075`, `mt080`, `mt085`, `mt090`, `mt095`, `regularized`) is run through screen and select at seeds 1 to 5 on `native_only`, half A, resolution 7. The design's figures stay the same for each seed: 10 trials, 20 bootstraps of 3 trials and 99 null maps.
- **Per k.** For each configuration and k, the stability is the median of the five seeds' stabilities. A seed that has no map at that k (its units above the size floor give fewer realms) counts as 0 there. The excess is the median over the seeds that have the map.
- **The choice.** The bar stays 0.8 and the tie stays 0.01. `select.select_level`'s rule is applied to the pooled table of every configuration and k: of the maps whose median stability is at least 0.8, the one with the highest median excess, the coarser map within 0.01. A tie in k goes to the higher median excess. If no map reaches 0.8, the most stable one is reported as falling short, and the choice goes back to Markus.
- **The headline map** is seed 1's map of the chosen configuration at the chosen k. Its V-measure against the other four seeds' maps at that k is reported beside it.
- **Unchanged.** The headline mode stays `native_only`. If it changes, the configuration has to be chosen again on that mode, as before. Cost: 27 more screen-and-select runs at about 17 minutes each, about 8 hours. Seed 2 of `mt080` is already done.

### Over five seeds, `mt090` at k = 7 (2026-11-03)

The rule of "Stability over five seeds" was applied as fixed above (`scripts/seeds.py`; `results/explore/res7/native_only/seeds/`: `per_seed.csv`, `across.csv`, `choice.json`, `agreement.csv`). Seeds 2 to 5 of the seven configurations ran on 2026-11-03, 27 screen-and-select runs in 7 h 44 min, none failing.

**Chosen: `mt090` (Markov time 0.90) at k = 7.** Its median stability is 0.809 and its median excess 0.148. It is the only map of the 46 whose median stability reaches 0.8.

| config | k | median stability | lowest-highest over seeds | median excess |
|---|---|---|---|---|
| mt090 | 7 | **0.809** | 0.736-0.880 | 0.148 |
| mt090 | 6 | 0.789 | 0.728-0.869 | 0.149 |
| default | 6 | 0.767 | 0.674-0.778 | 0.176 |
| mt080 | 9 | 0.765 | 0.737-0.814 | 0.146 |
| mt095 | 8 | 0.756 | 0.709-0.798 | 0.154 |
| mt080 | 7 | 0.752 | 0.739-0.806 | 0.150 |
| mt085 | 8 | 0.730 | 0.706-0.776 | 0.160 |
| mt075 | 11 | 0.701 | 0-0.722 | 0.172 |
| regularized | 12 | 0.630 | 0.617-0.673 | 0.149 |

The table gives each configuration's most stable map; `across.csv` has every k.

- **`mt090`'s map does not depend on the seed.** Infomap at 10 trials finds the same 122 modules at all five seeds; the partitions are identical, and so are the maps at k = 5, 6 and 7 (V 1.0). The seeds change only the bootstraps and the null maps. `mt080`'s partition changes with the seed: 853 modules at seed 1, 889 at seeds 2 and 3.
- **The margin is small again.** Over the five seeds the stability at k = 7 ranges from 0.736 to 0.880. The median clears the bar by 0.009, and that median is what the rule asks for. No other map comes within 0.02.
- **`mt090` has only k = 5 to 7.** Seven of its modules clear the size floor, so it has no map with more realms.
- **The map against `mt080`'s** at k = 7: V 0.892 by part. Two realms differ (named by the TDWG level 1 regions of their land parts' centres):
  - `mt080`'s Northern American realm (2) is inside `mt090`'s realm 0, with Eurasia and North Africa.
  - `mt080`'s South American realm (3) is split: 994 parts stay in the Neotropical realm (2), and the south of South America becomes its own realm 9 (102 parts, centres from 23 to 57°S).
  - The ocean realm, Africa, Australia and the tropical and eastern Asian realm are the same up to a few hundred parts.
- **Not yet redone on `mt090`.** The circularity gate (which settled `native_only`), the cross-checks and the review maps were all run on `mt080`. Until they are redone, the headline mode rests on `mt080`'s gate. Redoing them needs select on `mt090` for the four other cleaning runs and for the merged-names run.

### The gate, the cross-checks and the review on `mt090` (2026-11-03)

At Markus's request, everything that had been run on `mt080` was run again on `mt090`. First, select on `mt090` at seed 1 for the four other cleaning runs and for the merged-names run (`results/explore/res7/mt090-runs.log`, 1 h 29 min). Then `scripts/circularity.py`, `scripts/crosscheck.py` and `scripts/review_maps.py`, whose default configuration is now `mt090`. These scripts read the headline with `seeds.read_headline`, which replaces seed 1's own selection (k = 6, stability 0.770, fell short) with the k chosen over five seeds (k = 7). The other runs keep their seed-1 selection. The `mt080` outputs were kept as `circularity-mt080/`, `crosscheck-mt080/` and `review-mt080/`.

**The circularity gate passes for every mode, under both readings, and the headline mode stays `native_only`.** The base `off` has 0.773 of its land boundary near a TDWG edge, against 0.833 under `mt080`. The raises against it:

| run | k (seed 1) | stability (seed 1) | raise | null q95 (0.5 / none) |
|---|---|---|---|---|
| off | 7 | 0.806 | base | |
| introduced | 7 | 0.733 (short) | +0.010 | 0.064 / 0.072 |
| native_only | 7 (over seeds) | 0.761 (seed 1; median 0.809) | +0.003 | 0.059 / 0.056 |
| native_only, no extinct / no wild range | 6 | 0.806 | −0.026 | 0.057 / 0.077 |
| native_only, no GRIIS | 8 | 0.766 (short) | +0.051 | 0.063 / 0.065 |

- At k = 7, `native_only` and `off` agree at V 0.729, and `native_only` and `introduced` at V 0.945.
- The no-GRIIS run passes closest to its threshold (0.051 against 0.063). It is not a candidate for the headline mode.
- **Northern America is no longer an open point.** Under `mt090`, Northern America is in realm 0 in both runs. The land boundary that `native_only` adds there against `off` is 187 km (2 edges), and it removes 484 km. Under `mt080` it added 10,363 km. The one boundary left runs across northern Mexico.

**Cross-checks at k = 7.**
- The beta-sim consensus against Infomap gives V 0.554 to 0.571 with Ward and 0.455 with UPGMA, against 0.586 and 0.463 on `mt080`.
- The share of area on which they agree is 0.610 to 0.621 (Ward) and 0.765 (UPGMA), against 0.593 and 0.723.
- The consensus is stable across seeds: Ward 0.974 to 0.999, UPGMA 1.0.
- As before, it agrees on the split between ocean and land, not on the land realms.
- The merged-names run is further from the run on ids than on `mt080`: V 0.836 at k = 7, against 0.888. At seed 1 it selects k = 5, stable at 0.803.

**The review package** (`review/`) has k = 6 and 7, but not k = 8, because `mt090` has no map with more than seven realms. It also has k = 19 and 30 without the floor (now merged from 122 modules), a panel per realm and the Northern American boundary.
- `review_maps.py` does not remove the outputs of an earlier run, so `k8.png`, `k8.gpkg` and `realm-6-k7.png` from `mt080` were deleted by hand. Copies are in `review-mt080/`.
- The realm areas are on the units' basis; on the parts' basis, imputed parts included, the ocean realm is 68.9%.

| realm | mostly | unit Mkm² | unit share | endemic species at 95% |
|---|---|---|---|---|
| 1 | the oceans, with New Guinea | 175.21 | 57.8% | 86,375 of 141,850 |
| 0 | Eurasia, North Africa, North America to northern Mexico, Greenland, New Zealand, the Antarctic coast | 60.75 | 20.0% | 146,218 of 194,056 |
| 4 | Africa south of the Sahara, Madagascar, southern Arabia | 22.75 | 7.5% | 64,466 of 80,388 |
| 2 | Mexico, Central America, the Caribbean, South America north of the southern Andes | 18.35 | 6.1% | 106,206 of 136,584 |
| 3 | South and East Asia, the Philippines, Wallacea | 16.83 | 5.6% | 72,347 of 101,593 |
| 5 | Australia | 7.72 | 2.5% | 28,705 of 44,737 |
| 9 | the southern Andes and Patagonia | 1.43 | 0.5% | 4,534 of 9,365 |

- **The review's figures are seed 1's.** `summary.csv` gives seed 1's stability and excess (k = 7: 0.761, 0.136). The medians over the seeds are in `native_only/seeds/across.csv` (k = 7: 0.809, 0.148).

## 2026-11-04

### Is there a natural k? (2026-11-04)

**The question.** Markus asked whether the number of realms can be chosen from the data alone. The selection rule cannot answer that. On `mt090`, stability rises slowly (median 0.731, 0.789 and 0.809 at k = 5, 6 and 7) and excess falls slowly (0.161, 0.149 and 0.148). Neither peaks, and the steps between neighbouring k are smaller than the spread between seeds, so the 0.8 bar and the tie-break pick k = 7, not a feature of the data. Infomap does choose a number by itself, its modules (122 under `mt090`), but the step from modules to realms is a choice of scale: the Markov time, the 1M km2 floor and the merge sequence.

**Two diagnostics, which choose nothing** (`realms.breaks`, `scripts/breaks.py`; half A, resolution 7, release 315557, `native_only`, every configuration at seed 1, 199 null maps). The script rebuilds the steps of `candidates.candidate_maps` from each select result. It checks that merging the floored units gives `candidates.parquet` at every k, which holds for all seven configurations. It then runs both diagnostics on the units without the floor (down to one realm) and on the floored ones. `candidates.py`, the selection and the Costello gate are unchanged. Outputs: `results/explore/res7/native_only/breaks/` (`heights.csv`, `splits.csv`, `run.log`).

- **Merge heights** (`breaks.merge_heights`): the beta-sim of the pair merged at each step of `merge_sequence`. Under all six Markov-time configurations (`default` and `mt075` to `mt095`), the last four merges, from 5 units to 1, all happen at about the same beta-sim, 0.72 to 0.75. So the top of the hierarchy is a five-way tie, not a tree: five blocks are about equally distinct, and the data do not say which of them belong together. Under `mt090` the five are the Holarctic, the Neotropics, Africa south of the Sahara, South and East Asia, and the ocean block. Below five the merges get cheaper, steadily and without a clear gap:
  - Madagascar and Australia at 0.66 and 0.65;
  - New Zealand, the southern Andes and Patagonia, and New Guinea at 0.59, 0.54 and 0.52;
  - then islands, from 0.45 down.

  The largest step for k ≤ 10 is from Madagascar (0.657) to the five blocks (0.729), only a little larger than the steps of about 0.06 below it. `regularized` merges lower throughout (0.53 to 0.75) and has no such tie.

  Two caveats:
  - The tie may partly be a ceiling. The matrix holds genera and families, which unrelated continents share, so any two distinct biotas may sit near 0.72.
  - Beta-sim divides by the poorer unit's taxa, so a small unit looks nested in the vast ocean block and tends to merge into it. That is why Madagascar and Australia join the ocean block rather than Africa or Asia.

  Under `mt090` the floor keeps Patagonia (1.43M km2, 0.54) but absorbs Madagascar (0.82M km2, 0.66), which is more distinct by this measure, into Africa. The floor, an area rule, decides that, not the biota.

- **Split tests** (`breaks.split_tests`): each split from k − 1 to k realms against 199 random splits of the same realm into contiguous regions with the same areas and number of pieces (`nullmodel.null_partitions`). Endemism is `nullmodel`'s, scored two ways.
  - **Equal weight** (each side counts the same): every split beats every random one, with p = 0.005, the smallest possible. That holds at every k of every configuration down to k = 30 under `mt090`, including islands of a few thousand km2. A boundary the clustering drew from the data beats random boundaries almost by construction, so significance gives no place to stop. The gain over the null shrinks smoothly under `mt090`: 0.16 to 0.18 for k = 2 to 5, 0.156 at 6, 0.146 at 7, 0.115 at 8, about 0.10 at 9 and 10, then 0.03 to 0.10 for islands.
  - **Area weight** (each side counts by its area): the larger side's share is high under any boundary, so in a lopsided split this score cannot tell a real boundary from a random one (z mostly below 2). It is reported but not read.

**What the data say.** They give a floor and a nested order, not a single k:
- five blocks, about equally distinct, so fewer than five realms is arbitrary;
- then Madagascar and Australia;
- then New Zealand, Patagonia and New Guinea;
- then islands.

The headline at k = 7 is the five blocks plus Australia plus Patagonia, with Madagascar left out by the floor. Nothing about the selection changes here. Whether to report the nested levels instead of one k, and whether the floor should yield to distinctness, is for Markus. Any such rule is to be fixed before M3 opens half B.

## 2026-11-05

### Where the oceans separate (2026-11-05)

Markus asked at what point the oceans get separated. On half A at resolution 7 (release 315557, `native_only`, seed 1, before the layer 4 change below):
- **Under `mt080` to `mt095` they never do.** All of the sea is a single Infomap module of 169 to 173M km², every unit in it a sea part. Merging only ever joins units, so no level of the hierarchy divides it. It joins land at the top (k = 3, beta-sim 0.73). The land in the k = 7 ocean realm (New Guinea, the Bismarcks and the Solomons) is put there by the size floor, not by Infomap.
- **Under `mt075` the sea splits.** Its units are:
  - the Indo-West Pacific, 41.8M km², which merges back at beta-sim 0.734, as distinct as a continent;
  - the Caribbean and tropical West Atlantic, 7.2M km² (0.47);
  - the North Pacific, 4.4M km² (0.42);
  - New Zealand's waters, 2.6M km² (0.55);
  - the tropical East Pacific, 1.9M km² (0.49).

  The cold and temperate rest of the sea is one module with the Holarctic land: 62.6M km², half of it sea. `mt075` never reached a stability of 0.8, so these are hints only.

### Terrestrial unless COL says marine (2026-11-05)

**The leak.** The ocean realm acts as a sink in the merge sequence ("Is there a natural k?"). The cause is land species in the sea parts.
- **At the leaf cut of 0.5 under `mt090`,** the ocean realm's sea parts hold 49,020 species with no environment in COL. These include 13,518 flowering plants, 12,579 insects, 2,359 monocots, 1,996 arachnids, 1,274 reptiles, 812 ferns, 703 mosses and 526 mammals.
- **Of the 10,335 species Australia shares with the ocean realm,** 7,191 are only in its sea parts. Three-quarters of the shared species have no environment.
- **Why layer 4 missed them.** It makes a species land only if COL says exactly terrestrial or 95% of its records are on land. These species have a median of 80% on land, and a quarter have less than a third.

**What does not work.** Markus asked whether to query GBIF twice by shape (marine and terrestrial), or to ignore coastal cells.
- **Shape queries:** cells are already cut into land and sea parts, and a record's side is GBIF's country assignment (no GADM country means sea). Two shape queries would do the same thing, so the geometry is not the problem.
- **Ignoring coastal cells:** that would lose about 3,200 cells at resolution 7, and every island.
- **Inferring the environment from relatives** (the genus, then the family, where at least 95% of the species with an environment agree): this decides only 9,888 of the 49,020.
- **Inferring it from the class:** this fails too, because COL's environment is patchy by source. The only mammals with one are marine.

**Ruling (Markus, 2026-11-05).** Two changes to layer 4 (`clean.Settings`, commit 0ffa897):
- **`land_share` from 0.95 to two thirds.** Clearly more than half on one side is enough. Diatoms or red algae that are mostly at sea stay at sea.
- **Terrestrial unless COL says marine.** A species with no environment in COL is taken as terrestrial, as if COL said so, if it is of a class of `clean.TERRESTRIAL_CLASSES`:
  - insects, springtails, arachnids, centipedes, millipedes;
  - amphibians, reptiles, mammals;
  - the seed plants, ferns, lycophytes, mosses, liverworts and hornworts.

  COL's environment, where it has one, still decides. Freshwater, brackish and all of Aves stay exempt.

**The exceptions.** Markus expected WoRMS's marine flags to cover the marine members of these classes. They are complete in WoRMS but not in release 315557: COL's environment comes from the sector's source, and for these groups that source has none. No whale, seal, sea turtle, seagrass or sea skater is marked marine, and only 4 of 459 Elapidae are. So a short list, `clean.MARINE_EXCEPTIONS`, keeps the rule of their records (two thirds) for:
- **Mammals:** the orders Cetacea and Sirenia; the families Phocidae, Otariidae and Odobenidae; the genera Enhydra, Lontra and Ursus.
- **Reptiles:** the sea turtles (Cheloniidae, Dermochelyidae), the sea snakes (seven genera of Hydrophiinae, and Laticauda), and Amblyrhynchus.
- **Seagrasses:** Zosteraceae, Posidoniaceae, Cymodoceaceae, Halophila, Enhalus and Thalassia.
- **Marine insects and arachnids:** Halobates, Hermatobatidae, Pontomyia, Clunio, Halacaridae, Pontarachnidae and Desis.

Mangroves and salt-marsh plants are left as land. Folding their sea records into the land part of the same cell loses nothing. COL might want to carry WoRMS's flags into the environment of the merged checklist; that is a data question for COL, not for this analysis.

**The calibration** of layer 5 ran on layer 4 as it was. It is sealed and is not run again. `calibrate.DESIGN` now pins layer 4 to what it ran with (95%, no land classes), so its recorded results still describe the code. Setting 132's thresholds now act on layer 4's new output, which differs from the one they were chosen on. How much is in the numbers below.

**The Costello gate** is byte-identical after the change. The cleaning hash changes, so `native_only` is run again from the cleaning at seeds 1 to 5, followed by the choice over seeds and the breaks. The results before the change are kept in `results/explore/res7-before-layer4/`.

### WoRMS's environment (2026-11-05)

**Why.** Markus suggested filling COL's missing environment from WoRMS, mapped to COL ids, instead of the hand-made exceptions of "Terrestrial unless COL says marine". He will look into why release 315557 lacks the marine flags of those groups.

**Inputs:**
- **The match:** a ChecklistBank matching job of WoRMS (dataset 2011) against COL XR 315557, requested by Markus. Job `59fcfcc1-565c-470f-963e-57be8814a05d`, relaxed mode, synonyms on, 1,160,000 rows in 18 minutes.
- **WoRMS's environment:** its extended ColDP export (version 2026-10-01, export of 2026-10-27), which holds `col:environment`.
  - The plain export (the first redirect of `export.zip?format=ColDP`) has no environment column.
- `scripts/worms_environment.py` builds `data/reference/worms_environment.parquet` and enters it in the reference manifest with its sha256. That changes the cleaning hash of every run.

**The mapping:**
- **Which WoRMS usages count:** a WoRMS species usage takes its own environment if accepted, its accepted name's if a synonym.
- **Which matches count:** exact or variant matches to a COL name of rank species, accepted, provisionally accepted or a plain synonym. A synonym stands for its accepted species. Ambiguous synonyms and misapplied names are left out.
- **Several usages per species:** their environments are joined. That happens for 194,270 species, mostly through synonyms.
- **Result:** 470,752 COL species get a WoRMS environment.
  - **Where COL has its own,** the two agree exactly for 423,947 of 433,989 (98%). Most of the rest differ only in COL's "brackish".
  - **Where COL has none,** WoRMS fills 36,763 species: 14,017 marine, 13,798 terrestrial, 5,924 freshwater, and the rest mixed. In the land classes it marks 122 mammals, 87 reptiles, 402 insects and 144 arachnids as marine, the groups the exceptions were for.

**What WoRMS lacks here.** No seagrass is in this WoRMS export. Its *Posidonia* and *Halophila* are fossil molluscs and bryozoans of the same names, and it has no *Zostera* or *Thalassia*. So `clean.MARINE_EXCEPTIONS` stays, as a net for marine species that neither COL nor WoRMS marks. It acts only on species with no environment from either.

**The rule (`clean.Settings.worms_environment`, on):**
- WoRMS's environment stands in for COL's where COL has none; COL's own is never overruled. Layer 4 then applies its rules as before: exactly marine is a sea species, exactly terrestrial a land species, freshwater or brackish exempt, mixed values left to the records.
- Mangroves (*Rhizophora*, *Avicennia*), which WoRMS mostly calls marine, become sea species where COL has no environment.
- `Cleaned.sides` counts the species WoRMS gave an environment as `environment_from_worms`.
- The calibration's `DESIGN` pins the setting off, as it ran.

**The rerun.** The rerun of `native_only` started for "Terrestrial unless COL says marine" was stopped during seed 1, before any select result, and started again with WoRMS's environment. A new manifest would have cleaned again in the middle of it. The Costello gate is byte-identical.

## 2026-11-06

### Species-days (2026-11-06)

**Why.** Our record counts are our evidence in several places: layer 4's land share, layer 5's satellite record share, the under-reporting rates and the group minimum. Many GBIF records repeat a species on one day at one spot: repeat photos, many observers at one rare bird, a seabird colony. On 2.9 billion GBIF records, van Merriënboer et al. ([arXiv:2611.00284](https://arxiv.org/abs/2611.00284), 2026) found that below about 600 m and one day two records are no more likely to be different species. Keeping one record per species, day and S2 level 14 cell there dropped 28.8% of the records. On eBird hotspots, complete checklists grew only as records to the power 0.85. So a count of visits measures effort and evidence better than a count of records.

**The query.** `gbif.build_sql` takes an optional `visit_resolution` and then adds two columns: `n_visits` and `n_visits_open`, the distinct visits among the records counted by `n` and `n_open`.
- **A visit** is a species with one date (year, month and day, a missing part as 0) in one ISEA3H cell of resolution 17, about 0.4 km². Resolution 7 is about 23,000 km² and each step divides by 3; 17 is the step nearest the paper's S2 level 14, about 0.32 km².
- **Undated records** collapse to one visit per species and resolution 17 cell, and records with only a year to one per year there.
- **The query file** is `queries/occurrence_res7_visits.sql`, written and validated by `scripts/write_queries.py`. GBIF's validate endpoint accepts it. Without `visit_resolution` the query is byte-identical to the downloads so far, so `occurrence_res6.sql` and `occurrence_res7.sql` are unchanged.

**The first download failed.** Markus submitted it as 0003444-261104105008071 on 2026-11-06. It failed after three minutes with no reason given. That version built one text per record with `CAST(... AS VARCHAR)`, which GBIF's validator accepts but its engine does not: `DataType "VARCHAR" requires a length parameter` (the error Markus found). `PRINTF('%04d-%02d', "year", "month")` is GBIF's simpler way to build such a text. Resubmitted as 0003457-261104105008071. The query now counts `COUNT(DISTINCT year, month, day, cell)` over four columns, with no casts. The condition sits in the first column, whose null makes the row skipped.

**The download.** 0003457-261104105008071 (DOI [10.15468/dl.ud77we](https://doi.org/10.15468/dl.ud77we)) succeeded on 2026-11-06. It is ingested at `data/gbif/0003457-261104105008071/`, with its QA in `results/visits/res7/qa.json` so that the res 7 QA of 0020744 stays as it was. 28,142,522 rows (seven cell codes remapped, as for 0020744) and 3,280,964,672 clean records: GBIF has changed since 0020744 (28.05M rows), so its counts are not the same. **Visits are 76.2% of clean records** (2,501,637,023), 76.9% of open ones: about a quarter of the records repeat a species on one day in one resolution 17 cell (the paper found 28.8% at 600 m). No row has more visits than records. `occurrences.load_presence` sums `n_visits` and `n_visits_open` over countries when the ingest has them.

**What it is for.** Markus submits the download from his account. Nothing reads the new columns yet: `occurrences.COLUMNS` is exact, so ingesting the download needs the two columns added there first.
- **On half A only, as a sensitivity check:** the cleaning with visits in place of records where records are evidence, and whether the realms move.
- **Not as a change before M3:** the layer 5 thresholds were calibrated on record counts and that calibration is sealed (setting 132). A switch to visits would need a new calibration, so it belongs after M3, and is recorded here as a known limitation: records overstate effort and evidence most in the busiest cells.

### After WoRMS: mt075 at k = 5 over five seeds (2026-11-06)

The rerun of `native_only` with the new layer 4 (seeds 1-5, seven configurations, then `scripts/seeds.py` and `scripts/breaks.py --n-null 199`) finished on 2026-11-06 with no failures.
- **The five-seed rule now picks `mt075` at k = 5:** excess 0.340, within 0.01 of the highest, among 32 maps with median stability of at least 0.8. Before, it picked `mt090` at k = 7 (excess 0.148, stability 0.809), the only map at or above the bar.
- **The same map before and after:** `mt090` at k = 7 now has stability 0.971 (0.966-0.972) and excess 0.334, against 0.809 and 0.148 before. The likely cause is the sea leakage: land species recorded in sea parts made thousands of land species shared with the ocean realm rather than endemic. Checked against the backed-up cleaning (`results/explore/res7-before-layer4/native_only/cleaned/`): species with presences in both land and sea parts fell from 54,971 to 12,369, and land-only species rose from 457,225 to 498,578 (sea-only 76,398 to 77,453). About 42,600 species no longer reach the ocean realm.
- **Merge heights (`mt075`, Infomap's units):** k = 2-5 at 0.80-0.82, then 0.747, 0.710, 0.689, 0.623, 0.561, 0.537, 0.525, 0.525, 0.464, 0.412. All Markov-time configurations agree on the top five (0.79-0.82) and on 0.75 for the sixth.
- **The ocean** stays one realm at every cut. It no longer splits under `mt075`, as it did before ("Where the oceans separate").

### Nested levels: the cuts (2026-11-06, provisional until Markus confirms)

**Candidates.** `scripts/nested.py` on `mt075` (the five-seed choice), seed 1, leaf cut 0.4 and cuts 0.45 to 0.78 in steps of 0.05 (and 0.78), with 20 resamples of 3 trials (`results/explore/res7/native_only/nested-candidates/`, 39 min).
The nested merges (land with land, the ocean apart):

| Realms | Height of the merge that ends the level | What joins there |
|---|---|---|
| 15 → 14 | 0.412 | St Helena (one unit) into the Afrotropics |
| 14 → 13 | 0.464 | the Mascarenes into the Afrotropics |
| 13 → 11 | 0.525 | New Caledonia and Hawaii into Indo-Malaya with Australasia |
| 11 → 10 | 0.537 | Patagonia and Chile into the Neotropics |
| 10 → 9 | 0.561 | New Guinea into Indo-Malaya with Australasia |
| 9 → 8 | 0.623 | New Zealand into the Palearctic |
| 8 → 7 | 0.689 | Madagascar into the Afrotropics |
| 7 → 6 | 0.710 | the Nearctic into the Palearctic |
| 6 → 5 | 0.747 | Australia into Indo-Malaya |
| 5 → 1 | 0.80-0.82 | the five blocks (Holarctic, Neotropics, Afrotropics, Indo-Malaya with Australasia, the ocean), last the land with the sea |

**Recovery says little about the cut.** The median best Jaccard of every realm at every cut is at least 0.85 (most above 0.95), and each is recovered (Jaccard at least 0.5) in 70-100% of the resamples. The misses come mostly from the same resamples at every cut (6 of 20 at the top level, resamples 2, 8, 10, 13, 17 and 18): there Infomap's own modules join Indo-Malaya with the Holarctic, which no cut can undo. So no cut is more stable than its neighbours, and the cuts go by the merge heights.

**Rule.** A level is reported where the realms stay the same over a range of cuts at least 0.05 wide (a gap of at least 0.05 between successive nested merges), cut at the middle of the gap, rounded to 0.01:
- **0.77:** 5 realms
- **0.66:** 8, adding Australia, the Nearctic and Madagascar
- **0.59:** 9, adding New Zealand
- **0.49:** 13, adding New Guinea, Patagonia, Hawaii and New Caledonia
- **leaf cut 0.44:** 14, adding the Mascarenes

The full merge order (`merges.csv`) is the result; the levels are where it is drawn. Levels 6 and 7 (gaps of 0.037 and 0.021) are not drawn. Below 0.41 come single islands and units, which "distinct areas" would keep, but not as realms of their own.

**To look at.** New Zealand joins the Palearctic (0.623), not Australia. The merge compares the union of a realm's taxa, and the Palearctic block holds most genera and families, so a small realm tends to join the largest. That is a property of union beta-sim and is not changed here.

### The other runs and the circularity gate on `mt075` (2026-11-06)

**The runs.** `off`, `introduced`, `native_only-no-extinct-no-no_wild_range`, `native_only-no-griis-no-griis_subnational` and `native_only-merged-names` were cleaned again with the new layer 4 and selected on `mt075` at seed 1, about 18 minutes each (log `res7/mt075-runs.log`). Each selects k = 5. The stability at seed 1:

| Run | Stability |
|---|---|
| `off` | 0.813 |
| `introduced` | 0.874 |
| `native_only-no-extinct-no-no_wild_range` | 0.804 |
| `native_only-no-griis-no-griis_subnational` | 0.957 |
| `native_only` (median over five seeds) | 0.969 |

The `mt090` outputs from before the change are kept in `res7/{circularity,crosscheck,review}-mt090/`.

**The gate passes every mode, and the headline mode stays `native_only`** (`scripts/circularity.py --config mt075`, 99 rotations, `res7/circularity/`). The share of realm boundary along TDWG edges, against `off`'s 0.846, with the 95th percentile of the rotated boundaries:

| Mode | Share | Raise over `off` | Null q95 |
|---|---|---|---|
| `introduced` | 0.881 | +0.034 | +0.057 |
| `native_only` | 0.825 | −0.021 | +0.055 |
| `native_only-no-extinct-no-no_wild_range` | 0.845 | −0.001 | +0.053 |
| `native_only-no-griis-no-griis_subnational` | 0.872 | +0.026 | +0.053 |

The figures are at `land_tolerance` 0.5; with no tolerance the shares are the same and every mode still passes.
- **WCVP's native ranges do not pull the boundaries onto TDWG edges:** `native_only` lies further from them than `off`.
- **The maps of the cleaning modes agree:** V is 0.975-0.990 among them. Against `off` it is 0.86-0.87, with about 41,000-45,000 km of sea boundary moved.

### The cross-checks on `mt075` (2026-11-06)

`scripts/crosscheck.py --config mt075` (`res7/crosscheck/`; the `mt090` ones are in `res7/crosscheck-mt090/`).
- **Beta-sim clustering at k = 5:** Ward (20 shuffles) has V 0.66-0.69 against Infomap's map, with area agreement 0.61-0.63. Average linkage has V 0.63 and area agreement 0.77. Its seeds agree completely (V 1.000), Ward's at 0.95-1.00. On `mt090` at k = 7 the area agreement was 0.59-0.62 (Ward) and 0.72-0.77 (average), so the agreement with an independent method has not changed.
- **Duplicate names merged:**
  - V 0.984 against the map on ids, with the same selected k (5) and the same parts.
  - Excess differs by at most 0.007 at any k from 5 to 12.
  - 972 genus names (2,001 ids) and 238 family names (487 ids) are merged.
  - Repeated names do not shape the realms.

### Nested levels: the cut rules (2026-11-06, Markus's choice)

The gap rule above ("Nested levels: the cuts") rested on a width of 0.05 that nothing in the data set. Markus chose two rules instead, rules 2 and 6 of the alternatives discussed (fixed cuts, a measured gap, the inconsistency coefficient, explained dissimilarity, split significance, endemism against a null, bootstrap recovery). The gap rule's levels are kept in `res7/native_only/nested-gap-rule/`.

- **Rule 2, robust to the sampling of species** (`nested.robust_levels`, `nested.height_replicates`):
  - Infomap's units are merged freely to a base cut of 0.3 (`nested.leaf_realms`), then land with land and sea with sea (`nested.nested_merges`).
  - The merges are replayed, the same groups in the same order, on 199 resamples of the species.
  - A level, the realms after the first *d* merges, is robust where one cut separates the merges made from the merges to come in at least 95% of the resamples.
  - Its cut is the middle of the reference's gap, between the highest merge made and the lowest to come.
  - This measures how much the heights move with the sampling of species, for fixed groups. Infomap's own variation is what the per-realm recovery shows.
- **Rule 6, more endemic than chance** (`nested.realm_excess`):
  - At each robust level, each realm's endemic share (at least 95% of a taxon's cells inside, averaged over species, genera and families as the excess averages them) is set against the same realm in 99 random maps with the same realm areas and pieces (`nullmodel.null_partitions`).
  - A realm passes at p ≤ 0.05, and a level passes when all its realms do.

**Result** (`scripts/nested_cuts.py`, `res7/native_only/nested-cuts/`, `mt075`):
- **Rule 6 passes every realm at every robust level**, each at p = 0.01, the least that 99 null maps can show. At 13 realms the smallest excess is the Nearctic's, 0.19 (0.30 endemic against a null mean of 0.11), then New Guinea's 0.20 and Hawaii's 0.21; the islands' null means are 0.01-0.03, the large realms' 0.07-0.15, and the ocean's 0.35 (0.72 observed). Every realm the merge order produces holds far more endemics than a random region of its size. So rule 6 removes nothing here: rule 2 sets the levels.
- **Rule 2:** the robust levels and their cuts:

| Realms | Cut | Separated in | Added at this level |
|---|---|---|---|
| 5 | 0.773 | 100% | Holarctic, Neotropics, Afrotropics, Indo-Malaya with Australasia, the ocean |
| 6 | 0.728 | 100% | Australia |
| 7 | 0.699 | 97% | the Nearctic |
| 8 | 0.656 | 100% | Madagascar |
| 9 | 0.592 | 100% | New Zealand |
| 10 | 0.549 | 99.5% | New Guinea |
| 13 | 0.494 | 100% | Patagonia and Chile, Hawaii, New Caledonia |

  - **Not robust:** 11 and 12 realms (Patagonia, and Hawaii with New Caledonia, merge within 0.012 of each other: 14% and 31%); 14 (the Mascarenes, 94.5%, just under 95%); 15 (St Helena, 88%). Nothing below 0.376 comes near.
- **Against the gap rule:** the same top and bottom, with three levels more (Australia alone, the Nearctic alone, New Guinea) and the Mascarenes not a level. The finest level is 13 realms.

The levels are drawn with `scripts/nested.py --leaf-cut 0.494 --cuts 0.773 0.728 0.699 0.656 0.592 0.549` and `scripts/nested_maps.py`.

### Provisionally accepted species count as accepted (2026-11-06, Markus's ruling)

Markus: provisionally accepted species count as accepted. "Not ideal as there are conflicting or unresolved names/species, but ignoring them also seems wrong."

This is what the pipeline has done since step 3: `reference.ACCEPTED` holds both statuses, so they are taxa and nodes like accepted species. In 315557 they are 174,878 of 2,529,340 species, with 0.29% of the records directly.

Two known costs, already handled or reported:
- **A provisional species that shares its name with an accepted one** of the same kingdom (7,948 such pairs among WCVP's names, *Zea mays* `5CXX5` and `VJSFD`): the name index sends a source's record on the name to the accepted twin. `clean_report.provisional_twins` lists them.
- **Duplicate genus and family ids, many of them provisional:** the merged-names run shows that they do not shape the realms (V 0.984 on `mt075`).

Removed from the open points.

### Visits against records: the check (2026-11-06)

`native_only` was cleaned and selected on `mt075` (seed 1) from download 0003457 twice, with records and with visits as the evidence (`scripts/explore.py --evidence`, `results/visits/res7/`), and the two cut rules were applied to each (`scripts/nested_cuts.py`, `scripts/nested.py` without the bootstrap). A third column is the headline run on 0020744, so that the drift of GBIF between the two downloads (28.05M and 28.14M rows) shows apart from the evidence.

| | 0020744, records | 0003457, records | 0003457, visits |
|---|---|---|---|
| selected k | 5 | 5 | 5 |
| stability at k = 5 | 0.969 | 0.926 | 0.854 |
| excess at k = 5 | 0.341 | 0.341 | 0.340 |
| V at k = 5 | | 0.975 against 0020744 | 0.869 against records |
| levels by the two rules | 5, 6, 7, 8, 9, 10, 13 | 5, 6, 8, 9, 10, 13 | 5, 6, 7, 8, 9, 12 |

**Layer 5 with visits** (rows, against records on 0003457): `detached_satellite` 386,016 → 332,325, `thin_long_range` 55,451 → 78,554 (more rows fall under its 5), `under_reported` 154,719 → 124,055. Its thresholds were calibrated on records, so the shift is expected and the check stays a check.

**What all three share:**
- the coarsest level: the same five realms (the Holarctic, the Neotropics, the Afrotropics, Indo-Malaya with Australasia, the ocean);
- at the finest, the Palearctic, the Nearctic, the Neotropics, Patagonia, the Afrotropics, Madagascar, Indo-Malaya, Australia, New Zealand, Hawaii, New Caledonia and the ocean.

**What moves:**
- **New Guinea:** its own realm on 0020744; inside an Indo-Malayan realm on 0003457 with records, where Indo-Malaya splits in two instead; with Australia with visits.
- **New Zealand before it splits off:** with the Palearctic with records, with Australia with visits, together with the Antarctic coast and the subantarctic islands (247 parts; at k = 6 of the select stage as well). Its attachment to the Palearctic depends on how records are counted.
- **The top five** merge in another order with visits (South and East Asia joins the Holarctic first at k = 4), as their heights tie at 0.79-0.82.

**Reading:** the coarsest and the finest levels are robust to the evidence, apart from New Guinea; the intermediate attachments of the small realms are not. This supports reporting the coarsest and the finest levels that pass both rules, with the levels between as the merge order. Records stay the evidence, as the calibration was made on them; a switch to visits would need a new calibration, after M3.

### New Zealand: what the records say (2026-11-07, Markus's ruling)

Markus: "for NZ keep what the records tell us." Records stay the evidence, and New Zealand stays with the Palearctic until it splits off (at the 9-realm level on 0020744). Nothing changes.

**Why it is close** (download 0003457, records, the 10-realm level):
- **All ranks together:** New Zealand's beta-sim is 0.621 with the Palearctic and 0.631 with Australia. With visits it is 0.612 and 0.610, and with visits Australia's realm already holds New Guinea, which makes it richer.
- **By rank:**
  - species: 0.741 and 0.746;
  - genera: 0.387 and 0.406;
  - families: 0.110 and 0.154.

  The Palearctic wins at every rank, but only just.
- **Shared species:** of New Zealand's 12,477 species, 1,632 are shared with the Palearctic only, 1,571 with Australia only, 1,600 with both and 7,674 with neither.
- **What the Palearctic-only species are:** they have a median of 4 records in New Zealand, against 27 for the Australia-only ones. Most are microfungi (Sordariomycetes, Dothideomycetes, Leotiomycetes, Pezizomycetes), diatoms, bacteria (Actinomycetes, Cyanophyceae), Agaricomycetes and insects. These are widespread, thinly recorded organisms, often found where European workers sampled.
- **What the Australia-only species are:** insects, lichens, liverworts, mosses, flowering plants, ferns and birds, with far more records each.

So the attachment rests on a near tie that thinly recorded microorganisms tip. It is noted for the review, not acted on. Whether bacteria and diatoms belong in a realm analysis is a separate question, raised on the review page.

### Eukaryotes only (2026-11-07, Markus asked for the variant)

`scripts/explore.py --eukaryotes-only` (`explore.eukaryote_counts`, run `native_only-eukaryotes`, its own hash) keeps the species of COL's Animalia, Chromista, Fungi, Plantae and Protozoa, and leaves out bacteria, archaea, viruses and species with no kingdom: about 6,900 species, 1.2% of species and 0.6% of presences. It was selected on `mt075` at seed 1, and the two cut rules and the levels were run on it, all on 0020744.

| | headline | eukaryotes only |
|---|---|---|
| species after cleaning | 588,400 | 581,164 |
| units | 3,293 | 2,897 |
| selected k | 5 | 5 |
| stability at k = 5 (seed 1) | 0.969 | 0.805 |
| excess at k = 5 | 0.341 | 0.351 |
| levels by the two rules | 5, 6, 7, 8, 9, 10, 13 | 5, 6, 7, 8, 9, 14 |

- **The coarsest level** is the same five realms, and the order of the splits to nine realms is the same: Australia, the Nearctic, Madagascar, New Zealand.
- **New Zealand attaches to Australia** before it splits off, not to the Palearctic. The direct beta-sims on the headline's finest level, all species against eukaryotes only:
  - New Zealand with the Palearctic: 0.623 → 0.646;
  - New Zealand with Australia: 0.631 → 0.640.
- **The finest level, 14 realms,** holds the headline's 13, New Guinea among them, and a fourteenth: temperate East Asia (Japan, Korea and north-east China, 3.0M km², centre 45°N 128°E) splits from Indo-Malaya.
- **Stability at k = 5 is lower** (0.805, one seed), and excess is higher (0.351).

Whether prokaryotes and viruses stay is Markus's call on scope. The per-class signal test ("Which groups carry the realm signal") is the data-driven way to decide which groups to keep.

### Which groups carry the realm signal (2026-11-07)

Markus chose a data-driven test of which groups to keep, with the condition that it "must avoid removing all true, rarely sampled endemics which are good signals". The test (`nested.class_signal`, `scripts/group_signal.py`) works on the headline's finest level, 13 realms, against 199 random maps of the same realm areas and pieces (`nullmodel.null_partitions`).

**The protections for rare groups:**
- Only species in at least two cells are scored. A species in one cell is endemic to any map and would make a thinly recorded class look signal-free.
- A class is judged only if it has at least 30 scored species. A smaller class is kept unjudged: 245 classes, 19,706 species.
- A class is dropped only if it is not above the null (p > 0.05).

**Two designs failed on birds, and why.** On the whole map, with "endemic" meaning at least 95% of a species' cells in one realm, the test dropped birds (Aves 0.367 against a null of 0.373, p 0.59). With a continuous measure, the mean share of a species' cells in its main realm, it dropped them again (0.820 against 0.797, p 0.42), together with the mushrooms (Agaricomycetes, +0.108, p 0.075). Two things distort the whole map:
- **The ocean** is one realm at every level and 58% of the units. A coastal or sea-going species splits its cells between a land realm and the ocean, in every map alike.
- **The yes/no measure** counts range size more than boundaries. A vagrant cell in twenty fails a wide-ranging species in every map, and a small, compact range passes in every map.

**The test as adopted:** the continuous measure (`concordance`), on land only. The land realms' land units are used, and a species' cells are counted there (`totals_within`). The whole-map runs are kept as `group_signal-all-*.csv`. On land the test behaves as it should:
- birds +0.161, mammals +0.160, insects +0.162, flowering plants +0.143, amphibians +0.142, reptiles +0.110, mosses +0.143, mushrooms +0.153 and lichens +0.132, all at p 0.005;
- the microfungi as well: Dothideomycetes +0.129 (p 0.005), Leotiomycetes +0.154 (0.010), Sordariomycetes +0.105 (0.015).

**Result** (`res7/native_only/nested/group_signal-land-concordance.csv`): 95 classes testable, 26 without signal.
- **Prokaryotes: 10 of 12 testable classes have no signal**, which is 4,789 of their 6,078 species: Actinomycetes, Gamma-, Alpha- and Betaproteobacteria, Flavobacteriia, Bacilli, Clostridia, Cytophagia, Chitinophagia and Bacteroidia. Only the Cyanophyceae (+0.147, p 0.005) and one small class pass.
- **Eukaryotes: 16 of 83 have no signal**, which is 6,228 of 556,541 species (1.1%). Every one has a positive excess (0.05-0.18) at p 0.055-0.15:
  - slime moulds, yeasts, Kickxellomycetes, Mucoromycetes and Coniocybomycetes;
  - cryptophytes, golden algae, ciliates (Oligohymenophorea), diatoms (+0.101, p 0.09) and desmids;
  - tapeworms and two classes of thorny-headed worms;
  - rotifers;
  - Andreaeopsida;
  - sharks, which are marine and cannot be judged on land.

  These are small, thinly sampled groups with a weak but positive signal: the case the condition protects.

**Reading:** the test separates prokaryotes, which on land do not follow the realms, from eukaryotes, which do even where weakly. Dropping eukaryote classes on borderline p values would remove the rare true signal Markus wanted kept. Prokaryotes are the one group the data rule out. That points to the eukaryotes-only scope ("Eukaryotes only") rather than to dropping classes. Viruses are too few to test, and they are host records anyway. The decision is Markus's.

### The scope is eukaryotes (2026-11-07, Markus's ruling)

Markus: "yes, drop prokaryotes and start". The analysis keeps the kingdoms Animalia, Chromista, Fungi, Plantae and Protozoa (`explore --eukaryotes-only`, `explore.EUKARYOTES`) and drops COL's Bacteria, Archaea and Viruses. The ground is the class test above: on land, prokaryotes do not follow the realms (10 of 12 testable classes), while every eukaryote class shows a positive excess, even where it is too weak to pass. No eukaryote class is dropped. Rarely sampled groups keep their signal, as Markus required.

The headline run is now `native_only-eukaryotes` and replaces `native_only`. Everything downstream is redone on it, in this order, with the rules already fixed and nothing chosen by hand:
1. the five-seed choice: seeds 1-5 x the seven configurations, then `scripts/seeds.py --run native_only-eukaryotes`;
2. the cut rules (rules 2 and 6, `scripts/nested_cuts.py`) on the chosen configuration;
3. the nested levels with their bootstrap (`scripts/nested.py`, 20 resamples) and maps;
4. the other cleaning runs, each eukaryotes only: off, introduced, the two rule-off variants and merged names;
5. the circularity gate on them, then the cross-checks and the review maps.

The `native_only` outputs stay where they are, as the record of what the scope changed. The provisional eukaryotes-only run on mt075 at seed 1, whose `nested/` and `nested-cuts/` were drawn before the choice, moves to `nested-provisional/` and `nested-cuts-provisional/`.

### The reporting rule (2026-11-07, Markus's ruling)

Markus confirmed it. **The report gives the coarsest and the finest nested levels that pass both cut rules** (rule 2, robust to species resampling, and rule 6, every realm above its null; "Nested levels: the cut rules"). The levels between them are given as the nested merge order, not as further results.

The rule is what gets preregistered, not its numbers. On half A's `native_only` run it gave 5 and 13 realms, but those numbers are not carried forward:
- **Eukaryotes headline:** for `native_only-eukaryotes` (above), the rule's levels are whatever the rerun gives.
- **Half B:** in the confirmation at M3, its levels are whatever the rule gives there. Fixing "5 and 13" in advance would make the confirmation a comparison with what looked right on half A.

The rule adds no threshold of its own: the two cut rules already define the levels. It goes into `PREREGISTRATION.md` as written here, once that file's contents are settled. That is still an open point.

### The eukaryotes headline: results (2026-11-08)

This is the chain from "The scope is eukaryotes", run in full with no hand choices (log `res7/eukaryotes-headline.log`, 19:55 to 08:32). One step failed and was rerun: the gate's `--suffix -eukaryotes` was read by argparse as an option, so the gate ran again with `--suffix=-eukaryotes`.

**The five-seed choice** (`native_only-eukaryotes/seeds/`) is again mt075 at k = 5: excess 0.355 (was 0.340), with 35 stable maps (was 32). Median stability over the seeds is 0.961 (range 0.805-0.969).

**The cut rules** (`nested-cuts/`):
- Rule 2 finds robust levels at 5, 6, 7, 8, 9 and 14 realms, with cuts 0.796, 0.756, 0.709, 0.664, 0.596 and 0.484. Without the prokaryotes, `native_only`'s levels at 10 and 13 are no longer robust, and 14 is.
- Rule 6 passes every realm at every level (p = 0.01), with the smallest excess 0.136.
- **The reporting rule therefore gives 5 and 14 realms.**

**The nesting** (`nested/`, 20 Infomap resamples; `nested/maps/`):
- **5 realms:** Holarctic; Afrotropics with Madagascar; Neotropics with Patagonia; Indo-Malaya with Australasia; the ocean.
- **6:** Australia with New Zealand splits from Indo-Malaya.
- **7:** the Nearctic splits from the Palearctic.
- **8:** Madagascar.
- **9:** New Zealand.
- **14:** New Guinea, Patagonia, Hawaii, New Caledonia, and **temperate East Asia**: Manchuria, the Amur, Korea and Japan, 3.0M km², split from Indo-Malaya.

The order matches `native_only` except in two places:
- New Zealand now sits with Australia rather than the Palearctic before it splits off. This is the attachment the visits and eukaryotes-only checks had found.
- East Asia is new. In `native_only` it lay inside Indo-Malaya at every level.

**Recovery.** Every realm is recovered in the Infomap resamples at 75-100%, with median Jaccard 0.80-1.00, with two exceptions:
- **East Asia** is recovered in 5% of resamples (median Jaccard 0.17). In nearly every resample it is absorbed into a larger realm.
- **The Palearctic** falls to 55% (median Jaccard 0.63), from 80%. In those resamples it either takes in East Asia or is joined with Indo-Malaya, as before.

The cut rules do not use this recovery. Rule 2 replays the merges of fixed base realms on species resamples, and East Asia's merge height has a wide gap there. Recovery reruns Infomap itself, which moves the modules under it. Whether a realm needs both is a question for Markus (below).

**The gate** (`res7/circularity/`, all runs eukaryotes only, mt075, `--suffix=-eukaryotes`): every run selects k = 5, and every mode passes. The headline mode stays `native_only`.

| Run | Share of boundary along TDWG edges | Change from `off` | 95th percentile of rotations |
|---|---|---|---|
| `off` | 0.801 | | |
| `native_only` | 0.813 | +0.012 | 0.126 |
| `introduced` | 0.892 | +0.091 | 0.129; 0.102 without the land tolerance, the closest margin |
| no extinct / no-wild-range rules | 0.871 | +0.069 | |
| no GRIIS | 0.752 | -0.049 | |

**The cross-checks** (`res7/crosscheck/`):
- **Beta-sim area agreement with Infomap at k = 5:**
  - Ward 0.687, 0.842 and 0.687 over three seeds (was 0.61-0.63); seed 2 differs, with V 0.87 between seeds;
  - average linkage 0.769, as before.
- **Merged names:** V 0.983, and the cleaned presence and the units are identical.

**Outputs.**
- The review maps are in `res7/review/`.
- The `native_only` versions of the gate, cross-check and review outputs are in `res7/*-native_only/`.

**East Asia** (settled by "The reporting rule: a third level", 2026-11-08). It passes both preregistered rules but is not reproduced when Infomap reruns on resampled species. My recommendation is to keep the rules as they are and report 5 and 14, with each realm's recovery beside it and East Asia marked as the one realm the resamples do not reproduce. Two alternatives were considered:
- **A recovery bar on whole levels.** Its only effect here would be to drop level 14, so the finest level would be 9. That loses New Guinea, Patagonia, Hawaii and New Caledonia, which recover at 85-100%.
- **A bar on single realms.** It would fold East Asia back into Indo-Malaya and report a 13-realm level that no rule defines.

Either alternative would be a rule chosen after seeing which realm fails.

### The reporting rule: a third level (2026-11-08, Markus's ruling)

This entry replaces "The reporting rule" of 2026-11-07. Markus asked to report 5, 9 and 14 on the eukaryotes headline. The 9 comes from a rule, so the rule is what he confirmed and what gets preregistered.

**The report gives three of the levels that pass both cut rules** (`nested.reported_levels`, written as `reported` to `nested-cuts/cuts.json`):
1. **The coarsest.**
2. **The end of the unbroken run from the coarsest.** Each level in the run has exactly one realm more than the one before it, so every split down to it is robust on its own: the hierarchy is resolved one split at a time.
3. **The finest.** Each of its realms is given with its recovery in the Infomap resamples (`nested/realms.csv`).

The levels between them are given as the merge order. When the run reaches the finest level, two levels are reported; when only one level passes, one is.

**The rule on half A:**
- **`native_only-eukaryotes`** has robust levels 5, 6, 7, 8, 9 and 14, so the report gives **5, 9 and 14**. From 9 to 14, five realms split at once: New Guinea, Patagonia, Hawaii, New Caledonia and temperate East Asia. No level between 9 and 14 is robust, so their order is not resolved.
- **The earlier `native_only`** run has robust levels 5-10 and 13, so the rule would give 5, 10 and 13.
- **Half B:** at M3 the levels are whatever the rule gives there.

**Why not a recovery bar.** Recovery does not single out 9, and a bar would need a threshold chosen to land on it. Without the prokaryotes, the Palearctic (the Holarctic at 5 and 6) is weakly recovered at every level: 45-55% of resamples, median Jaccard 0.41-0.63, where `native_only` had 70-80%. In those resamples Infomap joins it with Indo-Malaya. The entry "The eukaryotes headline: results" gave the Palearctic's 55% for the finest level only; it holds at every level.

**East Asia** stays in the 14-realm level, which the rule reports. It is shown as the one realm the resamples rarely reproduce: 5%, median Jaccard 0.17.

### Why the Palearctic is hard to reproduce (2026-11-09)

Markus asked why the Palearctic's recovery is low (45-55% of resamples at every level, see "The eukaryotes headline: results"). The diagnostic is in the session scratchpad as `palearctic_boots.py`; its outputs are in `res7/native_only-eukaryotes/nested-diagnostics/`. It reruns the 20 resamples of `nested.bootstrap_recovery` with the same draws and the same Infomap seeds, and keeps each resample's maps and Infomap modules. It reproduces the recovery in `nested/recovery.csv` exactly.

**The Palearctic is never broken up.** In every resample, 95-100% of its area stays in one realm. When it is "not recovered", it has been joined to something else:
- **the ocean** in 7 of the 20 resamples;
- **the Afrotropics or the Nearctic** in 4.

The failures are all-or-nothing per resample: Jaccard 0.92-0.99 or 0.16-0.38, the same at every level.

**Infomap makes the join, not the merges after it.** In the failing resamples, one Infomap module already holds a large share of the Palearctic together with ocean, Nearctic or Afrotropical units. The ocean in that module is not the seas next to Eurasia. It is open-ocean cells from every basin and latitude, down to south of 30°S. These resamples also have many more modules: 3,200-5,400, against 2,200-2,900 in the resamples that recover the Palearctic. That is the mark of a poor local optimum: a catch-all module plus fragments.

**It is Infomap's search, not the data.** The recovery resamples run Infomap with 3 trials; the reference map is made with 10. The 11 failing resamples were rerun with 10 trials (`replicates-t10.parquet`):
- at 5 realms, the Palearctic comes back in all 11 (Jaccard 0.97-0.99);
- at 14 realms, it comes back in all 11 as well (0.69-0.98);
- the module counts fall to 2,100-3,000.

**Temperate East Asia stays weak with 10 trials.** It is recovered in 2 of the 11 resamples (Jaccard 0.92 and 0.53, the rest 0.15-0.28). Its weakness is in the data.

**Reading.** The low recovery of the Palearctic, and the misses in earlier runs (Infomap joining Indo-Malaya with the Holarctic), measure how often 3 trials of Infomap miss on a resample, not how stable the realm is. Without the prokaryotes this happens more often: 9 of 20 resamples against about 4 of 20 with all species. The reason is not established. The recovery check would measure the data if its resamples ran Infomap with the reference's 10 trials (`scripts/nested.py --boot-trials 10`). The cut rules and the selection do not use this recovery.

The selection's own bootstrap stability (20 bootstraps of 3 trials, the design's) may be lowered by the same misses. It stays as it is, because changing it changes the selection, the Costello gate and the calibrated design.

### The recovery check runs Infomap with 10 trials (2026-11-09, Markus's ruling)

Following "Why the Palearctic is hard to reproduce", each resample of the recovery check (`nested.bootstrap_recovery`) now runs Infomap with 10 trials, as many as the reference map has; it used to run 3. `scripts/nested.py --boot-trials` defaults to 10. Recovery then measures how stable a realm is on the data, not how often a short search misses.

Only the recovery column (`nested/realms.csv`, `recovery.csv`) changes. The levels, the cut rules, the selection and its 3-trial bootstrap stability stay as they were. The 3-trial outputs are kept in `native_only-eukaryotes/nested-boot3/`.

**Result (2026-11-09).** `scripts/nested.py` was rerun with the same cuts and 20 resamples of 10 trials. The levels, realms and areas are identical; only the recovery column changed.
- **Recovery:** every realm at every level is recovered in 85-100% of the resamples, with median Jaccard 0.68-1.00. The Palearctic is recovered in 100% at every level (Jaccard 0.97-0.98, up from 0.41-0.63). The Holarctic at 5 and 6 realms is at 100% (0.98, up from 0.53).
- **Temperate East Asia** is the one exception: 15% (Jaccard 0.17, up from 5%). It is the realm the report flags.
- **The ocean** at 14 realms has median Jaccard 0.70, down from 0.81, and is recovered in 95% of resamples.

The 3-trial outputs are in `nested-boot3/`. The review page and `nested/maps/tables.md` use the new recovery.

### The reported levels, the Antarctic coast and the preregistration (2026-11-09, Markus's rulings)

Markus answered the review page's three questions with "yes to all 3":
- **The reported levels** of the eukaryotes headline, 5, 9 and 14 realms, read right to him. Temperate East Asia is flagged by its recovery.
- **The Antarctic coast and Greenland stay Palearctic**, as the data and the imputation give them. Their imputed parts are hatched on the maps and flagged in the per-cell outputs.
- **`PREREGISTRATION.md` is drafted from DECISIONS** for Markus to edit. The pass marks for the confirmation on half B, which the design leaves to the freeze, are drafted as proposals for him to settle.

### Endemism by rank on its own (2026-11-09)

Rule 6 and the endemism excess of the configuration choice average species, genus and family equally. The proposed half B pass mark C1 asks each rank to pass separately, which had not been checked on half A. `scripts/rank_endemism.py` splits rule 6 by rank, at every level of `nested-cuts/cuts.json`. It rebuilds the levels and draws the null maps as `scripts/nested_cuts.py` does (99, seed 1), and checks that the rank average reproduces `nested-cuts/excess.csv` before writing anything; it does at all six levels. It changes no method. Outputs are in `native_only-eukaryotes/rank-endemism/`: `realms.csv`, `levels.csv`, `report.md`, and one map per reported level with a panel per rank (`level1.png`, `level5.png`, `level6.png`). The matrix has 580,158 species, 71,484 genera and 5,586 families.

**Per level.** At every level and rank the area-weighted share is above all 99 null maps (p = 0.01). At 14 realms: species 0.820 against a null mean of 0.436, genus 0.589 against 0.234, family 0.390 against 0.084. The excess falls only a little from species to family (0.38, 0.36, 0.31).

**The family score of a level is mostly the ocean.** Sea realms are 58% of the area. The ocean's family share is 61% (null 13%); the land realms' are 0.1-12%. Over land realms alone, weighted by area (the null as the area-weighted null means, so without a p):
- at 5 realms: species 0.84 against 0.33, genus 0.51 against 0.12, family 0.14 against 0.03;
- at 14 realms: species 0.76 against 0.25, genus 0.37 against 0.08, family 0.08 against 0.02.

The land signal at family rank is about five times the null, but small in absolute terms. A level-wide test at family rank would pass on the ocean alone.

**Per realm.** Every realm passes at every rank at every level (p = 0.01), with one exception: **New Guinea at family rank at 14 realms**, 1 endemic family of 746 (0.1%, null 0.1%, p = 0.35). The small realms pass at family rank on a handful of families against a null near zero: Hawaii 2 of 520, New Caledonia 4 of 442, Madagascar 12 of 692, Patagonia 13 of 855, New Zealand 15 of 1,039, East Asia 20 of 1,476. Counts that small can fall to zero or one on another sample.

**What this means for C1.** As proposed, C1 is met on half A at the level of whole maps, at each rank. But the family part of that is carried by the ocean. Per realm, the strict per-rank reading would already fail New Guinea on A, and on B would be at the mercy of a few families in each small realm. Markus decides the pass mark (PREREGISTRATION.md, section 9).

### Realms from one rank on its own (2026-11-09)

Markus asked for realms built from each rank on its own before the freeze: the robustness run "species only, species plus genus, and family only" of the design, here species, genus and family. `scripts/explore.py --rank R` makes a run `<run>-rank-<R>` (`explore.run_name`). It copies the cleaning and the units of `<run>` (`explore.copy_base_stages`); the rank enters neither stage's hash, so the copies are current, and the run never writes into its base run.

**The rank's taxa are the units of the analysis** (`explore.rank_alone`). The matrix holds that rank's columns only, marked as the species rank. Every step then treats them as a run of every rank treats species:
- the network weighs their links 1, not by `network.information_weight`, which is for genera and families beside species;
- the pieces and merges compare them;
- the bootstraps of the selection, the rule 2 replicates and the recovery resample them;
- the endemism excess, the configuration choice and rule 6 score them alone.

So a family run is what the pipeline makes of a world that has only families. The rules are the headline's and nothing is tuned. `network_matrix` takes the rank from the run's name, so the seeds, nested-cuts and nested scripts need no change, and the rank is in the configuration's hash (`config_inputs`, `read_selected`).

**First pass, seed 1 only.** The headline's five-seed choice took about 10 hours. Here each rank runs the seven configurations at seed 1. The choice is `scripts/seeds.py --seeds 1`, followed by the cut rules, the nested levels with a 10-trial recovery, and the maps. Five seeds can follow if a rank's map is worth it. The runs are compared with the headline by area-weighted V-measure at each reported realm count.

### Realms from one rank on its own: results (2026-11-09)

Half A, release 315557, download 0020744-261023074918325, seed 1 (log `results/explore/res7/rank-runs.log`).

**Family alone: no map.** All seven configurations are degenerate. default and regularized give one module. mt075 to mt095 give two big modules, about land (4,366 units) and sea (3,222), and 3,540 single-unit modules. With every family link weighted 1, the widespread families tie the whole land together.

**Genus alone: no stable map.** Six configurations are degenerate (2 to 4 units above the size floor). mt075 gives 5 realms but falls short: its stability is 0.69, below 0.8, beside 1,590 single-unit modules. The seeds script stops there, so there are no cut rules and no maps. Nothing was tuned to rescue either rank: a lower Markov time or other link weights would change `network.py` or `matrix.py`, which stay as they are.

**Species alone: close to the headline.**
- The choice is mt085 at k = 7 (excess 0.408). The headline's choice is mt075 at k = 5.
- Rules 2 and 6 pass at 3, 6, 7, 8 and 9 realms. The reporting rule gives 3 and 9: the run from 6 to 9 ends at the finest level, so there is no third.
- The 10-trial recovery finds every realm at every level in all 20 resamples.
- At 3 realms: the Afrotropics with Madagascar, the rest of the land, and the sea.
- At 9 realms, the headline's land realms are there with two differences:
  - New Zealand stays with Australia;
  - the sea splits in two, the Atlantic with the Arctic against the Indo-Pacific with the Southern Ocean.

`scripts/compare_runs.py` cuts each run's merge order at a realm count and compares it with the headline's map at that count, over the 11,130 units both map. It reports the V-measure over units and the area agreement (`crosscheck.area_agreement`). The headline against itself gives 1.000 at every count.

| realms | why | V-measure | area agreement |
|---|---|---|---|
| 3 | species-alone level | 0.944 | 0.997 |
| 5 | headline level | 0.837 | 0.901 |
| 9 | both | 0.854 | 0.794 |
| 14 | headline level | 0.783 | 0.578 |

At 3 realms the headline is cut from the same merge order, not from one of its own levels. At 14 the species-alone map is not one of its own levels: no level finer than 9 passes the rules.

**What this says.** The realms down to 9 come from species. Genera and families add little to where the boundaries fall, and they cannot make realms on their own at this resolution under these rules. In the headline they are weighted by information and sit beside the species. There they decide two things:
- New Zealand separates from Australia at 9;
- the sea stays one realm.

Both rest on genera and families, which matches the per-rank endemism (`Endemism by rank on its own`), where the ocean carries most of the family signal. The headline method stays as it is. Five seeds for the species-alone run would only change the result if its configuration choice moved; that is for Markus to ask for.

### Genus and family presence stays one species (2026-11-09)

A genus or family is present in a cell when one of its species is (`matrix._assemble`, an OR over the species). Markus asked whether a threshold that grows with the taxon's size should replace it, so that a single species of a family of hundreds counts for little.

On half A at the headline's 9-realm map:
- one-species presences are the rule: 68% of a genus's presences in a cell and 45% of a family's rest on one species, and 19% even in families of 1,000 species or more;
- a threshold of 5% of the family's species raises the share of families endemic to one realm from 0.49 to 0.57, and for families of 100 or more species from 0.21 to 0.55, by shrinking their median range from 953 cells to 29;
- a threshold of two species drops 55% of the genera.

Such a threshold measures sampling effort more than range. It shrinks large families to their best-sampled cells, and the threshold that gives them the most endemism would be chosen against the C1 pass mark. In the headline, genera and families already count for little: their links are weighted by information and they are about 1% of the columns of the merges. The realms come from species (`Realms from one rank on its own: results`). Markus: leave it. The roll-up stays as it is.

### The half B pass marks and the scaling of `group_min_records` (2026-11-09)

Markus settled the last two method points of the preregistration (PREREGISTRATION.md, sections 4 and 9).

**C1** passes when, at every reported level, half B's area-weighted endemism excess is above the null (p ≤ 0.05) at each of the three ranks, and the family rank is also above the null over A's land realms alone. On half A the sea realm carries most of a level's family score: at 9 realms the family share is 0.392 on the whole map but 0.083 over land. Without the land test the ocean alone could pass the families.

The land test is now scope `land` in `scripts/rank_endemism.py`. It weighs the land realms only (`nested.sea_realms`), and each null map's realm at the same position stands in for a land realm. Scope `all` is checked against `endemism.area_weighted`.

On half A, every rank passes at every level on both scopes, p = 0.01 (the smallest p 99 null maps allow):

| realms | family, whole map (null q95) | family, land (null q95) |
|---|---|---|
| 5 | 0.417 (0.142) | 0.141 (0.041) |
| 9 | 0.392 (0.124) | 0.083 (0.023) |
| 14 | 0.390 (0.114) | 0.077 (0.021) |

Per realm, the only failure is New Guinea at family rank at 14 realms. A realm that fails is reported as not confirmed and does not fail the test. C2 stays as proposed: V ≥ 0.8 at the coarsest reported level.

**`group_min_records`** scales with the records (open point I4, 2026-11-01). The calibration chose 300 on half A, which holds about half the records of all species. A run uses 300 times its records over half A's, which keeps the meaning the threshold had on A: about 600 on all species. The exact values, and the code that sets them, come at the freeze and in M2; half A's runs keep 300.

### The move to COL 316691, with a new download (2026-11-11)

GBIF production switched on 2026-11-10 from COL 315557 (COL26.6 XR) to **316691 (COL26.10 XR, issued 2026-10-16)**. Its match index was built on 2026-11-05. Markus: every COL release has issues, but the one GBIF used before had many duplicates and known problems, so the analysis moves to 316691 with a new full download, before the freeze rather than after it. The cleaning waits for the GRIIS re-fetch anyway (Markus is fixing ChecklistBank's interpretation of GRIIS-US's native flags), so both changes go into one reprocessing.

**What is new:**
- **The download.** 0001502-261110092616502, submitted by Markus on 2026-11-11 with `queries/occurrence_res7_visits.sql`. That query is the res 7 query plus `n_visits` and `n_visits_open`, so one download carries both records and visits.
- **The COL export.** 316691's extended export with classification is ChecklistBank job 3a5a79df-3a08-4877-a1c7-5d7ec869699d, requested by Markus with the settings of 315557's export. Its classification-only export, fd5642b7, already existed.
- **The WoRMS environment.** It comes from a new matching job of WoRMS (dataset 2011) against 316691, 6cb549e1-b347-4175-a7bc-f2e3661dbbb1, with the settings of 59fcfcc1 (relaxed, TSV, synonyms).
- **`fetch_reference.py --col-only`.** It fetches only the two COL exports and replaces just their manifest entries, so the COL release can change before the GRIIS lists are fetched again. The reference data of 315557, which every M1 result so far rests on, is kept in `data/reference-315557/`.

**What stays:**
- the calibration, setting 132; `calibrate.py` is not run again, and its held-out check on half B stays a record of 315557;
- every rule chosen on half A: the configuration, levels and cuts are re-derived by the same rules on the new data.

**The split is carried over by lineage (Markus, 2026-11-11).** Half A and half B were drawn from 315557's lineages: family names, or the next rank up. A seeded shuffle and a greedy balance over that whole list placed them. Run fresh on 316691, `write_split.py` would shuffle and balance a different list, so many families explored in half A would land in half B, and M3's test would no longer be independent of the exploration.

Instead, `split.carry_over` keeps every lineage of 315557 in its half. Only lineages new in 316691 are placed as before: within their group, largest first, into the half with fewer species, counting the carried-over lineages from the start. A species follows its family of 316691, so a species that changed family may change half. `split.crossed` lists those species, and `scripts/carry_split.py` reports them and writes them to `data/split/crossed-315557-<release>.parquet`. The carry-over reads taxonomy only.

On 315557 itself, `carry_over` reproduces the frozen split (hash 88f237b4…), with no species crossing. The refactored `assign_halves` also still gives that hash.

**The release review of 315557.** Markus shared ChecklistBank's review of COL26.6 XR (`download.checklistbank.org/releases/3/598/review.html`), dated 2026-11-11, with the verdict "do not publish". Its flags, and what they meant for M1:
- **Duplicate higher taxa: the one that touches us.** 22,848 duplicate provisionally accepted genera and 2,273 such families were injected by six merge sectors (GRIN, Dyntaxa, Artsnavnebasen, MadeiraSelvagens, SNSB 33082 and 51518). The duplicate-name issue went from 2,613 to 50,293, 26,672 of them species. Duplicate genus and family ids split a taxon's presences in the network's genus and family links. The merged-names run already measured this on half A: 1,030 of 73,398 genus ids and 249 of 5,964 family ids merge away, and the map stays at V 0.983 against the headline. Genera and families count for little in the realms, which come from species ("Realms from one rank on its own: results"). The split is unaffected, as it goes by family name.
- **Distributions of WoRMS and FishBase: not used.** All 493,624 were lost, entering through base release 315448. Our natives come from the WCVP, WRiMS and GRIIS exports, and the environment comes from COL's name usages and WoRMS's own export.
- **BOLD BIN and UNITE identifiers: not used.** They are unranked usages, and we use accepted species only.
- **Identifier churn at family and genus rank: matters for the carry-over.** In 315557, 79 accepted families and 1,371 genera changed ids, and 1,070 family ids came back. Carrying the split over by family name, not by id, makes it immune to this churn. Species ids churn too: 16,548 accepted species ids were deleted against 315192. A species that changed both its id and its family is counted by `crossed` as gone and new, not as crossed, so `carry_split.py` also reports the species matched by name.
- **GBIF impact:** no impact report existed for 315557, which download 0020744 was matched against.

**The split of 316691** (`scripts/carry_split.py`, `data/split/halves-316691.parquet`, canonical hash `73c940c3…`, kept as `split.SPLIT_SHA256_316691`; the pipeline's default stays 315557 until the reprocessing starts):
- 316691 has 2,504,814 accepted species in 27,471 lineages: 1,246,127 in half A and 1,258,687 in half B, balanced within every large group to between 0.48 and 0.52. 2,440,692 species are in both releases by id; 64,122 are new and 88,648 of 315557's are gone.
- **13,076 species change half (0.54% of those in both),** and 447 of the 23,410 that took a new id, paired by name. Every one of them moved into a family that 315557 already held in the other half, so whole-family halves can't avoid it. Examples: 1,019 Ectobiidae to Pseudophyllodromiidae, 729 to Chrysididae, 205 Margarodidae to Monophlebidae. The 729 are the cuckoo wasps *Chrysis*: 315557 had them under Vespidae (Masarinae), an error of that release that 316691 corrects.
- **A refinement of the carry-over (my ruling within Markus's approval, recorded here).** A lineage new in 316691 takes the half that most of its species had in 315557, when they were there and not tied. That covers a family split off or renamed (Monophlebidae from Margarodidae, Euglenaceae from Euglenidae) or moved to another kingdom (the Microsporidia families, Protozoa to Fungi, which renames their lineage). Only the rest are placed by balance. Without it, 15,221 species changed half; with it, 13,076. The approval's purpose was to keep the explored and held-out families apart, and this does so more closely. On 315557 itself, `carry_over` still reproduces the frozen split.

**The WoRMS environment of 316691.** Rebuilt from match job 6cb549e1: **547,068 COL species** get one, against 470,752 on 315557. Of WoRMS's 1,419,220 match rows, 1,055,073 are usable.

**The release review of 316691** (Markus shared it, dated 2026-10-19; verdict "publish with caveats"):
- The duplicate-name issue fell from 65,111 in COL26.9 XR to 1,367. No accepted name is duplicated at any rank down to family, and provisionally accepted usages fell by 2,354, among them GRIN's 3,031 placeholder names such as *Abroma spp*. That is the defect of 315557 that touched us.
- The caveats are far from our data:
  - 129 taxa have a synonym parent;
  - the Coleoptera type-genera sector was deleted, losing 686 family-group synonyms (the accepted classification stays);
  - 8,561 merge-source genus subtrees were skipped as homonyms, most of them already in COL from another source;
  - WoRMS's subphylum Endomyxa and class Nucleohelea were lost;
  - the IRMNG Protozoa synonyms were removed.
- Its GBIF impact reports add up four monthly steps from 315557. In the alltaxa report, 791 rows holding 6.8M occurrences change family, most of them protist reclassifications. Those changes reach our data through the new download.

**The occurrence status vocabulary changed, which emptied the first download.** Download 0001502-261110092616502 (DOI 10.15468/dl.x64e8g) succeeded with a header and no rows. GBIF's OccurrenceStatus vocabulary deprecated `Present` (concept 40) and `Absent` (41) on 2026-11-05 and replaced them with `Detected` (1750) and `NotDetected` (1749). Our query kept `occurrencestatus = 'PRESENT'`, which no record matches since. Two one-row diagnostics, which Markus submitted, found this:
- **0001523** (`queries/check_table.sql`): the table holds 3,892,311,392 records, and every filter of the query matches as before except the status, which matches 0.
- **0001583** (`queries/check_values.sql`): `occurrencestatus` holds `DETECTED` (3,754,511,859), `NOT_DETECTED` (137,799,513) and 29 empty values. TAXON_MATCH_FUZZY counts the same in `issue` and in the new `taxonomicissue` (8,996,859 among DETECTED), so the `issue` filter is unchanged. `taxonomicstatus`, `establishmentmeans.concept` and `degreeofestablishment.concept` keep their values.

`gbif.build_sql` takes `status`, `DETECTED` by default, and `status="PRESENT"` rebuilds the queries of 0020744 and 0003457 byte for byte (`test_sql_of_the_downloads_before_the_occurrence_status_vocabulary_changed`). The new download's query is `queries/occurrence_res7_visits_detected.sql`, which differs from 0003457's only in that word.

The new concepts' hidden labels, the verbatim strings that map to them, differ from the old ones. None of Absent's 9 (`0`, `ABSENT`, `Ausente`, `not observed` …) is on NotDetected, which has 2 labels, one of them garbled (`Non observÃ©`). Detected dropped 194 of Present's 204 and gained about 950 numbers and abundance classes, a few of them garbled too.

As far as the API shows, absences still reach NOT_DETECTED (137.8M, which leaves 3.65M records with an unparsable status at the default, DETECTED). Those come from shifted columns: kingdom names, UUIDs, dates and red-list codes in the status field, mostly Hatikka.fi (1.83M) and the Kalkalpen database (1.76M). No sign of absences read as presences was found. Markus has been told about the labels.

**The GRIIS lists after ChecklistBank's reimport (2026-11-11).** Markus reimported every GRIIS list with ChecklistBank's fix of how it reads their establishment means, and requested an extended export of each list's new import (`request_exports.py --griis-current`). `fetch_reference.py` took all 305 lists as extended exports of their current import, none missing. The old files are kept in `data/griis-before-reimport/` and `data/reference-315557/`. COL, WCVP and WRiMS were unchanged, and so was the WoRMS environment, which `worms_environment.py` wrote again byte for byte (the full fetch rewrites the manifest without it).

The old import read almost every GRIIS term as `native`: `Alien`, `Native|Alien`, `Cryptogenic|Uncertain` and the US registers' `introduced: assisted colonization`. That is why GRIIS-EC and GRIIS-SR applied nothing (the "What the lists say" entry above). The new import maps the sources' terms like this (checked against the source archives of AU, US-United, EC, JP, GB and SR):

| GRIIS source term | old import | new import |
|---|---|---|
| `Alien`, `introduced` | `introduced`, or `native` in some lists (all of EC and SR) | `introduced` |
| `Cryptogenic\|Uncertain` | `native` | `uncertain` |
| `introduced: assisted colonization` (US registers) | `native` | `introduced assisted colonisation` |
| `Native\|Alien`, `native\|introduced` | `native` | empty |

Species entries, by `griis_lists`: before, 127,338 introduced, 51 unstated, 5,077 native and 13 uncertain; after, 129,203 introduced, 477 unstated, 1 native and 2,894 uncertain. Of the 5,077 native entries, 2,752 are now uncertain, 1,767 introduced and 412 empty. *Passer domesticus* and *Columba livia* in GRIIS-US-United, the cases that started this, went from native to introduced.

- **`introduced assisted colonisation` counts as introduced (my ruling).** It is a kind of `introduced` in Darwin Core's establishmentMeans vocabulary. In the US registers it marks species released on purpose: 390 of the 431 rows in US-United are biocontrol insects, the rest include *Passer*, *Columba* and *Trapa natans*. `natives.GRIIS_MEANS_KIND` reads it as `introduced`. Any other new value still raises.
- **An empty value is left out, as `uncertain` is (Markus, 2026-11-11).** This replaces the rule of "GRIIS by its own establishment means" (2026-11-01), under which an empty value counted as introduced. That rule was made when only GRIIS-TR's 28 entries were empty. Now the 412 empty entries that were native before are GRIIS's `Native|Alien`, native in part of the country and alien in another, which ChecklistBank leaves empty. As aliens they would have taken native species out of their own country: *Trichosurus vulpecula*, *Melaleuca lanceolata* and *Eucalyptus cinerea* in Australia; *Oncorhynchus nerka* and *Popillia japonica* in Japan. The rule change was chosen over a ChecklistBank mapping, so the work goes on. The cost: the other 65 empty entries are left out too, 36 of them GRIIS-TR's. `natives.GRIIS_ALIEN` is `introduced` alone, and an introduced row beside an empty one still makes an entry alien. `griis_lists` still counts the empty ones as `unstated`.
- **`check_natives.py` on 316691 with the new lists** (`results/natives.json`; 315557's report is kept in `data/reference-315557/natives.json`): 116,792 GRIIS alien pairs of 27,063 species in 305 lists, against 115,371 of 26,734 in 303. Ecuador (690 pairs) and Suriname (172) apply for the first time. Names matched: GRIIS 92.3% after GBIF's fallback (92.0% before), WCVP 98.5% and WRiMS 98.8%; GBIF's match API was asked about 26 names (103 before) and placed 6. The same 62 subdivisions have no cells, KY-Brac in place of KY-CaymanBrac. Checked by hand: *Passer domesticus* is an alien of US, AU and EC, *Columba livia* of US, AU, EC and JP; *Popillia japonica* of US only, not of Japan; *Trichosurus vulpecula*, *Melaleuca lanceolata*, *Eucalyptus cinerea* and *Oncorhynchus nerka* of none. WCVP and WRiMS moved only with the release (357,706 WCVP species, 2,139 WRiMS).
- **GRIIS-AM** (dataset 24471) was not reimported with the others: its source archive had a data bug, and Markus fixed it by hand for one import (attempt 20, 2026-11-11), exported as 9cbcf7d8-edf4-4386-9d42-97ec1aa68fbb. Only its file and manifest entry were replaced; import 5's file is kept in `data/griis-before-reimport/am-import-5/`. Import 5 had broken encoding in 7 names (`Arctotis sp. Ã—hybrida`, `Larix kaempferi Ã¨re`). The list now has 386 species entries against 403: 334 kept, 75 gone, 52 new, with animals from 13 to 36 (*Salmo trutta*, *Carassius gibelio*, *Aedes albopictus*). Its one native entry, *Medicago arabica*, is now uncertain, so no list has a native entry left. `check_natives.py` again: **116,794 GRIIS alien pairs of 27,057 species in 305 lists**; nothing else changed.
- ChecklistBank now gives Cayman Brac its own alias, `GRIIS-KY-Brac` (it was `GRIIS-KY-CaymanBrac`, made up from the title). Natural Earth has no polygon for it under either name, as before.

### WCVP `native_only` stays the headline, `introduced` is reported beside it (2026-11-11, Markus)

Markus asked whether WCVP in mode `native_only` predetermines the result. There are two kinds of circularity.
- **Edges.** WCVP's ranges are coded in TDWG level 3 regions, so trimming plants to them could draw realm boundaries along those edges.
  The circularity gate tests this, and `native_only` passes it with room (+0.012 over `off`, against +0.126 by chance).
- **Knowledge.** `native_only` keeps a plant only where regional floras record it as native, so the botanists' judgement enters the data.
  No test against chance exists for that. The gate also has little power: 60% of land cells lie within one cell of a TDWG edge.

**The three runs side by side** (half A, eukaryotes, mt075, seed 1; `results/explore/res7/wcvp-modes/`, drawn by its `draw.py`):
- At 5 realms, `off` differs from `native_only` on 908 parts and 14.4M km² of land (V 0.908). Australia (13.3M km²), New Zealand and the
  Antarctic coast join the Holarctic.
- At 5 realms, `introduced` differs on 334 parts and 4.3M km² (V 0.946). Australia is with Indo-Malaya, as in `native_only`. New Zealand and
  the Antarctic coast (3.3M km²) stay Holarctic, and so do southern Arabia (0.5M km²) and specks of northeast Asia.
- At 9 realms (the select stage's candidates, not the nested cuts), `off` differs on 24.1M km² and `introduced` on 7.5M km².
  - `off`: the Holarctic does not split, and the extra realms are at sea.
  - `introduced`: the Nearctic, Patagonia and Australia as in `native_only`; New Zealand Palearctic; a North Atlantic–Arctic sea realm instead of
    temperate East Asia.
- **Why New Zealand moves.** `introduced` and `native_only` differ only in plant presences in regions WCVP says nothing about. In New Zealand's 70
  land parts, 5,995 such presences of 1,390 plant species are kept under `introduced`. They are led by European weeds and garden plants that WCVP
  does not record as introduced there (*Leontodon saxatilis*, *Cerastium fontanum*, *Mycelis muralis*, *Lapsana communis*, *Carduus nutans*,
  *Lunaria annua*, the potato). They also include a few natives that WCVP's native range misses and `native_only` wrongly removes (*Centipeda
  aotearoana*, *Coprosma pumila*, *Euchiton audax*, *Glossostigma diandrum*).
- `off` is no neutral baseline: alien plants pull whole continents into the Holarctic.

**The ruling.** `native_only` stays the headline, by the gate's rule. The `introduced` run is reported beside it at its own levels, with where
and why it differs (PREREGISTRATION sections 4, 8 and 10).

**My addition, for Markus to confirm or strike.** C1 and C2 on half B are also run with `introduced`, in the same single pass, as a reported
check that does not decide (section 9). It shows how much of the confirmation rests on the floras. It has to be fixed before B is opened, or it
would be a second look at B.

**The taxon lists go to Appendix A.** Section 4 said "domesticated species, major crops and humans". The list drops 16 species: humans, 13
domestic animals, bread wheat and rapeseed. Maize, rice, potato and seven others were looked at and kept, as COL files wild populations under
their names. Both lists are now in the preregistration with their reasons.

**Cultivated trees (open, Markus is thinking about it).** Markus suggested adding plantation trees such as eucalyptus and spruce to the list.
`native_only` already removes them outside their WCVP native regions:
- *Eucalyptus globulus* is native to TAS and VIC only, and introduced in 49 regions;
- *Pinus radiata* is native to CAL, MXI and MXN;
- *Picea abies* is native to 24 regions from France to the Urals.
Inside a native region a plantation can't be told from the wild at this resolution (spruce in lowland Germany). `NO_WILD_RANGE` is for species
with no wild population, so putting such a tree on it would delete its wild range as well.

**WCVP against World Plants and WFO (Appendix B).** In ChecklistBank:
- WCVP has 1,986,823 distribution rows, all coded as TDWG level 3 regions, native or introduced, for 371,856 species;
- World Plants (1141), COL's main source of vascular plants, has one free-text line per taxon (391,924 rows; introduced areas marked `[I]`);
- the WFO Plant List (2004) and POWO's checklist (2000) have none.

COL 316691 takes its coded plant distributions from WCVP: 599,462 of its 909,648 plant distribution rows are from the WCVP and WCVP-Fabaceae
sectors, and 275,333 are World Plants' text.

Coverage of COL 316691: 357,534 of 378,499 vascular species (94.5%) have a WCVP native range, from 85.5% (Pinopsida) to 96.7% (Liliopsida).
The 29,610 mosses, liverworts, hornworts and algae have none.

Coverage of half A's data (315557, WCVP off), as shares of plant presences (records in brackets):
- 86.9% (94.6%) belong to species with a native range;
- 9.5% (3.5%) are non-vascular;
- 3.4% (1.9%) are vascular plants WCVP doesn't match, mostly where the names differ. Some differences were 315557's and are gone in 316691
  (*Elytrigia repens*, now *Elymus repens*). Others remain: WCVP sinks *Taraxacum officinale* into an unresolved section, and keeps *Hydrogalium
  palustre* in *Galium*.

Unmatched species are never trimmed. The median vascular species is native to one level 3 region (57%), and the mean is 3.4.

### Crops: FAO's census list, dropped in full (2026-11-11, Markus's ruling)

Layer 2 drops every species on FAO's list of crops for agricultural censuses, whether or not it also grows wild. The list is
**World Programme for the Census of Agriculture 2030, Annex 7**: an alphabetical list of crops with botanical names and ICC 1.2 and CPC 3.0
codes (FAO 2026, Statistical Development Series 22, doi:10.4060/cd9437en). Markus builds the list as a ChecklistBank dataset, which gives the
mapping to COL ids. The pipeline fetches it like the other references.

**Why this source.** It is FAO's current, global list for national agricultural censuses, citable by crop code, and the successor to WCA 2010's
Appendix 4. The alternatives were weighed and set aside:
- the plant genetic resources treaty's Annex I lists genera only;
- FAOSTAT's items are commodity groups without botanical names;
- Mansfeld's database (~6,100 species) and GRIN's World Economic Plants (~12,000) are not official bodies, and are much broader.

**Why in full** (Markus). A per-species check of wild populations, as `NO_WILD_RANGE` did for 26 names, does not scale to some 250 species.
Dropping a few hundred species, whose records are mostly of cultivated plants, is easier to defend than a selection. The rule replaces the
reasoning of `NO_WILD_RANGE_DROPPED` for plants: maize, rice, barley and the rest go, wild relatives filed under their names included.

**What it costs on half A** (`results/explore/res7/fao-crops/`, `cost.py`). Annex 7 resolves to 253 COL 316691 species. 190 are named outright.
The other 63 are my reading of its genus entries as the species grown (*Coffea* → *arabica*, *canephora*, *liberica*; *Trifolium* → *pratense*,
*repens*, *alexandrinum*, *incarnatum*), for Markus to review in the dataset. 133 of the 253 occur in half A. Beyond what WCVP `native_only` already
removes, dropping them takes 20,474 presences (0.20%) and 8.87M records (0.49%). Most of those records are of plants growing wild in their own
range: hazel 1.36M, orchard grass 1.28M, bilberry 1.00M, creeping bent 1.00M, perennial ryegrass 0.99M, timothy 0.68M, wild carrot 0.66M.
Under `native_only` maize keeps 46 presences.

**Animals: open, asked of Markus.** WCA 2030's livestock classification (Annex 8) names categories only. FAO's DAD-IS names the species, about
38–40 mammals and birds plus managed bees. Several of them are wild species with large native ranges, and their records are mostly of wild animals.
Taken whole, they would remove 2.9% of half A's records (51.4M): mallard 27.0M, rock dove 9.8M, greylag goose 4.9M, wild turkey 2.8M, red fox 1.8M
(fur animals), honeybee 1.2M. Proposed instead: keep the current rule for animals, dropping a domestic animal where COL names it apart from its wild
ancestor (cattle and aurochs, pig and wild boar), and add the ones COL also names apart (yak, gayal, silkworm).
