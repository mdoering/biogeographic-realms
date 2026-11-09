# Preregistration: global biogeographic realms from GBIF occurrences

**Status: DRAFT for Markus to edit, 2026-11-09.** Drawn from `docs/design.md` and `DECISIONS.md` up to 2026-11-09, which hold the
reasons. Markus registers the final version on OSF before half B is unsealed (design, "Exploring and confirming", step 3). Items marked
**TO SETTLE** are proposals, not yet decided.

The rules are what is registered, not the numbers they gave on half A. On half B and on all species, the numbers are whatever the
rules give there.

## 1. Question and scope

- **Question.** Which biogeographic realms do GBIF occurrence records support, when realms are defined by shared restricted taxa, and
  do lineages not used to find them confirm them?
- **Taxonomy:** Catalogue of Life release [315557](https://www.checklistbank.org/dataset/315557) (COL26.6 XR), the release GBIF
  production uses. Accepted species only, provisionally accepted species included (DECISIONS, 2026-11-06). Taxa are identified by
  COL id, not by name.
- **Taxa: eukaryotes.** Animalia, Chromista, Fungi, Plantae and Protozoa. Bacteria, archaea and viruses are left out: on half A, 10 of
  12 testable prokaryote classes showed no realm signal on land (DECISIONS, "Which groups carry the realm signal" and "The scope is
  eukaryotes", 2026-11-07). No eukaryote class is left out.
- **Ranks:** species, genus and family.
- **Not in scope:** whether the realms resemble existing schemes (Costello et al. 2017, RESOLVE). That comparison comes after the
  confirmation (milestone M4) and never feeds back into a choice registered here.

## 2. Data

- **Occurrences:** GBIF SQL download `0020744-261023074918325`, DOI [10.15468/dl.baekdm](https://doi.org/10.15468/dl.baekdm):
  species × ISEA3H resolution 7 cell × land/sea part × country, with per-record count columns. The query is `queries/occurrence_res7.sql`.
  - Hard filters: coordinates present, no geospatial issue, `occurrencestatus = PRESENT`, species-level match, not fossil, coordinate
    uncertainty at most 100 km.
  - Grid: `GBIF_ISEA3HCode` at resolution 7 with uncertainty argument 0. A cell's land and sea parts are separate nodes.
- **Evidence is records**, not species-days. The visits check on half A changed no coarse realm (DECISIONS, "Visits against records:
  the check", 2026-11-06), and the cleaning was calibrated on records.
- **Reference data:**
  - the COLDP export of 315557: classification, extinct flags, `environment`;
  - WoRMS environment flags (ChecklistBank match of WoRMS dataset 2011 to 315557), filling where COL has none;
  - WCVP (CLB 2232): native regions of plants;
  - GRIIS country and sub-national lists, with the frozen country-code table;
  - WRiMS (CLB 26853);
  - TDWG level 3 regions;
  - Natural Earth land.

  `data/reference/manifest.json` records each file with its hash.

## 3. The split into halves

- Whole families (lineages) are assigned at random to halves A and B, balanced by major group and species count. The seed is `20261031`.
  The split file is `data/split/halves-315557.parquet`, content hash `SPLIT_SHA256 = 88f237b47d4964108c161e08acb00c9cf03601418dd2e0f16b49a43704845fa9`
  (DECISIONS, "Family split frozen", 2026-10-31).
- Everything in sections 4 to 7 was explored and fixed on half A only.
- Half B has been used once only. That one load scored the cleaning calibration on B's calibration species (3.19% of native presences
  lost, against 3.26% on A; DECISIONS, 2026-11-01). Nothing has been clustered, mapped or inspected per taxon on B.

## 4. Cleaning

The layers run in order (design, "Cleaning"; `src/realms/clean.py`):

1. **Record evidence:** drops living specimens, non-native establishment flags, fuzzy matches, kingdom clashes, matches to ambiguous
   synonyms and misapplied names, and the multilingual text patterns of the query.
2. **Taxon lists:** drops extinct species, domesticated species, major crops and humans.
3. **Authoritative ranges:**
   - WCVP in mode `native_only`: plants outside their native regions are dropped, and a cell touching a native region counts as native;
   - GRIIS by each entry's own establishment means, with sub-national lists on the cells they overlap;
   - WRiMS.
4. **Land/sea consistency.** A species is clearly terrestrial or marine by its COL environment, else WoRMS's, else when two thirds of its
   records lie on one side. Species of land classes (`clean.TERRESTRIAL_CLASSES`) with no environment count as terrestrial, except the
   known marine groups (`clean.MARINE_EXCEPTIONS`). A clear species' records on the wrong side fold into the other part of the same cell
   (DECISIONS, 2026-11-05).
5. **Rare-outlier heuristics:** at **setting 132** of the calibration:

   | Threshold | Value |
   |---|---|
   | satellite record share | 2% |
   | satellite cell share | 20% |
   | under-reported rate factor | 0.5% |
   | group minimum | 300 records |

   Everything else is as in the design: satellite gap 500 km, at least 2 datasets and 3 years, at most 5 records for the under-reported
   and thin long-range rules, 1,000 km, the core share 90%, no heuristic below 20 records. Setting 132 was chosen by the convex-hull
   walk at break-even, with at most 5% of native presences lost (DECISIONS, 2026-11-01/02).

**The WCVP mode is chosen by the circularity gate** (section 8). The headline mode is the first of `native_only` and `introduced` whose
gate passes, else WCVP off.

**TO SETTLE: how `group_min_records` (300 on half A) scales to all species.** This is an open point from 2026-11-01.

## 5. Units of the network

- A cell part with Good-Turing coverage under 0.5 or under 50 species merges with its neighbours into a resolution 5 hexagon.
- A part still sparse after that stays out of the network. So does a part where fewer than two major groups are sampled.
- A part left out takes the majority realm of its nearest network parts of the same kind, and is flagged as imputed.
- Units that Infomap puts in no module leave the scores (DECISIONS, 2026-11-02).

## 6. Network, clustering and the choice of configuration

- **Network.** Bipartite: cell-part units against taxa at species, genus and family rank. Links are presence only. Genus and family
  links are weighted by the spatial information they carry, as in Infomap Bioregions 2. Genus and family presence comes from the
  cleaned species only.
- **Configurations:** seven of Infomap's scale, as in `explore/configs.toml`:
  - `default`;
  - Markov time 0.95, 0.90, 0.85, 0.80 and 0.75;
  - `regularized`.
- **Candidates.** Infomap's finest modules are merged by beta-sim into one nested map per realm count, from 60 down to 5, with a size
  floor of 1M km² (a smaller unit becomes a sub-realm).
- **Per seed:** select at 10 trials, 20 bootstraps of 3 trials and 99 null maps. Stability is the median bootstrap V-measure. Excess is
  the area-weighted endemism excess over the null, averaged equally over species, genus and family. The null maps are random
  partitions with the same realm count, realm areas and number of pieces.
- **The five-seed choice** (DECISIONS, "Stability over five seeds", 2026-11-03):
  - Every configuration is run at seeds 1 to 5. A map's stability is the median of its five seeds, with a seed lacking the map counted
    as 0.
  - Of the maps with stability at least 0.8, the highest median excess is chosen. Within 0.01, the coarser map wins.
  - If none reaches 0.8, the most stable map is reported as falling short.
  - The headline map is seed 1's map of the chosen configuration.
- **On half A** this chose `mt075` (Markov time 0.75) at k = 5: excess 0.355, median stability 0.961.

## 7. Nested levels and what is reported

- **Merge order.** The chosen configuration's finest modules, with their pieces resolved, are merged by beta-sim on the union of their
  taxa. Land merges with land and sea with sea; the land-sea join comes last and is never cut (`realms.nested`).
- **Rule 2, robust to the sampling of species** (`nested.robust_levels`):
  - The units are merged freely to a base cut of 0.3.
  - The nested merges are then replayed on 199 resamples of the species.
  - A level is robust where one cut separates the merges made from those to come in at least 95% of the resamples.
  - Its cut is the middle of the gap.
- **Rule 6, more endemic than chance** (`nested.realm_excess`):
  - At each robust level, each realm's endemic share (at least 95% of a taxon's cells inside, averaged over the three ranks) is set
    against the same realm in 99 random maps of the same realm areas and pieces.
  - A level passes when every realm has p ≤ 0.05.
- **The reported levels** (`nested.reported_levels`; DECISIONS, "The reporting rule: a third level", 2026-11-08). Of the levels that pass
  both rules, the report gives:
  1. the coarsest;
  2. the end of the unbroken run from it, where each level adds exactly one realm;
  3. the finest.

  The levels between are given as the merge order. If the run reaches the finest level, two levels are reported; if one level passes,
  one is.
- **Recovery, reported and not a gate.** Each realm's recovery is computed over 20 resamples of the species: Infomap with 10 trials,
  then the same cuts (`nested.bootstrap_recovery`, `scripts/nested.py --boot-trials 10`). It is reported as the median Jaccard and the
  share of resamples with Jaccard at least 0.5. A realm recovered in under half of them is flagged.
- **On half A** the rules give 5, 9 and 14 realms. Every realm at every level passes rule 6 at p = 0.01. Temperate East Asia is the one
  flagged realm (recovered in 15%).
- **The ocean** is one realm wherever the data give one. No split is imposed (DECISIONS, 2026-11-02).

## 8. Gates and cross-checks

- **Circularity gate** (the TDWG-edge gate; DECISIONS, 2026-11-02; `scripts/circularity.py`):
  - The cleaning runs are: WCVP off (the base), `introduced`, `native_only`, and the layer variants without the extinct and
    no-wild-range rules and without GRIIS. Each runs on the chosen configuration at seed 1.
  - For each mode, the share of land boundary near TDWG level 3 edges is compared with the base's. The test uses 99 rotated region
    sets (tilt up to 30°, rotations kept within 0.5 of the regions' land boundary length).
  - A mode passes when its raise over the base is at most the 95th percentile of the rotated raises. The verdict at land tolerance 0.5
    decides; the reading without a tolerance is reported too.
  - On half A every mode passes, and the headline mode is `native_only`.
- **Border gate** (design, "Robustness runs"):
  - A boundary segment within one cell width of a national border, which also disappears without layer 3 or in the equal-coverage
    run, is flagged.
  - If flagged segments exceed 20% of the boundary length, the headline switches to the run without layer 3 when that run stays under
    20%. Otherwise the map is published as qualified.
- **Cross-checks, reported and not deciding:**
  - beta-sim over the three ranks with Ward and with average linkage, by area agreement with Infomap at the chosen k;
  - the run with duplicate genus and family names merged, by V-measure.

## 9. Confirmation on half B (once)

The frozen pipeline (sections 4 to 8) runs once on half B. Any failure is reported as a failure, and half B cannot be reused.

- **C1. Do B's unseen lineages show endemism inside A's realms?**
  - At each of A's reported levels, B's cleaned taxa are scored against A's realms. The measure is the area-weighted endemism excess,
    tolerant (95%), at species, genus and family rank separately.
  - The null is the same: 99 random maps with the same realm areas and pieces.
  - Each realm's excess on B against its null is reported too.
- **C2. Does a map built from B alone match A's?**
  - B goes through the whole frozen pipeline: the five-seed choice, the nested levels, the cut rules and the reporting rule.
  - B's merge order is cut at each of A's reported realm counts. The area-weighted V-measure against A's map at that count is taken
    over the parts both map.
  - B's own reported levels are reported beside A's.

**TO SETTLE: the pass marks.** Proposals:
- **C1 passes** when, at every reported level, B's excess is above the null (p ≤ 0.05) at each of the three ranks. A realm whose own
  excess on B is not above the null (p > 0.05) is reported as not confirmed. It does not fail the whole test.
- **C2 passes** when V ≥ 0.8 at the coarsest reported level. This is the bar the stability already uses. On half A, maps at different
  seeds and in different cleaning modes agree at V 0.87-0.99. V at the finer reported levels is reported against the same 0.8, but
  does not decide.
- **The overall result:**
  - confirmed if C1 and C2 both pass;
  - partly confirmed if one passes;
  - not confirmed if neither does.

## 10. Headline on all species

- The frozen pipeline runs on all species (halves A and B) with the thresholds above. The configuration, k and levels are chosen again
  by the rules (sections 6 and 7), not carried over from half A.
- The WCVP mode is chosen again by the circularity gate.
- Any change after registration is reported as a deviation, with its reason and its effect.

## 11. Robustness runs (reported, not tuned)

As listed in the design ("Robustness runs"):
- raw data, then each cleaning layer added in turn;
- equal-coverage subsampling;
- without eBird, and without iNaturalist;
- per taxon group, and without birds;
- whether a cell's sampled composition predicts its realm;
- a minimum of 5 species instead of 50 for a unit;
- species only, species plus genus, and family only;
- dated trees against ranks where trees exist;
- resolution 6 against 7;
- beta-sim with UPGMA against Infomap;
- lenient and strict heuristic settings;
- without the most widespread species;
- CC0 and CC BY records only;
- the border gate.

None of these changes the headline.

## 12. Still open before registration

- **The pass marks of section 9.**
- **The scaling of `group_min_records`** to all species.
- **Two source issues**, to be noted rather than fixed in the pipeline:
  - GRIIS-US flags *Passer domesticus* and *Columba livia* as native;
  - ChecklistBank dataset 28940 (Cayman Brac) carries the alias GRIIS-KY.
- **The exact commit and data hashes** that the registration cites. These are filled in at the freeze.
