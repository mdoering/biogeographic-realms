# Preregistration: global biogeographic realms from GBIF occurrences

**Status: DRAFT for Markus to edit, 2026-11-09.** Drawn from `docs/design.md` and `DECISIONS.md` up to 2026-11-09, which hold the
reasons. Markus registers the final version on OSF before half B is unsealed (design, "Exploring and confirming", step 3). Every
method choice is settled; section 12 lists what is still to fill in.

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
2. **Taxon lists:** drops the species COL flags as extinct, and the 16 species of `lists.NO_WILD_RANGE`: humans, 13 domestic animals and
   2 crops (bread wheat and rapeseed) with no wild population. Ten other domestic animals and crops (maize, rice, potato …) are kept, as
   COL files wild populations under their names. Both lists, with the reason for each entry, are in Appendix A.
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
gate passes, else WCVP off. The `introduced` run is always reported beside the headline (section 8). What WCVP covers, and why it is the
source of plant ranges rather than World Plants or WFO, is in Appendix B.

**`group_min_records` scales with the records** (Markus, 2026-11-09). The calibration chose 300 on half A, whose group totals per cell are
about half those of all species. A run uses 300 times its records over half A's records, so that the threshold keeps the meaning it had
on A: about 600 on all species, and B's own ratio on half B. The exact values are filled in at the freeze.

## 5. Units of the network

- A cell part with Good-Turing coverage under 0.5 or under 50 species merges with its neighbours into a resolution 5 hexagon.
- A part still sparse after that stays out of the network. So does a part where fewer than two major groups are sampled.
- A part left out takes the majority realm of its nearest network parts of the same kind, and is flagged as imputed.
- Units that Infomap puts in no module leave the scores (DECISIONS, 2026-11-02).

## 6. Network, clustering and the choice of configuration

- **Network.** Bipartite: cell-part units against taxa at species, genus and family rank. Links are presence only. Genus and family
  links are weighted by the spatial information they carry, as in Infomap Bioregions 2. Genus and family presence comes from the
  cleaned species only: a genus or family is present in a cell when any one of its species is (DECISIONS, 2026-11-09).
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
  - **What the gate does not test.** It tests one kind of circularity: realm boundaries drawn along the edges of the regions WCVP's
    ranges are coded in. It does not test the other: `native_only` keeps a plant only where regional floras record it as native, so the
    botanists' judgement behind those floras enters the data. That cannot be measured against chance, so it is shown by comparison
    (next point). The gate also has little power: at resolution 7, 60% of the cells with land lie within one cell of a TDWG edge.
- **`introduced`, reported beside the headline** (Markus, 2026-11-11). The run with WCVP in mode `introduced` uses WCVP only where it
  records an introduction, as GRIIS does, and imports no flora's native range. It is reported at its own levels by the same rules, with the
  parts and area where it differs from the headline and why. On half A, at 5 realms, it differs on 4.3M km² of land (V 0.946). New Zealand
  and the Antarctic coast join the Holarctic, held there by European weeds that WCVP does not record as introduced in New Zealand. WCVP
  off differs on 14.4M km² (V 0.908): Australia and New Zealand join the Holarctic (`results/explore/res7/wcvp-modes/`).
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

**The pass marks** (Markus, 2026-11-09):
- **C1 passes** when, at every reported level, B's area-weighted excess is above the null (p ≤ 0.05) at each of the three ranks, and the
  family rank is also above the null over A's land realms alone. The land test weighs the land realms only (`nested.sea_realms`: under half
  their area in sea parts), each null map's realm at the same position standing in for a land realm (`scripts/rank_endemism.py`, scope
  `land`). It is there because the sea realm carries most of a level's family score. A realm whose own excess on B is not above the null
  (p > 0.05) is reported as not confirmed. It does not fail the whole test. On half A every rank passes at every level on the whole map
  and on land (p = 0.01), and per realm only New Guinea fails, at family rank at 14 realms.
- **C2 passes** when V ≥ 0.8 at the coarsest reported level. This is the bar the stability already uses. On half A, maps at different
  seeds and in different cleaning modes agree at V 0.87-0.99. V at the finer reported levels is reported against the same 0.8, but
  does not decide.
- **The overall result:**
  - confirmed if C1 and C2 both pass;
  - partly confirmed if one passes;
  - not confirmed if neither does.
- **C1 and C2 with WCVP `introduced`, reported and not deciding.** In the same single pass, half B is also cleaned with WCVP in mode
  `introduced` and scored by C1 and C2 against A's `introduced` map. This shows how much of the confirmation rests on the floras' native
  ranges (section 8). It is fixed here, before B is opened, so that it is not a second look at B.

## 10. Headline on all species

- The frozen pipeline runs on all species (halves A and B) with the thresholds above. The configuration, k and levels are chosen again
  by the rules (sections 6 and 7), not carried over from half A.
- The WCVP mode is chosen again by the circularity gate, and the `introduced` run is reported beside the headline, as in section 8.
- Any change after registration is reported as a deviation, with its reason and its effect.

## 11. Robustness runs (reported, not tuned)

As listed in the design ("Robustness runs"):
- raw data, then each cleaning layer added in turn;
- equal-coverage subsampling;
- without eBird, and without iNaturalist;
- per taxon group, and without birds;
- whether a cell's sampled composition predicts its realm;
- a minimum of 5 species instead of 50 for a unit;
- species only, species plus genus, and family only. Each rank on its own was run on half A (DECISIONS, 2026-11-09): species alone
  comes close to the headline (V 0.84 at 5 realms, 0.85 at 9), genus and family alone give no stable map;
- dated trees against ranks where trees exist;
- resolution 6 against 7;
- beta-sim with UPGMA against Infomap;
- lenient and strict heuristic settings;
- without the most widespread species;
- CC0 and CC BY records only;
- the border gate.

None of these changes the headline.

## 12. Still open before registration

- **Two source issues, both settled on 2026-11-11** (DECISIONS, "The GRIIS lists after ChecklistBank's reimport"):
  - GRIIS-US flagged *Passer domesticus* and *Columba livia* as native. This was how checklistbank.org read GRIIS's terms. ChecklistBank
    now reads them correctly, and all 305 lists were fetched again. An establishment means ChecklistBank leaves empty, mostly GRIIS's
    `Native|Alien`, is left out, as `uncertain` is;
  - ChecklistBank dataset 28940 (Cayman Brac) carried the alias GRIIS-KY. It now has its own, GRIIS-KY-Brac.
- **The move to COL 316691** (DECISIONS, 2026-11-11). Sections 2 and 3 still describe release 315557 and download 0020744. They are
  brought up to date once the new download is cleaned.
- **The exact commit and data hashes** that the registration cites. These are filled in at the freeze.

## Appendix A. The taxon lists of layer 2

`src/realms/lists.py`. Each name was checked against the COL release: it must be an accepted species there, with nothing wild filed under
it (`lists.resolve_no_wild_range` checks this again at every run and drops an entry that fails). The extinct species are those COL flags.

**Dropped (`NO_WILD_RANGE`, 16 species):**

| Species | Why |
|---|---|
| *Homo sapiens* | humans: records are people, observed or photographed, wherever they were |
| *Bos taurus* | domestic cattle; the wild aurochs is *Bos primigenius*, a species of its own in COL |
| *Bos indicus* | zebu, domestic cattle with no wild population |
| *Bubalus bubalis* | domestic water buffalo; the wild water buffalo is *Bubalus arnee* |
| *Equus caballus* | domestic horse, whose feral herds have no natural range; the wild horse is *Equus ferus* |
| *Equus asinus* | domestic donkey; the wild ass is *Equus africanus* |
| *Canis familiaris* | domestic dog, held apart from the wolf; the dingo is filed here, a feral dog with no natural range |
| *Felis catus* | domestic cat; the wildcat is *Felis silvestris* |
| *Sus domesticus* | domestic pig; the wild boar is *Sus scrofa* |
| *Camelus dromedarius* | dromedary, domestic, and extinct in the wild |
| *Camelus bactrianus* | domestic Bactrian camel; the wild camel is *Camelus ferus* |
| *Lama glama* | llama, domestic; the guanaco is *Lama guanicoe* |
| *Vicugna pacos* | alpaca, domestic; the vicuña is *Vicugna vicugna* |
| *Cavia porcellus* | domestic guinea pig; wild cavies are *Cavia aperea* |
| *Triticum aestivum* | bread wheat, a cultigen with no wild population |
| *Brassica napus* | rapeseed, a cultigen of two other species with no wild population |

**Looked at and kept (`NO_WILD_RANGE_DROPPED`, 10 species)**, because COL files wild populations under the same accepted name, so
dropping the species would drop its wild records too:

| Species | Wild populations COL files under it |
|---|---|
| *Ovis aries* | the mouflon (*Ovis orientalis*, *O. orientalis ophion*, *O. aries musimon* are synonyms) |
| *Gallus gallus domesticus* | no such species: *Gallus domesticus* is a synonym of *Gallus gallus*, the red junglefowl |
| *Capra hircus* | the Chiltan markhor |
| *Zea mays* | the teosintes (subspecies *parviglumis* and *huehuetenangensis*) |
| *Hordeum vulgare* | wild barley (subspecies *spontaneum*) |
| *Glycine max* | wild soybean (subspecies *soja*) |
| *Solanum tuberosum* | wild potatoes, as 186 species-rank synonyms; WCVP gives it a native range in the Andes |
| *Oryza sativa* | wild perennial rice (*Oryza perennis*) |
| *Sorghum bicolor* | wild sorghum (subspecies *verticilliflorum*) |
| *Capsicum annuum* | the wild chiltepin (variety *glabriusculum*) |

Kept crops and plantation trees outside their native range are removed by WCVP in mode `native_only` (layer 3). Inside a native region
they stay: a planted spruce in lowland Germany can't be told from a wild one at the resolution of WCVP's regions.

## Appendix B. WCVP: what it covers, and why it is the source of plant ranges

**The global plant checklists in ChecklistBank** (2026-11):

| Dataset | Species | Distribution rows | Form |
|---|---|---|---|
| WCVP (2232; Govaerts et al. 2021, *Scientific Data* 8:215) | 371,856 | 1,986,823 | TDWG level 3 regions, each native or introduced |
| World Plants (1141), COL's main source of vascular plants | 355,419 | 391,924 | one free-text line per taxon, introduced areas marked `[I]` |
| WFO Plant List 2026-06 (2004) | 382,438 | none | — |
| POWO's checklist (2000) | 368,929 | none | — |

**Why WCVP:**
- **Only WCVP has native ranges coded for nearly every vascular plant species.** The method needs, for each plant, the areas where it is
  native, in units that map onto the grid. World Plants holds the same knowledge as text ("China (Guangdong, Hainan, S-Yunnan); India
  (Arunachal Pradesh …)"), which would first have to be parsed into areas, a project of its own. The WFO Plant List and POWO's checklist carry
  no distributions in ChecklistBank.
- **COL itself takes its coded plant distributions from WCVP.** Of the 909,648 plant distribution rows in COL 316691, 599,462 come from WCVP
  (its WCVP and WCVP-Fabaceae sectors, all TDWG level 3), 275,333 are World Plants' text, and the rest come from small sources.
- **Only ranges are taken from WCVP, never its taxonomy.** The species are COL's. WCVP's accepted species are matched by name to COL's
  (98.5% of names), and a COL species WCVP does not name is never trimmed.

**What it covers in COL 316691:**

| Group | Species | With a WCVP native range |
|---|---|---|
| Magnoliopsida | 275,331 | 258,949 (94.1%) |
| Liliopsida | 87,188 | 84,307 (96.7%) |
| Polypodiopsida | 12,759 | 11,711 (91.8%) |
| Lycopodiopsida | 1,536 | 1,422 (92.6%) |
| Pinopsida | 902 | 771 (85.5%) |
| Cycadopsida | 391 | 373 (95.4%) |
| All vascular plants | 378,499 | 357,534 (94.5%) |
| Mosses, liverworts, hornworts, algae | 29,610 | none: WCVP covers vascular plants only |

The vascular plants also include 390 Ginkgoopsida, of which 389 are fossils. A native range is a set of TDWG level 3 regions (botanical countries
and states): the median species is native to one region (57%), and the mean is 3.4.

**What it covers in the data** (half A, release 315557, WCVP off; shares of plant presences, with shares of records in brackets):
- 86.9% (94.6%) belong to species with a WCVP native range: these are what `native_only` can judge;
- 9.5% (3.5%) are mosses, liverworts and algae, which WCVP never touches;
- 3.4% (1.9%) are vascular plants WCVP does not match. Mostly COL and WCVP disagree on the name. Some of these differences were 315557's
  and are gone in 316691 (*Elytrigia repens*, now *Elymus repens* in both). Others remain: WCVP sinks *Taraxacum officinale* into an unresolved
  section, and keeps *Hydrogalium palustre* in *Galium*.

Where WCVP is silent, the presences are kept, so a gap in WCVP never cuts a range. In the other direction, `native_only` removes a native plant
where WCVP's native range misses a region. In New Zealand this happens to a few species (*Centipeda aotearoana*, *Coprosma pumila*), against
many more European weeds removed that WCVP does not record as introduced there (section 8).
