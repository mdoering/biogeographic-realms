# Global biogeographic realms derived from GBIF occurrences

## Context

The `realm` gazetteer in col-gazetteers tiles the globe with 19 realms: 8 RESOLVE land realms plus 11 marine realms grouped
by hand from Costello et al. 2017. None of it comes from GBIF, and land and sea come from different methods.

This project derives realms for the whole globe, land and sea together, from GBIF records alone, defined by endemism at
species, genus and family level. It then tests how well the existing schemes hold up. The result is a **reproducible
scientific analysis**, not a new gazetteer. Whether a gazetteer follows is decided afterwards.

**COL release.** Taxonomy is whichever release GBIF production interprets occurrences against on the day of the download,
because SQL downloads run on production. Read it from `/v2/species/match/metadata` and repeat the probe below.

Checked 2026-10-31:
- **Production is 315557 (COL26.6 XR).** The registry and matcher metadata say so, and the probe agrees: all 10 sampled ids
  first issued after 315557 got a different match on production, with no occurrences. The 5 rechecked individually matched
  only at genus level.
- **UAT and test hold a newer XR**, key 316700 on the *test* ChecklistBank ("COL XR for GBIF TEST"). On production ChecklistBank
  that key is an unrelated dataset. Three of the four sampled names match there as species, while one first issued in COL26.9
  XR doesn't.
- Classification, environment and extinct flags must come from the same release as the download, because ids change between
  releases.
- A newer release isn't automatically safer. *Alopecurus pratensis* L. (`C56N`, 616,505 GBIF records) is accepted in 315557 but
  absent from COL26.9 XR, COL26.10 and COL26.10 XR. Cause not investigated. If production moves, spot-check common species first.

## Goals, non-goals, anti-goals

**Goals**
- A full nested hierarchy of realms and one headline level, covering every point on Earth exactly once.
- A per-cell confidence value, so poorly sampled areas look uncertain.
- Evidence of how much data errors and non-native records move the boundaries.
- A quantitative comparison with `realm` (19) and `costello` (30).

**Non-goals**
- Changing the backend or the published gazetteers.
- Correcting GBIF records or matches.
- Fine-scale bioregions or ecoregions.
- Modelled ranges. Costello's trial with AquaMaps returned the models' environmental inputs rather than endemism.
- Depth layers. The sea is two-dimensional here, as in Costello.

**Anti-goals: these make it a failure even if it runs**
- Realms that trace where GBIF has data, or national borders.
- Tuning thresholds until the map resembles the traditional realms.
- Merging clusters or smoothing boundaries by hand, the main objection to Costello's realms.
- Cleaning that removes narrow endemics, since those are what define a realm.

## Decisions

**Project decisions**
- One network for the whole globe, with coastal cells split into a land part and a sea part.
- All species with usable records, as presence per cell, with reruns per major group.
- The data decides the number of realms, from a nested hierarchy.
- Errors, misidentifications, wrong matches and non-native records (cultivated, zoo) must be handled explicitly, using
  heuristics against rare outliers where possible.
- Units below the size floor merge into the realm they share the most species with, and stay visible as named sub-realms.
- A realm may consist of several separate pieces.
- The headline uses all records, including CC BY-NC. Open-only is a robustness run.
- A paper is possible, so the downloads carry Markus's name.
- The work lives in a new repo in the `gbif` org.
- Evolutionary depth matters: it counts if whole genera or families are endemic, not only species.

**Implementation choices**
- Everything else below, including all numeric thresholds. They are starting values, frozen as described under "Exploring
  and confirming".
- OSF registration happens when the pipeline is frozen, just before the confirmation run, not up front.
- The repo is `gbif/biogeographic-realms` (local `~/code/gbif/biogeographic-realms`), Apache-2.0, in Python (DuckDB/Polars,
  `infomap`, GeoPandas, in a venv). It starts private and goes public at the freeze, since it carries the GBIF name and
  private can become public but not the reverse. Approving this plan authorises creating it. Nothing goes into the backend repo.
- GBIF Engineering Handbook rules apply: `README.adoc`, an issue for each non-trivial piece of work, Conventional Commits
  referencing the issue (`feat(cleaning): ... #12`), docs updated with the code. Markus decided that all work happens directly
  on `main`, with no feature branches, because this is a one-person repo. The handbook makes the committer responsible for
  having read AI-written code and discourages unreviewed AI-written issues, so Markus reviews the tracking issue.

## What Costello et al. 2017 did with the same problems

Read from the paper, the critique (Castro-Insua et al. 2018) and the authors' reply. OBIS is GBIF's marine subset, so their
data had the same kinds of error, with about 200 times fewer records (18.1M named records in the 2009 download, against 3.7B).

| Problem | Costello et al. 2017 | This project |
|---|---|---|
| Wrong names and matches | Matched to WoRMS and inspected by hand: 110,000 names to 93,000, then 65,000 used | GBIF's COL match. Fuzzy matches and kingdom clashes dropped. Hand inspection doesn't scale to this size |
| Wrong coordinates | Records on land dropped; "questionable" cells omitted, with five c-squares excluded by hand | By rule, in layers 1, 4 and 5 |
| Misidentifications | Not addressed | Outlier heuristics, tolerant endemism |
| Non-native species | Not filtered; introductions mentioned only as an open question | Layers 1, 2, 3 and 5 |
| Rare outliers | Coarse 5° cells and presence counts. Isolated clusters merged by hand. Single-cell species excluded only from the list of characteristic species | Per-species heuristics |
| Endemism | Strict: unique to one realm. The 42% appears to include single-cell species | Tolerant and strict, all species and well-recorded species |
| Similarity measure | Jaccard as headline, which the critique says is biased by richness differences | Infomap, with beta-sim as cross-check. Both were among Costello's alternatives |
| Number of realms | Over 200 clusters merged by hand into 30, at varying similarity levels, with one boundary smoothed along a temperature front. The critique could not reproduce it and saw no reason for 30; its own test gave 7 | One rule-based map per count, with selection rules and a hold-out test |
| Grid | 5° latitude-longitude, 2,056 cells, boundaries up to 1,200 km wide | Equal-area hexagons |
| Evolutionary depth | Realms from species only. A genus-level rerun gave the same overall structure but less well defined, and they warn that genus concepts are arbitrary | Species, genus and family in one network, with ranks tested against real phylogenies |

**What this changes in the plan**
- **Sampling composition.** The reply names an effect easily missed: neighbouring cells can look different because one
  had only benthic samples and the other only pelagic. In GBIF, the equivalent is a birds-only cell beside a plants-only cell.
  Hence the composition gate.
- **Sparse cells look alike.** Distant cells with only widespread species cluster together. That directly threatens allowing
  separate pieces, hence the shared-species test for pieces.
- **Empty sea.** They rejected 1° cells because of empty squares, so sparse sea cells are coarsened before any are imputed.
- **Per-group runs may fail.** They dropped per-group results over poor spatial coverage, so ours report only where coverage allows.
- **Costello is a comparator, not truth.** The critique found few similarities between their 30 realms and a beta-sim map. The
  reply holds that the overall pattern is similar across methods. The marine half of `realm` rests on these realms.
- **A public test case.** Their species-by-cell matrix is public, so the pipeline can be tested on it.
- **Inference, not from the papers:** with 200 times more records, each species collects more stray records, so strict
  endemism and unfiltered presence should suffer more here than they did for them.

## Data

One GBIF SQL download gives a table of `species x cell x land/sea x country` with per-record counts. Counting each filter
as a column rather than applying it in `WHERE` means every cleaning layer and sensitivity run works from the same download.

- **Grid:** `GBIF_ISEA3HCode`, equal-area and global, uncertainty argument 0 so cells are deterministic. Headline at
  resolution 7, about 22k cells of about 23k km² each. Resolution 6 is both the pilot and the scale test. Sizes come from
  the aperture-3 formula, so confirm them against `download.gbif.org/grids/ISEA3H/`.
- **Land/sea part:** `level0gid IS NULL`. Check how lakes and the Caspian behave.
- **Hard filters only:** coordinate present, no geospatial issue, `occurrencestatus = PRESENT`, species-level match, not
  fossil, uncertainty at most 100 km.
- **Count columns:** `n_raw` counts everything past the hard filters. The rest count records passing the layer 1 flag checks
  (not `LIVING_SPECIMEN`, no non-native establishment flag, no fuzzy match, no kingdom clash): `n`, `n_datasets`, `n_years`,
  `n_nc` (CC BY-NC), `n_open` (CC0 and CC BY), `n_datasets_open`, `n_years_open`, `n_ebird`, `n_inat`, `n_human_text`. The `_open` columns give the
  open-only run its own support counts.

**Checked 2026-10-31**
- The query in exactly this shape passes `POST /occurrence/download/request/validate`, with a placeholder text pattern list.
  Validation doesn't run it, so size and runtime are unknown.
- 3.70 billion georeferenced records.
- 17.4% of them (646M) are CC BY-NC, led by iNaturalist (140M) and Observation.org (138M). The share is 31% in Europe, 19%
  in Oceania, 12-14% in Africa, Asia and South America, 6% in North America. eBird is CC BY.
- Only 0.37% are flagged `introduced`, and `cultivated` plus `captive` flags cover under 120k records. About 97% carry no
  establishment information. **Flags alone cannot do this job.**
- `LIVING_SPECIMEN` is 0.18%, `TAXON_MATCH_FUZZY` 0.15%.

**Traps**
- `specieskey` holds COL id strings, although `describe/sql` says INT.
- `taxonomicIssues` was null on the sample record while `issues` held the flag, and `ARRAY_CONTAINS` on a null array makes
  `NOT` drop the record. Use `issue`, wrapped in `COALESCE(..., FALSE)`. `issue` merges flags from both COL and the GBIF
  backbone, so a fuzzy flag from either drops the record. At 0.15% that loss is acceptable.
- `isinvasive` exists in the schema but its meaning is unverified (the API ignored it as a filter), so it isn't used.

**Reference data**
- COLDP export of the release GBIF production uses: extinct flag, `environment`, classification. In 315557 only about 36% of
  2.35M accepted species have an environment value, so it can't be the sole land/sea evidence. How many species have a
  family and genus is not yet measured. Genus and family presence is derived offline from this classification, so the
  download query doesn't change.
- WCVP (CLB 2232, CC BY): native and introduced TDWG regions per plant. Potato shows 10 unflagged (native) and 24
  `introduced` regions; confirm on a larger sample that an empty value means native.
- GRIIS (about 380 CLB datasets, per country): alien species for all groups.
- WRiMS (CLB 26853): introduced marine species.
- Natural Earth coastline (public domain) for the final geometry; `realm` and `costello` GeoJSON from col-gazetteers.
- Costello's species-by-cell matrix (figshare `10.17608/k6.auckland.5086654`, CC BY 4.0, 270 MB).
- Prefer existing exports and pace requests to the production API.

## Cleaning

This is the core of the project. Layers run in order, and each can be switched off so its effect on the map can be measured.

| Layer | Removes | Catches |
|---|---|---|
| 1. Record evidence | `LIVING_SPECIMEN`, non-native establishment flags, fuzzy matches, kingdom clashes, multilingual text patterns in `locality`, `habitat`, `occurrenceremarks` (botanic garden, zoo, cultivated, captive, aquarium...) | Flagged garden, zoo and wrong-match records |
| 2. Taxon lists | Extinct species, and a short list of domesticated species and major crops, plus humans | Species with no meaningful wild range |
| 3. Authoritative ranges | Plants outside WCVP native regions; species x country pairs listed in GRIIS; species flagged in WRiMS, at the level WRiMS allows | Naturalised species, the biggest group that no flag or outlier test can see |
| 4. Land/sea consistency | A clearly terrestrial species' sea records fold into the land part of that cell, and vice versa. Records far from any coast on the wrong side are dropped. "Clearly" means COL `environment`, or at least 95% of its records on one side. Seabirds stay on both. Freshwater and brackish species aren't folded, as they are what set the Baltic and Black Sea apart in Costello | Coordinate rounding across the shoreline |
| 5. Rare-outlier heuristics | See below | Unflagged zoo and garden records, vagrants, misidentifications, wrong coordinates |

Test the text patterns on a sample first for false positives; they are part of the query, so they are final before the
download. Layer 3 removes by country and by TDWG region, which can follow political lines, so it is softened: a cell touching
a native region counts as native, sub-national GRIIS lists are used where they exist (Hawaii, Alaska), and its effect on
boundaries has its own gate under "Robustness runs".

### Rare-outlier heuristics

A single wrong record is weak evidence on its own. The heuristics judge it against the species' own distribution and against
how heavily the cell is sampled.

- **Detached satellite.** Group a species' occupied cells into connected patches (gap of 500 km or less). Drop a patch holding
  under 1% of the species' records and under 5% of its cells, unless it has at least 2 datasets and at least 3 years. Real
  disjunct ranges survive because they have repeated, independent records.
- **Under-reported where well sampled.** Reporting rate = the species' share of its group's records in that cell. If a cell
  has at least 1,000 records of the group, the species has at most 5 records there, and its rate is under 1% of its median
  rate across its core cells, drop it. A tiger in Berlin is exactly this. Poorly sampled cells are never judged this way,
  because there a single record may be all anyone has seen.
- **Thin long-range link.** A species with under 5 records whose cells are over 1,000 km apart stays out of the network, as
  it would otherwise tie distant cells together on almost no evidence.
- **Suspect cells.** Cells where a far higher share of species were removed than in neighbouring cells point to a zoo,
  herbarium or institutional coordinate. Listed in the report for inspection.
- **Tolerant endemism.** A species counts as endemic to a realm with at least 95% of its cells inside. Strict 100% is reported
  too, because a single leftover error would otherwise erase it.

**Guardrails**
- A species' core cells are never removed: its top cells holding 90% of records, or all of them under 20 records.
- The record-share heuristics don't apply below 20 records, so rare endemics aren't erased.
- Removals are counted per rule, group and country, so over-cleaning is visible.

### Calibration

The thresholds are checked against independent knowledge rather than guessed.

1. On a random half of WCVP plant species, with layer 3 off, measure two rates per threshold setting: the share of
   out-of-native presences removed, and the share of native presences lost.
2. Choose settings that lose at most 5% of native presences and remove as many out-of-native presences as possible.
3. Evaluate on the other half of plants, and on GRIIS animal entries as a cross-check outside plants.
4. Report the remaining error as an estimate.

Limits to state plainly: WCVP is coarse (TDWG level 3) and incomplete, plants may not generalise to other groups, and layer
3 is far stronger for plants than for animals. The per-group runs will show any effect of that imbalance.

Misidentifications can't be detected directly. They are covered by the heuristics and tolerant endemism, and their effect is
measured through the out-of-native rate.

## Network and clustering

- **Nodes:** cell parts and taxa at three ranks: species, genus and family. Links are presence only, so bird record volume
  doesn't dominate.
- **Method:** Infomap on the bipartite network, run hierarchically with a fixed seed. It finds groups of cells sharing
  restricted taxa, which is closest to "defined by endemism".
- **Evolutionary depth:** this follows Infomap Bioregions 2, which links ancestors to the cells where their descendants occur.
  Here a genus or family links to every cell holding one of its cleaned species. Those links are weighted by the spatial
  information they carry, as in that paper: full weight for a taxon in one cell, falling to zero for one in every cell. An
  endemic family therefore counts heavily, and a cosmopolitan one hardly at all.
- **Why ranks, not a dated tree:** no dated phylogeny covers all groups, but COL's classification does. Ranks are only a proxy,
  as a family of birds and a family of beetles aren't the same age. The phylogenetic robustness run tests how much that matters.
- **Side benefits:** that paper reports more robust regions under lumping and splitting and under sparse data. A
  misidentification within a genus also doesn't change genus or family presence, though that is an inference, not a finding of the paper.
- **Cleaning comes first:** genus and family presence is derived only from cleaned species presence, as one stray species would
  otherwise give its whole family a presence.
- **Cross-check and fallback:** beta-sim, averaged across the three ranks, with Ward and with UPGMA. Ward is there because the critique found average linkage
  gives unbalanced trees of mostly one-cell clusters. Beta-sim with Ward becomes the headline only if Infomap gives fewer than
  5 modules above the size floor, and the report says so if that happens. If that fails too, the result is that no supported
  partition was found, and that gets reported rather than forced.
- **Sparse cells:** a cell part with Good-Turing coverage under 0.5 or under 50 species first merges with its neighbours into a
  coarser hexagon (resolution 5, about 210k km²), assigned by cell centroid so no extra download is needed. If still sparse,
  it stays out of the network, takes the majority label of its nearest network cells of the same land/sea type, and is flagged
  as imputed.
- **Composition gate:** a cell part also stays out unless at least two major groups are sampled there (for example plants,
  insects and vertebrates on land; fish and invertebrates at sea).
- **Separate pieces:** a piece stays with its realm only if it shares restricted species with the rest of that realm. Otherwise
  it goes to the adjacent realm it shares most with. This also clears stray patches, which Costello did by hand.

## Choosing the headline level

**Candidates.** Start from Infomap's finest modules and merge them one pair at a time, always the pair sharing the most
taxa by beta-sim, across every rank in the network. That gives exactly one nested map per realm count from 60 down to 5. It avoids a huge candidate set of
uneven cuts, and guarantees maps with 19 and 30 realms for the comparison. Infomap's own levels are reported alongside.

**Size floor.** In every candidate, a unit under 1M km² merges into the realm it shares the most species with and stays
visible as a named sub-realm. Size therefore never disqualifies a candidate. The floor applies only to the new realms, not to
the existing schemes being compared.

**Selection**
1. **Eligible:** median bootstrap V-measure at least 0.8.
2. **Pick:** the highest area-weighted endemism excess over a null model, averaged equally across species, genus and family.
   The null is random partitions matched on realm count, realm areas and number of pieces, so scattered maps gain no advantage.
3. **Tie:** within 0.01, take the coarser map.
4. **If none is eligible:** take the most stable map and state that it fell short of the threshold.

### Exploring and confirming

Registering only after the final map exists wouldn't show that the rules came first. Splitting the taxa allows free
exploration and still ends in a real test, so registration moves to the freeze instead of the start.

1. **Split.** Whole families are assigned at random to halves A and B, balanced by major group and species count, with a fixed
   seed. Splitting individual species would put relatives of one endemic family in both halves and inflate the test. Species
   without a family are assigned by their next rank up. Nothing is clustered on B's taxa, and they aren't inspected, until
   step 4. The only use before that is scoring the cleaning calibration on B's plants, against WCVP ranges, with no realms
   involved.
2. **Explore on A.** Methods and thresholds may change and maps may be viewed until we are happy. Every change is logged in
   `DECISIONS.md` with its reason. Changes are justified by endemism excess, stability and calibration, never by resemblance to
   an existing scheme.
3. **Freeze.** The pipeline, thresholds, selection rules, border gate and the pass marks for step 4 go into
   `PREREGISTRATION.md`. The repo goes public, and Markus registers the file on OSF (about 5 minutes with his login). An OSF
   registration can't be edited afterwards, which a Git history can.
4. **Confirm once on B.** Do B's unseen lineages show endemism excess inside A's realms, at each rank, against the null? How
   closely does a map built from B alone match A's (V-measure)? This asks whether unrelated lineages share the same boundaries,
   which is what makes a realm more than the history of one group.
5. **Headline.** The frozen pipeline runs on all species. Any later change is reported as a deviation.

**Limit:** B has the same providers and sampling gaps as A, so this guards against fitting noise in the species, not against
sampling artefacts. The robustness runs cover those.

## Comparison with existing schemes

On the same cleaned species and cells:
- Realm endemism at species, genus and family rank (tolerant and strict; all species and well-recorded species only) for
  `realm`, `costello`, and the hierarchy cut with the same number of realms, so scheme size doesn't drive the difference.
- Area-weighted V-measure between each scheme and the new realms, land and sea separately.
- Boundaries found in both, and those found in only one.
- Whether any realm crosses the shoreline.
- Endemism by Costello's definition as well (unique to one realm, no tolerance), so their 42% is compared like with like.

## Robustness runs

All come from the same download, except resolution 6, which is the pilot download.

| Question | Run |
|---|---|
| Do errors and non-natives matter at this scale? | Raw, then each cleaning layer added in turn |
| Is it driven by sampling? | Equal-coverage subsampling per cell |
| Is it driven by one provider? | Without eBird; without iNaturalist |
| Does it depend on the taxa? | Birds, plants, insects, marine fish, others, each only over cells that group samples |
| Do long-distance migrants matter? | Without birds |
| Is it driven by what was sampled? | Whether a cell's group composition predicts its realm once location is accounted for |
| Do the sparse-cell rules matter? | Minimum of 5 species, the cut-off the critique inferred for Costello, against 50 |
| Does evolutionary depth matter? | Species only, species plus genus, and family only, against the three-rank headline |
| Are ranks a fair proxy for evolutionary depth? | Infomap Bioregions 2 on dated trees where they exist (birds, mammals, amphibians, seed plants), against the rank-based run on the same groups |
| Does it depend on scale? | Resolution 6 against 7 |
| Does it depend on the method? | Infomap against beta-sim/UPGMA |
| Do the cleaning settings matter? | Lenient and strict heuristic settings |
| Do cosmopolitan species matter? | Without the most widespread species |
| Is the licence limiting? | CC0 and CC BY records only |
| Are boundaries political? | The border gate below |

Per-cell confidence combines bootstrap co-assignment, coverage and the imputed flag.

**Border gate.** Some overlap with national borders is expected, as borders often follow rivers, mountains and coasts. Overlap
alone therefore isn't the test.
- A boundary segment is flagged as a possible artefact if it lies within one cell width of a national border, and it also
  disappears in the run without layer 3 or in the equal-coverage run. Flagged segments are marked in the outputs.
- Total overlap is also compared against randomly rotated boundaries, and reported.
- **If flagged segments exceed 20% of total boundary length:** when the run without layer 3 stays under 20%, the headline
  switches to that run and says so. When it doesn't, the map is published as qualified, with no claim that it is a headline
  result, and Markus decides what happens next.

## Steps

1. Create the private `gbif/biogeographic-realms` repo and the tracking issue, fix the family split seed, and write the
   synthetic planted-realm test.
2. Run the clustering code on Costello's public matrix: reproduce the critique's beta-sim/Ward result, then see what Infomap
   and the selection rules give on the same data.
3. Read which COL release GBIF production uses. Fetch its export and the other reference data, measure family and genus
   coverage, and test the text patterns.
4. **Needs Markus:** provide the GBIF login, kept in the session scratchpad only. Run the resolution 6 download, which mints a
   public DOI. SQL downloads may need enabling per account. If Markus's account isn't enabled, stop and ask rather than
   falling back to a service account. Measure rows, memory and Infomap runtime on this 48 GB, 14-core Mac.
5. Run the resolution 7 download if the pilot shows it fits. Otherwise stay at resolution 6 and log it.
6. Build the cleaning layers. Calibrate on half A's plants and score on half B's.
7. Explore on half A: cluster, build the hierarchy, select a level, iterate. Review maps with Markus.
8. **Needs Markus:** agree the freeze. Write `PREREGISTRATION.md`, make the repo public, register on OSF.
9. Confirm once on half B, then run the headline on all species.
10. Run comparisons and robustness runs, including the border gate.
11. Build outputs:
    - realms in GeoJSON and GeoPackage with sub-realms;
    - a per-cell table with labels, confidence and flags;
    - a per-realm list of endemic species, genera and families;
    - maps on an equal-area projection;
    - an HTML report opening with the COL release used (key and alias) and the GBIF download DOIs, linking COL ids to
      checklistbank.org, published as a private artifact.

    Realm polygons and results stay out of the repo until Markus decides on publication.

## Verification

- **Synthetic:** a planted-realm network is recovered. Injected zoo-style and vagrant records are removed, while injected
  narrow endemics survive.
- **Published result:** beta-sim/Ward on Costello's matrix, with cells under 5 species removed, reproduces the critique's map.
- **Hold-out:** half B passes the marks set at the freeze. A failure is reported as a failure, and half B can't be reused.
- **Known cases, checked after cleaning:** no tiger presence in Europe; potato removed outside its native regions; starling
  in North America removed by GRIIS; a Hawaiian endemic bird keeps its range; Lemuridae stays endemic to Madagascar despite
  zoo records.
- **Calibration:** held-out rates are reported, and the 5% native-loss limit holds.
- **Geometry:** areas add up to 510.06M km², no overlaps, no gaps.
- **Reproducibility:** rerunning from the cached download gives identical labels.
- **Plausibility, reported but not tuned:** where the Wallace line, Sahara and Isthmus of Panama boundaries fall.

## Risks

- **Memory:** edge count is unknown until the pilot. Resolution 6 is the fallback.
- **Infomap may degenerate** into one giant module and many tiny ones at this density. The pilot would reveal this, and the
  fallback rule above handles it.
- **Uneven cleaning:** alien lists and range data are far better for plants than for animals, particularly marine invertebrates.
- **Ranks aren't ages:** families differ in age and concept across groups, and COL's classification comes from many sources. Some
  boundaries may reflect taxonomic tradition rather than evolution. The phylogenetic run only tests this for well-studied groups.
- **Release change:** if GBIF production moves to a newer COL release between the pilot and the final download, the two can't be
  joined by id. Rerun the pilot query rather than mixing them.
- **Open ocean:** even after coarsening, much of it may be imputed rather than observed. The confidence map shows where.
- **Over-strict gates:** the composition gate and sparse-cell rules could exclude much of the tropics, where they are most needed
  to work. Exploration on half A measures how many cells they remove, per continent, before the thresholds are frozen.
- **Licence:** the whole-index download inherits CC BY-NC, and whether that binds polygons derived from it is an open question.
  Settle it with GBIF before publishing the polygons. If a CC BY gazetteer is wanted later, the open-only run is the candidate,
  and it may be a different map from the headline. That run measures how far they differ.
