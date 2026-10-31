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
  stability threshold (0.785 against 0.8). The design's fallback doesn't trigger, because 5 units aren't "fewer than 5".
- **The usable scale lies somewhere else.** A scan of Infomap settings on this matrix, as (leaf modules, share of the largest):
  species only, default (26, 0.57); species and genus, default (27, 0.73); two-level (16, 0.92 with genus); Markov time 0.75
  (372, 0.33 species only; 354, 0.36 with genus); Markov time 0.5 (over 1,600, 0.01). Nothing here was tried between 0.75 and 1.
- **Adding genus makes modules coarser.** Infomap Bioregions 2 warns that ancestral nodes collapse modular structure, and that
  shows here. Evolutionary depth may need a lower weight than the paper's whole-tree formula, or integration at one rank.
- To settle while exploring half A, using only endemism excess and stability: a Markov time sweep, or recursing Infomap inside
  large modules. Runtime and memory aren't a concern at this size, with 110 s and 7.5 GB including 20 bootstraps.
- The area-weighted endemism score of a map dominated by one huge realm is high by construction (0.678 here), since nearly every
  species occurs in it and most occur nowhere else. Only the excess over the null, which has the same areas, is meaningful.
