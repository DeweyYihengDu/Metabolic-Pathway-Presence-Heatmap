# Changelog

## 3.1.0

### Scientific correctness
- **Module scorer**: nested module references (`M#####`) now resolve recursively
  (cycle-guarded) instead of always scoring 0; `--` placeholder steps are
  excluded from numerator and denominator; completeness scores only `Pathway`
  modules by default (`--all-modules` to include signature/reaction modules).
- **Presence mode** now restricts to the `Metabolism` BRITE top-level category
  by default (`--top-category`, `--all-categories`); the tool no longer mixes
  in Human-Diseases / Genetic-Information-Processing pathways silently.
- Distance metric now defaults per data type — Jaccard for binary presence,
  Euclidean for continuous completeness — and Jaccard/Dice are rejected on
  completeness data (added `braycurtis`, `cosine`).

### Reproducibility & QC
- Per-organism fetch errors are recorded rather than fatal; a `*_qc.csv` reports
  each organism's annotated-feature count and status.
- Clustered order is exported: `*_row_order.csv`, `*_col_order.csv`,
  `*_ordered_matrix.csv`, and `*_organism_tree.nwk` / `*_feature_tree.nwk`.
- Enriched `*_manifest.json`: full command, Python/platform, KEGG release,
  before/after feature counts, excluded organisms, and all filter settings.

### Usability
- `--match {word,prefix,exact,regex}` for taxon name matching (`--exact` kept as
  a deprecated alias for `prefix`).
- Input sources (taxon / `--codes` / `--user`) are now mutually exclusive.
- Two-object clustering works; row and column dendrograms are controlled
  independently. Feature ids are shown on the x-axis for small matrices.
- Stable category colours (fixed hues for common KEGG categories).
- Honest framing: figures/labels call binary results "pathway-map association".

### Housekeeping
- Removed the unused `seaborn` dependency; `requires-python >= 3.9`.
- Legacy prototype notebooks moved to `legacy/`.
- Added `CITATION.cff`, `DATA_SOURCES.md`, this changelog.

## 3.0.0
- Repackaged the single script into the installable `mpph` package with a
  console entry point; added KEGG module-completeness mode, user MAG/KO input,
  Newick export, on-disk caching, a pytest suite, and GitHub Actions CI.
