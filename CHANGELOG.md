# Changelog

## 3.3.1

Correctness and security hotfix. Every item below was reproduced with a failing
case first and now has a regression test.

### Missing data (NaN means *unknown*, never absent)
- `benjamini_hochberg`: a single NaN p-value no longer turns **every** q-value
  into NaN; only finite p-values are corrected and set the rank denominator.
- `differential_features`: NaN is dropped per feature and per group, so
  prevalences use a **known-only denominator**; reports
  `n_*_total/known/unknown` and a `warning`, and skips untestable features
  (`--min-known-per-group`) instead of silently treating unknown as absent.
- `pan_classify`: known-only prevalence + `n_present/n_known/n_unknown/
  known_fraction`; features below `min_known_fraction` are labelled
  `insufficient-data` rather than mis-classified as core/shell/cloud.

### Statistics
- `pcoa` now raises a clear error when no positive axis exists (all-zero
  distances) instead of returning an empty result that crashed `ordination`,
  and returns diagnostics (positive/negative eigenvalues, negative fraction).
- `permanova` excludes samples whose group label is missing (NaN was previously
  counted as its own group) and rejects groups with <2 samples; writes
  `*_permanova.csv`.
- `accumulation_curve` is documented as **permutation/rarefaction**, not
  bootstrap, and now reports SD and a 95% interval (plotted as a band).

### Security
- HTML report: the payload is embedded in a `<script type="application/json">`
  block with `<` escaped (no `</script>` breakout), NaN is serialised as null,
  and all user-controlled text is written via `textContent`/DOM construction
  instead of `innerHTML`. A CSP is set as defence in depth.

### Tree comparison
- `robinson_foulds` no longer loses a root-level singleton outgroup from the
  shared leaf set (identical `(a,(b,(c,d)))` reported 3 of 4 leaves); leaf sets
  are parsed from the tree and non-overlapping leaves are reported.

### Provenance
- Manifest records `prevalence_state`, `complete_threshold` and
  `score_semantics`, and uses the argv actually passed (correct for
  `main(argv=...)`).
- Examples regenerated at 3.3.1 with relative paths — the previous manifests
  leaked an absolute local build path.

## 3.3.0

- **Completeness states**: `--prevalence-state {any,complete}` +
  `--complete-threshold` — filter on how often a module is *fully complete*,
  not merely detectable.
- **More built-in trait panels**: `respiration` (terminal oxidases) and
  `carbon_fixation` (Calvin / rTCA / Wood-Ljungdahl / 3-HP), alongside
  `biogeochemistry`.
- **Analysis plots**: `compare` now writes a volcano plot; `pan` writes a
  pan-class prevalence histogram and a pan/core accumulation curve.
- **Trait figures** are labelled "trait completeness" / "trait category".
- **`docs/`**: methods, input formats, subcommands, and output schema.

## 3.2.0

Turns MPPH from a heatmap tool into a functional-profile analysis toolkit. The
CLI is now subcommand-based; `mpph <taxon> ...` still works as a shortcut for
`mpph run`.

### New analysis
- **Evidence-aware scoring** (`mpph.modules.evaluate_module`): per-step
  matched/missing KOs, complete/partial/absent/unknown state, unresolved
  references, parser status. Surfaced by **`mpph explain`**.
- **`mpph compare`** — differential features between two groups (Fisher exact /
  Mann-Whitney U, odds ratio / Cliff's delta, BH-FDR).
- **`mpph pan`** — core / soft-core / shell / cloud classes + bootstrap
  pan/core accumulation curve.
- **`mpph ordination`** — PCoA (classical MDS) with a scatter plot, and
  PERMANOVA when a grouping is supplied.
- **`mpph community`** — pairwise metabolic complementarity (modules completed
  by a union of organisms that neither completes alone).
- **`mpph traits`** — score JSON/YAML metabolic-trait panels; ships a built-in
  biogeochemistry marker panel (N/S/C/methane cycles, photosynthesis).
- **`mpph report`** — self-contained interactive HTML report (searchable colour
  heatmap, QC and manifest tabs).
- **`mpph validate`** — sample-sheet checks.

### New input / metadata
- `mpph.samplesheet`: reproducible sample-sheet input; CheckM2 and GTDB-Tk
  metadata import.
- `--input-format eggnog` reads the eggNOG-mapper `KEGG_ko` column.

### Public API
- The KEGG REST client and data-access functions are exported from the
  top-level package (see the README "Public Python API").

### Comparison / trees
- `mpph.treecompare`: Robinson-Foulds distance to a reference tree and
  feature-bootstrap clade support.

## 3.1.1
- Expose the KEGG connection interface as public API; derive the request
  User-Agent from the package version.

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
