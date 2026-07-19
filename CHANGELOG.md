# Changelog

## 3.5.1

Correctness hotfix for Newick tree export. Reproduced with a failing case
first, now has a regression test.

### Tree export
- `plot_matrix`'s returned linkage (`row_link`/`col_link`) indexes leaves in
  the *pre*-clustering row/column order, but the CLI's Newick export was
  passing it together with `row_labels`/`col_labels` -- the *post*-dendrogram,
  reordered display labels. The topology was correct but names could attach to
  the wrong leaves (e.g. two genuinely distant organisms shown as each
  other's closest relative, and vice versa). `plot_matrix` now also returns
  `row_link_labels`/`col_link_labels` (labels in the same order the linkage
  was computed on), and `organism_tree.nwk`/`feature_tree.nwk` export uses
  those instead.

## 3.5.0

New subcommand: **`mpph gsea`** — rank-based enrichment ("GSEAPreranked" style,
Subramanian et al. 2005), complementing `mpph enrich`'s hypergeometric ORA.
Every gene is scored and ranked (no threshold chosen to define a study set
first), and categories are tested for skewing toward either end of the ranking.

- `--ranked-list FILE`: bring your own `gene<TAB>score` ranking (e.g. a
  DESeq2/edgeR/limma statistic) -- the recommended, rigorous path.
- `--expression FILE --metadata FILE --group-column X --group-a A --group-b B`:
  rank genes automatically with `signal2noise` (the original GSEA paper's
  statistic) or `log2fc`. Explicitly documented as simple/transparent, **not a
  substitute for a dedicated DE tool**.
- `--ontology kegg-pathway` / `kegg-module` / `go`, same category-membership
  machinery as `enrich` (GO again needs a user-supplied gene-to-GO mapping).
- Weighted running-sum statistic (`--weight`, default 1.0 = standard GSEA
  weighting), significance by gene-set permutation batched per distinct
  category size (`--permutations`), NES, BH-FDR q-values, and leading-edge
  genes per category.
- Unlike `enrich`, **no background is required and no gene is dropped** --
  every ranked gene counts as a "miss" in the running sum (documented as a
  deliberate, important difference from ORA's annotated-universe restriction).
- Two new plots: `plot_gsea_running` (the classic running-enrichment-score
  plot with a hit rug and ranking-metric bars) and `plot_gsea_summary` (top
  hits by NES, coloured by enrichment direction and significance).
- New public API: `load_ranked_list`, `load_expression_matrix`,
  `rank_from_expression`, `enrichment_score`, `gsea_analysis`.
- 21 new tests: the running-sum statistic cross-checked on synthetic rankings
  with known top/bottom/scattered gene sets, leading-edge extraction, the
  expression-ranking statistics, and an offline CLI end-to-end test (KEGG
  pathway, KEGG module, GO, and expression-based ranking) served from
  committed cache fixtures.
- New example, computed live against real KEGG data: every KO across 25
  *Synechococcus* genomes ranked by (marine prevalence − freshwater
  prevalence), tested against KEGG modules. Top hit is "Sulfate-sulfur
  assimilation" (raw p≈0.001, correct ecological direction — marine water is
  sulfate-rich, freshwater is not) -- shown deliberately even though it does
  not clear q<0.05 across all 131 modules tested, as an honest example of
  reading GSEA output past a single threshold.
- docs/methods.md, docs/subcommands.md, docs/outputs.md, README updated.

## 3.4.0

New feature: **`mpph enrich`** — hypergeometric over-representation (ORA) of a
gene/KO study set against KEGG pathways, KEGG modules, or GO terms.

- `--ontology kegg-pathway` / `kegg-module`: category membership from the
  global (non-organism-specific) `link/pathway/ko` / `link/module/ko` KEGG
  endpoints, so results don't depend on which organism the study set came from.
  `--background-organism CODE` is a shortcut for "this organism's full KO
  complement" (built on the existing `organism_kos()`).
- `--ontology go`: GO enrichment via a **user-supplied** gene-to-GO mapping (a
  long table, or an eggNOG-mapper `.annotations` file's `GO_terms` column) --
  KEGG itself has no GO annotations, so this is bring-your-own-mapping, not a
  new external data source.
- One-sided hypergeometric test (over-representation only) with BH-FDR
  correction (reusing the existing `benjamini_hochberg`); both study and
  background are restricted to the annotated universe before testing, and
  dropped/used counts are reported.
- A ranked bar chart (`plot_enrichment`) of the top hits, coloured by
  significance.
- New public API: `fetch_ko_pathway_membership`, `fetch_ko_module_membership`,
  `fetch_pathway_names`, `hypergeometric_enrichment`, `invert_membership`,
  `load_gene_go_map`.
- 21 new tests: pure hypergeometric math (cross-checked directly against
  `scipy.stats.hypergeom`), GO/id-list file parsing, and an offline CLI
  end-to-end test served from committed KEGG cache fixtures. Also validated
  live against the real KEGG API (see the new example below).
- New example: KEGG pathway enrichment of the 190 genes unique to
  *Prochlorococcus* MIT 9313 (low-light ecotype) versus AS9601 (high-light,
  streamlined) — recovers "Photosynthesis - antenna proteins" as a top hit,
  the published low-light antenna-gene expansion in this genus.

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
