# Output schema

## `run`

| File | Contents |
|---|---|
| `<slug>_matrix.csv` | organism × feature matrix (0/1 presence, or 0–1 completeness) |
| `<slug>_features.csv` | feature id → name → category → top category |
| `<slug>_qc.csv` | per-organism annotated-feature count and status |
| `<slug>_heatmap.<fmt>` | the figure |
| `<slug>_ordered_matrix.csv` | matrix in clustered row/column order |
| `<slug>_row_order.csv`, `<slug>_col_order.csv` | clustered order (traceable columns) |
| `<slug>_organism_tree.nwk`, `<slug>_feature_tree.nwk` | Newick trees (`--newick`) |
| `<slug>_manifest.json` | command, versions, KEGG release, filters, excluded organisms, outputs |

## Post-processing

| Command | Extra files |
|---|---|
| `pan` | `<slug>_pan_classes.csv`, `<slug>_accumulation.csv`, `<slug>_prevalence.<fmt>`, `<slug>_accumulation.<fmt>` |
| `compare` | `<slug>_differential_<A>_vs_<B>.csv` + volcano `<fmt>` |
| `ordination` | `<slug>_pcoa.csv`, `<slug>_pcoa_diagnostics.json`, `<slug>_pcoa.<fmt>`, `<slug>_permanova.csv` (with `--metadata`/`--color`) |
| `community` | `community_complementarity.csv` |
| `report` | `<slug>_report.html` (self-contained) |
| `traits` | `<slug>_matrix.csv`, `<slug>_heatmap.<fmt>` |
| `enrich` | `<label>_enrichment.csv`, `<label>_enrichment.<fmt>` |
| `gsea` | `<label>_gsea.csv`, `<label>_gsea_summary.<fmt>`, `<label>_gsea_top.<fmt>` (+ `<label>_ranked_list.tsv` with `--expression`) |

## Manifest fields

`mpph_version`, `command`, `python`, `platform`, `source`, `mode`, `metric`,
`match`, `clustered`, `generated_utc`, `kegg_release`, `n_organisms_matched`,
`n_organisms`, `n_features_before_filtering`, `n_features`, `organisms`,
`excluded_organisms` (with reason), `filters`, `outputs`.

The manifest is the reproducibility record: it captures the exact command, the
software/KEGG versions, and every filter applied.
