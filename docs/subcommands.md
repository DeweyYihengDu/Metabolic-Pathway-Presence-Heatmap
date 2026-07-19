# Subcommands

`mpph <taxon> …` is a shortcut for `mpph run …`.

| Command | Input | Produces |
|---|---|---|
| `run` | taxon / codes / user | matrix, features, QC, figure, ordered order, trees, manifest |
| `traits` | taxon / codes / user | trait matrix + heatmap for a JSON/YAML panel |
| `explain` | taxon / codes / user | per-step matched/missing KO evidence for one module + organism |
| `compare` | a `run` output + metadata | differential features (Fisher / Mann-Whitney, FDR) + volcano plot |
| `pan` | a `run` output | core/soft-core/shell/cloud classes, prevalence + accumulation plots |
| `ordination` | a `run` output (+ metadata) | PCoA coords + scatter, PERMANOVA |
| `community` | taxon / codes / user | pairwise metabolic complementarity |
| `report` | a `run` output | self-contained interactive `report.html` |
| `validate` | a sample sheet | input checks |
| `enrich` | a study-set gene/KO list + background | hypergeometric enrichment table + bar chart |
| `gsea` | a ranked gene list, or an expression matrix + groups | ES/NES/leading-edge table + running-enrichment plot |

Examples:

```bash
mpph run Prochlorococcus --completeness --cluster --outdir results/
mpph pan        --results results/
mpph ordination --results results/ --metadata meta.tsv --color habitat
mpph compare    --results results/ --metadata meta.tsv \
                --group-column habitat --group-a surface --group-b deep
mpph explain    Prochlorococcus --module M00002 --organism pmt
mpph traits     Prochlorococcus --panel carbon_fixation --cluster
mpph community  --codes genomes.txt --max-combination-size 2
mpph report     --results results/
mpph enrich     --study unique.txt --background-organism pmt --ontology kegg-pathway
mpph gsea       --ranked-list ranked.tsv --ontology kegg-module
```

Built-in trait panels: `biogeochemistry`, `respiration`, `carbon_fixation`.
Supply your own with `--panel path/to/panel.json` (or `.yaml`); see the
`all_of` / `any_of` / `optional` step schema in `mpph/data/traits/`. A panel
file can declare `panel_version` (a version for *its own content*, not this
package) alongside the required `traits` map; `mpph traits` records the
panel's resolved path, content SHA-256 and any such metadata in
`<slug>_manifest.json`, so a later run can detect that a panel changed
underneath a previously-recorded analysis.

## `mpph enrich`

Hypergeometric over-representation (ORA): does a *study set* of genes/KOs
contain more members of a category than expected from the *background*
(universe) it was drawn from?

```bash
# KEGG pathway enrichment: study set vs. one organism's full KO complement
mpph enrich --study unique_genes.txt --background-organism pmt \
            --ontology kegg-pathway

# KEGG module enrichment: an explicit background file instead
mpph enrich --study degs.txt --background all_annotated_kos.txt \
            --ontology kegg-module

# GO enrichment: KEGG has no GO annotations, so supply your own mapping
# (a long table, or an eggNOG-mapper .annotations file's GO_terms column)
mpph enrich --study degs.txt --background all_genes.txt --ontology go \
            --gene-go-map sample.emapper.annotations --go-map-format eggnog \
            --go-names go_names.tsv
```

- `--study` / `--background`: plain id lists, one per line (`#` comments OK).
  `--background-organism CODE` is a shortcut for "this organism's full KO set"
  (kegg-pathway/kegg-module only — GO mapping is not organism-specific).
- Both study and background are restricted to genes *annotated in the category
  system being tested* before the test runs (an unannotated gene can't support
  or refute enrichment for any category); the printed summary reports how many
  inputs were dropped and why.
- One-sided (over-representation only): a zero-hit category is reported as
  non-significant, never as "significantly depleted".
- BH-FDR (`q_value`) is applied across every category actually tested (i.e.
  after `--min-category-size` filtering).
- `--gene-go-map` accepts a `gene<TAB>GO:#######` long table (repeat the gene
  for multiple terms, or comma/semicolon-separate them in one row) or an
  eggNOG-mapper `.annotations` file (`--go-map-format eggnog`, reads `GO_terms`).
  `--go-names` is an optional `GO:#######<TAB>name` table for readable labels;
  without it, categories are labelled by their raw GO id.

## `mpph gsea`

Rank-based enrichment (GSEAPreranked-style): every gene is scored and ranked,
then categories are tested for skewing toward either end of that ranking --
no need to first threshold a "study set" the way `enrich` does.

```bash
# Bring your own ranking (recommended): any gene<TAB>score file, e.g. a
# DESeq2/edgeR/limma statistic
mpph gsea --ranked-list deseq2_stat.tsv --ontology kegg-pathway

# Or let mpph rank genes from a raw expression matrix + two sample groups
mpph gsea --expression counts.tsv --metadata meta.tsv \
          --group-column condition --group-a treated --group-b control \
          --rank-metric signal2noise --ontology kegg-module

# GO, same bring-your-own-mapping rule as `enrich`
mpph gsea --ranked-list ranked.tsv --ontology go \
          --gene-go-map sample.emapper.annotations --go-map-format eggnog
```

- `--ranked-list`: a `gene<TAB>score` file (any statistic, higher = more
  associated with the "top" phenotype). `--expression` + `--metadata` +
  `--group-column` + `--group-a` + `--group-b` ranks genes automatically with
  `--rank-metric signal2noise` (the original GSEA paper's statistic) or
  `log2fc`; the computed ranking is also saved as `<label>_ranked_list.tsv`.
  These are simple, transparent statistics, **not a substitute for a dedicated
  DE tool** — prefer `--ranked-list` with your own DESeq2/edgeR/limma statistic
  for a rigorous analysis.
- Unlike `enrich`, **no background is needed and unannotated genes are not
  dropped** — every ranked gene counts as a "miss" in the running-sum statistic
  (dropping them would inflate the score by removing genuine background noise).
- `--min-size`/`--max-size` (default 15/500, standard GSEA defaults) filter
  categories by how many members actually appear in the ranking; KEGG modules
  are often small, so lower `--min-size` when using `--ontology kegg-module`.
- `--weight` is the score exponent in the weighted running-sum statistic (1.0 =
  standard GSEA weighting by |score|; 0 = unweighted KS statistic).
- Significance is by **gene-set permutation**: for each distinct category size,
  many random gene sets of that size are drawn from the same fixed ranking to
  build a null distribution (`--permutations`, default 1000). This is the
  standard "preranked" approach (as used by fgsea's simple mode) — weaker than
  phenotype permutation (which re-derives the ranking per permutation from the
  raw samples), which is not implemented.
- Output includes **leading-edge genes** per category — the genes actually
  driving the enrichment, up to (positive ES) or from (negative ES) the
  running-sum peak.
