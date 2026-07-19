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
```

Built-in trait panels: `biogeochemistry`, `respiration`, `carbon_fixation`.
Supply your own with `--panel path/to/panel.json` (or `.yaml`); see the
`all_of` / `any_of` / `optional` step schema in `mpph/data/traits/`.

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
