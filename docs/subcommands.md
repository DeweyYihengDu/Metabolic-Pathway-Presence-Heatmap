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
```

Built-in trait panels: `biogeochemistry`, `respiration`, `carbon_fixation`.
Supply your own with `--panel path/to/panel.json` (or `.yaml`); see the
`all_of` / `any_of` / `optional` step schema in `mpph/data/traits/`.
