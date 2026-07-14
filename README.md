# Metabolic Pathway Presence Heatmap (MPPH)

[![CI](https://github.com/DeweyYihengDu/Metabolic-Pathway-Presence-Heatmap/actions/workflows/ci.yml/badge.svg)](https://github.com/DeweyYihengDu/Metabolic-Pathway-Presence-Heatmap/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2023.06.27.546232-b31b1b.svg)](https://doi.org/10.1101/2023.06.27.546232)

MPPH profiles KEGG metabolism across a set of genomes and renders a
category-aware heatmap with UPGMA dendrograms. It works two ways:

- **Pathway presence** — for every sequenced genome in a taxon, which KEGG
  pathways are present/absent.
- **Module completeness** — the fraction of each KEGG *module*'s reaction steps
  covered by an organism's KO repertoire (0–1), a far more quantitative signal
  than binary presence.

It runs on genomes **KEGG already has**, on an explicit **list of organism
codes**, or on **your own MAGs** via their KO annotations.

<p align="center">
  <img src="examples/Prochlorococcus_presence_heatmap.png" width="90%"
       alt="Clustered pathway-presence heatmap for the genus Prochlorococcus"><br>
  <em>Pathway presence — cells coloured by KEGG functional category.</em>
</p>
<p align="center">
  <img src="examples/Prochlorococcus_completeness_heatmap.png" width="90%"
       alt="Module-completeness heatmap for the genus Prochlorococcus"><br>
  <em>Module completeness — sequential ramp, category strip on top.</em>
</p>

## Features

- **Two analyses** — pathway presence/absence *or* KEGG module completeness.
- **Three input sources** — a taxon name, a file of organism codes, or your own
  KO annotations (KofamScan / eggNOG / any `K#####` list).
- **Category-aware figures** — a functional-category colour strip + legend,
  gridded cells, UPGMA dendrograms, italic organism labels.
- **Newick export** (`--newick`) — the pathway/module tree for iTOL / FigTree.
- **Informative-feature filtering** — `--drop-core`, `--min/--max-prevalence`,
  and automatic removal of aggregate "Global and overview maps".
- **Robust & reproducible** — on-disk KEGG caching, retries, a run
  `manifest.json` recording the KEGG release, filters, and organism list.

## Installation

```bash
git clone https://github.com/DeweyYihengDu/Metabolic-Pathway-Presence-Heatmap.git
cd Metabolic-Pathway-Presence-Heatmap
pip install .          # or: pip install -e ".[dev]" for development
```

This installs the `mpph` command. Requires Python 3.8+.

## Usage

```bash
# Pathway presence across a genus (clustered, drop pan-conserved pathways)
mpph Prochlorococcus --cluster --drop-core

# KEGG module completeness, focused on variable modules, export the tree
mpph Prochlorococcus --completeness --cluster --min-prevalence 0.1 \
     --max-prevalence 0.9 --newick

# An explicit set of KEGG organism codes
mpph --codes my_codes.txt --cluster

# Your own MAGs: a directory with one KO list per genome (implies completeness)
mpph --user annotations/ --cluster --format png pdf
```

`mpph` and `python -m mpph` are equivalent.

### Key options

| Option | Description |
| --- | --- |
| `taxon` | Taxon name to query (whole-word match; `--exact` for a full `Genus species`). |
| `--codes FILE` / `--user PATH` | Use organism codes, or your own KO annotations, instead of a taxon. |
| `--completeness` | Score KEGG module completeness (0–1) instead of pathway presence. |
| `--cluster` / `--metric` | UPGMA-cluster both axes; distance `euclidean\|jaccard\|dice\|hamming`. |
| `--newick` | Export the organism tree as `.nwk`. |
| `--drop-core`, `--min/--max-prevalence` | Filter uninformative features by prevalence. |
| `--format {pdf,png,svg}`, `--outdir DIR` | Figure format(s) and output directory. |
| `--refresh`, `--no-cache`, `--cache-dir` | Control KEGG response caching. |

Run `mpph --help` for the full list.

## Input for your own MAGs

`--user` accepts either:

- **a directory** with one file per genome/MAG (the file stem is the sample
  name); each file is scanned for KO ids, so KofamScan `--format mapper` output
  (`gene<TAB>K#####`) or a bare KO list both work; or
- **a single `sample<TAB>KO` table** (long form).

For eggNOG-mapper output, extract the `KEGG_ko` column into per-sample lists
first.

## Output

For a run on `<name>`, MPPH writes to the output directory:

| File | Contents |
| --- | --- |
| `<name>_matrix.csv` | Organism × feature matrix (0/1 presence, or 0–1 completeness). |
| `<name>_features.csv` | Feature id → name → functional category. |
| `<name>_heatmap.<fmt>` | The figure. |
| `<name>_tree.nwk` | Newick organism tree (with `--newick`). |
| `<name>_manifest.json` | Run parameters, KEGG release, filters, organism list. |

## How it works

1. Resolve organisms from a taxon (`list/genome`), a code file, or user KO files.
2. **Presence:** fetch each organism's pathway list (`list/pathway/<code>`).
   **Completeness:** fetch module definitions (`get/md:` in batches of 10) and
   each organism's KO set (`link/ko/<code>`), then score every module's
   top-level steps against the KO set.
3. Assemble and filter the matrix; drop aggregate overview maps.
4. Render the heatmap, clustering both axes when `--cluster` is set, and colour
   by KEGG functional category (`br08901` / module `CLASS`).

## Citation

> Y.-H. Du and J.-H. Mu, "Metabolic-Pathway-Presence-Heatmap (MPPH):
> Constructing phylogenetic trees based on metabolic pathways," bioRxiv, Jun.
> 2023. doi:[10.1101/2023.06.27.546232](https://doi.org/10.1101/2023.06.27.546232)

## License

Released under the [MIT License](LICENSE).
