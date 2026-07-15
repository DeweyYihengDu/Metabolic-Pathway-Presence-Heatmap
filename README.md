# Metabolic Pathway Presence Heatmap (MPPH)

[![CI](https://github.com/DeweyYihengDu/Metabolic-Pathway-Presence-Heatmap/actions/workflows/ci.yml/badge.svg)](https://github.com/DeweyYihengDu/Metabolic-Pathway-Presence-Heatmap/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2023.06.27.546232-b31b1b.svg)](https://doi.org/10.1101/2023.06.27.546232)

MPPH profiles KEGG metabolism across a set of genomes and renders a
category-aware heatmap with UPGMA dendrograms. It works two ways:

- **Pathway-map association** — for each genome, whether KEGG provides an
  organism-specific version of a pathway map. This means the pathway *has some
  annotated KOs*, **not** that the complete pathway is present. By default only
  `Metabolism` top-level pathways are kept.
- **Module completeness** — the fraction of each KEGG *module*'s reaction steps
  covered by an organism's KO repertoire (0–1), a far more quantitative signal
  and the recommended mode for metabolic-capability analysis.

It runs on genomes **KEGG already has**, on an explicit **list of organism
codes**, or on **your own MAGs** via their KO annotations. The clustered
dendrogram reflects **functional-profile similarity**, and should not be read as
a sequence-based species phylogeny without external validation.

<p align="center">
  <img src="examples/Prochlorococcus_presence_heatmap.png" width="90%"
       alt="Clustered pathway-presence heatmap for the genus Prochlorococcus"><br>
  <em>Pathway-map association, <code>Prochlorococcus</code> (15 genomes) — cells
  coloured by KEGG functional category, UPGMA-clustered on Jaccard distance.</em>
</p>
<p align="center">
  <img src="examples/Synechococcus_completeness_heatmap.png" width="90%"
       alt="Module-completeness heatmap for the genus Synechococcus"><br>
  <em>Module completeness, <code>Synechococcus</code> (25 genomes) — the
  dendrogram separates several thermophilic and streamlined strains (top four
  rows) from the rest, visible as a near-empty completeness profile.</em>
</p>

<p align="center">
  <img src="examples/Prochlorococcus_traits_heatmap.png" width="90%"
       alt="Biogeochemical-trait completeness for the genus Prochlorococcus"><br>
  <em>Trait panel (<code>mpph traits</code>) — every <code>Prochlorococcus</code>
  scores complete carbon fixation (RuBisCO) and oxygenic photosynthesis, and
  (correctly) lacks nitrogen fixation, nitrification and methanogenesis.</em>
</p>

More examples, including the matching `Synechococcus` pathway-presence figure
and every underlying data/QC/tree file, are in [`examples/`](examples/).

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

This installs the `mpph` command. Requires Python 3.9+.

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

### Subcommands

`mpph <taxon> …` is a shortcut for `mpph run …`. The full toolkit:

```bash
# Score a metabolic-trait panel (built-in biogeochemistry, or your own JSON/YAML)
mpph traits Prochlorococcus --panel biogeochemistry --cluster

# Explain one module's evidence for one organism (matched / missing KOs per step)
mpph explain Prochlorococcus --module M00002 --organism pmt

# --- these post-process a previous `run --outdir results/` ---
mpph pan        --results results/                       # core/shell/cloud
mpph ordination --results results/ --metadata meta.tsv --color habitat   # PCoA + PERMANOVA
mpph compare    --results results/ --metadata meta.tsv \
                --group-column habitat --group-a surface --group-b deep   # differential
mpph report     --results results/                       # interactive report.html

# Metabolic complementarity across organisms (union completes a module)
mpph community  --codes genomes.txt --max-combination-size 2

# Validate a sample sheet before a user-data run
mpph validate   --samples samples.tsv
```

### Key options

| Option | Description |
| --- | --- |
| `taxon` | Taxon name to query. `--match {word,prefix,exact,regex}` (default `word`). |
| `--codes FILE` / `--user PATH` | Use organism codes, or your own KO annotations (mutually exclusive with a taxon). |
| `--completeness` / `--all-modules` | Score module completeness (0–1); by default only `Pathway` modules. |
| `--top-category` / `--all-categories` | Restrict presence mode to a BRITE top-level category (default `Metabolism`). |
| `--cluster` / `--metric` | UPGMA-cluster both axes. Metric defaults to `jaccard` (presence) / `euclidean` (completeness); also `dice`, `hamming`, `braycurtis`, `cosine`. |
| `--newick` | Export organism and feature trees as `.nwk`. |
| `--drop-core`, `--min/--max-prevalence`, `--present-threshold` | Filter uninformative features by prevalence. |
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
| `<name>_features.csv` | Feature id → name → functional category → top category. |
| `<name>_qc.csv` | Per-organism annotated-feature count and status. |
| `<name>_heatmap.<fmt>` | The figure. |
| `<name>_ordered_matrix.csv`, `<name>_row_order.csv`, `<name>_col_order.csv` | Clustered order, so every column in the figure is traceable. |
| `<name>_organism_tree.nwk`, `<name>_feature_tree.nwk` | Newick trees (with `--newick`). |
| `<name>_manifest.json` | Full command, KEGG release, filters, excluded organisms, output list. |

## How it works

1. Resolve organisms from a taxon (`list/genome`), a code file, or user KO files.
2. **Presence:** fetch each organism's pathway list (`list/pathway/<code>`).
   **Completeness:** fetch module definitions (`get/md:` in batches of 10) and
   each organism's KO set (`link/ko/<code>`), then score each module — top-level
   steps split on spaces, `,`=OR, `+`=AND, `-`=optional, `--`=gap (excluded),
   and nested `M#####` references resolved recursively.
3. Assemble and filter the matrix (top-level category, overview maps, empty
   genomes, prevalence).
4. Render the heatmap, clustering when `--cluster` is set, coloured by KEGG
   functional category (`br08901` / module `CLASS`).

## Public Python API (KEGG interface)

The KEGG REST client and the data-access functions are a **public, open part of
the package** — nothing is gated. They always connect directly to the official
KEGG endpoint (`https://rest.kegg.jp`), so you can reuse them to talk to KEGG
programmatically:

```python
import mpph

session = mpph.make_session()            # pooled, retrying HTTPS session

# Low-level: any KEGG REST endpoint, with optional on-disk caching
text = mpph.kegg_get(session, "list/genome", cache_dir=None)
print(mpph.KEGG_API_BASE)                # https://rest.kegg.jp

# Higher-level helpers
genomes = mpph.list_genomes(session, None)
vibrio  = mpph.select_by_taxon(genomes, "Vibrio", match="word")
paths   = mpph.get_pathways(session, "vch", None)     # {map_id: name}
kos     = mpph.organism_kos(session, "vch", None)     # {K#####, ...}
score   = mpph.module_completeness("K00844 K12407 K00845", kos)
```

`mpph.KEGG_API_BASE`, `make_session`, `kegg_get`, `kegg_release`,
`list_genomes`, `select_by_taxon`, `select_by_codes`, `get_pathways`,
`fetch_pathway_categories`, `list_modules`, `fetch_module_definitions`,
`organism_kos`, and the matrix/plot helpers are all exported from the top-level
`mpph` package. (Note the [KEGG terms](DATA_SOURCES.md) apply to the *data* you
retrieve; the client code is MIT-licensed.)

## Limitations

- **`pathway-map` ≠ complete pathway.** A `1` means KEGG lists an
  organism-specific map (some KOs annotated), not that the whole pathway is
  functional. Use `--completeness` for capability analysis.
- **Functional similarity ≠ phylogeny.** The dendrogram clusters KEGG functional
  profiles; horizontal transfer, niche convergence, gene loss, annotation
  coverage and genome completeness all affect it. Validate against a reference
  tree before making evolutionary claims.
- **Annotation depth varies** between genomes, and a MAG's missing features may
  reflect incomplete assembly rather than true absence.
- **Taxon matching is by organism name**, so it is reliable for genus/species
  but not for family/phylum. Results depend on the KEGG release (recorded in the
  manifest).

## Citation

> Y.-H. Du and J.-H. Mu, "Metabolic-Pathway-Presence-Heatmap (MPPH):
> Constructing phylogenetic trees based on metabolic pathways," bioRxiv, Jun.
> 2023. doi:[10.1101/2023.06.27.546232](https://doi.org/10.1101/2023.06.27.546232)

## License

Released under the [MIT License](LICENSE).
