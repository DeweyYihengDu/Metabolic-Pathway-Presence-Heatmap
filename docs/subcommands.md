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
| `pathmap` | a KEGG pathway/global-map id + two organism/user groups | KEGG-style network diagram, enzymes/reactions coloured by which group has them |
| `annotate` | a protein FASTA (+ a one-time KOfam database) | `gene<TAB>KO` mapper TSV + manifest, ready for `--user` |

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
mpph pathmap    --map ko00010 --codes-a marine.txt --label-a Marine \
                --codes-b freshwater.txt --label-b Freshwater
mpph annotate   --setup-db .mpph_kofam_db          # one-time, ~1.5 GB
mpph annotate   --fasta genome.faa --out genome_annotated.tsv
```

Built-in trait panels: `biogeochemistry`, `respiration`, `carbon_fixation`.
Supply your own with `--panel path/to/panel.json` (or `.yaml`); see the
`all_of` / `any_of` / `optional` step schema in `mpph/data/traits/`. A panel
file can declare `panel_version` (a version for *its own content*, not this
package) alongside the required `traits` map; `mpph traits` records the
panel's resolved path, content SHA-256 and any such metadata in
`<slug>_manifest.json`, so a later run can detect that a panel changed
underneath a previously-recorded analysis.

`mpph report --max-cells N` (default 50,000): above this many
organisms x features, the interactive grid is replaced with a summary --
a browser table that large (each cell its own styled DOM node) can hang the
page. The QC/manifest tabs and the underlying CSV/figure outputs are
unaffected either way.

## `mpph compare --tree`

Corrects the differential test for phylogenetic non-independence. Without it,
`compare` treats every genome as an independent sample; genomes that share a
recent ancestor share features because of that ancestry, which inflates
significance (measured Type I error 27.7% at a nominal 5% when the two groups
are two clades — see [methods](methods.md#phylogenetic-non-independence-compare---tree)).

```bash
mpph compare --results results/ --metadata meta.tsv \
             --group-column habitat --group-a surface --group-b deep \
             --tree gtdbtk_reference.nwk
```

- `--tree FILE` — a Newick reference phylogeny **with branch lengths**, whose
  tip labels match the matrix's organism names. It must be independent of the
  features being tested: a GTDB-Tk / 16S / marker-gene tree. mpph's own
  `*_organism_tree.nwk` is rejected (by topology, so a renamed copy is caught
  too) because it is derived from those very features.
- `--phylo-permutations N` (default 9999) — simulated replicates. A
  permutation q cannot fall below `n_features/(N+1)`, so 999 is too few for a
  few-hundred-feature matrix; `compare` warns if the floor still exceeds 0.05.
- `--seed N` (default 0) — same seed, same p-values.

Adds to the output CSV: `p_value_phylo` / `q_value_phylo`, `n_state_changes`
(minimum state changes on the tree — `1` means the feature arose once, so
association and shared ancestry cannot be separated), `phylo_confounded`,
`n_phylo_replicates_accepted`, `phylo_k_tolerance` and `phylo_warning`. The
uncorrected `p_value`/`q_value` stay in the CSV; the volcano plot switches to
the corrected q, with the axis labelled accordingly.

## `mpph enrich`

Hypergeometric over-representation (ORA): does a *study set* of genes/KOs
contain more members of a category than expected from the *background*
(universe) it was drawn from?

```bash
# KEGG pathway enrichment: study set vs. one organism's full KO complement
mpph enrich --study unique_genes.txt --background-organism pmt \
            --ontology kegg-pathway

# Restrict to metabolic pathways only (KEGG PATHWAY also covers genetic/
# environmental information processing, cellular processes, diseases, ...)
mpph enrich --study unique_genes.txt --background-organism pmt \
            --ontology kegg-pathway --top-category Metabolism

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
- `--top-category NAME` (`--ontology kegg-pathway` only, e.g. `Metabolism`):
  restricts the tested universe to one BRITE top-level category (the same
  `br08901` source `run --top-category` uses) — KEGG PATHWAY is not
  metabolism-only, it also covers Genetic Information Processing,
  Environmental Information Processing, Cellular Processes, Organismal
  Systems, Human Diseases and Drug Development, so a plain `--ontology
  kegg-pathway` run can surface any of these unless restricted. Unset by
  default (tests every category, unchanged prior behaviour); the output's
  `top_category` column is always present for `kegg-pathway`, so you can
  also just filter the CSV yourself after the fact without this flag.

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

# Metabolic pathways only, same --top-category as `enrich`
mpph gsea --ranked-list deseq2_stat.tsv --ontology kegg-pathway \
          --top-category Metabolism
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
- `--top-category NAME`: same meaning and same `top_category` output column
  as `enrich` (see above) — `--ontology kegg-pathway` only.

## `mpph pathmap`

Draws a KEGG pathway (e.g. `ko00010`, Glycolysis) or the global metabolic
network map (`ko01100` — every KEGG pathway map stitched into one diagram,
~3,800 reactions) using KEGG's own KGML layout, then colours each enzyme
(KO) node and the reaction(s) it catalyzes by whether it's present in
group A, group B, both, or neither — a comparative overlay, not KEGG's own
default static colouring.

```bash
# Two organism groups from KEGG codes files (one code per line each)
mpph pathmap --map ko00010 \
             --codes-a marine_codes.txt   --label-a Marine \
             --codes-b freshwater_codes.txt --label-b Freshwater

# The full global metabolic map, comparing your own MAG annotations
mpph pathmap --map ko01100 \
             --user-a marine_mags/ --label-a Marine \
             --user-b freshwater_mags/ --label-b Freshwater --format svg
```

- `--map` accepts `01100`, `ko01100`, `map01100`, or `path:ko01100` — all
  normalize to the KO-centric KGML KEGG actually serves (the bare `map#####`
  reference id has no KGML of its own).
- Each group is `--codes-X FILE` (organism codes, one per line; their KOs are
  unioned) or `--user-X PATH` (your own annotations, same formats as `--user`
  elsewhere) — exactly one of the two per group.
- An enzyme entry sometimes lists several alternative KOs (isozymes) under
  one node; it counts as present in a group if **any** of its KOs are, and a
  reaction catalyzed by more than one such entry takes the union of all of
  them — the same any-of-these-genes convention used elsewhere in this tool.
- Compound nodes and map-reference boxes (linking to other pathway maps) are
  drawn from KGML's own fixed layout for context; they are not coloured by
  your data. Only enzyme/reaction presence is comparative.
- `<slug>_pathmap_manifest.json` records both groups' sources/KO counts and
  the per-status enzyme node counts alongside the usual provenance fields.

## `mpph annotate`

Locally annotates a protein FASTA with KEGG KO ids -- KOfam HMM profiles
searched via `pyhmmer`, entirely offline after a one-time database download.
No external HMMER/KofamScan install needed. See
[Methods](methods.md#local-ko-annotation-annotate) for the scoring rule.

```bash
# One-time setup (~1.5 GB compressed, skipped if DIR already looks populated)
mpph annotate --setup-db .mpph_kofam_db

# Annotate a genome's proteins
mpph annotate --fasta genome.faa --kofam-db .mpph_kofam_db \
              --out genome_annotated.tsv

# Restrict to a known domain of life (much faster than the full database) --
# KEGG's own prokaryote.hal/eukaryote.hal work directly, or a bare KO list
mpph annotate --fasta genome.faa --ko-subset prokaryote.hal

# The output chains straight into --user, no adapter needed
mpph run --user genome_annotated.tsv --completeness
```

- Needs `pip install mpph[annotate]` (the optional `pyhmmer` dependency) --
  deliberately not part of the base install, so `pip install mpph` stays
  thin. Running `annotate` without it prints an actionable error instead of
  a traceback.
- `--fasta` and `--setup-db` are mutually exclusive: one call downloads the
  database, a separate call annotates.
- `--kofam-db` defaults to `.mpph_kofam_db` (mirroring `.mpph_cache`'s
  role for the KEGG REST cache); point it at a manually-arranged directory
  (`tar xzf profiles.tar.gz` + `gunzip ko_list.gz`) if you already have one.
- `--out` defaults to `<fasta stem>_annotated.tsv` in the current directory
  (not next to the input FASTA) -- a `gene<TAB>K#####` table, one row per
  significant assignment, plus a `<out stem>_manifest.json` alongside it.
- `--cpus` (default `0` = auto-detect all cores, matching pyhmmer's own
  default) controls search parallelism; searching the full ~20,000+ profile
  database against a genome's worth of proteins takes real time (minutes,
  not seconds) -- `--ko-subset` is the main lever for cutting that down.
- **`--cpus` is also the memory knob**, and it is the *only* thing peak
  memory really tracks -- not the size of your proteome. Measured on one
  bacterial genome: 1 thread 1.0 GB / 21 min, 4 threads 2.9 GB / 5 min,
  28 threads 11.8 GB / 2 min, with byte-identical output throughout. Parallel scaling is already well past linear
  by 28 (4 threads is a perfect 4.0x speed-up; 28 threads only 10x), so on a
  memory-constrained machine lowering `--cpus` costs far less time than it
  saves memory. This tool does not need a large-memory node.
- **Any organism works, not just microbes.** Measured F1 against KEGG's own
  assignments: 0.89-0.93 across bacteria and archaea, 0.911 for
  *S. cerevisiae*, 0.867 for *A. thaliana*. Scale is the only practical
  difference -- *Arabidopsis*' 48,265 proteins take ~35 min on 28 threads
  versus ~3 min for a bacterial genome. See
  [Methods](methods.md#local-ko-annotation-annotate) for what changes with a
  eukaryotic input (lower KO coverage, splice isoforms) and
  `benchmarks/annotation/` for the full comparison.
- `--sequence-loading {auto,prefetch,stream}` chooses whether the proteome is
  held in memory or streamed. Both give identical output; the difference is
  ~1 kB per protein, so this only matters for metagenome-scale protein
  catalogues in the millions. `auto` (default) prefetches up to 1,000,000
  proteins, which covers every single-organism proteome.
