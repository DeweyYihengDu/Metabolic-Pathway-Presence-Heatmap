# Benchmarks

Head-to-head validation of `mpph` against the established tools it overlaps
with, on real data, with every reference and parameter held constant so the
numbers measure the implementations rather than their inputs.

| Benchmark | Compares | Result |
|---|---|---|
| [`annotation/`](annotation/) | `mpph annotate` vs **KofamScan** vs **eggNOG-mapper**, 8 genomes across 3 difficulty tiers | `mpph` reproduces KofamScan at Jaccard **1.000** on 6 of 8 genomes (0.9997, 0.9994 on the other two); F1 vs KEGG's own KO assignments 0.842–0.934 |
| [`enrichment/`](enrichment/) | `mpph gsea`/`mpph enrich` vs **clusterProfiler** | GSEA ES Spearman ρ = **1.000**, NES ρ = 0.999; ORA p-value ρ = **1.000** |

Figures: [`figures/`](figures/) — regenerate with `python make_figures.py`
(reads only the committed result tables; no reference databases needed).

## What is and isn't committed here

Committed: the scripts, their tests, the summary result tables, and the
figures — everything needed to check the analysis or redraw the figures.

Not committed: the reference databases (KOfam ~8.7 GB extracted, eggNOG
~48 GB), the container images, and the per-gene raw tool outputs. Those are
reproducible from the documented accessions, container tags and commands in
each sub-README. The genomes used are identified by assembly accession
rather than vendored.

## Honest summary of what these benchmarks do and do not establish

**Do:** that `mpph annotate`'s pyhmmer reimplementation is call-for-call
equivalent to KofamScan, including on genomes KEGG never hand-curated and on
real MAGs with no reference at all; that `mpph`'s GSEA and ORA statistics
match clusterProfiler's on identical inputs.

**Do not:** establish that any of these tools' KO calls are biologically
correct — accuracy here is agreement with KEGG's own assignments, which are
themselves derived (by sequence-similarity best-hit search plus curation, a
method independent of profile-HMM scoring; see `annotation/README.md`).
Nor do they cover eukaryotic genomes, viral sequences, or ontologies other
than KEGG.

**Against `mpph`:** peak memory in `mpph annotate` is 10–17× KofamScan's and
scales erratically (2.3–12.3 GB vs 0.13–0.72 GB). Reported in full in
`annotation/README.md` rather than omitted.

## Running the tests

The comparison logic is tested against hand-constructed inputs with known
expected values, so a wrong table is caught before it is believed:

```bash
pytest benchmarks/ -q
```

These are deliberately outside the package's own `testpaths`, so `pytest -q`
at the repo root does not pick them up.
