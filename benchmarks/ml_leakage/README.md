# Phylogenetic leakage in genome-based module prediction

Tests whether reported accuracy for predicting KEGG module completeness from
a genome's KO vector survives evaluation on phylogenetically held-out clades,
or is inflated by pseudoreplication between train and test sets drawn from a
clade-redundant database.

## Why this is worth testing

MetaPathPredict (eLife 2024) reports mean F1 **0.96**, trained and tested on
~30,596 RefSeq+GTDB genomes split **at random**. Those databases are dominated
by many strains of one species and many species of one genus, so a random
split routinely puts near-identical relatives on both sides of it.

Yu et al. (PLOS Biology 2025) showed exactly this inflates accuracy badly for
antimicrobial-resistance prediction — models learn lineage markers rather than
causal determinants — and recommended clade-held-out validation. **That
critique has not been applied to metabolic-module prediction.**

## What this does *not* claim

It does not re-run MetaPathPredict, and it is not a measurement of that tool.
The claim under test is about an **evaluation protocol**, so a single fixed
model (a plain gradient-boosted tree) is evaluated under two splitting schemes
on the same data. Reimplementing someone else's architecture and reporting a
lower number would confound the protocol question with implementation
differences.

## Design

| | |
|---|---|
| Genomes | 3,000 prokaryotes with KEGG-curated KO sets, 535 genera |
| Features | KO presence vector, down-sampled to 30–90% to simulate an incomplete genome (MetaPathPredict's own protocol) |
| Label | Is the module complete in the **full** KO set |
| Splits | random vs **genus-blocked** (no genus on both sides) |

**The mandatory control.** Blocking by genus also reduces training diversity,
so part of any drop is expected from that alone. Every blocked split is
therefore compared against a random split *at the same training-set size*.
Without it the experiment cannot separate "leakage inflated the number" from
"a smaller training set performs worse", and the result would be
uninterpretable.

## Data

`fetch_ko_matrix.py` pulls KO sets and taxonomy from KEGG REST directly —
no proteomes, no HMM search, no GTDB, which is what makes this a network-bound
job that runs alongside CPU-bound ones. Requests are serialised with a delay
and cached on disk, so a rerun is free and an interrupted run resumes.

```bash
python fetch_ko_matrix.py --n-organisms 3000 --out-dir data/
python leakage_test.py --data-dir data/ --n-modules 25 --out results/leakage.csv
```

Results are written by the run and interpreted in
`../../paper/nmeth/` once complete.
