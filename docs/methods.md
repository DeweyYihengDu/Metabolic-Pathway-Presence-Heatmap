# Methods

## Pathway-map association (presence mode)

For each organism, `list/pathway/<code>` reports which KEGG pathway maps have an
organism-specific version. MPPH records this as 1/0. **A 1 means the pathway has
some annotated KOs — not that the complete pathway is present or functional.**
By default only pathways under the BRITE top-level category `Metabolism` are
kept (`--top-category`, `--all-categories`), and aggregate "Global and overview
maps" are dropped (`--keep-overview`).

## Module completeness (completeness mode)

A KEGG module `DEFINITION` is a boolean expression over KO ids:

| token | meaning |
|---|---|
| space | AND (sequential steps) |
| `,` | OR (alternatives) |
| `+` | essential complex component (AND) |
| `-` | optional component (ignored in the score) |
| `--` | placeholder step, no assigned KO (excluded from num & denom) |
| `M#####` | nested module reference (resolved recursively, cycle-guarded) |

Completeness = satisfied top-level steps ÷ total real steps (0–1). Only
`Pathway`-type modules are scored by default (`--all-modules` includes
signature/reaction modules). `mpph explain` shows, per step, which KOs matched
and which are missing.

**States.** `--prevalence-state complete` counts a module toward prevalence only
when its completeness ≥ `--complete-threshold` (default 1.0), rather than merely
detectable. This distinguishes "fully complete in X% of organisms" from
"detectable in X%".

## Clustering and the dendrogram

Rows/columns are UPGMA-clustered (`--cluster`). Distance defaults to Jaccard for
binary presence and Euclidean for continuous completeness (`--metric` also
offers `dice`, `hamming`, `braycurtis`, `cosine`; Jaccard/Dice are rejected on
continuous data).

**The dendrogram reflects functional-profile similarity, not a sequence-based
phylogeny.** Horizontal gene transfer, niche convergence, gene loss, annotation
coverage and genome completeness all affect it. Use `mpph`'s Newick export with
`treecompare` (Robinson–Foulds, bootstrap clade support) to compare against a
reference tree before making evolutionary claims.

## Enrichment (over-representation analysis)

`mpph enrich` tests whether a *study set* of genes/KOs contains more members of
a category (KEGG pathway, KEGG module, or a GO term) than expected by chance
given the *background* it was drawn from — the classic hypergeometric ORA test
(equivalent to a one-sided Fisher's exact test), the same design used by
clusterProfiler / DAVID / topGO.

For each category with at least `--min-category-size` background members:

```
p = P(X >= k)   where X ~ Hypergeom(N, K, n)
N = |background genes annotated in this category system|
K = |background genes in this category|
n = |study genes annotated in this category system|
k = |study genes in this category|
```

This is **one-sided** (over-representation only): a study set with zero hits in
a category is reported as non-significant (p→1), never as "significantly
depleted". P-values are BH-FDR corrected (`q_value`) across every category
actually tested. Both study and background are restricted to genes annotated in
the category system being tested first — an unannotated input gene can neither
support nor refute enrichment.

**Category membership**: KEGG pathway/module membership comes from the global,
non-organism-specific `link/pathway/ko` and `link/module/ko` endpoints, so the
same category set applies regardless of which organism the study/background
came from. **GO enrichment needs a gene-to-GO mapping you supply** (a long
table, or an eggNOG-mapper `.annotations` file) — KEGG itself carries no GO
annotations.

## Rank-based enrichment (GSEA)

`mpph gsea` scores *every* gene in a ranking (typically by a differential
expression statistic) rather than thresholding a discrete study set first —
the "GSEAPreranked" design from Subramanian et al. 2005 (PNAS).

For a ranking of `n` genes and a category with `k` members present in it, the
weighted running-sum statistic walks the ranking from top to bottom, stepping
up by `|score|^weight / (sum of |score|^weight over hits)` at each category
member and down by `1/(n-k)` at each non-member. The enrichment score (ES) is
the maximum deviation from zero (signed) of that walk. `weight=1` (default)
matches standard GSEA weighting; `weight=0` gives the unweighted
Kolmogorov-Smirnov statistic.

**Significance** comes from gene-set permutation: since ES under the null
depends only on a random gene set's *size* (not its identity) for a fixed
ranking, many random gene sets of each tested size are drawn once and reused
across every category of that size. The normalized ES (NES) divides the
observed ES by the mean absolute null ES of the same sign; the p-value is the
fraction of same-signed null draws at least as extreme, BH-FDR corrected
(`q_value`) across every category tested. This is the standard "preranked"
approach — weaker than *phenotype* permutation (which re-derives the ranking
from the raw samples on every permutation), which is not implemented here.

**Permutation depth caps p-value resolution**: with `--permutations N`, the
smallest achievable p-value is `1/(N+1)`, so with the default 1000 permutations
several genuinely different top categories can be tied at the same floor
p-value (and therefore the same `q_value`) — the reported "#1 hit" may just be
whichever tied category happened first. If several categories cluster at the
permutation floor, raise `--permutations` (10000–20000) to actually resolve
their order rather than trusting the tie-break.

**Unlike `enrich`, no background is chosen and no gene is dropped**: every
gene in the ranking — annotated or not — contributes to the running sum as a
potential "miss". Restricting to an annotated-only universe (as ORA does)
would inflate the statistic by removing genuine background noise.

**Ranking from expression**: `--expression` + two groups computes
`signal2noise` (`(mean_a - mean_b) / (std_a + std_b)`, the original GSEA
paper's statistic) or `log2fc` per gene. These are simple and dependency-free,
**not a replacement for a proper differential-expression tool** — for a
rigorous analysis, rank by DESeq2/edgeR/limma's own statistic and pass it via
`--ranked-list`.

**`--expression` must already be normalized for library size** (CPM/TPM/FPKM,
or DESeq2/edgeR size-factor-normalized counts) — MPPH does not normalize for
you. Raw read counts differ between samples simply because they were
sequenced to different depths, and computing `signal2noise`/`log2fc` directly
on raw counts confounds that depth difference with real biology. This is a
real mistake made (and caught) while preparing this project's own example: a
first pass ranked E. coli genes from *raw* per-sample counts that varied 1.75x
in library size between samples; switching to the same dataset's CPM values
shifted which pathway came out on top.

## Quality control

Organisms with zero retained features are dropped (`--keep-empty` to keep) and
recorded in `*_qc.csv`; per-organism fetch failures are recorded, not fatal. A
missing feature in a MAG may reflect incomplete assembly rather than true
absence — import CheckM2/GTDB-Tk metadata and interpret low-completeness genomes
with care.
