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

## Quality control

Organisms with zero retained features are dropped (`--keep-empty` to keep) and
recorded in `*_qc.csv`; per-organism fetch failures are recorded, not fatal. A
missing feature in a MAG may reflect incomplete assembly rather than true
absence — import CheckM2/GTDB-Tk metadata and interpret low-completeness genomes
with care.
