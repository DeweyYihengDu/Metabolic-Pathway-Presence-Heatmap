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

Completeness = satisfied top-level steps ÷ total *determinable* real steps
(0–1). Only `Pathway`-type modules are scored by default (`--all-modules`
includes signature/reaction modules) -- a `Pathway` module's definition can
reference a nested `M#####` that is itself signature/reaction-typed and so
isn't in the default scoring set. `mpph explain` shows, per step, which KOs
matched and which are missing.

**A step whose truth value can't be determined** -- a nested module reference
with no definition available (not fetched, filtered out by the
`Pathway`-only default, or a cyclic reference) or syntax this parser doesn't
handle -- is **unknown, not "confirmed absent"**: it is excluded from both the
numerator and denominator, the same treatment as a `--` placeholder step.
This uses two-valued (Kleene) evaluation: an undetermined token is
substituted with both 0 and 1; if the step's truth value doesn't change
either way (e.g. an OR already satisfied by a KO you have), it genuinely
doesn't depend on the undetermined part and the determined result is used
instead of discarding the step. If *every* real step in a module turns out
undetermined, `module_completeness` returns `NaN` rather than a fabricated
`0.0` -- NaN means unknown throughout this tool, never "confirmed absent".
`mpph explain` marks such steps `??` (not `--`) and lists them under
"unresolved module refs". A prevalence filter (`--min-prevalence` etc.)
computes prevalence among *known* organisms only, so a module that is
undetermined for some organisms isn't silently treated as absent in them.

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

**Newick labels needing special characters are quoted, not stripped**: a
label containing `(),:;[]'` is wrapped in single quotes (Newick's own quoting
rule) rather than having those characters deleted, so two different labels
(e.g. the microbiology convention `[Eubacterium] rectale` vs. a plain
`Eubacterium rectale`) can't collide into the same exported name. `treecompare`
reads quoted labels correctly, including a reference tree from another tool
that also uses standard Newick quoting. `robinson_foulds` compares clades as
**rooted** descendant-sets (matching a UPGMA dendrogram and a typically
outgroup-rooted reference tree) and reports the exact leaves unique to each
tree (`only_in_a`/`only_in_b`), not just their counts.

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
came from — this also means `--ontology kegg-pathway`/`mpph gsea`'s enrichment
already works for any species (any organism whose genes/proteins you can map to
KO ids, e.g. via `mpph annotate`), not only ones with a curated model-organism
database. **KEGG PATHWAY is not metabolism-only** — the same database also
covers Genetic Information Processing, Environmental Information Processing,
Cellular Processes, Organismal Systems, Human Diseases and Drug Development, so
an unrestricted `--ontology kegg-pathway` run can surface a significant hit
from any of these, not just metabolic pathways. `--top-category NAME` (e.g.
`Metabolism`) restricts the tested universe to one BRITE top-level category
(the same `br08901` source `run --top-category` already uses for the presence
heatmap); unset by default, matching this feature's behaviour before
`--top-category` existed. The output always includes a `top_category` column
when `--ontology kegg-pathway`, whether or not the flag is used, so you can
also just filter the CSV yourself afterward. **GO enrichment needs a
gene-to-GO mapping you supply** (a long
table, or an eggNOG-mapper `.annotations` file) — KEGG itself carries no GO
annotations. **This is flat GO term over-representation**: a gene annotated to
a specific child term counts only toward that term, not its ancestors, unless
your own mapping file already lists the ancestor terms explicitly. Tools that
do GO-DAG ancestor propagation (topGO, clusterProfiler with an
organism `.db`) will attribute more genes to a general parent term (e.g. "DNA
metabolic process") than this will from the same raw child-term annotations —
that is a difference in method, not a bug in either.

**Effect size**: `odds_ratio` (with a 95% CI, Haldane-Anscombe corrected for a
zero cell) is the standard 2x2-table effect size for this test and is more
comparable across categories of very different size than `fold_enrichment`
(a ratio of rates, which can look arbitrarily large for a tiny category with
a single lucky hit). A KO id with a KEGG `ko:` namespace prefix in your study
or background list (e.g. copied from `link/ko/<org>` or an eggNOG-mapper
`KEGG_ko` column) is normalized to the bare id automatically.

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
missing feature in a MAG may reflect incomplete assembly, a contig break, a
failed gene call, or missing KO annotation rather than true absence — not
just a caveat in this text, but something `run --qc-metadata` acts on:

- `--qc-metadata FILE`: a table keyed by `sample_id`/organism with a
  `completeness` column (optionally `contamination`, `taxonomy`) — e.g.
  `mpph.samplesheet.import_checkm2()` / `import_gtdbtk()` output saved to
  TSV. Joined into `*_qc.csv` (plus a `quality_tier` column, the MIMAG
  high/medium/low tiers from completeness+contamination only — the full
  MIMAG standard also needs rRNA/tRNA evidence this doesn't assess) and
  recorded in the manifest.
- `--min-genome-completeness N`: drops organisms below `N`% completeness
  before scoring (needs `--qc-metadata`); excluded organisms are reported the
  same way as `--keep-empty`'s zero-feature exclusions.
- Any *included* organism below 90% completeness gets a printed warning,
  whether or not `--min-genome-completeness` was set — a reminder that its
  apparent feature absences may be assembly gaps.

**This does not, and should not, correct scores for completeness** (e.g.
`observed / genome_completeness`) — different functional genes are not lost
uniformly at random as assembly quality drops, so a linear correction would
fabricate precision the data doesn't support. Filtering or flagging is the
honest option; silently "fixing" the number is not.

**`NaN` (unknown / not assessed) always gets its own colour on the heatmap,
distinct from both "confirmed absent" and any real value** — it is never
folded into the same visual as absent, and a legend entry ("unknown / not
assessed") appears whenever the plotted matrix actually contains one. This
applies to both `mode='presence'` (a `NaN` cell used to compare equal to
"absent" via `NaN > 0 == False`) and `mode='completeness'` (an undetermined
`module_completeness` result used to render as opaque black, since the
colormap's transparent "bad" colour lost its alpha channel when only RGB was
kept) — both are display bugs the value itself never had.

**`mpph compare`/`mpph pan` treat `NaN` the same way by default — excluded,
never counted as absent — but `--unknown-policy` makes this an explicit,
overridable choice** instead of a hardcoded one:

- `exclude` (default): `differential_features` tests, and `pan_classify`
  computes prevalence over, known values only — a known-only denominator.
- `absent`: counts `NaN` as absent instead of dropping it (the full sample/
  organism count becomes the denominator). This is never the default,
  because it may misrepresent an assembly or annotation gap as genuine
  absence — pick it only when there's a specific reason to believe missing
  really does mean absent here (e.g. a well-assembled, well-annotated
  genome where a true negative is far more likely than a missed call).
- `error`: refuses to run at all if any `NaN` is present among the compared
  groups (or, for `pan`, anywhere in the matrix) — forces missingness to be
  resolved first, e.g. with `--min-genome-completeness` QC filtering.

Whichever policy is active, `n_known`/`n_unknown`/`known_fraction` (and the
`--min-known-samples`/`--min-known-fraction` thresholds that gate
`insufficient_known_values`/`insufficient-data`) always reflect the *true*
evidence, never the policy's substitution — a feature with too little real
data is never waved through just because `absent` would otherwise happily
fill the gaps with zeros. `pan_classify`'s output also reports `n_absent`
alongside `n_present`, the complement under whichever denominator the
active policy selects.

## KGML pathway diagrams (`pathmap`)

`mpph pathmap` draws a KEGG pathway or the global metabolic map from KEGG's
own KGML data: compound and enzyme nodes at KEGG's own fixed layout
coordinates, connected by lines for each reaction's substrate(s)/product(s).
This is **not** KEGG's own default colouring (which shows a static reference
map) — it's a comparative overlay: each enzyme (KO) node, and every reaction
it catalyzes, is coloured by whether the KO is present in group A, group B,
both, or neither, the same "is this step present" question the rest of this
tool asks via heatmaps.

- **KGML is only served KO-centric** (`ko#####`), not the bare reference map
  (`map#####`, which 404s on KEGG's own `kgml` endpoint) — `--map` accepts
  any common spelling and normalizes it.
- **Isozymes**: KEGG sometimes lists several alternative KOs under one
  enzyme entry (`ko:K00001 ko:K00002 ko:K00003`). An entry counts as present
  in a group if *any* of its KOs are, matching the any-of-these-genes
  convention `mpph explain`/module scoring already use for alternative steps.
  A reaction catalyzed by more than one such entry (KEGG occasionally splits
  isozymes into separate entries pointing at the same reaction) takes the
  **union** of KOs across all of them.
- **Every substrate-product pair of a reaction is drawn**, not just a
  "primary chain" guess — KGML's own substrate/product order isn't a
  reliable primary-vs-cofactor signal, and some real steps (e.g. aldolase,
  one substrate splitting into two products) would lose a genuine branch
  under a "first substrate to first product only" simplification. This does
  mean a reaction with a cofactor pair (e.g. ATP/ADP) draws a few more lines
  than KEGG's own hand-tuned image shows for the same step.
- **Compound nodes and map-reference boxes are context, not data** — a
  compound is a shared intermediate, not an organism-specific feature, and a
  map-reference box just links to another pathway's own diagram.
- The global map (`ko01100`, ~7,600 nodes / ~3,800 reactions) renders in a
  few seconds; there is no size guard the way `report --max-cells` has one,
  since this is a single fixed-layout figure rather than a per-cell DOM table.

## Local KO annotation (`annotate`)

`mpph annotate` assigns KEGG KO ids to a protein FASTA using KEGG's own
KOfam HMM profile database (the same profiles KEGG uses to annotate genomes
itself) searched via [pyhmmer](https://pyhmmer.readthedocs.io/), reimplementing
KofamScan's own scoring logic rather than shelling out to compiled HMMER
binaries + KofamScan's Ruby wrapper. This generalizes across the tree of
life the way KOfam itself does — unlike tools that transfer annotation from
a fixed set of reference species (accuracy degrades with evolutionary
distance from whichever reference is closest), or tools of the KOBAS type
that ship curated backgrounds for human, mouse and a short list of other
model organisms, a profile built from many species' worth of a KO's
sequences has no such "nearest reference" problem.

**Scope**: KO assignment only, not GO — GO would need eggNOG's differently
structured per-clade databases, a separate addition. Input must already be
gene-called/translated protein sequences; `mpph annotate` does no ORF
prediction of its own (e.g. run Prodigal first for raw contigs).

**Domain of life is not a constraint.** Measured against KEGG's own KO
assignments (`benchmarks/annotation/`), F1 is 0.89–0.93 for bacteria and
archaea, 0.911 for *S. cerevisiae* and 0.867 for *A. thaliana*; the
48,265-protein *Arabidopsis* proteome runs in ~35 min at 9.7 GB peak RSS on
28 threads. Two things do differ when the input is a eukaryote, and neither
is a defect:

- **Coverage falls** — 76% of *E. coli* proteins receive a KO, 62% of yeast,
  24% of *Arabidopsis*, for every tool tested alike. KO describes metabolism
  and core cellular processes, not a plant's full protein complement. Recall
  against what KEGG *does* assign stays high (0.925 for *Arabidopsis*, the
  best of any genome benchmarked).
- **Isoforms inflate per-gene counts** — a eukaryotic proteome lists every
  splice variant. Per-genome KO *sets*, which is what MPPH's matrices are
  built from, are unaffected: isoforms of one gene collapse to the same KO.
  Per-gene call counts are not comparable to a prokaryotic genome's.

**Competing KO calls on one gene (`--multi-ko-policy`).** KOfam's thresholds
are fitted **per KO independently**, each maximising its own F-measure; nothing
in that procedure makes KOs compete with one another. A protein matching
several related profiles — paralogous subfamilies, different specificities of
one enzyme family — therefore clears all of their thresholds at once, and
KofamScan emits every one of them.

KEGG's own reference does not: across the seven benchmark genomes spanning all
three domains, **≥99.88% of genes carry exactly one KO** (max observed 2). The
extra calls are over-calls almost by construction, and they are where the
errors concentrate — on *E. coli* they are 13.5% of calls but **69% of all
false positives**, with precision 0.39 against 0.96 for single-KO genes.

`--multi-ko-policy best` keeps, per gene, the KO furthest above *its own*
threshold. Ranking is by that margin rather than raw score because thresholds
span roughly 30 to over 2000 bits, so raw scores are not comparable between
KOs — and empirically, ranking by *relative* margin (score/threshold) is
consistently worse than absolute margin. Measured against KEGG's assignments,
this raises mean precision from 0.874 to 0.903 and mean F1 from 0.885 to
0.893, improving **both on all seven genomes**; on *E. coli* it puts `mpph
annotate` ahead of the tool it reimplements (F1 0.941 vs KofamScan's 0.935).

`--min-ko-gap BITS` additionally requires the winner to beat the runner-up by
a margin, dropping the gene when the evidence does not separate them. It buys
precision (mean 0.909 at 40 bits) at a recall cost and leaves F1 flat, so it
is a preference about which error you would rather make, not a better setting.

**The default is `all`**, which reproduces KofamScan exactly — changing it
would silently alter existing users' results and invalidate that equivalence.
Full comparison of five candidate rules: `benchmarks/threshold_audit/`.

**Pushing precision further (`--min-margin BITS`).** Dropping calls close to
their threshold trades recall for precision on a measured frontier. With
arbitration on, mean precision runs 0.903 (0 bits) → 0.915 (10) → 0.924 (20) →
0.943 (50) → 0.956 (80), while recall falls 0.886 → 0.630 and **F1 declines
monotonically**. So this is worth using only when a false positive costs more
than a false negative — a property of your analysis, not of the tool. Default
is 0 (off).

It is not a substitute for arbitration and does not overlap with it:
arbitration alone (precision 0.903, recall 0.886) beats a raised margin alone
at 10 bits (0.894, 0.876) on *both* axes, because a wrong winner among
competing paralogous KOs and a weak hit are different errors.

**Memory is set by `--cpus`, not by your input.** Measured on one bacterial
genome, varying only the thread count: 1 thread 1.0 GB / 21 min, 4 threads
2.9 GB / 5 min, 28 threads 11.8 GB / 2 min — roughly 0.4 GB per thread, with
byte-identical output throughout. So a
48,265-protein plant proteome is not a memory problem and a large-memory node
is not required; if RAM is tight, lower `--cpus`. That costs little, because
parallel scaling is already well past linear by 28 threads (4 threads gives a
perfect 4.0x speed-up, 28 gives only 10x).

**Sequence loading.** `--sequence-loading {auto,prefetch,stream}` controls
whether the target proteome is held in memory or streamed. Both give
identical results; a prefetched block costs a measured ~1 kB per protein
(45 MB for all of *Arabidopsis*), so `auto` prefetches up to 1,000,000
proteins. This matters only for metagenome-scale protein catalogues, never
for a single organism.

**Significance rule**, ported from KofamScan's actual source (not its
README, which disagrees with its own code on one point):

- Every KOfam profile is searched against every input protein. A hit is
  *assigned* only if its bit score is **≥** the KO's own threshold from
  `ko_list` — KofamScan's README describes this as "higher than", but its
  real implementation (`result/hit.rb`) uses `>=`; this reimplements the
  actual behaviour, not the prose.
- Most KOs compare the **full-sequence** score against their threshold; a
  KO marked `score_type="domain"` in `ko_list` compares its **best single
  domain's** score instead (`result/parser.rb`) — a real branch, not every
  KO is scored the same way.
- **Some KOs have no threshold in `ko_list` at all** (too few reference
  sequences in KEGG GENES) — these can never be assigned, by KOfam's own
  design. This is expected, not a gap in this tool.

**Output** is `gene<TAB>K#####`, one row per significant assignment, no
header — the same shape as KofamScan `--format mapper`, so it's already
readable by `--user` with zero adapter code (see
[Input formats](input-formats.md)).

**The KOfam database is real KEGG data, not bundled with mpph.**
`mpph annotate --setup-db DIR` downloads it once (~1.5 GB compressed); this
is why `mpph`'s base install stays small — `pip install mpph[annotate]`
opts in to the extra `pyhmmer` dependency, and the database download is a
separate, explicit step. See `DATA_SOURCES.md` for KOfam's own terms.

## Phylogenetic non-independence (`compare --tree`)

`mpph compare`'s Fisher and Mann-Whitney tests assume every genome is an
independent observation. Genomes are not: close relatives share features by
common descent, so the same feature is counted many times over and p-values
are anti-conservative. Measured on a 32-tip balanced tree with the two groups
being the two clades descending from the root, and traits simulated with **no
group effect at all**, Fisher's exact test rejects at **27.7%** for a nominal
5% test. Across a wider grid of tree shapes, prevalences and rates it reaches
above 70%.

`--tree reference.nwk` replaces the exchangeable-samples null with one
simulated on the tree:

- **presence mode** — the feature is simulated under a symmetric 2-state Mk
  model whose rate is fitted by maximum likelihood (Felsenstein's pruning
  algorithm) on the observed data. Fitting a single rate with no group term
  *is* fitting the null, which makes this a parametric bootstrap.
- **completeness mode** — the feature is simulated under Brownian motion and
  scored with Cliff's delta. Because a rank statistic is invariant to
  positive scaling and to translation, the null distribution depends on
  neither the Brownian rate nor the ancestral state: **nothing is fitted**.
  The same null therefore serves every feature, which is why a large
  replicate count is cheap here.

**The presence null is conditioned on the observed number of present
genomes.** Without that condition roughly a third of simulated replicates
come out invariant, contribute a zero statistic, dilute the tail, and the
test then calls a purely clade-confounded feature significant (measured
p = 0.013 on a 16-tip tree, getting *worse* as the tree grows). Conditioning
gives p ≈ 0.8 on the same input. Fisher's exact test conditions on its
margins; this null has to as well. The accepted-replicate count and the
prevalence band actually used are reported per feature.

**The tree must be independent of the features being tested.** `mpph`'s own
dendrogram is built *from* those features, so using it here is circular; it
is rejected by comparing topologies (Robinson-Foulds distance 0), which
catches a renamed copy too. Supply a GTDB-Tk, 16S or concatenated
marker-gene tree. Branch lengths are required — a cladogram is rejected
rather than defaulted to 1, which would fabricate the distances the whole
correction is computed from.

**`n_state_changes` is the number that explains the result.** It is the
minimum number of state changes the tree implies (Fitch parsimony) — exact,
deterministic, no simulation. When it is 1 the feature arose once, and if
that single origin sits on the branch separating the groups then no method
can distinguish association from coincidence: the effective sample size is 1
(Maddison & FitzJohn 2015, *Syst. Biol.*). A large `p_value_phylo` there is
the correct answer, not a broken test, and `phylo_confounded` flags it.

**How much the correction moves a p-value depends on tree shape** — on how
much of the tree's total path length separates the two groups rather than
varying within them. Two deeply divergent clades are corrected strongly; two
interleaved sets of tips barely at all. This is a property of the question,
not a tuning knob.

**Resolution limit.** A permutation p-value cannot fall below
`1/(replicates+1)`, so after BH correction over `m` features the smallest
attainable q is `m/(replicates+1)`. With 999 replicates and 400 features that
floor is 0.40 — every q would look non-significant for purely arithmetic
reasons. `--phylo-permutations` therefore defaults to 9999, and `compare`
warns when the floor still exceeds 0.05.

Calibration and the uncorrected error rates above are reproducible with
`benchmarks/phylo/calibrate_type1.py`.
