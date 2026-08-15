# MPPH: an end-to-end, phylogeny-aware toolkit for comparative metabolic-pathway analysis across any genome

**Yiheng Du**

*Draft — Applications Note format (Bioinformatics). Every number below is
reproducible from `benchmarks/` in the repository; none is estimated.*

---

## Abstract

**Summary.** Comparative metabolic-profile analysis is normally assembled by
hand from three or four separate tools: one to assign KEGG Orthology (KO)
terms to genes, another to build a presence or completeness matrix, a third
for enrichment testing, a fourth for figures. Each hand-off invites the
silent errors that come from mismatched identifier namespaces and
incompatible defaults, and the group comparisons at the end of the chain
routinely treat genomes as independent samples when they are not. MPPH is a
single command-line toolkit covering that whole path — local KO annotation
from a protein FASTA, presence/completeness matrix construction, pan-genome
classification, ordination, over-representation and rank-based enrichment,
KEGG pathway-map rendering, and differential testing — with two properties
that are unusual in this space: its KO annotation reproduces KofamScan
essentially call for call while requiring no external HMMER installation,
and its differential test can be corrected for phylogenetic
non-independence against a user-supplied reference phylogeny. Uncorrected,
the same test's false-positive rate reaches 94.8% at a nominal 5% on
simulated data containing no group effect at all. The toolkit is not
restricted to any domain of life: annotation was validated against KEGG's own
assignments on bacterial, archaeal and eukaryotic genomes, including a
48,265-protein plant proteome, with F1 between 0.84 and 0.93 throughout.

**Availability and implementation.** MPPH is a pure-Python package
(Python ≥3.9, MIT licence) installable with `pip`. Source, documentation and
all validation code and results: https://github.com/DeweyYihengDu/Metabolic-Pathway-Presence-Heatmap

---

## 1 Introduction

Two gaps recur in comparative metabolic-profile work.

The first is upstream. Assigning KO terms to a newly sequenced genome,
proteome or metagenome-assembled genome (MAG) is a precondition for
everything downstream, but the accessible options are awkward: KEGG's own
KOALA servers are browser-and-email services with no programmatic interface;
KofamScan requires a separate HMMER installation and Ruby toolchain;
approaches that transfer annotation from a fixed panel of reference species
degrade with evolutionary distance, and tools of the KOBAS type ship
curated background sets for a short list of model organisms — human, mouse
and a handful of others — leaving everything else unsupported. Users
therefore either leave the annotation step outside their pipeline or accept a
narrowing of taxonomic scope.

The second is downstream and is a statistical error rather than an
inconvenience. Testing whether a metabolic feature differs between two
groups of genomes with Fisher's exact test or a Mann–Whitney test assumes
each genome is an independent observation. Genomes are not independent:
close relatives share features by common descent, so the same evolutionary
event is counted many times. The resulting p-values are anti-conservative,
and the effect is not marginal — see §3.3.

MPPH addresses both within one toolkit, and — because a tool's own claims
about itself are worth little — ships the validation that supports each
claim as runnable code with committed results.

---

## 2 Implementation

MPPH is organised as 13 subcommands over a shared KEGG REST layer with an
on-disk response cache; every run writes a manifest recording the command,
the KEGG release, dependency versions and the calling directory's git
commit, so a result can be traced to the exact software and reference data
that produced it.

**Local KO annotation (`mpph annotate`).** Protein sequences are scored
against KEGG's own KOfam profile HMM database using `pyhmmer` (Larralde and
Zeller, 2023), Cython bindings to HMMER3 that require no separate HMMER
binary. The assignment rule reimplements KofamScan's: a hit is accepted only
when its bit score meets that KO's own adaptive threshold from `ko_list`,
using the best single domain's score for KOs marked `score_type = domain`
and the full-sequence score otherwise. Because profile HMMs are built from
many species' sequences rather than transferred from a nearest reference,
the method's applicability does not depend on taxonomic proximity to a model
organism, nor on which domain of life the input belongs to — the same
profile set is searched for a plant proteome as for a bacterial genome
(§3.1). Output is `gene<TAB>KO`, identical in shape to KofamScan's mapper
format, so it feeds the rest of the toolkit with no conversion step.

**Matrix construction and NaN semantics.** Presence mode records whether
KEGG lists an organism-specific version of each pathway map; completeness
mode scores KEGG module DEFINITION expressions, evaluating them as a
three-valued logic so that a step depending on an unresolvable nested module
reference is reported as *unknown* rather than silently counted as absent.
Throughout the package NaN means "not assessed", never "confirmed absent" —
a distinction that matters for MAGs, where a missing feature is frequently
an assembly gap. `--unknown-policy` makes the handling of unknowns an
explicit choice (exclude, treat as absent, or refuse to run) rather than a
hidden default, and the reported evidence counts remain policy-independent.

**Phylogeny-aware differential testing (`mpph compare --tree`).** Given an
independent reference phylogeny with branch lengths, the exchangeable-samples
null is replaced by one simulated on the tree: a symmetric two-state Mk model
whose rate is fitted by maximum likelihood (Felsenstein's pruning algorithm)
for presence data, Brownian motion for completeness data. Two design points
are load-bearing:

*The presence null is conditioned on the observed number of genomes carrying
the feature.* Unconditioned, roughly a third of simulated replicates come out
invariant, contribute a zero statistic and dilute the tail; the test then
assigns p = 0.013 to a feature whose distribution is entirely explained by
clade membership, and the failure worsens as the tree grows. Conditioned, the
same input gives p ≈ 0.78. Fisher's exact test conditions on its margins;
this null must also.

*The completeness null estimates nothing.* Scored with Cliff's delta — a rank
statistic, invariant to positive scaling and to translation — the Brownian
null depends on neither the evolutionary rate nor the ancestral state. No
parameter is fitted, and one null distribution serves every feature.

Alongside the corrected p-value each feature reports `n_state_changes`, the
minimum number of state changes the tree implies (Fitch parsimony; exact and
simulation-free). When this equals 1 the feature arose once, and if that
single origin lies on the branch separating the two groups then no method can
distinguish association from coincidence — the effective sample size is one
(Maddison and FitzJohn, 2015). A large corrected p-value there is the correct
answer, and the reported integer is what makes it interpretable rather than
puzzling.

MPPH's own functional dendrogram is **rejected** as the reference tree. It is
built from the very features being tested, so correcting those tests with it
would be circular; the check compares topology rather than filenames, so a
renamed copy is caught as well.

---

## 3 Validation

All validation code, container digests, accessions and result tables are in
`benchmarks/`; figures regenerate from the committed tables alone.

### 3.1 KO annotation

`mpph annotate` was compared with KofamScan 1.3.0 (against the same KOfam
database) and eggNOG-mapper 2.1.15 on ten genomes spanning all three domains
of life, in tiers of increasing difficulty: five model organisms
(*Escherichia coli* K-12, *Bacillus subtilis* 168, *Methanocaldococcus
jannaschii*, *Saccharomyces cerevisiae* S288C, *Arabidopsis thaliana*
Col-0), two environmental genomes KEGG has never hand-curated
(*Verrucomicrobia* sp. S94, *Lentisphaerae* sp. WC36 — both PVC superphylum,
both carrying placeholder strain names), and three marine MAGs with no
reference annotation at all.

The reference is KEGG's own KO assignment for each genome. That reference is
methodologically independent of the methods being tested: KEGG's GENES
annotations are built on best-hit sequence-similarity tables, whereas both
MPPH and KofamScan score profile HMMs against per-KO thresholds.

Eukaryotic proteomes require one adjustment before they can be scored at all:
they list every splice isoform (*A. thaliana*: 48,265 proteins) while KEGG
names one representative protein per gene (27,562). Calls are therefore
restricted to the proteins the reference covers; without this, 10,025 calls
on isoforms the reference cannot represent would have been counted as false
positives, and a hypothetical perfectly-correct tool would score precision
0.5. The restriction is applied identically to all three tools and is
essentially a no-op for prokaryotes (*E. coli*: 4,288 of 4,300 covered).

**MPPH reproduces KofamScan almost exactly.** The Jaccard overlap of
(gene, KO) assignments is 1.000 on seven of the ten genomes, and 0.9997,
0.9999 and 0.9994 on the remaining three — one call in 3,517, one in 12,223
and one in 1,718. F1 against KEGG's assignments:

| Tier | Domain | Genome | MPPH | KofamScan | eggNOG-mapper |
|---|---|---|---|---|---|
| Model | Bacteria | *E. coli* K-12 | 0.934 | 0.934 | 0.851 |
| Model | Bacteria | *B. subtilis* 168 | 0.898 | 0.898 | 0.825 |
| Model | Archaea | *M. jannaschii* | 0.894 | 0.894 | 0.907 |
| Model | Eukaryota | *S. cerevisiae* | 0.911 | 0.911 | 0.862 |
| Model | Eukaryota | *A. thaliana* | 0.867 | 0.866 | 0.697 |
| Non-model | Bacteria | *Verrucomicrobia* sp. S94 | 0.848 | 0.848 | 0.755 |
| Non-model | Bacteria | *Lentisphaerae* sp. WC36 | 0.842 | 0.842 | 0.761 |

The objection this tier structure was built to test — that good accuracy on
long-curated model organisms says little about arbitrary genomes — is not
supported: the profile-HMM margin over eggNOG-mapper is *larger* on the two
non-model genomes (+0.09, +0.08) than on the model organisms (+0.08, +0.07).
On the three MAGs, where no reference exists to score against, MPPH and
KofamScan still agree at 1.000, 0.9994 and 1.000.

Accuracy does not degrade outside bacteria: yeast (0.911) scores above three
of the four prokaryotes, and *A. thaliana* (0.867) falls within the
prokaryotic range while carrying the benchmark's widest margin over
eggNOG-mapper (+0.17 F1) — the similarity-transfer approach loses most ground
precisely where the query is furthest from a close reference. What does fall
is coverage, for all three tools alike: the fraction of proteins receiving any
KO drops from 76% (*E. coli*) to 62% (yeast) to 24% (*A. thaliana*). That is
a property of KO's scope, which describes metabolism and core cellular
processes rather than a plant's full protein complement; recall against what
KEGG does assign remains highest of any genome tested (0.925).

The residual disagreements between MPPH and KofamScan have a single
explanation. All three differing calls across *A. thaliana* lie within 0.016
bits of the relevant KO threshold, and in both directions. HMMER prints bit
scores to one decimal place, and KofamScan parses that text output, so it
compares a rounded score against a threshold specified to two decimals
(`332.1 ≥ 332.13` is false; `164.1 ≥ 164.07` is true), whereas MPPH reads the
score from `pyhmmer`'s in-memory hit object at full precision. The two
implementations therefore agree on every call that is not a tie within the
reference database's own printing precision.

eggNOG-mapper is not thereby "worse": it annotates more genes at lower
precision, a different operating point, and produces far more than KO calls.

**Cost, including where MPPH is worse.** Wall-clock time is comparable —
MPPH is faster on seven of ten genomes. Peak memory is not: 2.3–11.8 GB against
KofamScan's 0.13–0.74 GB, 11–18× higher. It is, however, not a function of
input size, and the obvious explanation is the wrong one: a prefetched
sequence block costs a measured ~1 kB per protein (45 MB for the entire
*A. thaliana* proteome), and re-running *B. subtilis* with targets streamed
rather than prefetched changes peak RSS from 11.4 to 11.7 GB with
byte-identical output — less than the spread between two runs of the
identical prefetch command. The worst case in the benchmark is in fact the
4,237-protein *B. subtilis* (11.8 GB), not the 48,265-protein *A. thaliana*
(9.7 GB). What peak memory tracks is thread count: on one fixed input,
1/4/28 threads give 1.0/2.9/11.8 GB and 21:07/5:17/2:07, with identical
output. `--cpus` is
therefore the memory control, and lowering it is cheap — parallel scaling is
already well past linear by 28 threads (4 threads is a perfect 4.0× speed-up,
28 only 10×). The absolute level is reported in full rather than omitted, but
neither a large proteome nor a large-memory node is required.

### 3.2 Enrichment

`mpph gsea` and `mpph enrich` were compared with clusterProfiler 4.18.4 on
identical category definitions (MPPH's own KO-to-pathway snapshot exported as
a `TERM2GENE` table) and identical size filters, using a real ranked list
derived from *E. coli* RNA-seq (GEO GSE189154). Agreement on the 55 commonly
tested pathways: Spearman ρ = 1.000 for the enrichment score, 0.999 for NES,
0.994 for the p-value; for over-representation, ρ = 1.000 for the p-value.
The remaining difference in reported pathway counts is a reporting
convention — clusterProfiler's `enricher()` omits zero-hit categories from
its output while MPPH reports them as non-significant — not a difference in
any pathway's statistics.

### 3.3 Phylogenetic correction

Traits were simulated on trees under a symmetric Mk model with **no group
effect present**, with the two compared groups being the two clades
descending from the root; every rejection is a false positive by
construction. The grid covers balanced and unbalanced non-ultrametric trees,
16/32/64 tips, three prevalence strata and three evolutionary rates (54
cells, 50 evaluable).

| | uncorrected | `--tree` corrected |
|---|---|---|
| worst Type I error | **0.948** | **0.072** |
| median Type I error | 0.086 | **0.000** |
| cells above nominal 0.05 | 31 of 50 | **1 of 50** |

The single corrected cell above nominal rejected 18 of 250; under a true rate
of 0.05 that has P(≥18) = 0.079, and across 50 cells 3.9 such cells are
expected by chance, so one is fewer than chance expectation. The correction is
if anything conservative.

Stratifying by prevalence was essential rather than decorative: an earlier,
unconditioned formulation of the null averaged to an acceptable ~5% overall
while being roughly fourfold anti-conservative on exactly the
intermediate-prevalence features that produce reportable results.

---

## 4 Conclusion

MPPH covers the comparative metabolic-profiling path end to end, removing the
identifier hand-offs where errors enter silently, and closes a statistical gap
that is standard practice to document as a caveat rather than fix. Its KO
annotation is empirically interchangeable with KofamScan while dropping the
external toolchain requirement, and its enrichment statistics match the
reference R implementation to three decimal places. The correction for
phylogenetic non-independence is opt-in, requires a reference phylogeny
independent of the data, and reduces a 94.8% worst-case false-positive rate to
7.2%. None of this is restricted to a taxonomic scope: the same code and the
same reference database were validated on bacteria, archaea, a fungus, a
plant, and metagenome-assembled genomes with no reference annotation at all.

---

## Notes for revision (not for submission)

- **Author list and acknowledgements** are for the author to set. The 2023
  preprint (doi:10.1101/2023.06.27.546232) lists Y.-H. Du and J.-H. Mu; this
  manuscript describes a substantially different tool and its authorship
  should be decided deliberately, not inherited.
- **Framing decision outstanding**: this draft is written as a Bioinformatics
  Applications Note (tool-first, ~2 pages, no new biology). If a longer venue
  is preferred, the natural expansion is a real application section using the
  marine PVC MAGs already present in the benchmark, which would carry the
  biological result the Applications Note format has no room for.
- **Figures**: all three are submission-ready PDFs in
  `benchmarks/figures/` — `fig1_annotation_benchmark` (§3.1),
  `fig2_enrichment_benchmark` (§3.2), `fig3_phylo_calibration` (§3.3).
  An Applications Note normally allows one; Fig. 1 carries the most
  distinctive claim, Fig. 3 the most consequential one. Combining Fig. 1a
  with Fig. 3 into a two-panel figure is the likely compromise.
- **References** are cited in text but not yet formatted; the substantive ones
  are Kanehisa (KEGG), Aramaki et al. 2020 (KofamScan), Larralde and Zeller
  2023 (pyhmmer), Cantalapiedra et al. 2021 (eggNOG-mapper), Wu et al.
  (clusterProfiler), Felsenstein 1985/2004, Garland et al. 1993, and Maddison
  and FitzJohn 2015.
