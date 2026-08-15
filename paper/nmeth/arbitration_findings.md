# What the arbitration work actually established

Run 2026-08-15. Code and tables: `benchmarks/threshold_audit/`. Written after
the fact rather than before, so it records what the experiments returned
including where that contradicts how v3.20.0 was announced.

## The headline, and its boundary

`--multi-ko-policy best` raises per-gene KO precision from 0.874 to 0.903 and
F1 from 0.885 to 0.893, improving both on all seven benchmark genomes, with
every per-genome 95% CI excluding zero. That is real and it is not a tuning
artefact — F1 is flat (0.8927–0.8934) across every rule parameter tried, so
the gain comes from arbitrating at all.

**It does not reach the module level.** Measured against completeness computed
from KEGG's own KO sets over six genomes and 522 Pathway modules:

| level | precision | recall | F1 |
|---|---|---|---|
| per-gene (gene, KO) pairs | +2.9 pts | −1.3 | **+0.8** |
| KO set (deduplicated) | +1.6 pts | −1.4 | **+0.0004** |

| policy | module MAE | bias | modules broken | modules invented |
|---|---|---|---|---|
| `all` | 0.0272 | −0.0152 | **14.7** | 1.50 |
| `best` | 0.0267 | **−0.0184** | **16.0** | 1.17 |

MAE unchanged, bias more negative, and *more* modules that are complete in
truth get broken — worse on four of six genomes, equal on two, better on none.

**Mechanism, and it is not subtle once seen:** a genome has many genes but one
KO set. A false-positive KO call on gene X usually names a KO genuinely present
on another gene, so dropping it does not shrink the set. Dropping a true
positive that was a KO's only representative does. The precision gain is
roughly halved on the way up; the recall loss carries over intact.

## Where the lost recall goes

Not uniformly. Of 215 true calls dropped across six genomes, three pathways
survive BH correction against a 2.6% baseline loss rate:

| pathway | lost | enrichment | q |
|---|---|---|---|
| ABC transporters | 38/246 = 15.4% | 5.9× | 6×10⁻¹⁷ |
| Two-component system | 20/260 = 7.7% | 2.9× | 0.0024 |
| Bacterial chemotaxis | 5/24 = 20.8% | 7.9× | 0.035 |

The three most heavily paralogous families in bacterial genomes. Arbitration
exists to fix over-calling caused by paralogy and over-corrects hardest exactly
where paralogy is highest — **it trades one paralogy artefact for another.**
No single KO dominates (most-displaced is lost 4 times), so this is a property
of the families, not of a few bad profiles, and cannot be fixed by blacklisting.

## The comparison was unfair in this tool's favour

eggNOG-mapper also emits multi-KO genes (15.4% of *E. coli*, 9.9% of
*Arabidopsis*). Comparing arbitrated mpph against unarbitrated eggNOG measures
the arbitration step. eggNOG cannot be arbitrated the same way — no per-KO
score to rank by — which is a genuine advantage of the profile-HMM route and
is now claimed as one instead of quietly enjoyed. Under a one-KO-per-gene
bound eggNOG gains 5.7 points of precision (0.771 → 0.828), and the ranking
is unchanged: unarbitrated mpph is 0.883.

## Neighbour-based calibration: premise supported, test not possible here

The transfer result said a global confidence curve is not shippable and
calibration must be lineage-aware. Testing that properly is **not possible on
this benchmark** — seven genomes across three domains have no meaningful
"neighbour"; they are all maximally distant. Rather than fake it, only the
premise was tested: over all 42 ordered train→test pairs,

| training genome's relatedness | ECE of the margin model | Brier improvement | pairs |
|---|---|---|---|
| same phylum group | 0.046 | +0.0201 | 2 |
| same domain | 0.045 | +0.0163 | 12 |
| different domain | **0.073** | +0.0119 | 28 |

Calibrating across domains is ~60% worse, and Brier improvement declines
monotonically with distance. So relatedness does carry calibration
information — but only the domain-level contrast is resolvable here (tier 1
has n=2). The finer question the method depends on, whether genus-level
closeness beats family-level, **needs on the order of hundreds of genomes with
KEGG ground truth sampled across a range of pairwise distances**. That is the
next real experiment and it is a compute job, not an analysis.

## What this changes about the project's direction

The KO-precision work is a solid contribution to `mpph annotate` as a
gene-level annotator, and it is now correctly bounded in the documentation.
It is **not** the Nature Methods contribution — it makes the tool better at
reproducing KEGG's per-gene assignments, which is engineering, and its benefit
vanishes at the level this toolkit actually reports.

The live methods question remains the one from the calibration transfer test:
lineage-aware confidence for metabolic capability. This work sharpened it
rather than answered it, and supplied the premise check it needed.
