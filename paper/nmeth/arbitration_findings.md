# What the arbitration work actually established

Run 2026-08-15. Code and tables: `benchmarks/threshold_audit/`. Written after
the fact rather than before, so it records what the experiments returned
including where that contradicts how v3.20.0 was announced.

## The headline, and its boundary

`--multi-ko-policy best` raises per-gene KO precision from 0.874 to 0.903 and
F1 from 0.885 to 0.893, improving both on all seven benchmark genomes
(one-sided sign test p = 0.0078), with every per-genome 95% CI excluding zero
for both metrics. That is real and it is not a tuning
artefact — F1 is flat (0.8927–0.8934) across every rule parameter tried, so
the gain comes from arbitrating at all.

**It does not reach the module level.** Measured against completeness computed
from KEGG's own KO sets over seven genomes and 522 Pathway modules:

| level | precision | recall | F1 |
|---|---|---|---|
| per-gene (gene, KO) pairs | +2.9 pts | −1.3 | **+0.8** |
| KO set (deduplicated) | +1.6 pts | −1.4 | **+0.0006** |

| policy | module MAE | bias | modules broken | modules invented |
|---|---|---|---|---|
| `all` | 0.0268 | −0.0122 | **13.7** | 1.43 |
| `best` | 0.0266 | **−0.0176** | **15.7** | 1.00 |

MAE unchanged, bias 44% more negative, and **two more modules per genome**
that are complete in truth get broken.

**Mechanism, and it is not subtle once seen:** a genome has many genes but one
KO set. A false-positive KO call on gene X usually names a KO genuinely present
on another gene, so dropping it does not shrink the set. Dropping a true
positive that was a KO's only representative does. The precision gain is
roughly halved on the way up; the recall loss carries over intact.

## Where the lost recall goes

Not uniformly. Of 333 true calls dropped across seven genomes, 17 pathways
survive BH correction against a 2.78% baseline — but **those are not 17
independent findings.** Eight of them (Parkinson, Alzheimer, Huntington, prion
disease, neurodegeneration, diabetic cardiomyopathy, thermogenesis, chemical
carcinogenesis) each lost exactly 11 KOs, and those 11 are **100% the same
oxidative-phosphorylation subunits in every case** — KEGG's disease maps embed
the respiratory chain, so one signal is counted eight times. Checked directly
rather than inferred from the matching counts.

Collapsed to distinct signals:

| pathway | lost | enrichment | q |
|---|---|---|---|
| Glucosinolate biosynthesis | 5/16 = 31.3% | 11.2× | 0.0038 |
| Bacterial chemotaxis | 5/24 = 20.8% | 7.5× | 0.015 |
| β-Lactam resistance | 4/21 = 19.0% | 6.9× | 0.050 |
| ABC transporters | 38/249 = 15.3% | 5.5× | 1×10⁻¹⁵ |
| 2-Oxocarboxylic acid metabolism | 9/87 = 10.3% | 3.7× | 0.018 |
| Oxidative phosphorylation | 13/160 = 8.1% | 2.9× | 0.015 |
| Two-component system | 20/261 = 7.7% | 2.8× | 0.0033 |
| Biosynthesis of secondary metabolites | 41/859 = 4.8% | 1.7× | 0.015 |

Every one is a large paralogous family — transporters, two-component
kinases/regulators, chemotaxis proteins, respiratory-complex subunits (the
lost KOs are succinate dehydrogenase K00234/K00235, F-type ATPase
K02132/K02133/K02137, NADH dehydrogenase K03882/K03883), and in *Arabidopsis*
the P450/glucosyltransferase families behind glucosinolates. Arbitration
exists to fix over-calling caused by paralogy and over-corrects hardest
exactly where paralogy is highest — **it trades one paralogy artefact for
another.** No single KO dominates (most-displaced is lost 4 times), so this is
a property of the families, not of a few bad profiles, and cannot be fixed by
blacklisting.

## The comparison was unfair in this tool's favour

eggNOG-mapper also emits multi-KO genes (15.4% of *E. coli*, 9.9% of
*Arabidopsis*). Comparing arbitrated mpph against unarbitrated eggNOG measures
the arbitration step. eggNOG cannot be arbitrated the same way — no per-KO
score to rank by — which is a genuine advantage of the profile-HMM route and
is now claimed as one instead of quietly enjoyed. Under a one-KO-per-gene
bound eggNOG gains 5.8 points of precision (0.749 → 0.807), and the ranking
is unchanged: unarbitrated mpph is 0.874, arbitrated 0.903.

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

## On MAGs, the target use case, this is formally untested

MAGs carry no KEGG ground truth, so precision there cannot be measured at all
and none of the above is verified on the genomes the tool is most used for.
The measurable part is how much arbitration changes, and it is markedly less:
**2.07% of calls resolved on the three MAGs (1.75–2.26%) against 4.70% on the
reference genomes (1.91–8.29%)** — under half.

The likely cause is mundane: an incomplete assembly recovers fewer members of
each paralogous family, so fewer KO profiles compete on one protein. The
reading is bounded rather than reassuring — arbitration's effect on MAGs is
about half as large, so whatever its untested benefit or harm there, the
exposure is proportionally smaller. It is not evidence that it helps.

## What this changes about the project's direction

The KO-precision work is a solid contribution to `mpph annotate` as a
gene-level annotator, and it is now correctly bounded in the documentation.
It is **not** the Nature Methods contribution — it makes the tool better at
reproducing KEGG's per-gene assignments, which is engineering, and its benefit
vanishes at the level this toolkit actually reports.

The live methods question remains the one from the calibration transfer test:
lineage-aware confidence for metabolic capability. This work sharpened it
rather than answered it, and supplied the premise check it needed.
