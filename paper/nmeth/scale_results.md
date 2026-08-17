# Two experiments at scale: neighbour calibration, and phylogenetic leakage

Run 2026-08-16. Both were previously reported as premise-checks that the
seven-genome benchmark could not resolve. Both are now properly powered.

---

## 1. Calibration is lineage-dependent, and closeness is what determines it

**Panel:** 95 prokaryotic genomes, 24 genera, all with KEGG-curated ground
truth; 322,637 proteins searched against KOfam in one pass. Every ordered
train→test pair evaluated: **8,930 pairs**, of which **282 are same-genus** —
against 2 in the seven-genome version, which is why that one could only report
a premise rather than an answer.

| training genome's relatedness | pairs | **ECE of the margin model** | ECE of a constant | Brier, margin model | Brier, constant |
|---|---|---|---|---|---|
| same genus | 282 | **0.0173** | 0.0062 | **0.0843** | 0.0979 |
| same clade | 560 | **0.0380** | 0.0373 | **0.0896** | 0.1030 |
| same domain | 8,088 | **0.0464** | 0.0427 | **0.0898** | 0.1003 |

**The premise holds and is now measured, not inferred.** A confidence curve
fitted on a congeneric genome is calibrated **2.7× better** than one fitted on
a distant relative (ECE 0.0173 vs 0.0464), and the gradient is monotone across
all three tiers. This is the quantitative basis the neighbour-recalibration
idea needed.

**One result cuts against the simple story and is reported because of that.**
At the genus tier the *constant* baseline is better calibrated than the margin
model (ECE 0.0062 vs 0.0173) while being clearly worse at discrimination
(Brier 0.0979 vs 0.0843). Congeneric genomes have such similar overall
precision that simply borrowing a close relative's average precision is nearly
perfectly calibrated already — the margin buys sharpness on top of that, not
calibration. The margin model wins on Brier at **every** tier; it wins on ECE
at none.

The design consequence is specific rather than general: a shippable confidence
score should take its *level* from close relatives and its *shape* from the
margin, rather than fitting one global curve and hoping. That is a concrete,
testable next method — and it is not what the earlier "fit a global logistic"
attempt did, which is why that attempt failed its transfer test.

**Limits.** All 95 genomes are prokaryotes from one KEGG snapshot, so
"same domain" here means Bacteria-vs-Bacteria and the cross-domain tier from
the seven-genome run is not reproduced at this scale. Genus is the closest
tier available; strain-level pairs, where the effect should be strongest, are
untested.

---

## 2. Random splits inflate module-prediction accuracy, modestly and consistently

**Design:** 3,000 KEGG prokaryotes across 535 genera. Full KO set gives the
label (is this module complete), a 30–90% down-sample gives the features,
following MetaPathPredict's own protocol. One fixed gradient-boosted model,
evaluated under a random split and a genus-blocked split **at matched
training-set size** — the control without which leakage cannot be separated
from reduced training diversity.

| | mean F1 |
|---|---|
| random split | **0.9262** |
| genus-blocked split | **0.8926** |
| **inflation** | **+0.0337** |

Inflated on **22 of 25 modules** (median +0.033, max +0.075).

**Read this as real, consistent, and smaller than the rhetoric around it might
suggest.** A 3.4-point F1 inflation is not the collapse that Yu et al. (2025)
documented for AMR prediction, and it does not by itself invalidate anything.
What it does establish is that the effect exists for metabolic-module
prediction, in the expected direction, on nearly every module tested.

**Two reasons this is probably a lower bound on the published setting.**
First, redundancy: 3,000 genomes across 535 genera is far less clade-redundant
than the ~30,596 RefSeq+GTDB genomes MetaPathPredict trained on, and leakage
scales with redundancy. Second, blocking rank: genus-level blocking still
allows same-family relatives across the split, so this measures only the
leakage that genus blocking removes.

**What this is not.** It is not a measurement of MetaPathPredict. No attempt
was made to reproduce that model, and reporting a lower number from a different
architecture would confound the evaluation-protocol question with
implementation differences. The claim is about the protocol: random splits over
a clade-redundant database report an optimistic number, and the size-matched
control shows that is not merely an artefact of training-set size.

---

## Where this leaves the programme

The margin–reliability relationship is established (it survives an independent
reference), its lineage dependence is now quantified at proper scale, and the
evaluation-protocol critique is confirmed in the expected direction. The
open methods question is unchanged in shape but sharper in detail: build a
confidence score that takes its level from phylogenetic neighbours and its
shape from the margin, then test whether *that* is calibrated on a held-out
clade. The failed global-curve attempt and the genus-tier result above
together say exactly why that split of responsibilities is the one to try.
