# Pilot result — what the binary KO call throws away

Run 2026-08-15 on Hermes. Code: `benchmarks/threshold_audit/`. Raw tables:
`threshold_audit/summary.csv`, `threshold_audit/curve.csv`.

The question idea_003 posed was whether KOfam's adaptive threshold discards
recoverable true positives. **It largely does not** — and the more interesting
result is on the other side of the threshold.

## What was measured

Every KOfam profile was searched against each benchmark genome retaining all
hits scoring above zero, not only significant ones, and each (gene, KO) pair
was scored against KEGG's own curated `link/ko/<org>` assignments as a
function of `delta = score − that KO's adaptive threshold`. `delta ≥ 0` is
what every tool in the field keeps; `delta < 0` is what every tool discards.

## Result 1 — the discarded band is mostly not recoverable

Precision just below the threshold is ~0.10 in the first 5 bits (*E. coli*)
and decays to the background rate by about −20 bits. A label-shuffled control
gives 0.000 in essentially every bin, so the residual signal is real but
small: 9–15 KEGG-confirmed assignments per genome within 5 bits of the cut.

**KOfam's threshold is well placed.** The "rescue the false absences" framing
in idea_002's first channel is not supported and should be dropped.

## Result 2 — the *accepted* calls are not homogeneous, and that is the finding

Precision of accepted calls rises monotonically with margin above threshold,
in every genome tested:

| margin above threshold (bits) | *E. coli* | *B. subtilis* | *M. jannaschii* | *Verruco.* S94 | *Lentisph.* WC36 | *S. cerevisiae* | *A. thaliana* |
|---|---|---|---|---|---|---|---|
| 0–5 | 0.397 | 0.281 | 0.513 | 0.411 | 0.375 | 0.588 | **0.269** |
| 5–10 | 0.622 | 0.440 | 0.545 | 0.453 | 0.667 | 0.760 | 0.311 |
| 10–20 | 0.653 | 0.560 | 0.794 | 0.430 | 0.679 | 0.823 | 0.356 |
| 20–40 | 0.787 | 0.674 | 0.843 | 0.715 | 0.704 | 0.909 | 0.558 |
| 40–60 | 0.843 | 0.810 | 0.896 | 0.840 | 0.795 | 0.921 | 0.653 |
| 60–100 | 0.882 | 0.840 | 0.940 | 0.889 | 0.853 | 0.946 | 0.736 |
| 100–200 | 0.939 | 0.938 | 0.977 | 0.947 | 0.931 | 0.957 | 0.876 |
| ≥200 | 0.966 | 0.978 | 0.989 | 0.984 | 0.954 | 0.975 | 0.943 |

**A KO call is anywhere from 27% to 99% likely to be right, and every tool in
the field writes all of them as the same `1`.** The quantity that predicts
which is computed by the HMM search and then thrown away. The relationship is
monotone in all seven genomes, across all three domains of life — the
*Arabidopsis* column is the steepest, spanning 0.269 to 0.943.

## Result 3 — the errors concentrate in a small, identifiable minority

| genome | domain | precision, all accepted calls | % of calls within 20 bits | precision of those | **% of all false positives they contribute** | enrichment |
|---|---|---|---|---|---|---|
| *E. coli* K-12 | Bacteria | 0.900 | 5.7% | 0.572 | **24.4%** | 4.26× |
| *B. subtilis* 168 | Bacteria | 0.858 | 7.9% | 0.449 | **30.8%** | 3.89× |
| *M. jannaschii* | Archaea | 0.903 | 12.9% | 0.652 | **46.5%** | 3.61× |
| *Verrucomicrobia* sp. S94 | Bacteria (PVC) | 0.845 | 12.7% | 0.431 | **46.9%** | 3.68× |
| *Lentisphaerae* sp. WC36 | Bacteria (PVC) | 0.848 | 11.3% | 0.597 | **29.9%** | 2.66× |
| *S. cerevisiae* | Eukaryota | 0.945 | 5.6% | 0.752 | **25.4%** | 4.51× |
| *A. thaliana* | Eukaryota | 0.815 | 7.6% | **0.319** | **28.0%** | 3.68× |

Between 6% and 13% of calls carry between a quarter and a half of every
error a genome makes — a **2.7–4.5× enrichment in every one of seven genomes
spanning all three domains** — and membership of that set is known at
annotation time, for free.

*Arabidopsis* is the sharpest case: a low-margin call in the plant proteome
is correct **less than a third of the time** (0.319), against 0.815 for its
calls overall. 51 KEGG-confirmed assignments also sit within 5 bits *below*
its threshold, the largest sub-threshold residue of any genome — as expected
for the genome furthest from the sequences KOfam's profiles were built on.

## Result 4 — it gets worse exactly where the tools are used

The low-confidence fraction roughly **doubles** away from the model bacteria:

- model bacteria (*E. coli*, *B. subtilis*): 5.7%, 7.9%
- archaeon (*M. jannaschii*): 12.9%
- environmental PVC bacteria (*Verrucomicrobia* S94, *Lentisphaerae* WC36): 12.7%, 11.3%

The sub-threshold band shows the same gradient: precision below the cut rises
from 0.007 (*E. coli*) to 0.035 (*Lentisphaerae*), a 5× increase. Both effects
point the same way — the further a genome is from the sequences KOfam's
profiles were built on, the more of its annotation sits in the uncertain zone.
That is the regime environmental genomics and MAG-based work operate in, and
it is where every one of these tools is actually deployed.

## Why this is self-consistent with an earlier, independent result

v3.19.0 established that all three disagreements between `mpph annotate` and
KofamScan across *Arabidopsis*' 12,223 calls lie **within 0.016 bits of the
threshold**. That was reported as a rounding artefact, which it is. It is also
the same phenomenon from another angle: tool disagreement concentrates exactly
where precision is ~0.4. Two independent observations, one explanation.

## What this changes

The idea_002 framing survives but its emphasis inverts. The recoverable
information is not in rescuing absences — it is in **calibrating the
presences**, which are currently reported with a confidence they do not have.
Concretely:

1. **`mpph annotate` should emit the margin**, and a per-KO calibrated
   P(correct) derived from a curve like the one above. Cheap, immediate,
   and no other tool does it.
2. **Module completeness should propagate it.** A module scored 0.75 from four
   KOs each at P≈0.45 is a different claim from one scored 0.75 from four KOs
   at P≈0.97. These are currently indistinguishable in every tool, including
   this one.
3. **The differential test should weight by it.** `compare --tree` already
   corrects for phylogenetic non-independence; it still treats every feature
   value as equally certain.

## Result 5 — the decisive test: the signal transfers, the calibration does not

A curve measured on the same genomes it describes proves nothing. The gate is
whether a curve fitted on one lineage is *calibrated* on another. Fitted on
the two model bacteria (`P = sigmoid(−2.074 + 0.938·log1p(margin))`, n=6,209,
training precision 0.882) and applied to three held-out genomes, against the
only honest baseline — a constant confidence equal to the training precision,
which is what a tool implicitly asserts today when it emits a bare `1`.
95% CIs are paired bootstrap over calls, n=2,000.

| held-out genome | ECE improvement [95% CI] | Brier improvement [95% CI] |
|---|---|---|
| *M. jannaschii* (Archaea) | **−0.045 [−0.059, −0.037]** — significantly *worse* | +0.008 [−0.001, +0.018] — n.s. |
| *Verrucomicrobia* S94 (PVC) | +0.001 [−0.024, +0.017] — n.s. | **+0.030 [+0.022, +0.038]** |
| *Lentisphaerae* WC36 (PVC) | +0.002 [−0.023, +0.014] — n.s. | **+0.015 [+0.007, +0.025]** |

Read plainly, this is one positive and one negative result, and both matter:

- **The margin carries genuinely transferable information.** Brier improves
  significantly on both environmental PVC genomes and trends positive on the
  archaeon. A curve fitted on *E. coli* and *B. subtilis* makes materially
  better probability estimates on a Verrucomicrobium than the field's implicit
  constant does.
- **The absolute probabilities do not transfer.** ECE is never significantly
  improved and is significantly *worse* on the archaeon, where the model is
  systematically under-confident (stated 0.65 → observed 0.82; stated 0.86 →
  0.92). The constant baseline wins there only because *M. jannaschii*'s true
  precision (0.903) happens to sit near the training rate (0.882) — luck about
  that genome, not a property of the baseline.

**So a single global confidence curve is not shippable.** That is the finding,
and it is more useful than a clean win would have been, because it says
precisely what the method has to do: calibration must be **lineage-aware**.

The natural mechanism is already in this repo. `mpph.phylo` reads a reference
phylogeny; the genomes that *do* have KEGG-curated ground truth can be located
on that tree; and a query genome can be recalibrated against its phylogenetic
neighbours rather than against a global average. That is a concrete,
falsifiable next experiment — does neighbour-based recalibration close the ECE
gap on a held-out clade? — and it is not something any tool in the field does.

## Honest limits of this pilot

- **Precision measured here is a lower bound.** KEGG's reference is itself
  similarity-derived and incomplete, so a hit counted as a false positive may
  be a genuine function KEGG has not recorded. This biases the numbers down,
  which is the conservative direction, but they are bounds and must not be
  quoted later as point estimates.
- **Seven genomes.** Enough to show the effect holds across all three domains,
  not enough to characterise how the curve varies *between* lineages — which
  is exactly what Result 5 says the method needs. The three MAGs in the
  annotation benchmark have no ground truth and so cannot be audited this way
  at all; whether the curve behaves the same on genuinely fragmentary genomes
  is untested and is a real gap, since MAGs are the target use case.
- **Transfer was tested on three held-out genomes from a two-genome training
  set.** That is a small basis for a claim about calibration in general. The
  direction of the result (ranking transfers, level does not) is consistent
  across all three, but the CIs on the ECE differences are wide enough that
  "no improvement" and "small improvement" are not distinguishable for the two
  PVC genomes.
- **No part of this is shippable yet.** Result 5 is a negative result on the
  obvious implementation. Nothing here licenses adding a confidence score to
  `mpph annotate` until neighbour-based recalibration is tested.
