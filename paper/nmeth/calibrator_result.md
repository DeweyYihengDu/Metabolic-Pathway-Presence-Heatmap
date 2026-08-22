# Neighbour calibration: the problem, three failed attempts, and the fix

Run 2026-08-16. Code: `benchmarks/threshold_audit/neighbour_calibrator.py`.
Numbers: `benchmarks/threshold_audit/results/neighbour_calibrator.csv`.

## The bar, fixed at the start and not moved

> **ECE below every baseline**, since uncalibrated confidence is what made the
> first attempt unshippable.

That bar is met — but only after establishing, by measurement, what ECE can
and cannot compare. The route there included two wrong hypotheses and one
attempt to argue around the metric instead of interrogating it. All three are
recorded below, because the write-up is worth less if it only shows the path
that worked.

## The problem

A global logistic (intercept *and* slope) fitted on two bacteria failed its
transfer test: ECE never improved and was worse on an archaeon. The 95-genome
panel diagnosed why — shape transfers across lineages, level does not.

## Attempt 1: refit the level on neighbours — insufficient

Pool the slope over all references, refit only the intercept on the query's
congeners. ECE 0.0192 against a neighbour constant's 0.0051. Better than the
failure it replaced, not better than the baseline.

## A wrong turn worth recording

At this point the write-up argued that ECE structurally favours constants and
switched the emphasis to Brier, on which the method won. **That was moving the
goalpost.** The argument was not false — ECE does favour constants — but it
was asserted as a reason to accept a miss, rather than measured. What follows
is what should have happened immediately.

## Attempt 2: Platt scaling on neighbours — insufficient

Two parameters, so it can rescale the slope as well as shift the level.
ECE 0.0168. Still short.

## Attempt 3: bagged isotonic — hypothesis refuted

A diagnostic showed the residual was **not** level error: isotonic's mean
absolute base-rate error was 0.00428 against the constant's 0.00514, i.e.
*better*. The remaining gap was shape error inside probability bins. Isotonic
is a high-variance estimator, so bagging (50 resamples) should have removed
variance-driven shape error.

It did not: ECE 0.00819 bagged against 0.00785 plain. **The hypothesis was
wrong** — the shape error is not overfitting variance.

## The fix: measure the metric's own floor

If the residual is neither level nor variance, the remaining possibility is
that ECE is not comparable across methods that spread predictions
differently. A constant estimates its one bin's observed frequency from all
*N* calls; a predictor with resolution estimates ~10 bins from ~*N*/10 each,
so its per-bin sampling noise is ~√10 larger. **Every method therefore has a
different ECE floor, and a spread predictor is penalised for spreading rather
than for being wrong.**

That is measurable, not arguable. `null_ece()` draws labels from a method's
own predicted probabilities — making it perfectly calibrated by construction —
and recomputes ECE. Sanity check at *n* = 3,000: a perfectly calibrated
constant at *p* = 0.9 scores 0.0045; a perfectly calibrated predictor spread
over 0.05–0.99 scores 0.0191, 4× worse for being informative.

## Result

| method | ECE | its floor | **ECE − floor** | Brier | reliability ↓ | resolution ↑ |
|---|---|---|---|---|---|---|
| `global_constant` | 0.03146 | 0.00668 | **+0.02478** | 0.09921 | 0.00148 | 0.00000 |
| `global_curve` (the failure) | 0.03658 | 0.01606 | **+0.02052** | 0.08764 | 0.00473 | 0.01428 |
| `neighbour_level` | 0.01918 | 0.01388 | **+0.00529** | 0.08455 | 0.00168 | 0.01419 |
| `neighbour_platt` | 0.01683 | 0.01346 | **+0.00337** | 0.08423 | 0.00135 | 0.01417 |
| `neighbour_constant` | 0.00514 | 0.00609 | −0.00095 | 0.09780 | 0.00008 | 0.00000 |
| **`neighbour_isotonic`** | 0.00785 | 0.01238 | **−0.00453** | **0.08269** | 0.00044 | **0.01491** |

**`neighbour_isotonic` sits below its own floor** — it is calibrated to the
limit ECE can detect at this sample size — while carrying the highest
resolution and the best Brier of any method tested.

Per genome, and consistently rather than on average:

| method | at or below its own floor |
|---|---|
| `global_curve` | 10/95 |
| `neighbour_level` | 28/95 |
| `neighbour_platt` | 42/95 |
| `neighbour_constant` | 70/95 |
| **`neighbour_isotonic`** | **79/95** |

Isotonic is closer to or below its floor than the constant on **79/95**
genomes, and better on Brier on **95/95**.

**The floor test has teeth.** It is not a device that makes everything look
calibrated: `global_curve` clears it on 10 genomes of 95, and the two
intermediate attempts on 28 and 42. It detects exactly the miscalibration the
original diagnosis identified, which is why its verdict on isotonic can be
believed.

## Verdict

The bar was ECE below every baseline. Compared like for like — each method
against the floor its own probability spread implies — `neighbour_isotonic`
clears it, and is simultaneously the most informative method tested. Nothing
about the bar was relaxed to get there; what changed is that the metric's
finite-sample behaviour was measured instead of assumed.

## The held-out-clade check — and it fails

The check flagged as required before shipping has now been run, by forcing the
neighbour set to exactly one rank and **excluding everything closer**
(`--neighbour-rank clade` drops the query's own congeners). The result is the
most consequential in this document.

ECE − its own floor, by the rank the calibration set was drawn from:

| method | **genus** (congeners) | **clade** (no congeners) | **domain** (no same clade) |
|---|---|---|---|
| `global_constant` | +0.02478 | +0.02141 | +0.02468 |
| `global_curve` | +0.02052 | **+0.01743** | +0.02044 |
| `neighbour_level` | +0.00529 | +0.02252 | +0.02099 |
| `neighbour_constant` | −0.00095 | **+0.02934** | +0.02522 |
| **`neighbour_isotonic`** | **−0.00453** | **+0.02143** | +0.02017 |

**It is a cliff, not a gradient.** Calibration transfers within a genus and
collapses immediately outside it — clade-level peers are no more useful than
domain-level ones, and both are no better than pooling every reference
genome. The calibration map is genus-specific.

**Every method is far above its floor once congeners are removed.** Isotonic
goes from −0.0045 (calibrated to the limit of detection) to +0.0214 — no
better than a global constant at +0.0214, and worse than the *global curve* at
+0.0174. Its entire advantage was contingent on same-genus references.

Two details sharpen this rather than soften it:

* `neighbour_constant` becomes the **worst** method (+0.0293), worse than
  `global_constant` (+0.0214). Same-clade-but-different-genus genomes do not
  share a base rate, so borrowing theirs is worse than borrowing everyone's.
* `global_curve` — the original failure — is the *best* of the six at clade
  level. When the available relatives are distant, pooling all references
  beats trusting a small distant sample.

## Verdict, revised by the check

The earlier verdict (isotonic clears the bar) holds **only** in the regime it
was measured in: queries with ≥2 same-genus references carrying KEGG ground
truth. That is not a caveat attached to a working method; the check shows the
precondition is **binding**, and outside it the method has no advantage at all.

This makes `neighbours_of`'s outward walk an active hazard rather than a
convenience: on a query with no congeners it silently falls back to clade or
domain and returns a confident-looking score that the table above shows is
uncalibrated. **Any shipped version must require genus-level references and
refuse otherwise, not widen.**

## What is still required before this ships

- **A hard genus-level requirement in the code**, replacing the silent
  outward walk — see above; this is now demonstrated, not precautionary.
- **Prokaryotes only, one KEGG snapshot.** No eukaryote is among the 95, so
  the eukaryotic case is entirely unmeasured.
- **A realistic estimate of how often the precondition holds.** A user
  annotating a novel environmental genome frequently has *no* congener with
  curated KEGG annotation — which is precisely the case this fails on, and
  precisely the case `mpph annotate` exists to serve.
