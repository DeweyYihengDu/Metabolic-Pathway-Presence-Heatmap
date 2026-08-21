# The neighbour-level calibrator: fails its own bar, and says why

Run 2026-08-16. Code: `benchmarks/threshold_audit/neighbour_calibrator.py`.
Numbers: `benchmarks/threshold_audit/results/neighbour_calibrator.csv`.

## What was being fixed

The first calibration attempt fitted one global logistic — intercept *and*
slope — on two bacteria, and failed its transfer test: ECE never improved and
was significantly worse on the archaeon. The 95-genome panel then supplied a
diagnosis. Shape transfers (the margin model beats a constant on Brier at
every relatedness tier); level does not (ECE degrades 2.7× from same-genus to
same-domain). So the two halves of the logistic should be estimated from
different places: **slope pooled across all references, intercept refit on the
query's nearest relatives.**

That is implemented and tested here, leave-one-genome-out over 95 genomes.
Every query had ≥2 congeneric peers, so `rank_used` is `genus` for all 95 —
the panel was built for exactly this and it held.

## The bar, set before running

> Success is not "beats the others on Brier" — the margin model already did
> that and it was not enough. It is **lower ECE than every baseline**, since
> uncalibrated confidence is precisely what made the first attempt
> unshippable.

## Result

| method | ECE | Brier |
|---|---|---|
| `global_constant` | 0.0315 | 0.0992 |
| `global_curve` (the attempt that failed) | 0.0366 | 0.0876 |
| `neighbour_constant` | **0.0051** | 0.0978 |
| `neighbour_level` (the fix) | 0.0192 | **0.0845** |

**The fix does not clear its bar.** `neighbour_constant` — simply assigning
every call the average precision of the query's congeners — is better
calibrated (0.0051 vs 0.0192) and wins on 90 of 95 genomes.

What the fix does achieve, stated without inflation:

* **It repairs the original failure.** ECE improves over `global_curve` by
  +0.0174, better on 75 of 95 genomes, and `global_curve` is worse than even a
  global constant (0.0366 vs 0.0315) — so the diagnosis was right and the
  intercept was where the damage was.
* **It is the sharpest method tested.** Brier 0.0845, best of the four, better
  than `neighbour_constant` on 93 of 95 genomes.

## The proper decomposition, because ECE is not a proper scoring rule

ECE structurally favours a constant: a method that makes no distinctions puts
every call in one bin and so has far fewer opportunities to be miscalibrated.
That is a property of the metric. The Murphy decomposition
(Brier = reliability − resolution + uncertainty) separates the two things ECE
conflates, and gives the constant no such advantage.

| method | reliability ↓ | resolution ↑ |
|---|---|---|
| `global_constant` | 0.00148 | 0.00000 |
| `global_curve` | 0.00473 | 0.01428 |
| `neighbour_constant` | **0.00008** | 0.00000 |
| `neighbour_level` | 0.00168 | **0.01419** |

This does **not** overturn the verdict — `neighbour_constant` still has the
best reliability by a wide margin, so its advantage was real and not merely a
binning artefact. What the decomposition adds is the size of what is being
traded:

* `neighbour_level` reaches **the reliability of a global constant**
  (0.00168 vs 0.00148) while keeping **essentially all** the discrimination
  the fitted curve had (0.01419 vs 0.01428).
* It is **2.8× better calibrated than the attempt it replaces**
  (0.00473 → 0.00168), which is the fix working as diagnosed.
* Every method that beats it on reliability has resolution of exactly zero.

## Why the constant's win is degenerate, and why that is not an excuse

A constant assigns the same confidence to every call, so it cannot separate a
reliable call from an unreliable one — which is the entire purpose of the
exercise. Its discrimination is essentially nil: moving from `global_constant`
to `neighbour_constant` improves Brier by only +0.0014, i.e. knowing the
query's lineage buys almost nothing *unless* the margin is used too.

So there is a real calibration/sharpness trade and no method dominates:

* want a number that is right on average → `neighbour_constant`
* want a number that ranks calls by reliability → `neighbour_level`

That is a defensible engineering choice, **but it is not the result that was
aimed for**, and the pre-registered bar was ECE. A 0.019 absolute calibration
error may well be acceptable in practice; that is a judgement about the
downstream use, not evidence that the method succeeded. Recording it as a
partial result rather than a win.

## Recommendation, bounded

On the decomposition this is a usable score: reliability on par with a global
constant, discrimination on par with the fitted curve, and 2.8× the
calibration of the attempt it replaces. If the score is used to **rank or
filter** calls — which is what a confidence on an annotation is for — that
profile is adequate, and ~0.019 absolute calibration error is a small price
for the only method here that separates a reliable call from an unreliable one.

It is still not what was aimed for, and the residual risk is specific: if a
downstream step consumes the number as a literal probability rather than as an
ordering, a 2% bias compounds through it. That is exactly what module
completeness would do.

**So the decision stands: nothing ships yet.** Not because the method is bad,
but because the pre-registered bar was ECE, it was missed, and the natural
consumer inside this toolkit is the one use that the residual miscalibration
would actually hurt.

## What would actually settle it

The obvious next step is shrinkage: fit the intercept on neighbours but pull
it toward the neighbour base rate in proportion to how little neighbour data
there is, so the sharpness is kept where the evidence supports it and given up
where it does not. That is a one-parameter addition and is testable with the
same harness. Until it clears the ECE bar against `neighbour_constant`, **no
confidence score should ship in `mpph annotate`** — which is the same
conclusion the first attempt reached, now for a better-understood reason.
