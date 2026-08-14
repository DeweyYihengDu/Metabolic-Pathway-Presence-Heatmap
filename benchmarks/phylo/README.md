# Type I error calibration for `mpph compare --tree`

This is the gate the phylogenetic correction had to pass before shipping, not
a showcase. `calibrate_type1.py` simulates traits on a tree under a symmetric
Mk model with **no group effect whatsoever**, with the two compared groups
being the two clades descending from the root — the worst case for
phylogenetic non-independence. Every rejection is therefore a false positive
by construction, and a correctly calibrated test must reject at most 5% of
the time at a nominal 5% threshold.

## Why the grid is stratified by prevalence

An early version of the null was **not** conditioned on the observed number
of present genomes. Averaged over all prevalences it looked fine (~5%), and
that average was hiding the defect: it was absurdly conservative on
extreme-prevalence features (which nobody reports) and roughly 4x
anti-conservative on the intermediate-prevalence features that produce every
reportable result. Stratifying by `k/n` is what exposed it. The shipped
version rejection-samples replicates whose presence count matches the
observed one, exactly as Fisher's exact test conditions on its margins.

Cells report `nan` where the requested prevalence stratum is unreachable at
that evolutionary rate — e.g. `k/n = 0.10` under `q = 0.5`, where the trait
equilibrates near 50% and a 10%-prevalence draw essentially never occurs.
That is a property of the simulation, not a failure.

## Reproducing

```bash
python calibrate_type1.py --replicates 250 --n-sims 999
```

54 cells: {balanced, unbalanced non-ultrametric} x {16, 32, 64 tips} x
{10%, 25%, 50% prevalence} x {q = 0.02, 0.1, 0.5}. Results land in
`calibration.csv`; `calibration.log` is the run transcript.

Runtime is a few tens of minutes — it is a validation script, deliberately
kept out of CI (`pytest benchmarks/` does not pick it up). The regression
guards that *do* run in CI are `tests/test_phylo.py`, in particular
`test_clade_confounded_binary_feature_is_not_significant_after_correction`,
which fails if the conditioning ever regresses.

## Result (250 replicates/cell, 999 simulated replicates per test)

Gate: **passed.**

| | uncorrected Fisher | `--tree` corrected |
|---|---|---|
| worst Type I across cells | **0.948** | **0.072** |
| median Type I across cells | 0.086 | **0.000** |
| cells above nominal 0.05 | 41 of 50 | **1 of 50** |

50 of the 54 cells were evaluable; 4 requested prevalence strata are
unreachable at their rate (see above).

The single cell above nominal (unbalanced tree, 64 tips, k/n = 0.25,
q = 0.02) rejected 18/250 = 0.072. Under a true rate of 0.05 that has
P(>= 18) = 0.079, and across 50 cells you would expect **3.9** cells to
exceed 0.05 by chance — one is fewer than chance expectation, not evidence
of miscalibration. The correction is if anything conservative: its median
Type I across the grid is 0.

Raw numbers: `calibration.csv`.
