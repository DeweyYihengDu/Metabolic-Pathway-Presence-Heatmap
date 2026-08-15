# Power of `mpph compare --tree` — the other half of the gate

`calibrate_type1.py` showed the corrected test does not reject when there is
nothing to find. That is necessary and not sufficient: **a test that never
rejects is perfectly calibrated and useless.** Power was not measured before
the feature shipped, which was an omission — the correction cannot be free,
and how much sensitivity it costs decides whether it is usable.

Run: `power_curve.py --replicates 150 --n-sims 499`, 96 cells
(2 tree shapes × 2 sizes × 2 groupings × 2 rates × 6 effect sizes).
Raw numbers: `power.csv`.

## Two groupings, because either alone misleads

Traits are simulated on the tree under a symmetric Mk model, then a genuine
group effect is added: each tip in group A flips toward presence with
probability `effect`, each tip in B toward absence. `effect = 0` is the null.

- **`clade`** — groups are the two clades descending from the root, so a real
  difference is maximally confounded with phylogeny.
- **`scattered`** — group membership is random across the tree, same group
  sizes, so a real difference is not confounded.

## Result

Mean rejection rate at nominal α = 0.05:

| effect | Fisher, clade | Fisher, scattered | **`--tree`, clade** | **`--tree`, scattered** |
|---|---|---|---|---|
| 0.0 (null) | **0.217** | 0.023 | **0.024** | **0.024** |
| 0.2 | 0.303 | 0.180 | 0.058 | 0.170 |
| 0.4 | 0.622 | 0.637 | 0.209 | 0.542 |
| 0.6 | 0.903 | 0.944 | 0.352 | 0.752 |
| 0.8 | 0.998 | 1.000 | 0.326 | 0.792 |
| 1.0 | 1.000 | 1.000 | 0.205 | 0.935 |

**Type I error is 0.024 under both groupings** — the correction is calibrated
whether or not the groups follow the tree, while Fisher is fine when they do
not (0.023) and badly anti-conservative when they do (0.217). That localises
the problem precisely: it is phylogenetic confounding, not group testing.

**The correction costs about 14 points of power where power was real.** Under
`scattered`, at a large effect, `--tree` rejects 86.3% of the time against
Fisher's 100%. That is the price, it is modest, and it is paid in the regime
where Fisher was already correct.

## The non-monotonic column is the most informative part

Under `clade`, power *falls* as the effect grows: 0.352 at effect 0.6, 0.326
at 0.8, 0.205 at 1.0. That looks wrong and is not.

A trait that separates the two root clades perfectly is exactly what **one**
evolutionary event on the root branch produces. The stronger the effect, the
more the pattern resembles a single origin — and the conditioned null, which
simulates clade-structured traits on the same tree, then generates equally
extreme patterns more often. So the evidence *against* the null weakens as the
alignment becomes perfect.

This is Maddison and FitzJohn (2015) appearing as a measured curve rather than
a citation: when a trait arose once on the branch separating the groups, the
effective sample size is one and no method can distinguish association from
coincidence. `mpph compare --tree` reports `n_state_changes` alongside the
p-value for exactly this reason — a large corrected p with
`n_state_changes == 1` is the correct answer, and the integer is what makes it
interpretable rather than puzzling.

## What this means for using the tool

- If your groups are **not** aligned with the phylogeny, expect to lose a
  little sensitivity and gain nothing much — Fisher was already calibrated
  there. The correction is cheap insurance, not a fix.
- If your groups **are** clade-structured (habitat that tracks taxonomy,
  host-associated vs free-living, any comparison of two lineages), the
  uncorrected test rejects 22% of the time on pure noise and the correction is
  doing essential work. Low power there is not a defect to engineer around;
  it is the honest amount of evidence such a design contains.

## Reproducing

```bash
python power_curve.py --replicates 150 --n-sims 499
```

Runs on a laptop in about an hour; deliberately outside CI, like the Type I
sweep. The regression guards that run in CI are `tests/test_phylo.py`,
including `test_real_effect_scattered_across_the_tree_stays_significant`,
which fails if the corrected test ever stops detecting an unconfounded
difference.
