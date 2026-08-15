"""Statistical power of `mpph compare --tree`, the other half of the gate.

`calibrate_type1.py` establishes that the corrected test does not reject when
there is nothing to find. That is necessary and not sufficient: **a test that
never rejects is perfectly calibrated and useless.** Power was never measured
before the feature shipped, which is a real omission -- the correction cannot
be free, and how much sensitivity it costs decides whether it is usable.

Procedure, deliberately the mirror image of the Type I script. Traits are
simulated on the tree under a symmetric Mk model as before, but now a genuine
group effect is *added*: each tip in group A independently flips toward
presence with probability `effect`, each tip in B toward absence. `effect = 0`
reproduces the null and should return the Type I numbers; larger values give a
real difference of increasing size. Any failure to reject is now a false
negative.

**Two groupings are run, and reporting only one would be misleading.**

* `clade` -- the groups are the two clades descending from the root. A real
  difference is then maximally confounded with phylogeny, and a trait that
  separates the groups perfectly is exactly what *one* evolutionary event on
  the root branch produces. Low power here is the correct answer, not a
  defect: the effective sample size is one (Maddison and FitzJohn, 2015), and
  a test that rejected confidently would be wrong to.
* `scattered` -- group membership is assigned at random across the tree, so a
  real difference is *not* confounded with phylogeny. This is where power
  should be retained, and where losing it would be a genuine defect.

The contrast between the two is the result. Quoting the `clade` number alone
would understate the method badly; quoting `scattered` alone would hide what
the correction costs in the case it exists for.

Reported per cell: rejection rate of the uncorrected Fisher test and of the
corrected test. Fisher's rate is *not* a target to match -- under `clade` at
effect 0 it is mostly false positives -- it is context for what is being
given up.

    python power_curve.py --replicates 200 --n-sims 999
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).parent))

from calibrate_type1 import balanced_newick, root_split_groups, unbalanced_newick

from mpph.phylo import binary_phylo_test, parse_newick, simulate_mk


def apply_effect(states: np.ndarray, idx_a: np.ndarray, idx_b: np.ndarray,
                 effect: float, rng) -> np.ndarray:
    """Push group A toward presence and B toward absence, independently per tip.

    Applied *after* the phylogenetic simulation rather than by simulating a
    different rate per group, so the phylogenetic signal and the group effect
    are separately controlled: at `effect = 0` this is exactly the null the
    Type I script uses, and raising it adds group signal without altering the
    tree-driven correlation structure.
    """
    out = states.copy()
    flip_a = rng.uniform(size=len(idx_a)) < effect
    flip_b = rng.uniform(size=len(idx_b)) < effect
    out[idx_a[flip_a]] = 1.0
    out[idx_b[flip_b]] = 0.0
    return out


def run_cell(tree, idx_a, idx_b, q_true, effect, replicates, n_sims, seed):
    """Return (fisher_reject_rate, phylo_reject_rate, n_used)."""
    rng = np.random.default_rng(seed)
    fisher_hits = phylo_hits = used = 0
    for _ in range(replicates):
        states = simulate_mk(tree, q_true, 1, rng)[:, 0].astype(float)
        states = apply_effect(states, idx_a, idx_b, effect, rng)
        # An invariant feature carries no information for either test; both
        # correctly fail to reject, and counting it would understate power
        # for reasons that have nothing to do with the correction.
        if states.min() == states.max():
            continue
        used += 1
        a_pres = int(states[idx_a].sum())
        b_pres = int(states[idx_b].sum())
        _, fp = stats.fisher_exact([[a_pres, len(idx_a) - a_pres],
                                    [b_pres, len(idx_b) - b_pres]])
        if fp < 0.05:
            fisher_hits += 1
        res = binary_phylo_test(tree, states, idx_a, idx_b, n_sims=n_sims,
                                rng=rng, min_accepted=50)
        p = res["p_value_phylo"]
        if np.isfinite(p) and p < 0.05:
            phylo_hits += 1
    if used == 0:
        return float("nan"), float("nan"), 0
    return fisher_hits / used, phylo_hits / used, used


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--replicates", type=int, default=200)
    ap.add_argument("--n-sims", type=int, default=999)
    ap.add_argument("--out", default=str(Path(__file__).parent / "power.csv"))
    args = ap.parse_args()

    rows = []
    for shape in ("balanced", "unbalanced"):
        for levels, n_tips in ((5, 32), (6, 64)):
            newick = (balanced_newick(levels)[0] if shape == "balanced"
                      else unbalanced_newick(n_tips, seed=levels)[0])
            tree = parse_newick(newick)
            clade_a, clade_b = root_split_groups(tree)
            # Scattered groups: same sizes as the clade split so the two
            # groupings differ only in whether membership follows the tree,
            # not in sample size.
            perm = np.random.default_rng(levels).permutation(n_tips)
            scat_a = np.sort(perm[:len(clade_a)])
            scat_b = np.sort(perm[len(clade_a):])
            groupings = {"clade": (clade_a, clade_b),
                         "scattered": (scat_a, scat_b)}
            for grouping, (idx_a, idx_b) in groupings.items():
                for q_true in (0.05, 0.2):
                    for effect in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
                        fisher, phylo, used = run_cell(
                            tree, idx_a, idx_b, q_true, effect,
                            args.replicates, args.n_sims,
                            seed=hash((shape, levels, grouping, q_true,
                                       effect)) % (2 ** 31))
                        rows.append({"shape": shape, "n_tips": n_tips,
                                     "grouping": grouping,
                                     "q_true": q_true, "effect": effect,
                                     "fisher_reject": fisher,
                                     "phylo_reject": phylo,
                                     "n_replicates": used})
                        print(f"{shape:11s} n={n_tips:3d} {grouping:9s} "
                              f"q={q_true:<5} effect={effect:.1f} -> "
                              f"Fisher {fisher:.3f}  phylo {phylo:.3f}  "
                              f"(n={used})", flush=True)

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)

    print("\n=== mean rejection rate by grouping and effect ===")
    pivot = df.pivot_table(index="effect", columns="grouping",
                           values=["phylo_reject", "fisher_reject"])
    print(pivot.round(3).to_string())
    for g in ("scattered", "clade"):
        sub = df[df["grouping"] == g]
        null = sub[sub["effect"] == 0.0]["phylo_reject"].mean()
        strong = sub[sub["effect"] >= 0.8]["phylo_reject"].mean()
        print(f"\n{g:9s}: Type I {null:.3f}, power at effect >= 0.8 {strong:.3f}")
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
