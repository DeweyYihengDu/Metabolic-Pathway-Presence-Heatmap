"""Type I error calibration for `mpph compare --tree`.

This is a gate, not a demo. The conditioned binary null was designed against
one operating point; before it can be trusted it has to hold its nominal
level across tree sizes, feature prevalences and evolutionary rates, on an
unbalanced non-ultrametric tree as well as a tidy balanced one.

Procedure: simulate traits on the tree under a symmetric Mk model with **no
group effect**, with the two groups being the two clades descending from the
root (the worst case for phylogenetic non-independence). Any rejection is a
false positive by construction. A correctly calibrated test rejects at most
alpha of the time; the uncorrected Fisher test is reported alongside to show
the size of the problem being fixed.

    python calibrate_type1.py --replicates 400 --n-sims 1999
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from mpph.phylo import binary_phylo_test, parse_newick, simulate_mk


def balanced_newick(levels: int, blen: float = 1.0):
    labels = [f"t{i}" for i in range(2 ** levels)]
    cur = [f"{lab}:{blen}" for lab in labels]
    while len(cur) > 1:
        cur = [f"({cur[i]},{cur[i + 1]}):{blen}" for i in range(0, len(cur), 2)]
    return cur[0][:cur[0].rfind(":")] + ";", labels


def unbalanced_newick(n_tips: int, seed: int = 0):
    """A random-topology, non-ultrametric tree: the realistic shape. Built by
    repeatedly joining two random subtrees with random branch lengths."""
    rng = np.random.default_rng(seed)
    parts = [f"t{i}" for i in range(n_tips)]
    labels = list(parts)
    while len(parts) > 1:
        i, j = rng.choice(len(parts), size=2, replace=False)
        a, b = parts[int(i)], parts[int(j)]
        la, lb = rng.uniform(0.2, 2.0), rng.uniform(0.2, 2.0)
        merged = f"({a}:{la:.4f},{b}:{lb:.4f})"
        parts = [p for k, p in enumerate(parts) if k not in (int(i), int(j))]
        parts.append(merged)
    return parts[0] + ";", labels


def root_split_groups(tree):
    """The two clades descending from the root, as tip-position arrays."""
    pos = {lab: i for i, lab in enumerate(tree.tip_labels)}

    def tips_below(node):
        if tree.is_tip[node]:
            return [tree.labels[node]]
        out = []
        for child in tree.children[node]:
            out.extend(tips_below(child))
        return out

    kids = tree.children[0]
    left = tips_below(kids[0])
    right = [lab for lab in tree.tip_labels if lab not in set(left)]
    return (np.array(sorted(pos[x] for x in left)),
            np.array(sorted(pos[x] for x in right)))


def run_cell(tree, idx_a, idx_b, q_true, target_k, replicates, n_sims, seed):
    """Return (fisher_reject_rate, phylo_reject_rate, n_used)."""
    rng = np.random.default_rng(seed)
    n_used = len(idx_a) + len(idx_b)
    fisher_hits = phylo_hits = used = 0
    attempts = 0
    while used < replicates and attempts < replicates * 60:
        attempts += 1
        states = simulate_mk(tree, q_true, 1, rng)[:, 0].astype(float)
        k = int(states[np.concatenate([idx_a, idx_b])].sum())
        # Only score replicates near the prevalence stratum of interest --
        # this is what exposed the defect that a marginal average hides.
        if abs(k - target_k) > max(1, round(0.05 * n_used)):
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
    ap.add_argument("--replicates", type=int, default=300)
    ap.add_argument("--n-sims", type=int, default=1999)
    ap.add_argument("--out", default=str(Path(__file__).parent / "calibration.csv"))
    args = ap.parse_args()

    rows = []
    for shape in ("balanced", "unbalanced"):
        for levels, n_tips in ((4, 16), (5, 32), (6, 64)):
            newick = (balanced_newick(levels)[0] if shape == "balanced"
                      else unbalanced_newick(n_tips, seed=levels)[0])
            tree = parse_newick(newick)
            idx_a, idx_b = root_split_groups(tree)
            n_used = len(idx_a) + len(idx_b)
            for frac in (0.10, 0.25, 0.50):
                for q_true in (0.02, 0.1, 0.5):
                    target_k = max(1, round(frac * n_used))
                    fisher, phylo, used = run_cell(
                        tree, idx_a, idx_b, q_true, target_k,
                        args.replicates, args.n_sims, seed=hash(
                            (shape, levels, frac, q_true)) % (2 ** 31))
                    rows.append({"shape": shape, "n_tips": n_tips,
                                 "k_over_n": frac, "q_true": q_true,
                                 "fisher_type1": fisher, "phylo_type1": phylo,
                                 "n_replicates": used})
                    print(f"{shape:11s} n={n_tips:3d} k/n={frac:.2f} "
                          f"q={q_true:<5} -> Fisher {fisher:.3f}  "
                          f"phylo {phylo:.3f}  (n={used})", flush=True)

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    worst = df["phylo_type1"].max()
    n_fail = int((df["phylo_type1"] > 0.05).sum())
    print(f"\nworst phylo Type I: {worst:.3f} across {len(df)} cells; "
          f"{n_fail} cell(s) above the nominal 0.05")
    print(f"uncorrected Fisher reaches {df['fisher_type1'].max():.3f}")
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
