"""Precision of KOfam hits as a function of distance from the threshold.

Joins `collect_subthreshold.py`'s output against KEGG's own curated
`link/ko/<org>` assignments and reports, per bin of (score - threshold), what
fraction of (gene, KO) pairs KEGG also assigns.

Reads for the interpretation, stated before looking at any number:

* `delta >= 0` is the band every tool already keeps. Its precision is the
  reference level -- this repo measures ~0.90 for E. coli.
* `delta < 0` is the band every tool discards. If precision there falls
  immediately to the background rate, KOfam's threshold is well placed and
  there is nothing to recover. If it decays slowly, real assignments are being
  thrown away.

**The measured precision is a lower bound, not a point estimate.** KEGG's
reference is itself similarity-derived and incomplete, so a sub-threshold hit
counted here as a false positive may be a genuine function KEGG has not
recorded. That biases this number downward, which is the conservative
direction for the argument but must be reported as a bound.

Two controls are computed rather than assumed:

* `--shuffle-control` recomputes precision after permuting the KO labels
  across hits, giving the background rate a bin would show under no signal at
  all. Comparing against 0 is wrong; comparing against this is right.
* Eukaryotic proteomes are restricted to the proteins KEGG's reference can
  actually speak to (`--protein-id-map`), reusing the same correction the
  annotation benchmark applies, since isoforms would otherwise manufacture
  false positives.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "annotation"))
from compare_annotation import load_ground_truth, load_protein_id_map


def load_hits(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # Match the annotation benchmark's id handling: the FASTA carries RefSeq
    # version suffixes, KEGG's conv table does not.
    df["gene"] = df["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    return df


def annotate_truth(df: pd.DataFrame, truth: dict[str, set[str]]) -> pd.DataFrame:
    true_pairs = {(g, k) for g, kos in truth.items() for k in kos}
    df = df.copy()
    df["is_true"] = [
        (g, k) in true_pairs for g, k in zip(df["gene"], df["ko"])
    ]
    return df


def bin_precision(df: pd.DataFrame, *, width: float, lo: float,
                  hi: float) -> pd.DataFrame:
    edges = np.arange(lo, hi + width, width)
    idx = np.digitize(df["delta"].to_numpy(), edges) - 1
    rows = []
    for b in range(len(edges) - 1):
        sel = df[idx == b]
        if len(sel) == 0:
            continue
        n_true = int(sel["is_true"].sum())
        rows.append({
            "delta_lo": edges[b], "delta_hi": edges[b + 1],
            "n_hits": len(sel), "n_true": n_true,
            "precision": n_true / len(sel),
        })
    return pd.DataFrame(rows)


def shuffle_control(df: pd.DataFrame, truth: dict[str, set[str]], *,
                    seed: int, width: float, lo: float, hi: float) -> pd.DataFrame:
    """Precision after permuting KO labels -- the rate a bin shows with no
    real signal. Preserves each bin's size and the KO frequency distribution,
    so it is a like-for-like background rather than an assumed zero."""
    rng = np.random.default_rng(seed)
    shuffled = df.copy()
    shuffled["ko"] = rng.permutation(shuffled["ko"].to_numpy())
    return bin_precision(annotate_truth(shuffled, truth), width=width, lo=lo, hi=hi)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hits", required=True, help="collect_subthreshold.py output")
    p.add_argument("--org", required=True)
    p.add_argument("--ground-truth", required=True)
    p.add_argument("--protein-id-map", default=None)
    p.add_argument("--bin-width", type=float, default=5.0)
    p.add_argument("--delta-lo", type=float, default=-60.0)
    p.add_argument("--delta-hi", type=float, default=60.0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    df = load_hits(args.hits)
    id_map = (load_protein_id_map(args.protein_id_map, args.org)
              if args.protein_id_map else None)
    truth = load_ground_truth(args.ground_truth, args.org, id_map)

    if id_map is not None:
        evaluable = set(id_map.values())
        before = len(df)
        df = df[df["gene"].isin(evaluable)]
        if len(df) < before:
            print(f"Restricted to the {len(evaluable)} proteins KEGG's reference "
                  f"covers: {before - len(df)} hit(s) on unrepresented isoforms "
                  f"excluded.")

    scored = annotate_truth(df, truth)
    obs = bin_precision(scored, width=args.bin_width,
                        lo=args.delta_lo, hi=args.delta_hi)
    ctl = shuffle_control(scored, truth, seed=args.seed, width=args.bin_width,
                          lo=args.delta_lo, hi=args.delta_hi)
    merged = obs.merge(ctl[["delta_lo", "precision"]], on="delta_lo",
                       how="left", suffixes=("", "_shuffled"))
    merged.insert(0, "organism", args.org)
    merged["excess_over_shuffled"] = (merged["precision"]
                                      - merged["precision_shuffled"])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, index=False)

    kept = scored[scored["delta"] >= 0]
    disc = scored[scored["delta"] < 0]
    print(f"\n=== {args.org} ===")
    print(f"kept band    (delta >= 0): {len(kept):7d} hits, "
          f"precision {kept['is_true'].mean():.4f}")
    print(f"discarded    (delta <  0): {len(disc):7d} hits, "
          f"precision {disc['is_true'].mean():.4f}, "
          f"{int(disc['is_true'].sum())} KEGG-confirmed assignment(s) discarded")
    for band in (5.0, 10.0, 20.0):
        sel = scored[(scored["delta"] < 0) & (scored["delta"] >= -band)]
        if len(sel):
            print(f"  within {band:4.0f} bits below: {len(sel):6d} hits, "
                  f"precision {sel['is_true'].mean():.4f}, "
                  f"{int(sel['is_true'].sum())} recoverable")
    print()
    print(merged.to_string(index=False))
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
