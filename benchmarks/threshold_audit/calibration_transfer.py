"""Does the confidence curve transfer to a clade it was not fitted on?

The audit (`summarize_audit.py`) shows that a KO call's precision rises
monotonically with its margin above the KOfam threshold. That alone does not
license shipping a confidence score: the curve was measured on the same
genomes it describes. The question that decides whether this is a usable
method or just an observation is whether a curve fitted on one set of genomes
is *calibrated* on a genome from a different lineage.

Protocol: fit P(call is correct | margin) on training genomes, predict on a
held-out genome, and compare against the only honest baseline -- a constant
predictor set to the training genomes' overall precision, which is what a tool
implicitly asserts today when it emits a bare `1`.

Metrics are calibration metrics, not accuracy metrics:

* **ECE** (expected calibration error): mean |stated confidence - observed
  frequency| over equal-width bins. This is the quantity no incumbent reports.
* **Brier score**: mean squared error of the probability. Rewards sharpness as
  well as calibration, so a constant predictor cannot win on it by being vague.

A margin-aware model that does not beat the constant baseline on a held-out
clade means the margin does not generalise, and no confidence score should be
shipped. That is a real possible outcome of this script.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize

sys.path.insert(0, str(Path(__file__).parent))
from summarize_audit import load_scored


def _features(delta: np.ndarray) -> np.ndarray:
    """log1p compresses a margin spanning 0 to several thousand bits into a
    range where a two-parameter logistic is a reasonable shape. Margins are
    clipped at 0 because only accepted calls are being calibrated."""
    return np.log1p(np.clip(delta, 0.0, None))


def fit_logistic(delta: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x = _features(delta)

    def nll(params):
        a, b = params
        z = a + b * x
        # log(1 + exp(-|z|)) + max(-z, 0) -- the stable form of log(1+exp(-z))
        return float(np.sum(np.logaddexp(0.0, -z) + np.where(y > 0, 0.0, z)))

    res = optimize.minimize(nll, x0=np.array([0.0, 1.0]), method="Nelder-Mead",
                            options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 5000})
    return float(res.x[0]), float(res.x[1])


def predict(delta: np.ndarray, a: float, b: float) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-(a + b * _features(delta))))


def ece(p: np.ndarray, y: np.ndarray, n_bins: int = 10) -> float:
    """Expected calibration error over equal-width probability bins."""
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    idx = np.clip(np.digitize(p, edges) - 1, 0, n_bins - 1)
    total = 0.0
    for b in range(n_bins):
        sel = idx == b
        if not sel.any():
            continue
        total += sel.mean() * abs(p[sel].mean() - y[sel].mean())
    return float(total)


def brier(p: np.ndarray, y: np.ndarray) -> float:
    return float(np.mean((p - y) ** 2))


def reliability(p: np.ndarray, y: np.ndarray, n_bins: int = 10) -> pd.DataFrame:
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    idx = np.clip(np.digitize(p, edges) - 1, 0, n_bins - 1)
    rows = []
    for b in range(n_bins):
        sel = idx == b
        if not sel.any():
            continue
        rows.append({"bin_lo": edges[b], "bin_hi": edges[b + 1],
                     "n": int(sel.sum()), "mean_confidence": float(p[sel].mean()),
                     "observed_frequency": float(y[sel].mean())})
    return pd.DataFrame(rows)


def accepted(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    acc = df[df["delta"] >= 0]
    return acc["delta"].to_numpy(float), acc["is_true"].to_numpy(float)


def bootstrap_ci(p_model: np.ndarray, p_const: np.ndarray, y: np.ndarray, *,
                 n_boot: int = 2000, seed: int = 0) -> dict:
    """Percentile CIs for the *paired* improvement in ECE and Brier.

    Needed because the raw differences here span two orders of magnitude
    (~0.001 to ~0.03) and reading a 0.001 difference as a real effect would be
    a mistake. Resampling calls with replacement and recomputing both metrics
    on the same resample keeps the comparison paired, so the CI is on the
    difference rather than on each metric separately.
    """
    rng = np.random.default_rng(seed)
    n = len(y)
    d_ece = np.empty(n_boot)
    d_brier = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        ys = y[idx]
        d_ece[i] = ece(p_const[idx], ys) - ece(p_model[idx], ys)
        d_brier[i] = brier(p_const[idx], ys) - brier(p_model[idx], ys)
    return {
        "ece_improvement_lo": float(np.percentile(d_ece, 2.5)),
        "ece_improvement_hi": float(np.percentile(d_ece, 97.5)),
        "brier_improvement_lo": float(np.percentile(d_brier, 2.5)),
        "brier_improvement_hi": float(np.percentile(d_brier, 97.5)),
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit-dir", required=True)
    p.add_argument("--ground-truth-dir", required=True)
    p.add_argument("--train", nargs="+", required=True)
    p.add_argument("--test", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--out-reliability", default=None)
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    audit, gt = Path(args.audit_dir), Path(args.ground_truth_dir)

    def load(org):
        id_map = gt / f"protein_id_map_{org}.tsv"
        return load_scored(audit / f"subthreshold_{org}.tsv", org,
                           gt / f"link_ko_{org}.tsv",
                           id_map if id_map.exists() else None)

    tr = [accepted(load(o)) for o in args.train]
    tr_delta = np.concatenate([d for d, _ in tr])
    tr_y = np.concatenate([y for _, y in tr])
    a, b = fit_logistic(tr_delta, tr_y)
    base_rate = float(tr_y.mean())
    print(f"fitted on {args.train}: n={len(tr_y)}, "
          f"P = sigmoid({a:.4f} + {b:.4f}*log1p(margin)), "
          f"training precision {base_rate:.4f}")

    rows, rel_frames = [], []
    for org in args.test:
        d, y = accepted(load(org))
        p_model = predict(d, a, b)
        p_const = np.full_like(y, base_rate)
        rows.append({
            "test_organism": org, "n_calls": len(y),
            "observed_precision": float(y.mean()),
            "ece_margin_model": ece(p_model, y),
            "ece_constant_baseline": ece(p_const, y),
            "brier_margin_model": brier(p_model, y),
            "brier_constant_baseline": brier(p_const, y),
            "train": "+".join(args.train),
            **bootstrap_ci(p_model, p_const, y, n_boot=args.n_boot,
                           seed=args.seed),
        })
        rel = reliability(p_model, y)
        rel.insert(0, "test_organism", org)
        rel_frames.append(rel)

    out = pd.DataFrame(rows)
    out["ece_improvement"] = out["ece_constant_baseline"] - out["ece_margin_model"]
    out["brier_improvement"] = (out["brier_constant_baseline"]
                                - out["brier_margin_model"])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    if args.out_reliability:
        pd.concat(rel_frames, ignore_index=True).to_csv(args.out_reliability,
                                                        index=False)

    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print()
        print(out.to_string(index=False))
        print()
        print(pd.concat(rel_frames, ignore_index=True).to_string(index=False))
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
