"""Calibrate a genome from its own score distribution, needing no relatives.

The held-out check killed the neighbour approach for its intended use: the
calibration map is genus-specific, and a novel environmental genome usually
has no congener carrying curated KEGG annotation -- exactly the case
`mpph annotate` exists to serve.

But look at *what* congeners were supplying. `neighbour_constant` -- borrowing
a relative's mean precision, i.e. the **level** -- was calibrated to the floor
at genus rank and useless beyond it. The shape was already transferable from
the pooled fit. So the missing ingredient is one number per genome: its base
rate.

And a genome's base rate is largely a consequence of how well its proteins
match KOfam's profiles, which shows up directly in **its own margin
distribution** -- a quantity observable with no labels and no relatives at
all. If base rate is predictable from those observables, congeners were never
the requirement; they were a proxy for something the query already carries.

This fits that regression across reference genomes and evaluates it
leave-one-genome-out, against what a novel genome would actually get: the
pooled global curve, and the neighbour route restricted to clade-level peers.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "annotation"))
from calibration_transfer import _features, brier, ece, fit_logistic, predict
from compare_annotation import load_ground_truth, load_protein_id_map
from neighbour_calibrator import (
    isotonic_on_neighbours,
    murphy_decomposition,
    neighbours_at_rank,
    null_ece,
)

FEATURES = ("mean_margin", "median_margin", "q25_margin", "q75_margin",
            "frac_low_margin", "log_n_calls", "mean_log1p_margin")


def observable_features(delta: np.ndarray) -> dict:
    """Summary statistics of a genome's accepted-call margins.

    Deliberately label-free: every one is computable from the search output
    alone, so they exist for a genome nobody has ever annotated. That is the
    entire point -- if they predict the base rate, calibration no longer
    requires a curated relative.
    """
    d = np.asarray(delta, float)
    return {
        "mean_margin": float(d.mean()),
        "median_margin": float(np.median(d)),
        "q25_margin": float(np.percentile(d, 25)),
        "q75_margin": float(np.percentile(d, 75)),
        "frac_low_margin": float((d < 20).mean()),
        "log_n_calls": float(np.log1p(len(d))),
        "mean_log1p_margin": float(np.log1p(d).mean()),
    }


def fit_base_rate_model(rows_x, rows_y):
    """Ridge regression of base rate on observable features.

    Ridge rather than ordinary least squares because the margin summaries are
    strongly collinear -- mean, median and both quartiles move together -- and
    with ~94 training genomes an unregularised fit would chase that
    collinearity. Implemented directly so this needs no scikit-learn, which is
    not installed everywhere this benchmark runs.
    """
    x = np.asarray(rows_x, float)
    y = np.asarray(rows_y, float)
    mu, sigma = x.mean(0), x.std(0)
    sigma[sigma == 0] = 1.0
    xs = np.hstack([np.ones((len(x), 1)), (x - mu) / sigma])
    penalty = np.eye(xs.shape[1])
    penalty[0, 0] = 0.0                    # never penalise the intercept
    coef = np.linalg.solve(xs.T @ xs + penalty, xs.T @ y)
    return {"coef": coef, "mu": mu, "sigma": sigma}


def predict_base_rate(model, feats: dict) -> float:
    x = np.array([feats[f] for f in FEATURES], float)
    xs = np.concatenate([[1.0], (x - model["mu"]) / model["sigma"]])
    return float(np.clip(xs @ model["coef"], 0.01, 0.99))


def intercept_for_target_rate(delta, slope, target, lo=-30.0, hi=30.0):
    """Intercept making the mean predicted probability equal `target`.

    Not a likelihood fit. The earlier version reused `fit_intercept_only` with
    a constant array in place of labels, which that function reads as "every
    call is a true positive" -- it duly drove every probability to 1 and
    produced a degenerate constant (resolution 0.00000, Brier == ECE). This
    solves the intended equation directly: mean(sigmoid(a + slope*x)) =
    target, by bisection on the monotone left-hand side over the *query's own*
    margins, so the predicted level is imposed while the pooled shape is kept.
    """
    x = _features(np.asarray(delta, float))

    def mean_p(a):
        return float(np.mean(1.0 / (1.0 + np.exp(-(a + slope * x)))))

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mean_p(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def load_panel(root: Path, tax: dict) -> dict:
    data = {}
    for org in sorted(tax):
        h = root / "subthreshold" / f"subthreshold_{org}.tsv"
        g = root / "ground_truth" / f"link_ko_{org}.tsv"
        m = root / "ground_truth" / f"protein_id_map_{org}.tsv"
        if not (h.exists() and g.exists() and m.exists()):
            continue
        df = pd.read_csv(h, sep="\t")
        df["gene"] = df["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        idm = load_protein_id_map(m, org)
        truth = load_ground_truth(g, org, idm)
        df = df[df["gene"].isin(set(idm.values()))]
        acc = df[df["delta"] >= 0]
        if len(acc) < 200:
            continue
        tp = {(a, b) for a, ks in truth.items() for b in ks}
        y = np.array([(a, b) in tp for a, b in zip(acc["gene"], acc["ko"])], float)
        if y.mean() in (0.0, 1.0):
            continue
        data[org] = (acc["delta"].to_numpy(float), y)
    return data


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panel-dir", required=True)
    p.add_argument("--taxonomy", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    with open(args.taxonomy, encoding="utf-8") as fh:
        tax = {r["org"]: r for r in csv.DictReader(fh, delimiter="\t")}
    data = load_panel(Path(args.panel_dir), tax)
    orgs = sorted(data)
    feats = {o: observable_features(data[o][0]) for o in orgs}
    print(f"{len(orgs)} genomes usable", flush=True)

    rows = []
    for q in orgs:
        d_q, y_q = data[q]
        others = [o for o in orgs if o != q]
        d_all = np.concatenate([data[o][0] for o in others])
        y_all = np.concatenate([data[o][1] for o in others])
        a_g, b_g = fit_logistic(d_all, y_all)

        # Base-rate model trained on the other genomes only.
        model = fit_base_rate_model(
            [[feats[o][f] for f in FEATURES] for o in others],
            [data[o][1].mean() for o in others])
        rate_hat = predict_base_rate(model, feats[q])

        preds = {"global_curve": predict(d_q, a_g, b_g),
                 "global_constant": np.full_like(y_q, float(y_all.mean()))}
        # Level from the genome's own observable score distribution, shape
        # from the pooled fit. No relatives involved at any point.
        a_self = intercept_for_target_rate(d_q, b_g, rate_hat)
        preds["self_level"] = predict(d_q, a_self, b_g)
        preds["self_constant"] = np.full_like(y_q, rate_hat)

        # What a novel genome would actually get from the neighbour route.
        peers, _ = neighbours_at_rank(q, tax, orgs, "clade", 2)
        if peers:
            d_nb = np.concatenate([data[o][0] for o in peers])
            y_nb = np.concatenate([data[o][1] for o in peers])
            preds["clade_isotonic"] = isotonic_on_neighbours(
                predict(d_nb, a_g, b_g), y_nb, preds["global_curve"])

        for name, p_hat in preds.items():
            rows.append({"query": q, "method": name,
                         "true_base_rate": float(y_q.mean()),
                         "pred_base_rate": rate_hat,
                         "ece": ece(p_hat, y_q), "brier": brier(p_hat, y_q),
                         "null_ece": null_ece(p_hat, n_draws=60,
                                              seed=abs(hash(q)) % 10000),
                         **murphy_decomposition(p_hat, y_q)})

    out = pd.DataFrame(rows)
    out["ece_above_floor"] = out["ece"] - out["null_ece"]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)

    print("\n=== leave-one-genome-out ===")
    print(out.groupby("method")[
        ["ece", "null_ece", "ece_above_floor", "brier", "reliability",
         "resolution"]].mean().round(5).to_string())

    per = out[out["method"] == "self_level"]
    err = (per["pred_base_rate"] - per["true_base_rate"]).abs()
    print(f"\nbase rate predicted from the genome's OWN scores: "
          f"mean abs error {err.mean():.5f}")
    print("(neighbour_constant needed congeners and achieved 0.00514 there; "
          "at clade rank it was far worse)")
    for m, g in out.groupby("method"):
        print(f"  {m:18s} at/below its floor on "
              f"{int((g['ece_above_floor'] <= 0).sum())}/{len(g)}")
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
