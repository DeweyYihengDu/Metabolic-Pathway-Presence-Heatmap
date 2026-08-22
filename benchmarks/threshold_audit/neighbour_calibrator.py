"""The fix the scale result pointed to: neighbour level + global shape.

The first calibration attempt fitted one global logistic (both intercept and
slope) on two bacteria and failed its transfer test -- ECE never improved and
was significantly worse on the archaeon. The 95-genome panel then showed why,
in two measurements that only make sense together:

* **Shape transfers.** The margin model beats a constant on Brier at every
  relatedness tier, so the *slope* -- how fast reliability rises with margin --
  is a property of the HMM score and generalises.
* **Level does not.** ECE degrades 2.7x from same-genus to same-domain, and at
  the genus tier a plain constant set to a relative's average precision is
  *better* calibrated than the fitted curve (0.0062 vs 0.0173). Close relatives
  carry the right base rate; a global fit carries the wrong one.

So the two halves of a logistic should be estimated from different places.
This fits the slope by pooling every reference genome, and the intercept from
the query's nearest relatives only. Neither the failed global curve nor the
neighbour constant does that.

Four methods are compared under leave-one-genome-out, each predicting on a
genome that contributed nothing to its own fit:

* `global_curve`      -- intercept and slope both pooled (the attempt that failed)
* `global_constant`   -- every call gets the pooled mean precision
* `neighbour_constant`-- every call gets the neighbours' mean precision
* `neighbour_level`   -- pooled slope, intercept refit on neighbours  **(the fix)**

Success is not "beats the others on Brier" -- the margin model already did
that and it was not enough. It is **lower ECE than every baseline**, since
uncalibrated confidence is precisely what made the first attempt unshippable.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "annotation"))
from calibration_transfer import _features, brier, ece, fit_logistic, predict
from compare_annotation import load_ground_truth, load_protein_id_map

RANKS = ("genus", "clade", "domain")


def murphy_decomposition(p, y, n_bins: int = 10):
    """Brier = reliability - resolution + uncertainty (Murphy 1973).

    Needed because ECE structurally favours a constant predictor: a method
    that makes no distinctions puts every call in one bin and so has far
    fewer opportunities to be miscalibrated. That is a property of the metric,
    not evidence the constant is better calibrated in any useful sense.

    The decomposition separates the two things ECE conflates:

    * **reliability** -- do stated confidences match observed frequencies.
      This is the calibration comparison, and it is the fair one because a
      constant gets no structural advantage from having only one bin.
    * **resolution** -- how far bin frequencies depart from the base rate,
      i.e. how much the method actually discriminates. Exactly 0 for any
      constant predictor, by construction.

    Lower reliability is better; higher resolution is better.
    """
    p = np.asarray(p, float)
    y = np.asarray(y, float)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    idx = np.clip(np.digitize(p, edges) - 1, 0, n_bins - 1)
    base = y.mean()
    reliability = resolution = 0.0
    for b in range(n_bins):
        sel = idx == b
        if not sel.any():
            continue
        w = sel.mean()
        reliability += w * (p[sel].mean() - y[sel].mean()) ** 2
        resolution += w * (y[sel].mean() - base) ** 2
    return {"reliability": float(reliability), "resolution": float(resolution),
            "uncertainty": float(base * (1 - base))}


def platt_on_neighbours(p_raw_nb, y_nb):
    """Fit a 2-parameter Platt correction on the neighbour set.

    `fit_intercept_only` moves the curve up or down but cannot change how
    steeply confidence rises, so if the pooled slope is slightly wrong for
    this lineage the error shows up *within* probability bins -- which is
    exactly the residual the intercept-only version could not remove
    (reliability 0.00168 against a constant's 0.00008).

    Platt scaling fits `sigmoid(A * logit(p) + B)` on held-out data, so it can
    rescale the slope as well as shift the level, while still spending only
    two parameters on the (pooled, and therefore reasonably large) neighbour
    sample.
    """
    eps = 1e-6
    z = np.log(np.clip(p_raw_nb, eps, 1 - eps) / (1 - np.clip(p_raw_nb, eps, 1 - eps)))

    def nll(params):
        a, b = params
        t = a * z + b
        return float(np.sum(np.logaddexp(0.0, -t) + np.where(y_nb > 0, 0.0, t)))

    res = optimize.minimize(nll, x0=np.array([1.0, 0.0]), method="Nelder-Mead",
                            options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 5000})
    return float(res.x[0]), float(res.x[1])


def apply_platt(p_raw, a, b):
    eps = 1e-6
    q = np.clip(p_raw, eps, 1 - eps)
    return 1.0 / (1.0 + np.exp(-(a * np.log(q / (1 - q)) + b)))


def isotonic_on_neighbours(p_raw_nb, y_nb, p_raw_query):
    """Non-parametric monotone recalibration fitted on neighbours.

    Isotonic regression attacks reliability directly -- it maps predicted to
    observed frequency with no functional form at all -- at the cost of being
    the easiest of these to overfit. Included precisely because it is the
    strongest available attempt: if even isotonic recalibration on close
    relatives cannot beat a constant's reliability, the shortfall is not a
    matter of choosing a better correction.

    Implemented here (pool-adjacent-violators) rather than imported, so this
    benchmark does not require installing scikit-learn into whatever
    environment it runs in. PAVA is exact, not an approximation.
    """
    order = np.argsort(p_raw_nb, kind="mergesort")
    xs = np.asarray(p_raw_nb, float)[order]
    ys = np.asarray(y_nb, float)[order]

    # Each block holds (weighted sum, weight); merge left while the running
    # means are non-monotone.
    vals: list[float] = []
    wts: list[float] = []
    for v in ys:
        vals.append(v)
        wts.append(1.0)
        while len(vals) > 1 and vals[-2] / wts[-2] > vals[-1] / wts[-1]:
            v2, w2 = vals.pop(), wts.pop()
            vals[-1] += v2
            wts[-1] += w2
    fitted = np.concatenate([np.full(int(w), v / w) for v, w in zip(vals, wts)])

    # Step function evaluated at the query's raw scores, clipped at the ends.
    return np.clip(np.interp(np.asarray(p_raw_query, float), xs, fitted),
                   0.0, 1.0)


def null_ece(p_hat, n_draws=200, seed=0, n_bins=10):
    """The ECE a *perfectly calibrated* predictor with these probabilities gets.

    ECE is not comparable across methods that spread their predictions
    differently, and the reason is finite-sample noise rather than anything
    about the methods. A constant puts all N calls in one bin, so its observed
    frequency is estimated from N samples; a predictor with resolution spreads
    them over ~10 bins and estimates each from ~N/10. The expected |stated -
    observed| gap therefore scales like sqrt(bin occupancy), so a spread
    predictor is penalised for spreading, not for being wrong.

    This measures that floor directly instead of arguing about it: labels are
    drawn from the model's own probabilities, which makes it perfectly
    calibrated by construction, and ECE is recomputed. Any observed ECE at or
    near this value is as calibrated as ECE can detect at this sample size.
    """
    rng = np.random.default_rng(seed)
    p_hat = np.asarray(p_hat, float)
    return float(np.mean([
        ece((p_hat), (rng.uniform(size=p_hat.size) < p_hat).astype(float), n_bins)
        for _ in range(n_draws)]))


def bagged_isotonic(p_raw_nb, y_nb, p_raw_query, *, n_bags=50, seed=0):
    """Isotonic recalibration, bagged to cut its variance.

    The diagnosis that motivates this: plain isotonic's *level* is already
    better than a neighbour constant's (mean absolute base-rate error 0.00428
    vs 0.00514), so its remaining calibration gap is entirely **shape** error
    inside probability bins. Isotonic is a high-variance estimator -- an
    unconstrained step function fitted to whatever noise the neighbour sample
    happens to carry -- and shape error from variance is exactly what bagging
    removes.

    Averaging monotone functions preserves monotonicity, so the bagged fit is
    still a valid calibration map; it is simply smoother, with steps supported
    by evidence that recurs across resamples rather than by single points.
    """
    rng = np.random.default_rng(seed)
    p_raw_nb = np.asarray(p_raw_nb, float)
    y_nb = np.asarray(y_nb, float)
    n = len(y_nb)
    acc = np.zeros(len(p_raw_query), float)
    for _ in range(n_bags):
        idx = rng.integers(0, n, n)
        acc += isotonic_on_neighbours(p_raw_nb[idx], y_nb[idx], p_raw_query)
    return acc / n_bags


def fit_intercept_only(delta, y, slope):
    """Refit only the intercept, holding the pooled slope fixed.

    This is the whole idea in one function: the slope is the part that
    transfers, so it is imported rather than re-estimated, and only the level
    is learned from the (often small) neighbour sample. Fitting both on a
    handful of relatives would just reintroduce the variance the pooled slope
    exists to avoid.
    """
    x = _features(delta)

    def nll(a):
        z = a[0] + slope * x
        return float(np.sum(np.logaddexp(0.0, -z) + np.where(y > 0, 0.0, z)))

    res = optimize.minimize(nll, x0=np.array([0.0]), method="Nelder-Mead",
                            options={"xatol": 1e-8, "fatol": 1e-8})
    return float(res.x[0])


def neighbours_at_rank(query, tax, orgs, rank, min_n):
    """Neighbours at exactly `rank`, **excluding everything closer**.

    The held-out-clade check the shipping decision depends on. Isotonic is the
    most overfit-prone method tried here, so its advantage has to be
    reconfirmed when its calibration set is not the query's own congeners --
    otherwise the measured edge could be memorised lineage detail rather than
    a transferable calibration map.

    `clade` therefore excludes same-genus genomes and `domain` excludes
    same-clade ones, rather than the outward-walking `neighbours_of`, which
    would quietly hand back congeners and answer a different question.
    """
    closer = {"genus": [], "clade": ["genus"], "domain": ["genus", "clade"]}[rank]
    peers = []
    for o in orgs:
        if o == query or not tax[o][rank] or tax[o][rank] != tax[query][rank]:
            continue
        if any(tax[o][c] and tax[o][c] == tax[query][c] for c in closer):
            continue          # too close: excluded by construction
        peers.append(o)
    return (peers, rank) if len(peers) >= min_n else (None, rank)


def neighbours_of(query, tax, orgs, min_n, *, allow_widening=False):
    """Congeners, or nothing.

    This used to walk outward to clade and then domain when a genus had too
    few members. The held-out check makes that indefensible: calibration
    transfers within a genus and collapses immediately outside it, so widening
    returns a confident-looking score that is measurably uncalibrated --
    at clade rank every method sits ~0.02 above its own ECE floor, no better
    than pooling all references.

    So the default is now to refuse: a query with too few congeners returns
    ``(None, rank_reached)`` and the caller must handle it rather than receive
    a silently degraded answer. ``allow_widening=True`` restores the old
    behaviour for experiments that deliberately want distant peers, and is
    never the default.
    """
    peers = [o for o in orgs
             if o != query and tax[o]["genus"]
             and tax[o]["genus"] == tax[query]["genus"]]
    if len(peers) >= min_n:
        return peers, "genus"
    if not allow_widening:
        return None, "genus"
    for rank in RANKS[1:]:
        wider = [o for o in orgs
                 if o != query and tax[o][rank] and tax[o][rank] == tax[query][rank]]
        if len(wider) >= min_n:
            return wider, rank
    return [o for o in orgs if o != query], "all"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panel-dir", required=True)
    p.add_argument("--taxonomy", required=True)
    p.add_argument("--min-neighbours", type=int, default=2)
    p.add_argument("--neighbour-rank", default="auto",
                   choices=["auto", "genus", "clade", "domain"],
                   help="auto walks outward from genus until enough peers are "
                        "found. Naming a rank forces exactly that rank and "
                        "EXCLUDES everything closer -- the held-out-clade "
                        "check: does the edge survive without congeners?")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    root = Path(args.panel_dir)
    with open(args.taxonomy, encoding="utf-8") as fh:
        tax = {r["org"]: r for r in csv.DictReader(fh, delimiter="\t")}

    data = {}
    for org in sorted(tax):
        hits = root / "subthreshold" / f"subthreshold_{org}.tsv"
        gt = root / "ground_truth" / f"link_ko_{org}.tsv"
        idm = root / "ground_truth" / f"protein_id_map_{org}.tsv"
        if not (hits.exists() and gt.exists() and idm.exists()):
            continue
        df = pd.read_csv(hits, sep="\t")
        df["gene"] = df["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        m = load_protein_id_map(idm, org)
        truth = load_ground_truth(gt, org, m)
        df = df[df["gene"].isin(set(m.values()))]
        acc = df[df["delta"] >= 0]
        if len(acc) < 200:
            continue
        tp = {(g, k) for g, kos in truth.items() for k in kos}
        y = np.array([(g, k) in tp for g, k in zip(acc["gene"], acc["ko"])], float)
        if y.mean() in (0.0, 1.0):
            continue
        data[org] = (acc["delta"].to_numpy(float), y)
    orgs = sorted(data)
    print(f"{len(orgs)} genomes usable", flush=True)

    rows = []
    for query in orgs:
        d_q, y_q = data[query]
        others = [o for o in orgs if o != query]
        # Pooled fit excludes the query entirely -- leave-one-genome-out.
        d_all = np.concatenate([data[o][0] for o in others])
        y_all = np.concatenate([data[o][1] for o in others])
        a_glob, b_glob = fit_logistic(d_all, y_all)

        if args.neighbour_rank == "auto":
            peers, rank = neighbours_of(query, tax, orgs, args.min_neighbours)
            if peers is None:
                # Refusing is the point; see neighbours_of.
                print(f"  ({query}: fewer than {args.min_neighbours} congeners "
                      f"-- skipped rather than widened)")
                continue
        else:
            peers, rank = neighbours_at_rank(query, tax, orgs,
                                             args.neighbour_rank,
                                             args.min_neighbours)
            if peers is None:
                continue      # no peers at this rank; skipped, not widened

        d_nb = np.concatenate([data[o][0] for o in peers])
        y_nb = np.concatenate([data[o][1] for o in peers])
        a_nb = fit_intercept_only(d_nb, y_nb, b_glob)

        # Raw pooled-curve scores, on neighbours and on the query, are the
        # input both recalibrators correct.
        p_raw_nb = predict(d_nb, a_glob, b_glob)
        p_raw_q = predict(d_q, a_glob, b_glob)
        pa, pb = platt_on_neighbours(p_raw_nb, y_nb)

        preds = {
            "global_curve": p_raw_q,
            "global_constant": np.full_like(y_q, float(y_all.mean())),
            "neighbour_constant": np.full_like(y_q, float(y_nb.mean())),
            "neighbour_level": predict(d_q, a_nb, b_glob),
            "neighbour_platt": apply_platt(p_raw_q, pa, pb),
            "neighbour_isotonic": isotonic_on_neighbours(p_raw_nb, y_nb, p_raw_q),
            "neighbour_isotonic_bagged": bagged_isotonic(p_raw_nb, y_nb, p_raw_q),
        }
        for name, p_hat in preds.items():
            rows.append({"query": query, "rank_used": rank, "n_peers": len(peers),
                         "method": name, "ece": ece(p_hat, y_q),
                         "brier": brier(p_hat, y_q),
                         "null_ece": null_ece(p_hat, n_draws=60,
                                              seed=abs(hash(query)) % 10000),
                         **murphy_decomposition(p_hat, y_q)})

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)

    print("\n=== leave-one-genome-out, mean over queries ===")
    summary = out.groupby("method")[
        ["ece", "null_ece", "brier", "reliability", "resolution"]].mean()
    summary["ece_above_floor"] = summary["ece"] - summary["null_ece"]
    print(summary.round(5).to_string())
    print("(reliability: lower is better -- the fair calibration comparison.)")
    print("(resolution:  higher is better, and is 0 for any constant.)")
    print("\n=== by the rank the neighbours came from ===")
    print(out.groupby(["rank_used", "method"])[["ece", "brier"]]
          .agg(["mean", "count"]).round(4).to_string())

    # Paired comparison: the same query genome under two methods.
    wide = out.pivot(index="query", columns="method", values="ece")
    fix, base = "neighbour_level", "global_curve"
    d = wide[base] - wide[fix]
    print(f"\nECE improvement of {fix} over {base}: mean {d.mean():+.4f}, "
          f"better on {(d > 0).sum()}/{len(d)} genomes")
    d2 = wide["neighbour_constant"] - wide[fix]
    print(f"ECE improvement of {fix} over neighbour_constant: mean {d2.mean():+.4f}, "
          f"better on {(d2 > 0).sum()}/{len(d2)} genomes")
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
