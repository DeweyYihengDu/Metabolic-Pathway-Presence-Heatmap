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


def neighbours_of(query, tax, orgs, min_n):
    """Closest available relatives, walking outward until enough are found.

    Returns (list_of_orgs, rank_used). A query whose genus has no other member
    falls back to clade, then domain -- reporting which rank was actually used
    matters, because a "neighbour" calibration that silently fell back to
    domain is just the global fit wearing a different name.
    """
    for rank in RANKS:
        peers = [o for o in orgs
                 if o != query and tax[o][rank] and tax[o][rank] == tax[query][rank]]
        if len(peers) >= min_n:
            return peers, rank
    return [o for o in orgs if o != query], "all"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panel-dir", required=True)
    p.add_argument("--taxonomy", required=True)
    p.add_argument("--min-neighbours", type=int, default=2)
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

        peers, rank = neighbours_of(query, tax, orgs, args.min_neighbours)
        d_nb = np.concatenate([data[o][0] for o in peers])
        y_nb = np.concatenate([data[o][1] for o in peers])
        a_nb = fit_intercept_only(d_nb, y_nb, b_glob)

        preds = {
            "global_curve": predict(d_q, a_glob, b_glob),
            "global_constant": np.full_like(y_q, float(y_all.mean())),
            "neighbour_constant": np.full_like(y_q, float(y_nb.mean())),
            "neighbour_level": predict(d_q, a_nb, b_glob),
        }
        for name, p_hat in preds.items():
            rows.append({"query": query, "rank_used": rank, "n_peers": len(peers),
                         "method": name, "ece": ece(p_hat, y_q),
                         "brier": brier(p_hat, y_q)})

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)

    print("\n=== leave-one-genome-out, mean over queries ===")
    print(out.groupby("method")[["ece", "brier"]].mean().round(4).to_string())
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
