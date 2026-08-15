"""Does the KO-level precision gain survive into pathway-level accuracy?

Higher precision on (gene, KO) pairs is not automatically an improvement to
what this toolkit actually reports. Arbitration removes ~1.3 points of recall,
and a removed *true* KO can break a KEGG module outright -- so module
completeness could get **worse** even as KO precision gets better. That is a
real possible outcome and this script is written to detect it, not to confirm
the improvement.

Design: for each benchmark genome, compute KEGG module completeness three ways
from the same module definitions --

* **truth**  -- from KEGG's own curated KO set for that organism
* **all**    -- from `mpph annotate` with `--multi-ko-policy all` (= KofamScan)
* **best**   -- from `mpph annotate` with `--multi-ko-policy best`

-- and score `all` and `best` against `truth` over the modules KEGG can speak
to. Reported per genome:

* **MAE** of completeness against truth (lower is better)
* **signed bias** (mean of estimate - truth), which separates "closer" from
  "systematically lower": arbitration removing true KOs would show as a
  negative bias even if MAE improved
* **modules broken** -- modules complete in truth but not in the estimate, the
  failure mode that matters most for a metabolic-capability claim

A NaN completeness means "not assessed" throughout this codebase and is
excluded pairwise rather than filled, so the comparison is over modules both
sides actually scored.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "annotation"))

from compare_annotation import load_ground_truth, load_protein_id_map

from mpph.kegg import DEFAULT_CACHE, make_session
from mpph.modules import fetch_module_definitions, list_modules, module_completeness

COMPLETE = 1.0  # a module counts as "present" only when fully complete


def kos_from_mapper(path: Path) -> set[str]:
    out = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        if len(f) >= 2 and f[1].startswith("K"):
            out.add(f[1])
    return out


def completeness_vector(ko_set: set[str],
                        defs: dict[str, tuple[str, str, str]]) -> dict[str, float]:
    return {mid: module_completeness(d[0], ko_set, defs)
            for mid, d in defs.items()}


def ko_set_scores(est_kos: set[str], truth_kos: set[str]) -> dict:
    """Precision/recall/F1 of the *deduplicated KO set*, not of (gene, KO) pairs.

    This is the level module completeness actually consumes, and it behaves
    differently from the gene level in a way that turns out to explain the
    whole downstream result. A false-positive KO call on one gene frequently
    names a KO that is genuinely present on some *other* gene in the same
    genome, so dropping it does not remove the KO from the set -- but dropping
    a true positive that was a KO's only representative does. Arbitration's
    precision gain is therefore roughly halved at this level while its recall
    loss carries over intact.
    """
    tp = len(est_kos & truth_kos)
    precision = tp / len(est_kos) if est_kos else float("nan")
    recall = tp / len(truth_kos) if truth_kos else float("nan")
    f1 = (2 * tp / (len(est_kos) + len(truth_kos))
          if (est_kos or truth_kos) else float("nan"))
    return {"koset_precision": precision, "koset_recall": recall,
            "koset_f1": f1}


def compare(est: dict[str, float], truth: dict[str, float]) -> dict:
    mids = [m for m in truth
            if not np.isnan(truth[m]) and not np.isnan(est.get(m, np.nan))]
    if not mids:
        return {}
    e = np.array([est[m] for m in mids])
    t = np.array([truth[m] for m in mids])
    broken = int(np.sum((t >= COMPLETE) & (e < COMPLETE)))
    invented = int(np.sum((t < COMPLETE) & (e >= COMPLETE)))
    return {
        "n_modules_scored": len(mids),
        "mae": float(np.mean(np.abs(e - t))),
        "bias": float(np.mean(e - t)),
        "rmse": float(np.sqrt(np.mean((e - t) ** 2))),
        "n_exact": int(np.sum(np.isclose(e, t))),
        "modules_complete_in_truth": int(np.sum(t >= COMPLETE)),
        "modules_broken": broken,
        "modules_invented": invented,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir", required=True,
                   help="holds mpph_<org>.tsv (all) and <org>_best.tsv (best)")
    p.add_argument("--ground-truth-dir", required=True)
    p.add_argument("--organisms", nargs="+", required=True)
    p.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    p.add_argument("--all-modules", action="store_true",
                   help="include signature/reaction modules, not only Pathway")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    session = make_session()
    cache = Path(args.cache_dir)
    mods = list_modules(session, cache)
    print(f"{len(mods)} module(s) listed; fetching definitions ...", flush=True)
    defs = fetch_module_definitions(session, sorted(mods), cache)
    if not args.all_modules:
        defs = {m: d for m, d in defs.items() if d[2] == "Pathway"}
    print(f"{len(defs)} module definition(s) in scope")

    results, gt = Path(args.results_dir), Path(args.ground_truth_dir)
    rows = []
    for org in args.organisms:
        id_map_path = gt / f"protein_id_map_{org}.tsv"
        id_map = (load_protein_id_map(id_map_path, org)
                  if id_map_path.exists() else None)
        truth_kos = set()
        for kos in load_ground_truth(gt / f"link_ko_{org}.tsv", org, id_map).values():
            truth_kos |= kos
        truth_vec = completeness_vector(truth_kos, defs)

        for label, path in (("all", results / f"mpph_{org}.tsv"),
                            ("best", results / f"{org}_best.tsv")):
            if not path.exists():
                print(f"  (skipping {org}/{label}: {path} missing)")
                continue
            est_kos = kos_from_mapper(path)
            est = completeness_vector(est_kos, defs)
            row = compare(est, truth_vec)
            if row:
                rows.append({"organism": org, "policy": label,
                             "n_kos": len(est_kos),
                             "n_kos_truth": len(truth_kos),
                             **ko_set_scores(est_kos, truth_kos), **row})

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(out.to_string(index=False))
        if not out.empty:
            print("\n=== mean by policy ===")
            print(out.groupby("policy")[
                ["koset_precision", "koset_recall", "koset_f1",
                 "mae", "bias", "rmse", "modules_broken", "modules_invented"]
            ].mean().round(4).to_string())
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
