"""Does calibrating on a *closer* genome help? The minimum viable test.

The transfer experiment showed a confidence curve fitted on two bacteria is
not calibrated on an archaeon: the ranking transfers, the absolute
probabilities do not. The proposed fix was to recalibrate each query against
its phylogenetic neighbours rather than a global average.

**Scoping honesty first.** That idea cannot be properly tested on this
benchmark. Seven genomes spanning three domains have no meaningful notion of
"neighbour" -- they are all maximally distant from one another. A real test
needs on the order of hundreds of genomes with KEGG ground truth, sampled to
span a *range* of pairwise distances, which is a substantial compute job and
is not what this script does.

What *is* testable now is the premise the idea rests on: **if relatedness
matters at all, a curve fitted on a same-group genome should calibrate a query
better than one fitted on a distant genome.** If that ordering does not appear
even between "same phylum-ish" and "different domain", the neighbour idea has
no foundation and should be dropped rather than scaled up.

Every ordered pair (train -> test) is evaluated, and each pair is labelled by
a coarse relatedness tier from the benchmark's own taxonomy. No tree is
required for that, which is the point -- it keeps this a cheap falsification
test rather than a small version of the expensive experiment.
"""
from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from calibration_transfer import accepted, brier, ece, fit_logistic, predict
from summarize_audit import load_scored

# Coarse relatedness, from the benchmark's own composition. Deliberately not a
# tree: with seven maximally-separated genomes a tree would add precision the
# design cannot use.
GROUP = {
    "eco": ("Bacteria", "Proteobacteria"),
    "bsu": ("Bacteria", "Firmicutes"),
    "vbs": ("Bacteria", "PVC"),
    "lbac": ("Bacteria", "PVC"),
    "mja": ("Archaea", "Euryarchaeota"),
    "sce": ("Eukaryota", "Fungi"),
    "ath": ("Eukaryota", "Viridiplantae"),
}


def tier(a: str, b: str) -> str:
    da, pa = GROUP[a]
    db, pb = GROUP[b]
    if pa == pb:
        return "1_same_phylum_group"
    if da == db:
        return "2_same_domain"
    return "3_different_domain"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit-dir", required=True)
    p.add_argument("--ground-truth-dir", required=True)
    p.add_argument("--organisms", nargs="+", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    audit, gt = Path(args.audit_dir), Path(args.ground_truth_dir)
    data = {}
    for org in args.organisms:
        id_map = gt / f"protein_id_map_{org}.tsv"
        data[org] = accepted(load_scored(
            audit / f"subthreshold_{org}.tsv", org,
            gt / f"link_ko_{org}.tsv", id_map if id_map.exists() else None))

    rows = []
    for train, test in itertools.permutations(args.organisms, 2):
        d_tr, y_tr = data[train]
        d_te, y_te = data[test]
        a, b = fit_logistic(d_tr, y_tr)
        p_model = predict(d_te, a, b)
        p_const = np.full_like(y_te, float(y_tr.mean()))
        rows.append({
            "train": train, "test": test, "tier": tier(train, test),
            "ece_model": ece(p_model, y_te),
            "ece_constant": ece(p_const, y_te),
            "brier_model": brier(p_model, y_te),
            "brier_constant": brier(p_const, y_te),
        })

    df = pd.DataFrame(rows)
    df["ece_improvement"] = df["ece_constant"] - df["ece_model"]
    df["brier_improvement"] = df["brier_constant"] - df["brier_model"]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(df.round(4).to_string(index=False))
        print("\n=== mean by relatedness tier ===")
        summary = df.groupby("tier")[
            ["ece_model", "ece_improvement", "brier_model", "brier_improvement"]
        ].agg(["mean", "count"]).round(4)
        print(summary.to_string())
        print("\nThe premise predicts ece_model rises monotonically from tier 1 "
              "to tier 3.\nIf it does not, closeness is not carrying calibration "
              "information here and\nthe neighbour-recalibration idea has no "
              "foundation to scale up.")
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
