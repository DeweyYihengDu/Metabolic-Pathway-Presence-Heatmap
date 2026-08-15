"""Arbitration between competing KO calls on the same gene.

KOfam sets each KO's threshold independently, by maximising that KO's own
F-measure in isolation. Nothing in the procedure makes KOs compete. So when a
protein matches several related profiles -- paralogous subfamilies, different
specificities of one enzyme family -- every one of them can clear its own
threshold, and both KofamScan and `mpph annotate` emit all of them.

KEGG's own reference does not: across all seven benchmark genomes and all
three domains, >=99.88% of genes carry exactly one KO (max 2). Multi-KO calls
are therefore over-calls almost by construction, and on *E. coli* they account
for 13.5% of calls but 69% of all false positives.

This evaluates candidate arbitration rules against KEGG's assignments. Scores
are not comparable across KOs -- thresholds span roughly 30 to >2000 bits -- so
rules rank by a *normalised* quantity rather than by raw score:

* `all`          keep everything above threshold (what the field does today)
* `top_delta`    per gene, keep the largest absolute margin (score - threshold)
* `top_relative` per gene, keep the largest relative margin (score / threshold)
* `top_delta_tie` as `top_delta`, but keep any call within `--tie-window` bits
                  of the winner, since a near-tie is not evidence for either
* `gap_gated`    keep the winner only when it beats the runner-up by
                  `--min-gap` bits; otherwise keep nothing for that gene,
                  trading recall for precision where the evidence is ambiguous

Reported per rule: precision, recall and F1 over (gene, KO) pairs, so a rule
that buys precision by discarding true positives is visibly penalised rather
than flattered.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "annotation"))
from compare_annotation import load_ground_truth, load_protein_id_map

RULES = ("all", "top_delta", "top_relative", "top_delta_tie", "gap_gated")


def apply_rule(acc: pd.DataFrame, rule: str, *, tie_window: float,
               min_gap: float) -> pd.DataFrame:
    if rule == "all":
        return acc
    if rule == "top_relative":
        key = "relative"
    else:
        key = "delta"

    if rule in ("top_delta", "top_relative"):
        return acc.loc[acc.groupby("gene")[key].idxmax()]

    if rule == "top_delta_tie":
        # transform() on `acc` itself, so the result carries acc's own index.
        # Computing it on a sorted copy and comparing back against `acc` is an
        # index-alignment bug, not a cosmetic one.
        best = acc.groupby("gene")[key].transform("max")
        return acc[acc[key] >= best - tie_window]

    if rule == "gap_gated":
        # Work entirely within the sorted frame so winner and runner-up are
        # always drawn from the same row ordering.
        ordered = acc.sort_values(["gene", key], ascending=[True, False])
        best = ordered.groupby("gene")[key].transform("first")
        second = ordered.groupby("gene")[key].transform(
            lambda s: s.iloc[1] if len(s) > 1 else float("-inf"))
        # A gene whose top two are within min_gap contributes nothing: the
        # evidence does not distinguish them, and guessing costs precision.
        winners = ordered[(ordered[key] >= best) & (ordered[key] - second >= min_gap)]
        return winners.drop_duplicates(subset="gene", keep="first")

    raise ValueError(f"unknown rule {rule!r}")


def evaluate(kept: pd.DataFrame, n_true_pairs: int) -> dict:
    tp = int(kept["is_true"].sum())
    fp = len(kept) - tp
    fn = n_true_pairs - tp
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    f1 = (2 * precision * recall / (precision + recall)
          if precision and recall and (precision + recall) else float("nan"))
    return {"n_calls": len(kept), "tp": tp, "fp": fp, "fn": fn,
            "precision": precision, "recall": recall, "f1": f1}


def load_accepted(hits: Path, org: str, truth_path: Path,
                  id_map_path: Path | None) -> tuple[pd.DataFrame, int]:
    df = pd.read_csv(hits, sep="\t")
    df["gene"] = df["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    id_map = load_protein_id_map(id_map_path, org) if id_map_path else None
    truth = load_ground_truth(truth_path, org, id_map)
    if id_map is not None:
        df = df[df["gene"].isin(set(id_map.values()))]
    acc = df[df["delta"] >= 0].copy()
    true_pairs = {(g, k) for g, kos in truth.items() for k in kos}
    acc["is_true"] = [(g, k) in true_pairs for g, k in zip(acc["gene"], acc["ko"])]
    # score/threshold: thresholds span ~30 to >2000 bits, so an absolute
    # margin of 20 means very different things for different KOs.
    acc["relative"] = acc["score"] / acc["threshold"]
    return acc, len(true_pairs)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit-dir", required=True)
    p.add_argument("--ground-truth-dir", required=True)
    p.add_argument("--organisms", nargs="+", required=True)
    p.add_argument("--tie-window", type=float, default=5.0)
    p.add_argument("--min-gap", type=float, default=10.0)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    audit, gt = Path(args.audit_dir), Path(args.ground_truth_dir)
    rows = []
    for org in args.organisms:
        hits = audit / f"subthreshold_{org}.tsv"
        if not hits.exists():
            print(f"  (skipping {org})")
            continue
        id_map = gt / f"protein_id_map_{org}.tsv"
        acc, n_true = load_accepted(hits, org, gt / f"link_ko_{org}.tsv",
                                    id_map if id_map.exists() else None)
        for rule in RULES:
            kept = apply_rule(acc, rule, tie_window=args.tie_window,
                              min_gap=args.min_gap)
            rows.append({"organism": org, "rule": rule, **evaluate(kept, n_true)})

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)

    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(out.to_string(index=False))
        print("\n=== mean across genomes ===")
        print(out.groupby("rule")[["precision", "recall", "f1"]].mean()
              .reindex(RULES).round(4).to_string())
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
