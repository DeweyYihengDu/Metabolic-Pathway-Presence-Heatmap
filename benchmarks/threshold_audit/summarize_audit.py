"""Cross-genome summary of the sub-threshold precision audit.

Answers two questions per genome, both against KEGG's own curated assignments:

1. **Is the discarded band recoverable?** Precision of hits below each KO's
   threshold, i.e. what the field throws away.
2. **Are the accepted calls homogeneous?** Precision of accepted calls as a
   function of margin above threshold. Every tool writes these as the same
   `1`; if precision varies strongly with margin, that binary is a lossy
   compression of a well-behaved curve, and the loss is measurable.

The headline statistic is `frac_fp_from_low_margin` -- the share of a genome's
false positives contributed by the small, *identifiable* set of calls sitting
near the threshold. If a few percent of calls produce a quarter of the errors,
the margin is actionable information rather than a curiosity.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "annotation"))
from compare_annotation import load_ground_truth, load_protein_id_map

# Margin bands for the accepted-call calibration curve, in bits above the KO's
# own adaptive threshold. Open-ended at the top: most calls sit far above.
BANDS = [(0, 5), (5, 10), (10, 20), (20, 40), (40, 60),
         (60, 100), (100, 200), (200, float("inf"))]
LOW_MARGIN = 20.0  # "near the threshold" for the headline statistic


def load_scored(hits_path: Path, org: str, truth_path: Path,
                id_map_path: Path | None) -> pd.DataFrame:
    df = pd.read_csv(hits_path, sep="\t")
    df["gene"] = df["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    id_map = load_protein_id_map(id_map_path, org) if id_map_path else None
    truth = load_ground_truth(truth_path, org, id_map)
    if id_map is not None:
        df = df[df["gene"].isin(set(id_map.values()))]
    true_pairs = {(g, k) for g, kos in truth.items() for k in kos}
    df = df.copy()
    df["is_true"] = [(g, k) in true_pairs for g, k in zip(df["gene"], df["ko"])]
    return df


def summarize(df: pd.DataFrame, org: str) -> tuple[dict, pd.DataFrame]:
    acc = df[df["delta"] >= 0]
    disc = df[df["delta"] < 0]
    low = acc[acc["delta"] < LOW_MARGIN]
    n_fp = int((~acc["is_true"]).sum())
    row = {
        "organism": org,
        "n_accepted": len(acc),
        "precision_accepted": acc["is_true"].mean() if len(acc) else float("nan"),
        "n_discarded_within_60bits": len(disc),
        "precision_discarded": disc["is_true"].mean() if len(disc) else float("nan"),
        "n_recoverable_within_5bits": int(
            disc[disc["delta"] >= -5]["is_true"].sum()),
        "frac_calls_low_margin": len(low) / len(acc) if len(acc) else float("nan"),
        "precision_low_margin": low["is_true"].mean() if len(low) else float("nan"),
        "frac_fp_from_low_margin": (
            int((~low["is_true"]).sum()) / n_fp if n_fp else float("nan")),
    }
    # Enrichment of errors in the low-margin band: >1 means the band carries
    # more than its share of the genome's mistakes.
    row["fp_enrichment_low_margin"] = (
        row["frac_fp_from_low_margin"] / row["frac_calls_low_margin"]
        if row["frac_calls_low_margin"] else float("nan"))

    curve = []
    for lo, hi in BANDS:
        sel = acc[(acc["delta"] >= lo) & (acc["delta"] < hi)]
        if not len(sel):
            continue
        curve.append({"organism": org, "margin_lo": lo, "margin_hi": hi,
                      "n_calls": len(sel),
                      "frac_of_calls": len(sel) / len(acc),
                      "precision": sel["is_true"].mean()})
    return row, pd.DataFrame(curve)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit-dir", required=True,
                   help="directory holding subthreshold_<org>.tsv")
    p.add_argument("--ground-truth-dir", required=True)
    p.add_argument("--organisms", nargs="+", required=True)
    p.add_argument("--out-summary", required=True)
    p.add_argument("--out-curve", required=True)
    args = p.parse_args()

    audit, gt = Path(args.audit_dir), Path(args.ground_truth_dir)
    rows, curves = [], []
    for org in args.organisms:
        hits = audit / f"subthreshold_{org}.tsv"
        if not hits.exists():
            print(f"  (skipping {org}: {hits} not found)")
            continue
        id_map = gt / f"protein_id_map_{org}.tsv"
        df = load_scored(hits, org, gt / f"link_ko_{org}.tsv",
                         id_map if id_map.exists() else None)
        row, curve = summarize(df, org)
        rows.append(row)
        curves.append(curve)

    summary = pd.DataFrame(rows)
    curve = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame()
    for path, frame in ((args.out_summary, summary), (args.out_curve, curve)):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False)

    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(summary.to_string(index=False))
        print()
        print(curve.to_string(index=False))
    print(f"\nWrote {args.out_summary} and {args.out_curve}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
