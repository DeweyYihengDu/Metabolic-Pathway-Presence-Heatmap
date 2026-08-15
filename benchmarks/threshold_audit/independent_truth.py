"""Is the margin-precision curve a property of the method, or of KEGG?

The whole audit is scored against KEGG's `link/ko/<org>`. That reference is
itself similarity-derived and curated, so a sceptical reading of the result is
available: perhaps precision rises with margin because high-margin hits are
the ones KEGG's own similarity pipeline also finds, and the curve measures
agreement between two views of the same evidence rather than reliability.

That reading makes a testable prediction. Scored against a *different*
reference the curve should flatten -- if it does not, the relationship is a
property of the HMM score itself.

eggNOG-mapper is the available independent reference: it is
similarity-based like KEGG's pipeline but built from a different orthology
resource, with different clade-specific models, and it never sees KOfam's
profiles or thresholds. It is a *noisier* reference than KEGG, not a better
one, so absolute precision against it will be lower everywhere. Only the
**shape** of the curve is being compared, and only the shape is interpreted.

Because both references are similarity-flavoured, agreement between them
cannot rule out every shared bias. What this does rule out is the specific
and most likely alternative -- that the curve is an artefact of KEGG's
particular curation.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "annotation"))
sys.path.insert(0, str(Path(__file__).parent))
from compare_annotation import (
    load_emapper_annotations,
    load_ground_truth,
    load_protein_id_map,
)
from summarize_audit import BANDS


def curve_against(df: pd.DataFrame, truth_pairs: set, label: str) -> pd.DataFrame:
    acc = df[df["delta"] >= 0].copy()
    acc["is_true"] = [(g, k) in truth_pairs
                      for g, k in zip(acc["gene"], acc["ko"])]
    rows = []
    for lo, hi in BANDS:
        sel = acc[(acc["delta"] >= lo) & (acc["delta"] < hi)]
        if not len(sel):
            continue
        rows.append({"reference": label, "margin_lo": lo, "margin_hi": hi,
                     "n_calls": len(sel), "precision": sel["is_true"].mean()})
    return pd.DataFrame(rows)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit-dir", required=True)
    p.add_argument("--results-dir", required=True, help="holds emapper outputs")
    p.add_argument("--ground-truth-dir", required=True)
    p.add_argument("--organisms", nargs="+", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    audit = Path(args.audit_dir)
    results = Path(args.results_dir)
    gt = Path(args.ground_truth_dir)

    frames = []
    for org in args.organisms:
        hits = audit / f"subthreshold_{org}.tsv"
        emap = results / f"emapper_{org}.emapper.annotations"
        if not (hits.exists() and emap.exists()):
            print(f"  (skipping {org})")
            continue
        df = pd.read_csv(hits, sep="\t")
        df["gene"] = df["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)

        id_map_path = gt / f"protein_id_map_{org}.tsv"
        id_map = (load_protein_id_map(id_map_path, org)
                  if id_map_path.exists() else None)
        if id_map is not None:
            df = df[df["gene"].isin(set(id_map.values()))]

        kegg = load_ground_truth(gt / f"link_ko_{org}.tsv", org, id_map)
        kegg_pairs = {(g, k) for g, kos in kegg.items() for k in kos}
        egg = load_emapper_annotations(emap)
        egg_pairs = {(g, k) for g, kos in egg.items() for k in kos}

        for label, pairs in (("KEGG", kegg_pairs), ("eggNOG", egg_pairs)):
            c = curve_against(df, pairs, label)
            c.insert(0, "organism", org)
            frames.append(c)

    out = pd.concat(frames, ignore_index=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)

    pivot = out.pivot_table(index=["margin_lo", "margin_hi"],
                            columns="reference", values="precision")
    print(pivot.round(3).to_string())

    # The claim is about shape, so compare the *span* each reference shows
    # and the rank correlation with margin, not the absolute levels.
    print()
    for ref in ("KEGG", "eggNOG"):
        sub = out[out["reference"] == ref]
        by_band = sub.groupby("margin_lo")["precision"].mean()
        span = by_band.max() - by_band.min()
        rho = np.corrcoef(by_band.index.values, by_band.to_numpy())[0, 1]
        print(f"{ref:7s}: precision {by_band.min():.3f} -> {by_band.max():.3f} "
              f"(span {span:.3f}), correlation with margin {rho:+.3f}")
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
