"""Is the recall arbitration gives up biased toward particular pathways?

Arbitration costs ~1.3 points of per-gene recall. Whether that matters depends
entirely on *which* true calls are lost. Losing them uniformly at random is a
tolerable price; losing them systematically from one part of metabolism would
skew every downstream metabolic profile in the same direction, which would be
worse than the precision gain is good.

Test: take the true (gene, KO) pairs that `--multi-ko-policy all` finds and
`best` drops, and ask whether their KOs are enriched in particular BRITE
top-level categories or KEGG pathways relative to the true calls that survive.
A hypergeometric test per pathway with BH correction, pooled across genomes so
there is enough signal to see.

Reports the top enriched pathways and, crucially, whether *anything* survives
correction. "No pathway is significantly enriched" is a real and reassuring
result and is reported as such rather than buried.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "annotation"))

from compare_annotation import load_ground_truth, load_protein_id_map

from mpph.analysis import benjamini_hochberg
from mpph.enrichment import (
    fetch_ko_pathway_membership,
    fetch_pathway_names,
    invert_membership,
)
from mpph.kegg import DEFAULT_CACHE, make_session


def pairs_from_mapper(path: Path) -> set[tuple[str, str]]:
    out = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        if len(f) >= 2 and f[1].startswith("K"):
            gene = f[0].rsplit(".", 1)[0] if f[0].rsplit(".", 1)[-1].isdigit() else f[0]
            out.add((gene, f[1]))
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--ground-truth-dir", required=True)
    p.add_argument("--organisms", nargs="+", required=True)
    p.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    p.add_argument("--out", required=True)
    args = p.parse_args()

    results, gt = Path(args.results_dir), Path(args.ground_truth_dir)
    lost: list[str] = []
    kept: list[str] = []
    per_org = []
    for org in args.organisms:
        id_map_path = gt / f"protein_id_map_{org}.tsv"
        id_map = (load_protein_id_map(id_map_path, org)
                  if id_map_path.exists() else None)
        truth = load_ground_truth(gt / f"link_ko_{org}.tsv", org, id_map)
        true_pairs = {(g, k) for g, kos in truth.items() for k in kos}
        all_p = pairs_from_mapper(results / f"mpph_{org}.tsv") & true_pairs
        best_p = pairs_from_mapper(results / f"{org}_best.tsv") & true_pairs
        dropped = all_p - best_p
        lost.extend(k for _, k in dropped)
        kept.extend(k for _, k in best_p)
        per_org.append({"organism": org, "true_calls_all": len(all_p),
                        "true_calls_best": len(best_p),
                        "true_calls_lost": len(dropped),
                        "frac_lost": len(dropped) / len(all_p) if all_p else float("nan")})

    print(pd.DataFrame(per_org).to_string(index=False))
    print(f"\npooled: {len(lost)} true call(s) lost, {len(kept)} retained")
    if not lost:
        print("nothing lost -- no bias to test")
        return 0

    session = make_session()
    cache = Path(args.cache_dir)
    ko_to_pathways = fetch_ko_pathway_membership(session, cache)
    pathway_to_kos = invert_membership(ko_to_pathways)
    names = fetch_pathway_names(session, cache)

    lost_set, kept_set = set(lost), set(kept)
    universe = lost_set | kept_set
    rows = []
    for pw, kos in pathway_to_kos.items():
        in_pw = kos & universe
        if len(in_pw) < 5:
            continue                      # too small to say anything
        k = len(in_pw & lost_set)
        # P(X >= k) for X ~ Hypergeom(N=universe, K=pathway members, n=lost)
        pval = stats.hypergeom.sf(k - 1, len(universe), len(in_pw),
                                  len(lost_set))
        rows.append({"pathway": pw, "name": names.get(pw, ""),
                     "n_in_pathway": len(in_pw), "n_lost": k,
                     "frac_lost": k / len(in_pw), "p_value": pval})

    df = pd.DataFrame(rows)
    if df.empty:
        print("no pathway large enough to test")
        return 0
    df["q_value"] = benjamini_hochberg(df["p_value"].to_numpy())
    df = df.sort_values("p_value")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    n_sig = int((df["q_value"] < 0.05).sum())
    baseline = len(lost_set) / len(universe)
    print(f"\nbaseline loss rate across all KOs: {baseline:.4f}")
    print(f"{len(df)} pathway(s) tested, {n_sig} with q < 0.05")
    if n_sig == 0:
        print("=> the recall loss is NOT concentrated in any pathway; it looks "
              "uniform across metabolism, which is the reassuring outcome")
    print("\ntop 10 by p-value:")
    print(df.head(10).to_string(index=False))
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
