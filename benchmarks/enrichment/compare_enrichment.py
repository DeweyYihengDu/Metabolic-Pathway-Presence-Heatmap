"""Head-to-head comparison: mpph gsea/enrich vs. clusterProfiler GSEA()/enricher().

Both tools are run against the *identical* TERM2GENE category definitions
(see export_term2gene.py) and the identical ranked list / study set, so this
measures agreement in the statistical implementation, not database
differences. Not part of the installed `mpph` package.

Usage:
    python compare_enrichment.py \
        --mpph-gsea results/mpph_gsea.csv \
        --clusterprofiler-gsea results/clusterprofiler_gsea.csv \
        --out report/enrichment_comparison.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy import stats


def load_mpph_gsea(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"category_id": str})
    return df[["category_id", "ES", "NES", "p_value", "q_value"]].rename(
        columns={"ES": "ES_mpph", "NES": "NES_mpph",
                 "p_value": "p_mpph", "q_value": "q_mpph"})


def load_clusterprofiler_gsea(path: str | Path) -> pd.DataFrame:
    """clusterProfiler's `as.data.frame(GSEA(...))` export: ID, enrichmentScore,
    NES, pvalue, qvalue (or p.adjust, depending on version -- both handled)."""
    df = pd.read_csv(path, dtype={"ID": str})
    q_col = "qvalue" if "qvalue" in df.columns else "p.adjust"
    return df[["ID", "enrichmentScore", "NES", "pvalue", q_col]].rename(
        columns={"ID": "category_id", "enrichmentScore": "ES_cp", "NES": "NES_cp",
                 "pvalue": "p_cp", q_col: "q_cp"})


def load_mpph_ora(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"category_id": str})
    return df[["category_id", "p_value", "q_value"]].rename(
        columns={"p_value": "p_mpph", "q_value": "q_mpph"})


def load_clusterprofiler_ora(path: str | Path) -> pd.DataFrame:
    """clusterProfiler's `as.data.frame(enricher(...))` export: ID, pvalue,
    qvalue/p.adjust. `GeneRatio`/`BgRatio` are string fractions ("45/606"),
    not directly comparable to mpph's numeric `fold_enrichment`/`odds_ratio`
    without parsing -- out of scope here; p-value/significance agreement
    (the part every ORA implementation shares regardless of which effect-size
    statistic it additionally reports) is what's compared."""
    df = pd.read_csv(path, dtype={"ID": str})
    q_col = "qvalue" if "qvalue" in df.columns else "p.adjust"
    return df[["ID", "pvalue", q_col]].rename(
        columns={"ID": "category_id", "pvalue": "p_cp", q_col: "q_cp"})


def _compare_categories(mpph_df: pd.DataFrame, cp_df: pd.DataFrame,
                        score_cols: list[tuple[str, str, str]],
                        alpha: float = 0.05) -> dict:
    """Shared merge/correlate/significant-set-agreement logic for both GSEA
    and ORA: join on category_id, Spearman-correlate whichever numeric
    columns both sides have, compare which categories each calls significant
    (q<alpha). A category one tool reports and the other doesn't (n_only_*)
    is not necessarily a bug -- e.g. clusterProfiler's enricher() omits
    zero-hit categories from its output entirely, while mpph enrich reports
    them (correctly, as non-significant) -- see this module's docstring/the
    accompanying README for that specific, confirmed, intentional difference.
    """
    merged = mpph_df.merge(cp_df, on="category_id", how="outer", indicator=True)
    both = merged[merged["_merge"] == "both"]

    correlations = {}
    for metric, col_a, col_b in score_cols:
        rho, pval = stats.spearmanr(both[col_a], both[col_b])
        correlations[f"spearman_{metric}"] = rho
        correlations[f"spearman_{metric}_pvalue"] = pval

    sig_mpph = set(both.loc[both["q_mpph"] < alpha, "category_id"])
    sig_cp = set(both.loc[both["q_cp"] < alpha, "category_id"])
    union = sig_mpph | sig_cp
    jaccard = len(sig_mpph & sig_cp) / len(union) if union else float("nan")

    return {
        "n_tested_mpph": len(mpph_df), "n_tested_clusterprofiler": len(cp_df),
        "n_tested_both": len(both),
        "n_only_mpph": int((merged["_merge"] == "left_only").sum()),
        "n_only_clusterprofiler": int((merged["_merge"] == "right_only").sum()),
        "n_significant_mpph": len(sig_mpph),
        "n_significant_clusterprofiler": len(sig_cp),
        "n_significant_both": len(sig_mpph & sig_cp),
        "significant_set_jaccard": jaccard,
        **correlations,
    }


def compare_gsea(mpph_df: pd.DataFrame, cp_df: pd.DataFrame, alpha: float = 0.05) -> dict:
    return _compare_categories(mpph_df, cp_df, [
        ("ES", "ES_mpph", "ES_cp"), ("NES", "NES_mpph", "NES_cp"),
        ("p_value", "p_mpph", "p_cp"),
    ], alpha=alpha)


def compare_ora(mpph_df: pd.DataFrame, cp_df: pd.DataFrame, alpha: float = 0.05) -> dict:
    return _compare_categories(mpph_df, cp_df, [("p_value", "p_mpph", "p_cp")], alpha=alpha)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=["gsea", "ora"], default="gsea")
    p.add_argument("--mpph", dest="mpph_path", required=True,
                   help="mpph gsea's *_gsea.csv, or mpph enrich's *_enrichment.csv "
                        "for --mode ora.")
    p.add_argument("--clusterprofiler", dest="cp_path", required=True,
                   help="as.data.frame(GSEA(...)) or as.data.frame(enricher(...)) export.")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    if args.mode == "gsea":
        mpph_df = load_mpph_gsea(args.mpph_path)
        cp_df = load_clusterprofiler_gsea(args.cp_path)
        result = compare_gsea(mpph_df, cp_df, alpha=args.alpha)
    else:
        mpph_df = load_mpph_ora(args.mpph_path)
        cp_df = load_clusterprofiler_ora(args.cp_path)
        result = compare_ora(mpph_df, cp_df, alpha=args.alpha)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([result]).to_csv(out, index=False)
    for k, v in result.items():
        print(f"{k}: {v}")
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
