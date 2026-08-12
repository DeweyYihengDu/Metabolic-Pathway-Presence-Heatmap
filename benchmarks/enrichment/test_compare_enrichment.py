"""Verifies compare_enrichment.py's merge/correlation/agreement logic against
hand-constructed synthetic GSEA outputs with known expected correlation and
significant-set overlap, before trusting it on real mpph/clusterProfiler
output. Not part of the installed mpph package -- run explicitly via
`pytest benchmarks/`, not picked up by the main `pytest -q`.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from compare_enrichment import (
    compare_gsea,
    compare_ora,
    load_clusterprofiler_gsea,
    load_clusterprofiler_ora,
    load_mpph_gsea,
    load_mpph_ora,
)


def test_load_mpph_gsea_selects_and_renames_columns(tmp_path):
    p = tmp_path / "mpph.csv"
    pd.DataFrame({
        "category_id": ["00010"], "category_name": ["Glycolysis"],
        "size": [50], "ES": [0.9], "NES": [2.1], "p_value": [0.001],
        "leading_edge_size": [10], "leading_edge_genes": ["K1,K2"],
        "q_value": [0.004],
    }).to_csv(p, index=False)
    result = load_mpph_gsea(p)
    assert list(result.columns) == ["category_id", "ES_mpph", "NES_mpph", "p_mpph", "q_mpph"]
    assert result.iloc[0]["ES_mpph"] == 0.9


def test_load_clusterprofiler_gsea_handles_qvalue_or_p_adjust(tmp_path):
    p = tmp_path / "cp.csv"
    pd.DataFrame({
        "ID": ["00010"], "enrichmentScore": [0.88], "NES": [2.0],
        "pvalue": [0.0015], "qvalue": [0.005],
    }).to_csv(p, index=False)
    result = load_clusterprofiler_gsea(p)
    assert list(result.columns) == ["category_id", "ES_cp", "NES_cp", "p_cp", "q_cp"]

    p2 = tmp_path / "cp_no_qvalue.csv"
    pd.DataFrame({
        "ID": ["00010"], "enrichmentScore": [0.88], "NES": [2.0],
        "pvalue": [0.0015], "p.adjust": [0.005],
    }).to_csv(p2, index=False)
    result2 = load_clusterprofiler_gsea(p2)
    assert result2.iloc[0]["q_cp"] == 0.005


def test_compare_gsea_perfect_rank_agreement_gives_correlation_one():
    # 00010/00020/00030 present in both, in perfectly consistent rank order;
    # 00040 only in mpph, 00050 only in clusterprofiler.
    mpph_df = pd.DataFrame({
        "category_id": ["00010", "00020", "00030", "00040"],
        "ES_mpph": [0.9, -0.8, 0.1, 0.05], "NES_mpph": [2.1, -1.9, 0.3, 0.1],
        "p_mpph": [0.001, 0.002, 0.4, 0.9], "q_mpph": [0.004, 0.004, 0.5, 0.9],
    })
    cp_df = pd.DataFrame({
        "category_id": ["00010", "00020", "00030", "00050"],
        "ES_cp": [0.88, -0.79, 0.12, 0.5], "NES_cp": [2.0, -1.85, 0.25, 1.1],
        "p_cp": [0.0015, 0.0025, 0.35, 0.03], "q_cp": [0.005, 0.005, 0.45, 0.06],
    })
    result = compare_gsea(mpph_df, cp_df, alpha=0.05)

    assert result["n_tested_both"] == 3
    assert result["n_only_mpph"] == 1
    assert result["n_only_clusterprofiler"] == 1
    # 00010 and 00020 pass q<0.05 in both; 00030 in neither; 00050's q=0.06 fails.
    assert result["n_significant_mpph"] == 2
    assert result["n_significant_clusterprofiler"] == 2
    assert result["significant_set_jaccard"] == pytest.approx(1.0)
    assert result["spearman_ES"] == pytest.approx(1.0)
    assert result["spearman_NES"] == pytest.approx(1.0)
    assert result["spearman_p_value"] == pytest.approx(1.0)


def test_compare_gsea_disagreeing_significance_lowers_jaccard():
    mpph_df = pd.DataFrame({
        "category_id": ["00010", "00020"],
        "ES_mpph": [0.9, 0.1], "NES_mpph": [2.0, 0.2],
        "p_mpph": [0.001, 0.5], "q_mpph": [0.004, 0.6],
    })
    cp_df = pd.DataFrame({
        "category_id": ["00010", "00020"],
        "ES_cp": [0.88, 0.15], "NES_cp": [1.9, 0.3],
        "p_cp": [0.002, 0.01], "q_cp": [0.01, 0.02],  # 00020 significant here, not in mpph
    })
    result = compare_gsea(mpph_df, cp_df, alpha=0.05)
    assert result["n_significant_mpph"] == 1
    assert result["n_significant_clusterprofiler"] == 2
    assert result["significant_set_jaccard"] == pytest.approx(0.5)  # 1 shared / 2 union


def test_load_mpph_ora_selects_and_renames_columns(tmp_path):
    p = tmp_path / "mpph_ora.csv"
    pd.DataFrame({
        "category_id": ["00010"], "category_name": ["Glycolysis"],
        "k_study_hits": [10], "n_study_total": [100], "K_background_hits": [50],
        "N_background_total": [1000], "fold_enrichment": [2.0], "odds_ratio": [2.2],
        "p_value": [0.001], "q_value": [0.004],
    }).to_csv(p, index=False)
    result = load_mpph_ora(p)
    assert list(result.columns) == ["category_id", "p_mpph", "q_mpph"]


def test_load_clusterprofiler_ora_handles_qvalue_or_p_adjust(tmp_path):
    p = tmp_path / "cp_ora.csv"
    pd.DataFrame({
        "ID": ["00010"], "Description": ["Glycolysis"], "GeneRatio": ["10/100"],
        "BgRatio": ["50/1000"], "pvalue": [0.001], "qvalue": [0.004],
    }).to_csv(p, index=False)
    result = load_clusterprofiler_ora(p)
    assert list(result.columns) == ["category_id", "p_cp", "q_cp"]


def test_compare_ora_only_correlates_p_value_not_es_nes():
    # ORA has no ES/NES concept -- compare_ora must not require those columns.
    mpph_df = pd.DataFrame({
        "category_id": ["00010", "00020", "00030"],
        "p_mpph": [0.001, 0.5, 0.9], "q_mpph": [0.003, 0.6, 0.9],
    })
    cp_df = pd.DataFrame({
        "category_id": ["00010", "00020", "00030"],
        "p_cp": [0.002, 0.4, 0.8], "q_cp": [0.006, 0.5, 0.8],
    })
    result = compare_ora(mpph_df, cp_df, alpha=0.05)
    assert result["n_tested_both"] == 3
    assert result["n_significant_mpph"] == 1
    assert result["n_significant_clusterprofiler"] == 1
    assert "spearman_p_value" in result
    assert "spearman_ES" not in result


def test_compare_ora_reports_categories_only_one_tool_lists():
    # Real, confirmed, intentional difference: clusterProfiler's enricher()
    # omits zero-hit categories from its output; mpph enrich reports them
    # (as non-significant). n_only_mpph should surface this, not hide it.
    mpph_df = pd.DataFrame({
        "category_id": ["00010", "00020"],  # 00020 has zero study-set hits
        "p_mpph": [0.001, 1.0], "q_mpph": [0.003, 1.0],
    })
    cp_df = pd.DataFrame({
        "category_id": ["00010"],  # 00020 never appears here
        "p_cp": [0.002], "q_cp": [0.006],
    })
    result = compare_ora(mpph_df, cp_df, alpha=0.05)
    assert result["n_only_mpph"] == 1
    assert result["n_only_clusterprofiler"] == 0
    assert result["n_tested_both"] == 1
