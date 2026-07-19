"""Tests for the rank-based (GSEA-style) enrichment engine."""
import numpy as np
import pandas as pd
import pytest

from mpph.gsea import (
    enrichment_score,
    gsea_analysis,
    leading_edge,
    load_expression_matrix,
    load_ranked_list,
    rank_from_expression,
)


def _ranking(n=100):
    genes = [f"g{i}" for i in range(n)]
    scores = np.array(list(range(n, 0, -1)), dtype=float)
    return genes, scores


def test_es_top_concentrated_is_strongly_positive():
    genes, scores = _ranking()
    es, running = enrichment_score(genes, scores, set(genes[:10]))
    assert es > 0.9
    assert running[9] == pytest.approx(es)  # peak right after the last hit


def test_es_bottom_concentrated_is_strongly_negative():
    genes, scores = _ranking()
    es, _ = enrichment_score(genes, scores, set(genes[-10:]))
    assert es < -0.9


def test_es_empty_gene_set_is_zero():
    genes, scores = _ranking()
    es, running = enrichment_score(genes, scores, set())
    assert es == 0.0
    assert not running.any()


def test_leading_edge_positive_es_is_prefix_hits():
    genes = ["a", "b", "c", "d", "e"]
    hits = np.array([True, False, True, False, False])
    running = np.array([0.5, 0.3, 0.9, 0.6, 0.3])  # peak at index 2
    edge = leading_edge(genes, hits, running, es=0.9)
    assert edge == ["a", "c"]  # hits at/after start, up to and incl. the peak


def test_leading_edge_negative_es_is_suffix_hits():
    genes = ["a", "b", "c", "d", "e"]
    hits = np.array([False, False, True, False, True])
    running = np.array([0.1, 0.0, -0.9, -0.5, -0.3])  # trough at index 2, not last
    edge = leading_edge(genes, hits, running, es=-0.9)
    assert edge == ["c", "e"]  # hits at/after the trough (indices 2..4)


def test_gsea_analysis_recovers_known_signal():
    genes, scores = _ranking()
    ranked = pd.Series(scores, index=genes)
    categories = {"top_set": set(genes[:10]), "random_set": set(genes[40:50])}
    results, running_sums, order = gsea_analysis(
        ranked, categories, min_size=5, max_size=50, permutations=300, seed=1)
    assert order == genes
    top = results.set_index("category_id").loc["top_set"]
    assert top["ES"] > 0.9
    assert top["q_value"] < 0.05
    assert "top_set" in running_sums
    assert len(running_sums["top_set"]) == len(genes)


def test_null_es_batch_size_does_not_change_results():
    # The null distribution is now computed in chunks to bound peak memory
    # for a large gene universe -- this must not change which random draws
    # are consumed or the resulting statistics, regardless of chunk size.
    from mpph.gsea import _null_es_batch
    abs_scores = np.abs(np.random.default_rng(1).normal(size=200))
    n, size, weight, permutations = 200, 15, 1.0, 137  # not divisible by 17
    r1 = _null_es_batch(abs_scores, n, size, weight, permutations,
                        np.random.default_rng(42), batch_size=200)
    r2 = _null_es_batch(abs_scores, n, size, weight, permutations,
                        np.random.default_rng(42), batch_size=17)
    assert np.array_equal(r1, r2)


def test_gsea_size_filter_excludes_out_of_range_categories():
    genes, scores = _ranking()
    ranked = pd.Series(scores, index=genes)
    categories = {"tiny": set(genes[:2]), "right_size": set(genes[:20])}
    results, _, _ = gsea_analysis(ranked, categories, min_size=10, max_size=30,
                                  permutations=100, seed=0)
    assert "tiny" not in set(results["category_id"])
    assert "right_size" in set(results["category_id"])


def test_gsea_unannotated_genes_still_count_as_misses():
    # A gene set that covers HALF the ranked list but only among annotated
    # genes should not silently shrink N -- unlike ORA, GSEA keeps every
    # ranked gene (annotated or not) in the running-sum denominator.
    genes, scores = _ranking(20)
    gene_set = set(genes[:5])  # only 5 of 20 genes are "annotated" anywhere
    es, running = enrichment_score(genes, scores, gene_set)
    assert len(running) == 20  # not reduced to 5


def test_rank_from_expression_signal2noise():
    expr = pd.DataFrame({
        "s1": [10.0, 1.0], "s2": [12.0, 1.0],
        "s3": [1.0, 10.0], "s4": [1.0, 12.0],
    }, index=["up_gene", "down_gene"])
    groups = pd.Series({"s1": "A", "s2": "A", "s3": "B", "s4": "B"})
    ranked = rank_from_expression(expr, groups, "A", "B", metric="signal2noise")
    assert ranked.index[0] == "up_gene"     # higher in A -> ranked first
    assert ranked.index[-1] == "down_gene"  # higher in B -> ranked last
    assert ranked["up_gene"] > 0
    assert ranked["down_gene"] < 0


def test_rank_from_expression_log2fc():
    expr = pd.DataFrame({"s1": [8.0], "s2": [8.0], "s3": [2.0], "s4": [2.0]},
                        index=["gene1"])
    groups = pd.Series({"s1": "A", "s2": "A", "s3": "B", "s4": "B"})
    ranked = rank_from_expression(expr, groups, "A", "B", metric="log2fc")
    assert ranked["gene1"] == pytest.approx(2.0, abs=1e-6)  # log2(8/2) = 2


def test_rank_from_expression_requires_two_samples_per_group():
    expr = pd.DataFrame({"s1": [1.0], "s2": [2.0]}, index=["gene1"])
    groups = pd.Series({"s1": "A", "s2": "B"})
    with pytest.raises(ValueError, match=">=2 samples"):
        rank_from_expression(expr, groups, "A", "B")


def test_load_ranked_list(tmp_path):
    f = tmp_path / "ranked.tsv"
    f.write_text("geneA\t2.5\ngeneB\t-1.0\n# comment\ngeneC\t0.1\n")
    ranked = load_ranked_list(f)
    assert list(ranked.index) == ["geneA", "geneC", "geneB"]  # sorted descending


def test_load_ranked_list_rejects_duplicates(tmp_path):
    f = tmp_path / "ranked.tsv"
    f.write_text("geneA\t1.0\ngeneA\t2.0\n")
    with pytest.raises(ValueError, match="duplicate"):
        load_ranked_list(f)


@pytest.mark.parametrize("bad_value", ["nan", "inf", "-inf", "Infinity"])
def test_load_ranked_list_rejects_non_finite_scores(tmp_path, bad_value):
    # float() accepts these without raising -- a common real source is a DE
    # tool (e.g. DESeq2) reporting Inf/-Inf log-fold-change or NA for genes
    # with zero counts in one group. They must not silently enter the ranking
    # (an Inf score would sort to one end and dominate the running sum).
    f = tmp_path / "ranked.tsv"
    f.write_text(f"geneA\t2.5\ngeneB\t{bad_value}\ngeneC\t-1.0\n")
    with pytest.raises(ValueError, match="non-finite"):
        load_ranked_list(f)


def test_rank_from_expression_rejects_infinite_result(tmp_path):
    # An Inf already present in the input expression matrix survives the
    # epsilon-padding meant to avoid literal division by zero (log2fc's mean
    # ratio propagates a single Inf cell straight through to log2(Inf)=Inf).
    import pandas as pd
    expr = pd.DataFrame({
        "s1": [1.0, np.inf], "s2": [1.0, 5.0],
        "s3": [5.0, 2.0], "s4": [5.0, 2.0],
    }, index=["geneX", "geneY"])
    groups = pd.Series({"s1": "A", "s2": "A", "s3": "B", "s4": "B"})
    with pytest.raises(ValueError, match="non-finite"):
        rank_from_expression(expr, groups, "A", "B", metric="log2fc")


def test_load_expression_matrix(tmp_path):
    f = tmp_path / "expr.tsv"
    f.write_text("gene\ts1\ts2\ngeneA\t1.0\t2.0\ngeneB\t3.0\t4.0\n")
    expr = load_expression_matrix(f)
    assert expr.loc["geneA", "s1"] == 1.0
    assert list(expr.columns) == ["s1", "s2"]
