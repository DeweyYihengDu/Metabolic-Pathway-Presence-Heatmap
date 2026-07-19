"""Tests for the hypergeometric over-representation (ORA) engine."""
import numpy as np
from scipy.stats import hypergeom

from mpph.enrichment import hypergeometric_enrichment, invert_membership


def _synthetic_universe():
    """100 background genes: cat1 = {bg0..bg9}, catX (everything else) = {bg10..bg99}.

    Every background gene is annotated in exactly one category, so the study
    set's "used" size only shrinks when a study gene is outside the background.
    """
    cat1 = {f"bg{i}" for i in range(10)}
    catx = {f"bg{i}" for i in range(10, 100)}
    return {"cat1": cat1, "catX": catx}


def test_invert_membership():
    item_to_cats = {"a": {"cat1", "cat2"}, "b": {"cat1"}}
    inv = invert_membership(item_to_cats)
    assert inv == {"cat1": {"a", "b"}, "cat2": {"a"}}


def test_strong_enrichment_matches_scipy_directly():
    category_to_items = _synthetic_universe()
    background = {f"bg{i}" for i in range(100)}
    # study: 8 of cat1's 10 members, plus 12 catX members (n=20, k_cat1=8)
    study = {f"bg{i}" for i in range(8)} | {f"bg{i}" for i in range(50, 62)}

    results, stats = hypergeometric_enrichment(study, background, category_to_items)
    row = results.set_index("category_id").loc["cat1"]

    assert stats["n_study_used"] == 20
    assert row["K_background_hits"] == 10
    assert row["N_background_total"] == 100
    assert row["k_study_hits"] == 8
    assert row["fold_enrichment"] == (8 / 20) / (10 / 100)
    # Cross-check against scipy directly (guards off-by-one in the sf() call).
    expected_p = hypergeom.sf(8 - 1, 100, 10, 20)
    assert np.isclose(row["p_value"], expected_p)
    assert row["q_value"] < 0.01  # this one should survive BH-FDR easily


def test_zero_overlap_is_not_flagged_significant():
    # ORA is one-sided (over-representation only): a study set with NO hits in
    # a category must not be reported as "significantly depleted".
    category_to_items = _synthetic_universe()
    background = {f"bg{i}" for i in range(100)}
    study = {f"bg{i}" for i in range(50, 70)}  # all catX, none in cat1

    results, _ = hypergeometric_enrichment(study, background, category_to_items)
    row = results.set_index("category_id").loc["cat1"]
    assert row["k_study_hits"] == 0
    assert row["p_value"] == 1.0
    assert row["fold_enrichment"] == 0.0


def test_min_category_size_excludes_small_categories():
    category_to_items = {"tiny": {"bg0"}, "big": {f"bg{i}" for i in range(1, 20)}}
    background = {f"bg{i}" for i in range(20)}
    study = {"bg0", "bg1", "bg2"}

    results, stats = hypergeometric_enrichment(
        study, background, category_to_items, min_category_size=2)
    assert "tiny" not in set(results["category_id"])
    assert "big" in set(results["category_id"])
    assert stats["n_categories_tested"] == 1


def test_study_genes_outside_background_are_dropped_and_reported():
    category_to_items = _synthetic_universe()
    background = {f"bg{i}" for i in range(100)}
    study = {"bg0", "bg1", "not_in_background_1", "not_in_background_2"}

    _, stats = hypergeometric_enrichment(study, background, category_to_items)
    assert stats["n_study_input"] == 4
    assert stats["n_study_not_in_background"] == 2
    assert stats["n_study_used"] == 2


def test_category_names_applied_when_provided():
    category_to_items = {"map00010": {"bg0", "bg1", "bg2"}}
    background = {f"bg{i}" for i in range(10)}
    study = {"bg0", "bg1"}
    results, _ = hypergeometric_enrichment(
        study, background, category_to_items,
        category_names={"map00010": "Glycolysis / Gluconeogenesis"})
    assert results.iloc[0]["category_name"] == "Glycolysis / Gluconeogenesis"


def test_empty_results_when_nothing_meets_min_size():
    results, stats = hypergeometric_enrichment(
        {"a"}, {"a", "b"}, {"only": {"a", "b"}}, min_category_size=5)
    assert results.empty
    assert stats["n_categories_tested"] == 0
    assert list(results.columns) == [
        "category_id", "category_name", "k_study_hits", "n_study_total",
        "K_background_hits", "N_background_total", "gene_ratio", "bg_ratio",
        "fold_enrichment", "p_value", "study_items", "q_value",
    ]
