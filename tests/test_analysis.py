"""Tests for the comparative and community analyses."""
import numpy as np
import pandas as pd
import pytest

from mpph.analysis import (
    accumulation_curve,
    benjamini_hochberg,
    cliffs_delta,
    differential_features,
    pairwise_complementarity,
    pan_classify,
    pcoa,
    permanova,
)


def test_benjamini_hochberg_monotone():
    q = benjamini_hochberg(np.array([0.01, 0.02, 0.03, 0.04]))
    assert np.all((q >= 0) & (q <= 1))
    assert q[0] <= q[-1]


def test_bh_nan_does_not_destroy_finite_qvalues():
    # A single NaN p-value must not turn every other q-value into NaN.
    q = benjamini_hochberg(np.array([0.01, np.nan, 0.20]))
    assert np.isfinite(q[0]) and np.isfinite(q[2])
    assert np.isnan(q[1])
    assert q[0] <= q[2]


def test_differential_unknown_not_counted_as_absent():
    # NaN is unknown: group A is 1/1 known-present, not 1/2.
    matrix = pd.DataFrame({"f": [1.0, np.nan, 0.0, 0.0]},
                          index=["a1", "a2", "b1", "b2"])
    groups = pd.Series(["A", "A", "B", "B"], index=matrix.index)
    row = differential_features(matrix, groups, "A", "B").iloc[0]
    assert row["prevalence_a"] == 1.0
    assert row["n_a_known"] == 1 and row["n_a_unknown"] == 1
    assert row["prevalence_b"] == 0.0


def test_pan_prevalence_uses_known_denominator():
    matrix = pd.DataFrame({"f": [1.0, 0.0, np.nan, np.nan]}, index=list("abcd"))
    row = pan_classify(matrix).iloc[0]
    assert row["prevalence"] == 0.5          # 1 present / 2 known, not 1/4
    assert row["n_known"] == 2 and row["n_unknown"] == 2
    assert row["known_fraction"] == 0.5


def test_pan_insufficient_known_fraction_flagged():
    matrix = pd.DataFrame({"f": [1.0, np.nan, np.nan, np.nan]}, index=list("abcd"))
    row = pan_classify(matrix, min_known_fraction=0.8).iloc[0]
    assert row["pan_class"] == "insufficient-data"


def test_cliffs_delta_extremes():
    assert cliffs_delta([3, 4, 5], [0, 1, 2]) == 1.0
    assert cliffs_delta([0, 1, 2], [3, 4, 5]) == -1.0


def test_differential_binary_detects_difference():
    matrix = pd.DataFrame(
        {"f1": [1, 1, 1, 0, 0, 0], "f2": [1, 1, 1, 1, 1, 1]},
        index=[f"o{i}" for i in range(6)],
    )
    groups = pd.Series(["A", "A", "A", "B", "B", "B"], index=matrix.index)
    out = differential_features(matrix, groups, "A", "B")
    f1 = out.set_index("feature_id").loc["f1"]
    assert f1["prevalence_a"] == 1.0
    assert f1["prevalence_b"] == 0.0
    assert f1["p_value"] <= out.set_index("feature_id").loc["f2", "p_value"]


def test_pan_classify_buckets():
    # 8 organisms so the "cloud" bucket (prevalence < 0.15) is reachable.
    matrix = pd.DataFrame({
        "core": [1] * 8,
        "shell": [1, 1, 1, 0, 0, 0, 0, 0],   # 0.375
        "cloud": [1, 0, 0, 0, 0, 0, 0, 0],    # 0.125
    }, index=list("abcdefgh"))
    classes = pan_classify(matrix).set_index("feature_id")["pan_class"].to_dict()
    assert classes["core"] == "core"
    assert classes["shell"] == "shell"
    assert classes["cloud"] == "cloud"


def test_accumulation_curve_shape():
    matrix = pd.DataFrame(np.random.default_rng(0).integers(0, 2, (5, 8)),
                          index=[f"o{i}" for i in range(5)])
    curve = accumulation_curve(matrix, permutations=10, seed=1)
    assert list(curve["n_genomes"]) == [1, 2, 3, 4, 5]
    assert curve["pan_mean"].is_monotonic_increasing


def test_pcoa_and_permanova_separate_groups():
    # Two well-separated clusters in feature space.
    rng = np.random.default_rng(0)
    a = rng.normal(0, 0.1, (5, 6))
    b = rng.normal(5, 0.1, (5, 6))
    matrix = pd.DataFrame(np.vstack([a, b]), index=[f"o{i}" for i in range(10)])
    coords, explained, diag = pcoa(matrix, metric="euclidean")
    assert coords.shape[0] == 10
    assert explained[0] >= explained[-1]
    assert diag["n_positive_axes"] >= 1

    groups = pd.Series(["A"] * 5 + ["B"] * 5, index=matrix.index)
    res = permanova(matrix, groups, metric="euclidean", permutations=99, seed=1)
    assert res["pseudo_F"] > 1
    assert res["p_value"] <= 0.05


def test_pcoa_all_zero_matrix_raises_clearly():
    zeros = pd.DataFrame(np.zeros((3, 2)), index=list("abc"))
    with pytest.raises(ValueError, match="No positive PCoA axes"):
        pcoa(zeros, metric="euclidean")


def test_permanova_missing_group_label_is_excluded():
    rng = np.random.default_rng(0)
    matrix = pd.DataFrame(np.vstack([rng.normal(0, 0.1, (2, 4)),
                                     rng.normal(5, 0.1, (2, 4)),
                                     rng.normal(9, 0.1, (1, 4))]),
                          index=list("abcde"))
    groups = pd.Series(["A", "A", "B", "B", np.nan], index=matrix.index)
    res = permanova(matrix, groups, metric="euclidean", permutations=49, seed=1)
    assert res["n_samples"] == 4          # the NaN-labelled sample is dropped
    assert res["n_groups"] == 2           # NaN is not a group
    assert res["n_excluded_missing_label"] == 1


def test_permanova_rejects_singleton_group():
    matrix = pd.DataFrame(np.random.default_rng(0).normal(size=(3, 4)),
                          index=list("abc"))
    groups = pd.Series(["A", "A", "B"], index=matrix.index)
    with pytest.raises(ValueError, match=">=2 samples per group"):
        permanova(matrix, groups, metric="euclidean", permutations=9, seed=1)


def test_pcoa_zero_distance_pairs_counts_only_true_ties():
    # np.triu(d, 1) zeroes the diagonal AND the whole lower triangle too, so
    # counting "== 0" on that directly overcounts by construction, not by
    # measuring real ties. 3 pairwise-distinct samples must report 0 ties.
    distinct = pd.DataFrame({"f1": [1.0, 0.0, 0.5], "f2": [0.0, 1.0, 0.5],
                             "f3": [0.3, 0.7, 0.9]}, index=["s1", "s2", "s3"])
    _, _, diag = pcoa(distinct, metric="euclidean")
    assert diag["zero_distance_pairs"] == 0

    # s1 and s2 are identical; s3 differs -- exactly one true zero-distance pair.
    one_tie = pd.DataFrame({"f1": [1.0, 1.0, 0.5], "f2": [0.0, 0.0, 0.5]},
                           index=["s1", "s2", "s3"])
    _, _, diag2 = pcoa(one_tie, metric="euclidean")
    assert diag2["zero_distance_pairs"] == 1


def test_permanova_rejects_all_zero_distances():
    identical = pd.DataFrame({"f1": [1.0] * 4, "f2": [0.0] * 4},
                             index=["s1", "s2", "s3", "s4"])
    groups = pd.Series(["A", "A", "B", "B"], index=identical.index)
    with pytest.raises(ValueError, match="every pairwise distance is zero"):
        permanova(identical, groups, metric="euclidean", permutations=99)


@pytest.mark.parametrize("permutations", [0, -5])
def test_permanova_rejects_invalid_permutation_count(permutations):
    matrix = pd.DataFrame({"f1": [1.0, 0.0, 0.5, 0.2], "f2": [0.0, 1.0, 0.5, 0.8]},
                          index=["s1", "s2", "s3", "s4"])
    groups = pd.Series(["A", "A", "B", "B"], index=matrix.index)
    with pytest.raises(ValueError, match="permutations"):
        permanova(matrix, groups, metric="euclidean", permutations=permutations)


def test_pairwise_complementarity_finds_completion():
    # Neither org completes M1 alone; their union does.
    org_kos = {"A": {"K00001"}, "B": {"K00002"}}
    module_defs = {"M00001": ("K00001 K00002", "cat", "Pathway")}
    out = pairwise_complementarity(org_kos, module_defs)
    assert len(out) == 1
    row = out.iloc[0]
    assert row["individual_max"] == 0.5
    assert row["combined_score"] == 1.0
    assert bool(row["newly_completed"]) is True
