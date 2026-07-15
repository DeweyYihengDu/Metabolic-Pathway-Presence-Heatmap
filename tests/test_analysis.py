"""Tests for the comparative and community analyses."""
import numpy as np
import pandas as pd

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
    coords, explained = pcoa(matrix, metric="euclidean")
    assert coords.shape[0] == 10
    assert explained[0] >= explained[-1]

    groups = pd.Series(["A"] * 5 + ["B"] * 5, index=matrix.index)
    res = permanova(matrix, groups, metric="euclidean", permutations=99, seed=1)
    assert res["pseudo_F"] > 1
    assert res["p_value"] <= 0.05


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
