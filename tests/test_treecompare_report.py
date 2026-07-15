"""Tests for tree comparison, bootstrap support, and the HTML report."""
import json

import numpy as np
import pandas as pd

from mpph.report import build_report
from mpph.treecompare import bootstrap_support, clades, robinson_foulds


def test_clades_extraction():
    c = clades("((a,b),(c,d));")
    assert frozenset({"a", "b"}) in c
    assert frozenset({"c", "d"}) in c


def test_rf_identical_and_different():
    same = robinson_foulds("((a,b),(c,d));", "((a,b),(c,d));")
    assert same["rf_distance"] == 0
    diff = robinson_foulds("((a,b),(c,d));", "((a,c),(b,d));")
    assert diff["rf_distance"] > 0
    assert diff["n_shared_leaves"] == 4


def test_bootstrap_support_high_for_real_clade():
    rng = np.random.default_rng(0)
    a = rng.normal(0, 0.05, (4, 30))
    b = rng.normal(8, 0.05, (4, 30))
    matrix = pd.DataFrame(np.vstack([a, b]),
                          index=["a1", "a2", "a3", "a4", "b1", "b2", "b3", "b4"])
    sup = bootstrap_support(matrix, metric="euclidean", n_boot=50, seed=1)
    top = sup.iloc[0]
    assert top["support"] >= 0.8


def test_build_report(tmp_path):
    slug = "Demo"
    pd.DataFrame({"00010": [1, 0], "00020": [1, 1]},
                 index=["Org one", "Org two"]).to_csv(
        tmp_path / f"{slug}_matrix.csv", index_label="organism")
    pd.DataFrame({"feature_id": ["00010", "00020"],
                  "name": ["Glycolysis", "TCA"],
                  "category": ["Carbohydrate metabolism"] * 2}).to_csv(
        tmp_path / f"{slug}_features.csv", index=False)
    pd.DataFrame({"organism": ["Org one", "Org two"],
                  "n_annotated_features": [2, 1],
                  "status": ["included", "included"]}).to_csv(
        tmp_path / f"{slug}_qc.csv", index=False)
    (tmp_path / f"{slug}_manifest.json").write_text(
        json.dumps({"mode": "presence", "mpph_version": "test"}))

    out = build_report(tmp_path, slug)
    assert out.exists()
    html = out.read_text(encoding="utf-8")
    assert "Org one" in html and "Glycolysis" in html
    assert "<script>" in html  # self-contained interactive report
