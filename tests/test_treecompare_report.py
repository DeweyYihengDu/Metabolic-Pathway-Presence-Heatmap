"""Tests for tree comparison, bootstrap support, and the HTML report."""
import json

import numpy as np
import pandas as pd

from mpph.report import build_report
from mpph.treecompare import bootstrap_support, clades, leaves, robinson_foulds


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


def test_rf_root_singleton_outgroup_is_not_lost():
    # 'a' sits alone at the root and belongs to no non-trivial clade; it must
    # still be counted in the shared leaf set.
    res = robinson_foulds("(a,(b,(c,d)));", "(a,(b,(c,d)));")
    assert res["n_shared_leaves"] == 4
    assert res["rf_distance"] == 0
    assert res["n_leaves_a"] == 4 and res["n_leaves_b"] == 4


def test_rf_reports_non_overlapping_leaves():
    res = robinson_foulds("(a,(b,(c,d)));", "(a,(b,(c,e)));")
    assert res["n_shared_leaves"] == 3
    assert res["n_only_in_a"] == 1 and res["n_only_in_b"] == 1
    assert res["only_in_a"] == ["d"]
    assert res["only_in_b"] == ["e"]
    assert res["rooted"] is True


def test_parse_newick_handles_quoted_labels_with_special_characters():
    # mpph.tree.linkage_to_newick quotes a label containing Newick-structural
    # characters (brackets, spaces are fine unquoted) instead of stripping
    # them -- the parser here must round-trip that, including the doubled
    # single-quote escape for a literal quote inside the label.
    nwk = "(('[Eubacterium] rectale':0.5,'O''Brien strain':0.5):1.0,plain_c:1.5);"
    ls = leaves(nwk)
    assert ls == {"[Eubacterium] rectale", "O'Brien strain", "plain_c"}
    assert frozenset({"[Eubacterium] rectale", "O'Brien strain"}) in clades(nwk)


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


def _report_with_names(tmp_path, organism, feature_name, value=1):
    slug = "Evil"
    pd.DataFrame({"00010": [value]}, index=[organism]).to_csv(
        tmp_path / f"{slug}_matrix.csv", index_label="organism")
    pd.DataFrame({"feature_id": ["00010"], "name": [feature_name],
                  "category": ["c"]}).to_csv(
        tmp_path / f"{slug}_features.csv", index=False)
    (tmp_path / f"{slug}_manifest.json").write_text(
        json.dumps({"mode": "presence"}))
    return build_report(tmp_path, slug).read_text(encoding="utf-8")


def test_report_escapes_script_closing_tag(tmp_path):
    html = _report_with_names(tmp_path, "</script><script>alert(1)</script>",
                              "<img src=x onerror=alert(2)>")
    # No raw tag may survive: '<' is escaped in the JSON payload, and names are
    # written with textContent rather than innerHTML.
    assert "</script><script>alert(1)</script>" not in html
    assert "<img src=x onerror=alert(2)>" not in html
    assert "\\u003c" in html  # payload escaped


def test_report_handles_nan_and_unicode(tmp_path):
    import numpy as np
    html = _report_with_names(tmp_path, "Prochlorococcus α β 中文",
                              "naïve", value=np.nan)
    assert "NaN" not in html          # invalid JSON literal must not appear
    assert "null" in html             # NaN became null (unknown)
    assert "中文" in html
