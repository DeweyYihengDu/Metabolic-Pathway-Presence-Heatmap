"""Tests for figure labelling (rendered to SVG so text is inspectable)."""
import re

import pandas as pd

from mpph.plot import (
    plot_accumulation,
    plot_enrichment,
    plot_gsea_running,
    plot_gsea_summary,
    plot_matrix,
    plot_ordination,
    plot_prevalence,
    plot_volcano,
)
from mpph.tree import linkage_to_newick


def _matrix():
    return pd.DataFrame(
        {"M1": [1.0, 0.5, 0.0], "M2": [0.2, 1.0, 0.8], "M3": [0.0, 0.3, 1.0]},
        index=["orgA", "orgB", "orgC"],
    )


def test_custom_value_and_strip_labels(tmp_path):
    out = tmp_path / "fig.svg"
    plot_matrix(_matrix(), {"M1": "Nitrogen cycle", "M2": "Sulfur cycle",
                            "M3": "Carbon cycle"},
                out, "Traits", "sub", cluster=True, mode="completeness",
                value_label="trait completeness", strip_label="trait category")
    svg = out.read_text(encoding="utf-8")
    assert "trait completeness" in svg
    assert "trait category" in svg
    assert "module completeness" not in svg


def test_default_presence_label(tmp_path):
    out = tmp_path / "pres.svg"
    plot_matrix(pd.DataFrame({"00010": [1, 0], "00020": [1, 1]},
                             index=["a", "b"]),
                {"00010": "Carbohydrate metabolism",
                 "00020": "Energy metabolism"},
                out, "Presence", "sub", mode="presence")
    svg = out.read_text(encoding="utf-8")
    assert "KEGG functional category" in svg


def test_analysis_plots_produce_files(tmp_path):
    diff = pd.DataFrame({"feature_id": ["a", "b", "c"],
                         "prevalence_diff": [0.8, -0.2, 0.1],
                         "q_value": [0.001, 0.9, 0.4]})
    plot_volcano(diff, tmp_path / "v.png", "A vs B")
    assert (tmp_path / "v.png").exists()

    pan = pd.DataFrame({"feature_id": ["a", "b", "c"],
                        "prevalence": [1.0, 0.5, 0.1],
                        "pan_class": ["core", "shell", "cloud"]})
    plot_prevalence(pan, tmp_path / "p.png", "pan")
    assert (tmp_path / "p.png").exists()

    acc = pd.DataFrame({"n_genomes": [1, 2, 3], "pan_mean": [3, 5, 6],
                        "core_mean": [3, 2, 1]})
    plot_accumulation(acc, tmp_path / "a.png", "acc")
    assert (tmp_path / "a.png").exists()


def test_plot_enrichment_produces_file(tmp_path):
    results = pd.DataFrame({
        "category_id": ["map00010", "map00020"],
        "category_name": ["Glycolysis", "TCA cycle"],
        "gene_ratio": ["8/20", "1/20"],
        "q_value": [0.001, 0.8],
    })
    out = tmp_path / "e.png"
    plot_enrichment(results, out, "Enrichment test")
    assert out.exists()


def test_plot_enrichment_empty_raises(tmp_path):
    import pytest
    with pytest.raises(ValueError, match="Nothing to plot"):
        plot_enrichment(pd.DataFrame(), tmp_path / "e.png", "empty")


def test_plot_gsea_running_produces_file(tmp_path):
    import numpy as np
    n = 30
    running = np.concatenate([np.linspace(0, 1, 10), np.linspace(1, 0, 20)])
    hits = np.zeros(n, dtype=bool)
    hits[:10] = True
    scores = np.linspace(2, -2, n)
    out = tmp_path / "gsea_run.png"
    plot_gsea_running(scores, running, hits, out, "GSEA test", "subtitle")
    assert out.exists()


def test_plot_gsea_summary_produces_file(tmp_path):
    results = pd.DataFrame({
        "category_id": ["A", "B"], "category_name": ["catA", "catB"],
        "NES": [2.0, -1.5], "q_value": [0.01, 0.2],
    })
    out = tmp_path / "gsea_summary.png"
    plot_gsea_summary(results, out, "GSEA summary test")
    assert out.exists()


def test_plot_gsea_summary_empty_raises(tmp_path):
    import pytest
    with pytest.raises(ValueError, match="Nothing to plot"):
        plot_gsea_summary(pd.DataFrame(), tmp_path / "e.png", "empty")


def test_row_link_labels_match_linkage_leaf_order(tmp_path):
    # C and A sit almost on top of each other; B is a distant outlier.
    # Fed in as C, A, B so a linkage-index/display-order mix-up would
    # visibly attach the wrong name to the wrong branch.
    df = pd.DataFrame({"f1": [0.0, 0.05, 10.0], "f2": [0.0, 0.05, 10.0]},
                      index=["C", "A", "B"])
    layout = plot_matrix(df, {"f1": "cat", "f2": "cat"}, tmp_path / "t.svg",
                        "title", "sub", cluster=True, mode="completeness",
                        metric="euclidean")
    assert layout["row_link_labels"] == ["C", "A", "B"]

    nwk = linkage_to_newick(layout["row_link"], layout["row_link_labels"])
    # The innermost clade (no nested parens) must be exactly the two truly
    # close leaves {A, C}; the distant outlier B must sit outside it.
    inner = re.search(r"\(([^()]+)\)", nwk).group(1)
    inner_leaves = set(re.findall(r"[A-Z](?=:)", inner))
    assert inner_leaves == {"A", "C"}


def test_plot_ordination_rank_one_solution_does_not_crash(tmp_path):
    # Exactly 2 samples (or other degenerate distance structure) yields a
    # PCoA result with only one positive axis -- a 1-column coords frame.
    coords = pd.DataFrame({"PCo1": [-0.5, 0.5]}, index=["s1", "s2"])
    out = tmp_path / "ord.png"
    plot_ordination(coords, [1.0], out, "rank-one test")
    assert out.exists()
