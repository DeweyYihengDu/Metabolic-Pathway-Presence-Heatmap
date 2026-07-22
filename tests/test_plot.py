"""Tests for figure labelling (rendered to SVG so text is inspectable)."""
import re

import pandas as pd

import numpy as np

from mpph.kgml import parse_kgml
from mpph.plot import (
    plot_accumulation,
    plot_enrichment,
    plot_gsea_running,
    plot_gsea_summary,
    plot_kgml_map,
    plot_matrix,
    plot_ordination,
    plot_prevalence,
    plot_volcano,
)
from mpph.tree import linkage_to_newick

_SAMPLE_KGML = """<?xml version="1.0"?>
<pathway name="path:ko99999" org="ko" number="99999" title="Test Pathway">
    <entry id="1" name="cpd:C00001" type="compound">
        <graphics name="C00001" type="circle" x="100" y="100" width="8" height="8"/>
    </entry>
    <entry id="2" name="cpd:C00002" type="compound">
        <graphics name="C00002" type="circle" x="200" y="100" width="8" height="8"/>
    </entry>
    <entry id="3" name="cpd:C00003" type="compound">
        <graphics name="C00003" type="circle" x="300" y="100" width="8" height="8"/>
    </entry>
    <entry id="10" name="ko:K00001" type="ortholog" reaction="rn:R00001">
        <graphics name="K00001" type="rectangle" x="150" y="100" width="46" height="17"/>
    </entry>
    <entry id="11" name="ko:K00002 ko:K00003" type="ortholog" reaction="rn:R00002">
        <graphics name="K00002..." type="rectangle" x="250" y="100" width="46" height="17"/>
    </entry>
    <entry id="20" name="path:ko00020" type="map">
        <graphics name="Other pathway" type="roundrectangle" x="300" y="300" width="100" height="30"/>
    </entry>
    <reaction id="1" name="rn:R00001" type="irreversible">
        <substrate id="1" name="cpd:C00001"/>
        <product id="2" name="cpd:C00002"/>
    </reaction>
    <reaction id="2" name="rn:R00002" type="reversible">
        <substrate id="2" name="cpd:C00002"/>
        <product id="3" name="cpd:C00003"/>
    </reaction>
</pathway>
"""


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


def test_completeness_heatmap_no_nan_has_no_unknown_legend(tmp_path):
    out = tmp_path / "no_nan.svg"
    plot_matrix(_matrix(), {}, out, "t", "s", mode="completeness")
    svg = out.read_text(encoding="utf-8")
    assert "unknown" not in svg.lower()


def test_completeness_heatmap_nan_is_not_rendered_black(tmp_path):
    # An undetermined module score (mpph.modules.module_completeness can
    # return NaN) must not silently render as black -- NaN maps to a fully
    # transparent RGBA via the colormap's default "bad" color, so naively
    # dropping the alpha channel (rgba[:, :, :3]) leaves opaque black.
    df = pd.DataFrame({"M1": [1.0, 0.5, np.nan], "M2": [0.2, 1.0, 0.8]},
                      index=["a", "b", "c"])
    out = tmp_path / "nan.svg"
    plot_matrix(df, {}, out, "t", "s", mode="completeness")
    svg = out.read_text(encoding="utf-8")
    assert "unknown / not assessed" in svg
    assert not re.search(r"#000000|rgb\(0%?,\s*0%?,\s*0%?\)", svg, re.IGNORECASE)


def test_presence_heatmap_nan_differs_from_absent(tmp_path):
    df = pd.DataFrame({"00010": [1.0, 0.0, np.nan], "00020": [1.0, 1.0, 0.0]},
                      index=["a", "b", "c"])
    out = tmp_path / "pres_nan.svg"
    plot_matrix(df, {"00010": "Carbohydrate metabolism",
                     "00020": "Energy metabolism"},
                out, "t", "s", mode="presence")
    svg = out.read_text(encoding="utf-8")
    assert "unknown / not assessed" in svg


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


def test_plot_kgml_map_categorizes_by_group_and_unions_multi_ko_entries(tmp_path):
    pw = parse_kgml(_SAMPLE_KGML)
    out = tmp_path / "kgml.png"
    # K00001 (entry 10) in group A only; K00003 (one of entry 11's two KOs,
    # a multi-isozyme entry) in group B -- entry 11 must count as "both"
    # since only ONE of its KOs needs to match either group.
    stats = plot_kgml_map(pw, {"K00001"}, {"K00003"}, out, "kgml test",
                          group_a_label="Marine", group_b_label="Freshwater")
    assert out.exists()
    assert stats["n_orthologs"] == 2
    assert stats["n_orthologs_a_only"] == 1   # entry 10 (K00001)
    assert stats["n_orthologs_both"] == 0
    assert stats["n_orthologs_b_only"] == 1    # entry 11 (has K00003, not K00001)
    assert stats["n_orthologs_neither"] == 0


def test_plot_kgml_map_both_when_group_a_and_b_overlap_in_one_entry(tmp_path):
    pw = parse_kgml(_SAMPLE_KGML)
    out = tmp_path / "kgml_both.png"
    # entry 11 has KOs {K00002, K00003}; group A brings K00002, group B
    # brings K00003 -- the union means this single entry is "both".
    stats = plot_kgml_map(pw, {"K00002"}, {"K00003"}, out, "kgml both test")
    assert stats["n_orthologs_both"] == 1
    assert stats["n_orthologs_a_only"] == 0
    assert stats["n_orthologs_b_only"] == 0


def test_plot_kgml_map_all_absent_reports_neither(tmp_path):
    pw = parse_kgml(_SAMPLE_KGML)
    out = tmp_path / "kgml_neither.png"
    stats = plot_kgml_map(pw, {"K99999"}, {"K88888"}, out, "kgml neither test")
    assert stats["n_orthologs_neither"] == 2
    assert stats["n_orthologs_both"] == stats["n_orthologs_a_only"] == 0
    assert stats["n_orthologs_b_only"] == 0
