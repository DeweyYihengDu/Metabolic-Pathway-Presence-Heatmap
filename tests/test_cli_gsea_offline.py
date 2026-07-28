"""Offline tests for `mpph gsea`, served from the same committed KEGG cache
fixtures used by test_cli_enrich_offline.py:

  pathway map00010 <- K00001,K00002,K00003   (3 KOs)
  pathway map00020 <- K00004..K00010         (7 KOs)
  pathway map00030 <- K00020..K00025         (6 KOs)
"""
import shutil
from pathlib import Path

import pandas as pd
import pytest

from mpph import cli

FIXTURES = Path(__file__).parent / "fixtures" / "cache"


class _NoNetworkSession:
    def __init__(self):
        self.headers = {}

    def get(self, url, *a, **k):
        raise AssertionError(f"unexpected network call in offline test: {url}")

    def mount(self, *a, **k):
        pass


@pytest.fixture
def offline_cache(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    shutil.copytree(FIXTURES, cache)
    monkeypatch.setattr(cli, "make_session", lambda: _NoNetworkSession())
    return cache


def _ranked_list_favoring_map00010(tmp_path) -> Path:
    # map00010's 3 KOs at the very top of the ranking; map00020's 7 KOs in the
    # middle; map00030's 6 KOs at the very bottom.
    rows = (
        [("K00001", 100), ("K00002", 99), ("K00003", 98)]
        + [(f"K0000{i}", 50 - i) for i in range(4, 10)] + [("K00010", 40)]
        + [(f"K000{20 + i}", -50 - i) for i in range(6)]
    )
    f = tmp_path / "ranked.tsv"
    f.write_text("".join(f"{g}\t{s}\n" for g, s in rows))
    return f


def test_gsea_kegg_pathway_offline(tmp_path, offline_cache):
    ranked = _ranked_list_favoring_map00010(tmp_path)
    out = tmp_path / "out"
    rc = cli.main([
        "gsea", "--ranked-list", str(ranked), "--ontology", "kegg-pathway",
        "--min-size", "2", "--max-size", "20", "--permutations", "200",
        "--seed", "1", "--cache-dir", str(offline_cache), "--outdir", str(out),
        "--format", "svg",
    ])
    assert rc == 0

    results = pd.read_csv(out / "gsea_gsea.csv", dtype={"category_id": str})
    row = results.set_index("category_id").loc["00010"]
    assert row["category_name"] == "Glycolysis / Gluconeogenesis"
    assert row["ES"] > 0.9          # concentrated at the very top
    assert row["NES"] > 0
    assert row["q_value"] < 0.05
    assert set(row["leading_edge_genes"].split(",")) == {
        "K00001", "K00002", "K00003"}

    bottom = results.set_index("category_id").loc["00030"]
    assert bottom["ES"] < -0.9      # concentrated at the very bottom
    assert bottom["q_value"] < 0.05

    assert (out / "gsea_gsea_summary.svg").exists()
    assert (out / "gsea_gsea_top.svg").exists()


def test_gsea_kegg_pathway_top_category_filters_to_metabolism(tmp_path, offline_cache):
    # map00030 (the bottom-ranked pathway) is classified under "Genetic
    # Information Processing" in the br08901 fixture -- restricting to
    # Metabolism should drop it from the tested universe entirely.
    ranked = _ranked_list_favoring_map00010(tmp_path)
    out = tmp_path / "out"
    rc = cli.main([
        "gsea", "--ranked-list", str(ranked), "--ontology", "kegg-pathway",
        "--top-category", "Metabolism", "--min-size", "2", "--max-size", "20",
        "--permutations", "200", "--seed", "1", "--cache-dir", str(offline_cache),
        "--outdir", str(out),
    ])
    assert rc == 0
    results = pd.read_csv(out / "gsea_gsea.csv", dtype={"category_id": str})
    assert set(results["category_id"]) == {"00010", "00020"}
    assert (results["top_category"] == "Metabolism").all()


def test_gsea_kegg_pathway_default_has_top_category_column_without_filtering(
        tmp_path, offline_cache):
    ranked = _ranked_list_favoring_map00010(tmp_path)
    out = tmp_path / "out"
    rc = cli.main([
        "gsea", "--ranked-list", str(ranked), "--ontology", "kegg-pathway",
        "--min-size", "2", "--max-size", "20", "--permutations", "200",
        "--seed", "1", "--cache-dir", str(offline_cache), "--outdir", str(out),
    ])
    assert rc == 0
    results = pd.read_csv(out / "gsea_gsea.csv", dtype={"category_id": str})
    assert set(results["category_id"]) == {"00010", "00020", "00030"}
    by_id = results.set_index("category_id")["top_category"]
    assert by_id["00010"] == "Metabolism"
    assert by_id["00030"] == "Genetic Information Processing"


def test_gsea_top_category_warns_for_irrelevant_ontology(
        tmp_path, offline_cache, capsys):
    ranked = _ranked_list_favoring_map00010(tmp_path)
    rc = cli.main([
        "gsea", "--ranked-list", str(ranked), "--ontology", "kegg-module",
        "--top-category", "Metabolism", "--min-size", "2", "--max-size", "20",
        "--permutations", "200", "--seed", "1", "--cache-dir", str(offline_cache),
        "--outdir", str(tmp_path / "out"), "--label", "modtest",
    ])
    assert rc == 0
    assert "has no meaning for --ontology kegg-module" in capsys.readouterr().err


def test_gsea_kegg_module_offline(tmp_path, offline_cache):
    ranked = _ranked_list_favoring_map00010(tmp_path)
    out = tmp_path / "out"
    rc = cli.main([
        "gsea", "--ranked-list", str(ranked), "--ontology", "kegg-module",
        "--min-size", "2", "--max-size", "20", "--permutations", "200",
        "--seed", "1", "--cache-dir", str(offline_cache), "--outdir", str(out),
        "--label", "modtest",
    ])
    assert rc == 0
    results = pd.read_csv(out / "modtest_gsea.csv")
    row = results.set_index("category_id").loc["M00010"]  # K00001,K00002
    assert row["category_name"] == "Test module ten"
    assert row["ES"] > 0.9


def test_gsea_go_ontology_offline_no_kegg_access(tmp_path, offline_cache):
    go_map = tmp_path / "gene_go.tsv"
    go_map.write_text("K00001\tGO:0006096\nK00002\tGO:0006096\nK00003\tGO:0006096\n"
                      + "".join(f"K0000{i}\tGO:0005975\n" for i in range(4, 10)))
    ranked = _ranked_list_favoring_map00010(tmp_path)
    out = tmp_path / "out"
    rc = cli.main([
        "gsea", "--ranked-list", str(ranked), "--ontology", "go",
        "--gene-go-map", str(go_map), "--min-size", "2", "--max-size", "20",
        "--permutations", "200", "--seed", "1", "--outdir", str(out),
    ])
    assert rc == 0
    results = pd.read_csv(out / "gsea_gsea.csv")
    row = results.set_index("category_id").loc["GO:0006096"]
    assert row["ES"] > 0.9


def test_gsea_expression_ranking_offline(tmp_path, offline_cache):
    # 3 genes up in group A (map00010), 3 genes up in group B (map00030);
    # signal2noise ranking should reproduce the same enrichment story.
    genes = ["K00001", "K00002", "K00003", "K00020", "K00021", "K00022"]
    expr = pd.DataFrame({
        "a1": [10, 10, 10, 1, 1, 1], "a2": [11, 9, 10, 1, 2, 1],
        "b1": [1, 1, 1, 10, 10, 10], "b2": [2, 1, 1, 11, 9, 10],
    }, index=genes)
    expr_file = tmp_path / "expr.tsv"
    expr.to_csv(expr_file, sep="\t", index_label="gene")
    meta_file = tmp_path / "meta.tsv"
    meta_file.write_text("sample_id\tgroup\na1\tA\na2\tA\nb1\tB\nb2\tB\n")
    out = tmp_path / "out"

    rc = cli.main([
        "gsea", "--expression", str(expr_file), "--metadata", str(meta_file),
        "--group-column", "group", "--group-a", "A", "--group-b", "B",
        "--ontology", "kegg-pathway", "--min-size", "2", "--max-size", "20",
        "--permutations", "200", "--seed", "1", "--cache-dir", str(offline_cache),
        "--outdir", str(out),
    ])
    assert rc == 0
    assert (out / "gsea_ranked_list.tsv").exists()
    results = pd.read_csv(out / "gsea_gsea.csv", dtype={"category_id": str})
    row = results.set_index("category_id").loc["00010"]
    assert row["ES"] > 0.9  # higher in group A -> ranked to the top -> map00010 hit


def test_gsea_requires_exactly_one_rank_source(tmp_path):
    # argparse's own required mutually-exclusive-group check fires at parse
    # time (before our code runs), so it exits directly rather than returning.
    with pytest.raises(SystemExit) as exc:
        cli.main(["gsea", "--ontology", "kegg-pathway"])
    assert exc.value.code == 2


def test_gsea_expression_requires_group_args(tmp_path):
    expr_file = tmp_path / "expr.tsv"
    expr_file.write_text("gene\ta1\tb1\ng1\t1\t2\n")
    rc = cli.main(["gsea", "--expression", str(expr_file), "--ontology",
                   "kegg-pathway"])
    assert rc == 2


def test_gsea_go_without_gene_go_map_errors(tmp_path):
    ranked = tmp_path / "ranked.tsv"
    ranked.write_text("g1\t1.0\n")
    rc = cli.main(["gsea", "--ranked-list", str(ranked), "--ontology", "go"])
    assert rc == 2
