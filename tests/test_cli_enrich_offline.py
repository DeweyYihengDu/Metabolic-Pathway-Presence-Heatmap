"""Offline tests for `mpph enrich`, served from committed KEGG cache fixtures.

Synthetic fixture universe (see tests/fixtures/cache/):
  pathway map00010 <- K00001,K00002,K00003   (3 KOs; study set hits all 3)
  pathway map00020 <- K00004..K00010         (7 KOs; no overlap with study)
  pathway map00030 <- K00020..K00025         (6 KOs; no overlap with study)
  module  M00010   <- K00001,K00002          (2 KOs)
  module  M00020   <- K00003..K00010         (8 KOs)
  module  M00030   <- K00020..K00025         (6 KOs)
  background-organism "testorg" -> all 16 KOs above (link/ko/testorg)
"""
import shutil
from pathlib import Path

import pandas as pd
import pytest

import mpph.cli as cli

FIXTURES = Path(__file__).parent / "fixtures" / "cache"


class _NoNetworkSession:
    headers: dict = {}

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


def test_enrich_kegg_pathway_offline(tmp_path, offline_cache):
    study = tmp_path / "study.txt"
    study.write_text("K00001\nK00002\nK00003\nK99999\n")  # K99999 not in background
    out = tmp_path / "out"

    rc = cli.main([
        "enrich", "--study", str(study), "--background-organism", "testorg",
        "--ontology", "kegg-pathway", "--cache-dir", str(offline_cache),
        "--outdir", str(out), "--format", "svg",
    ])
    assert rc == 0

    # category_id "00010" is all-digits; force str or pandas strips the
    # leading zeros on read-back (bit us once before, in report.py).
    results = pd.read_csv(out / "enrichment_enrichment.csv",
                          dtype={"category_id": str})
    row = results.set_index("category_id").loc["00010"]
    assert row["category_name"] == "Glycolysis / Gluconeogenesis"
    assert row["k_study_hits"] == 3
    assert row["K_background_hits"] == 3
    assert row["N_background_total"] == 16
    assert row["n_study_total"] == 3  # K99999 dropped (not in background)
    assert row["q_value"] < 0.05
    # map00020/map00030 have zero overlap with the study set -> not significant.
    other = results.set_index("category_id").loc[["00020", "00030"]]
    assert (other["q_value"] > 0.05).all()
    assert (out / "enrichment_enrichment.svg").exists()


def test_enrich_kegg_module_offline(tmp_path, offline_cache):
    study = tmp_path / "study.txt"
    study.write_text("K00001\nK00002\n")  # both members of M00010
    out = tmp_path / "out"

    rc = cli.main([
        "enrich", "--study", str(study), "--background-organism", "testorg",
        "--ontology", "kegg-module", "--cache-dir", str(offline_cache),
        "--outdir", str(out), "--label", "modtest",
    ])
    assert rc == 0
    results = pd.read_csv(out / "modtest_enrichment.csv")
    row = results.set_index("category_id").loc["M00010"]
    assert row["category_name"] == "Test module ten"
    assert row["k_study_hits"] == 2
    assert row["K_background_hits"] == 2
    assert row["q_value"] < 0.05


def test_enrich_go_ontology_offline_no_kegg_access(tmp_path, offline_cache):
    # The GO path must not touch KEGG at all; offline_cache's no-network
    # session guards this, but the point is it isn't even invoked.
    # Universe: 20 genes; GO:0006096 has 6 members, all 6 fall in a 10-gene
    # study set (N=20,K=6,n=10,k=6 -> p ~= 0.0054, clearly significant).
    hits = [f"gene{i}" for i in range(1, 7)]
    rest = [f"gene{i}" for i in range(7, 21)]
    go_map = tmp_path / "gene_go.tsv"
    go_map.write_text(
        "".join(f"{g}\tGO:0006096\n" for g in hits)
        + "".join(f"{g}\tGO:0005975\n" for g in rest[:4])  # a second, unrelated term
        + "".join(f"{g}\tGO:0000001\n" for g in rest[4:])  # filler so all 20 are annotated
    )
    study = tmp_path / "study.txt"
    study.write_text("\n".join(hits + rest[:4]) + "\n")  # 6 hits + 4 filler = n=10
    background = tmp_path / "background.txt"
    background.write_text("\n".join(hits + rest) + "\n")  # N=20
    names = tmp_path / "names.tsv"
    names.write_text("GO:0006096\tglycolytic process\n")
    out = tmp_path / "out"

    rc = cli.main([
        "enrich", "--study", str(study), "--background", str(background),
        "--ontology", "go", "--gene-go-map", str(go_map),
        "--go-names", str(names), "--outdir", str(out),
    ])
    assert rc == 0
    results = pd.read_csv(out / "enrichment_enrichment.csv")
    row = results.set_index("category_id").loc["GO:0006096"]
    assert row["category_name"] == "glycolytic process"
    assert row["k_study_hits"] == 6
    assert row["K_background_hits"] == 6
    assert row["N_background_total"] == 20
    assert row["q_value"] < 0.05


def test_enrich_go_without_gene_go_map_errors(tmp_path, offline_cache):
    study = tmp_path / "study.txt"
    study.write_text("gene1\n")
    rc = cli.main(["enrich", "--study", str(study), "--ontology", "go",
                   "--background-organism", "testorg"])
    assert rc == 2


def test_enrich_requires_a_background(tmp_path, offline_cache):
    study = tmp_path / "study.txt"
    study.write_text("K00001\n")
    rc = cli.main(["enrich", "--study", str(study), "--ontology", "kegg-pathway",
                   "--cache-dir", str(offline_cache)])
    assert rc == 2
