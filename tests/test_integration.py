"""End-to-end pipeline test served entirely from committed cache fixtures.

No live KEGG access: the session is replaced with one that raises on any GET,
so the test passes only if every endpoint is satisfied from the fixture cache.
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


def test_presence_pipeline_offline(tmp_path, offline_cache):
    out = tmp_path / "out"
    rc = cli.main([
        "Footaxon", "--match", "word", "--cluster", "--newick",
        "--cache-dir", str(offline_cache), "--outdir", str(out), "--format", "svg",
    ])
    assert rc == 0

    matrix = pd.read_csv(out / "Footaxon_matrix.csv", index_col=0)
    # Barorganism excluded by whole-word match; 3 Footaxon genomes remain.
    assert matrix.shape[0] == 3
    # Overview map (01100) and the non-Metabolism pathway (03430) are dropped;
    # the four Metabolism pathways remain.
    assert set(matrix.columns) == {"00010", "00020", "00190", "00230"}
    assert "03430" not in matrix.columns
    assert "01100" not in matrix.columns

    foo = [i for i in matrix.index if i.startswith("Footaxon alpha")][0]
    fo2 = [i for i in matrix.index if i.startswith("Footaxon beta")][0]
    assert matrix.loc[foo, "00020"] == 1
    assert matrix.loc[fo2, "00020"] == 0

    features = pd.read_csv(out / "Footaxon_features.csv")
    assert set(features["top_category"]) == {"Metabolism"}

    for name in ("Footaxon_qc.csv", "Footaxon_manifest.json",
                 "Footaxon_organism_tree.nwk", "Footaxon_row_order.csv",
                 "Footaxon_col_order.csv", "Footaxon_heatmap.svg"):
        assert (out / name).exists(), name


def test_input_source_mutually_exclusive(tmp_path, offline_cache, capsys):
    codes = tmp_path / "codes.txt"
    codes.write_text("foo\n")
    rc = cli.main(["Footaxon", "--codes", str(codes)])
    assert rc == 2
    assert "one input source" in capsys.readouterr().err
