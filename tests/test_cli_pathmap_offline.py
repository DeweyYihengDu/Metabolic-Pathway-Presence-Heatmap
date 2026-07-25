"""Offline tests for `mpph pathmap`, served from committed KEGG cache fixtures.

Synthetic pathway ko99999 (see tests/fixtures/cache/get_ko99999_kgml_*.tsv):
  entry 10: ko:K00001              (single-KO enzyme, catalyzes rn:R00001)
  entry 11: ko:K00002 ko:K00003    (two-isozyme enzyme, catalyzes rn:R00002)
  reactions: C00001 --[10]--> C00002 --[11]--> C00003
"""
import json
import shutil
from pathlib import Path

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


def test_pathmap_offline_user_groups(tmp_path, offline_cache):
    group_a = tmp_path / "groupA"
    group_a.mkdir()
    (group_a / "mag1.txt").write_text("K00001\n")  # entry 10 only

    group_b = tmp_path / "groupB"
    group_b.mkdir()
    (group_b / "mag2.txt").write_text("K00003\n")  # one of entry 11's two KOs

    out = tmp_path / "out"
    rc = cli.main([
        "pathmap", "--map", "ko99999",
        "--user-a", str(group_a), "--label-a", "Marine",
        "--user-b", str(group_b), "--label-b", "Freshwater",
        "--cache-dir", str(offline_cache), "--outdir", str(out), "--format", "svg",
    ])
    assert rc == 0
    assert (out / "ko99999_pathmap.svg").exists()

    manifest = json.loads((out / "ko99999_pathmap_manifest.json").read_text())
    assert manifest["map_id"] == "path:ko99999"
    assert manifest["title"] == "Test Pathway"
    assert manifest["group_a"] == {"label": "Marine", "n_kos": 1,
                                   "source": "user:" + str(group_a)}
    assert manifest["group_b"] == {"label": "Freshwater", "n_kos": 1,
                                   "source": "user:" + str(group_b)}
    assert manifest["n_orthologs"] == 2
    assert manifest["n_orthologs_a_only"] == 1  # entry 10 (K00001)
    assert manifest["n_orthologs_b_only"] == 1  # entry 11 (has K00003)
    assert manifest["n_orthologs_both"] == 0
    assert "git_commit" in manifest  # provenance wired in like run/traits


def test_pathmap_requires_exactly_one_source_per_group(tmp_path, offline_cache, capsys):
    rc = cli.main([
        "pathmap", "--map", "ko99999",
        "--cache-dir", str(offline_cache), "--outdir", str(tmp_path / "out"),
    ])
    assert rc == 2
    assert "--codes-a / --user-a" in capsys.readouterr().err


def test_pathmap_rejects_both_sources_for_one_group(tmp_path, offline_cache, capsys):
    group_a = tmp_path / "groupA"
    group_a.mkdir()
    (group_a / "mag1.txt").write_text("K00001\n")
    rc = cli.main([
        "pathmap", "--map", "ko99999",
        "--user-a", str(group_a), "--codes-a", "codes.txt",
        "--user-b", str(group_a),
        "--cache-dir", str(offline_cache), "--outdir", str(tmp_path / "out"),
    ])
    assert rc == 2
    assert "--codes-a / --user-a" in capsys.readouterr().err
