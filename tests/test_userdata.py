"""Tests for loading user-supplied KO annotations."""
import pytest

from mpph.userdata import load_user_kos


def test_directory_of_ko_lists(tmp_path):
    (tmp_path / "mag1.txt").write_text("gene1\tK00001\ngene2\tK00002\n")
    (tmp_path / "mag2.txt").write_text("K00001\nK99999\n")
    kos = load_user_kos(tmp_path)
    assert kos["mag1"] == {"K00001", "K00002"}
    assert kos["mag2"] == {"K00001", "K99999"}


def test_long_table(tmp_path):
    f = tmp_path / "table.tsv"
    f.write_text("sampleA\tK00001\nsampleA\tK00002\nsampleB\tK00001\n")
    kos = load_user_kos(f)
    assert kos["sampleA"] == {"K00001", "K00002"}
    assert kos["sampleB"] == {"K00001"}


def test_single_sample_file(tmp_path):
    f = tmp_path / "onegenome.txt"
    f.write_text("K00001 K00002 K00003 some free text\n")
    kos = load_user_kos(f)
    assert kos["onegenome"] == {"K00001", "K00002", "K00003"}


def test_no_ko_raises(tmp_path):
    f = tmp_path / "empty.txt"
    f.write_text("no kegg orthologs here\n")
    with pytest.raises(ValueError):
        load_user_kos(f)
