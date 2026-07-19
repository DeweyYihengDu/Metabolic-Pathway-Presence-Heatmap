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


def test_ambiguous_two_column_file_is_not_split_into_fake_samples(tmp_path):
    # A single-genome KofamScan --format mapper file: one row per gene, every
    # gene id unique. Grouping by column 1 (the old heuristic: >1 distinct
    # key => "long table") would wrongly turn each gene into its own sample.
    f = tmp_path / "my_mag.txt"
    f.write_text("geneA\tK00001\ngeneB\tK00002\ngeneC\tK00003\n")
    kos = load_user_kos(f, fmt="auto")
    assert kos == {"my_mag": {"K00001", "K00002", "K00003"}}


def test_input_format_long_forces_long_table_reading(tmp_path):
    # Escape hatch for a genuine multi-sample table where every sample
    # happens to contribute exactly one row (indistinguishable from the
    # per-gene case above by structure alone -- the caller must say so).
    f = tmp_path / "my_mag.txt"
    f.write_text("geneA\tK00001\ngeneB\tK00002\ngeneC\tK00003\n")
    kos = load_user_kos(f, fmt="long")
    assert kos == {"geneA": {"K00001"}, "geneB": {"K00002"}, "geneC": {"K00003"}}


def test_directory_duplicate_stem_rejected(tmp_path):
    (tmp_path / "sample.txt").write_text("K00001\nK00002\n")
    (tmp_path / "sample.tsv").write_text("K00099\n")
    with pytest.raises(ValueError, match="duplicate sample name"):
        load_user_kos(tmp_path)


def test_eggnog_format(tmp_path):
    f = tmp_path / "MAG1.emapper.annotations"
    f.write_text(
        "##  eggNOG-mapper\n"
        "#query\tseed_ortholog\tKEGG_ko\tKEGG_Pathway\n"
        "gene1\tx\tko:K00844,ko:K12407\tmap00010\n"
        "gene2\tx\t-\t-\n"
        "gene3\tx\tko:K01810\tmap00010\n"
    )
    kos = load_user_kos(f, fmt="eggnog")
    assert kos["MAG1.emapper"] == {"K00844", "K12407", "K01810"}
