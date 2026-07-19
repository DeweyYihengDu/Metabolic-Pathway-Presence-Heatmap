"""Tests for enrichment I/O: id lists, gene-GO mapping, GO name tables."""
import pytest

from mpph.enrichment import load_gene_go_map, load_go_names, read_id_list


def test_read_id_list_strips_comments_and_blanks(tmp_path):
    f = tmp_path / "ids.txt"
    f.write_text("K00001\n# a comment\n\nK00002  # inline comment\nK00001\n")
    assert read_id_list(f) == {"K00001", "K00002"}


def test_gene_go_map_long_table_one_row_per_pair(tmp_path):
    f = tmp_path / "map.tsv"
    f.write_text("gene1\tGO:0006096\ngene1\tGO:0005975\ngene2\tGO:0006096\n")
    m = load_gene_go_map(f)
    assert m["gene1"] == {"GO:0006096", "GO:0005975"}
    assert m["gene2"] == {"GO:0006096"}


def test_gene_go_map_multi_go_per_row(tmp_path):
    f = tmp_path / "map.tsv"
    f.write_text("gene1\tGO:0006096;GO:0005975\ngene2\tGO:0006096,GO:0016301\n")
    m = load_gene_go_map(f)
    assert m["gene1"] == {"GO:0006096", "GO:0005975"}
    assert m["gene2"] == {"GO:0006096", "GO:0016301"}


def test_gene_go_map_eggnog_format(tmp_path):
    f = tmp_path / "out.emapper.annotations"
    f.write_text(
        "##  eggNOG-mapper\n"
        "#query\tseed_ortholog\tGO_terms\tKEGG_ko\n"
        "gene1\tx\tGO:0006096,GO:0005975\tko:K00844\n"
        "gene2\tx\t-\tko:K00845\n"
        "gene3\tx\tGO:0006096\tko:K00846\n"
    )
    m = load_gene_go_map(f, fmt="eggnog")
    assert m["gene1"] == {"GO:0006096", "GO:0005975"}
    assert "gene2" not in m  # no GO terms -> not included
    assert m["gene3"] == {"GO:0006096"}


def test_gene_go_map_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_gene_go_map(tmp_path / "nope.tsv")


def test_gene_go_map_no_pairs_raises(tmp_path):
    f = tmp_path / "empty.tsv"
    f.write_text("gene1\tno_go_here\n")
    with pytest.raises(ValueError):
        load_gene_go_map(f)


def test_load_go_names(tmp_path):
    f = tmp_path / "names.tsv"
    f.write_text("GO:0006096\tglycolytic process\nGO:0005975\tcarbohydrate "
                 "metabolic process\nnot_a_go_id\tignored\n")
    names = load_go_names(f)
    assert names["GO:0006096"] == "glycolytic process"
    assert "not_a_go_id" not in names
