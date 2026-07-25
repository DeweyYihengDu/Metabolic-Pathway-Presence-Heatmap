"""Tests for the pure logic in mpph.kofam -- no pyhmmer needed (the module
never imports it at module scope, only inside annotate_fasta()), so these
run regardless of whether the optional `annotate` extra is installed.
"""
import pytest

from mpph.kofam import (
    Assignment,
    discover_profiles,
    is_significant_hit,
    load_ko_thresholds,
    parse_ko_list,
    read_ko_subset,
    write_mapper_tsv,
)
from mpph.userdata import load_user_kos

# One "full"-type row, one "domain"-type row, one unassignable (no reference
# threshold) row, and a "-" in an unrelated numeric column to confirm that
# doesn't break parsing.
_SAMPLE_KO_LIST = """knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen\tmlen\teff_nseq\tre/pos\tdefinition
K00001\t329.57\tfull\tfull\t0.95\t100\t95\t350\t340\t10.5\t0.59\talcohol dehydrogenase
K00002\t123.45\tdomain\tdomain\t0.80\t50\t45\t200\t190\t8.2\t0.55\tsome domain enzyme
K99999\t-\t-\tfull\t-\t2\t2\t100\t95\t1.0\t0.40\ttoo few reference sequences
"""


def test_parse_ko_list_full_and_domain_and_unassignable_rows():
    entries = parse_ko_list(_SAMPLE_KO_LIST)
    assert set(entries) == {"K00001", "K00002", "K99999"}

    full = entries["K00001"]
    assert full.threshold == 329.57 and full.score_type == "full"
    assert full.assignable is True

    domain = entries["K00002"]
    assert domain.threshold == 123.45 and domain.score_type == "domain"
    assert domain.assignable is True

    unassignable = entries["K99999"]
    assert unassignable.threshold is None and unassignable.score_type is None
    assert unassignable.assignable is False
    assert unassignable.definition == "too few reference sequences"


@pytest.mark.parametrize("full_score,expected", [
    (329.58, True),   # above threshold
    (329.57, True),   # exactly at threshold -- pins down >= vs > (kofam_scan's
                      # own README says "higher than" but its actual source
                      # (result/hit.rb) uses >=; this implements the source)
    (329.56, False),  # just below threshold
])
def test_is_significant_hit_full_score_type_boundary(full_score, expected):
    entries = parse_ko_list(_SAMPLE_KO_LIST)
    assert is_significant_hit(entries["K00001"], full_score, None) is expected


def test_is_significant_hit_domain_score_type_uses_domain_not_full_score():
    entries = parse_ko_list(_SAMPLE_KO_LIST)
    domain_entry = entries["K00002"]  # threshold 123.45, score_type "domain"
    # Full-sequence score is high but irrelevant for a "domain"-type KO;
    # only the domain score is compared against the threshold.
    assert is_significant_hit(domain_entry, 999.0, 100.0) is False
    assert is_significant_hit(domain_entry, 1.0, 123.45) is True
    assert is_significant_hit(domain_entry, 999.0, None) is False  # no domain at all


def test_is_significant_hit_unassignable_ko_is_always_false():
    entries = parse_ko_list(_SAMPLE_KO_LIST)
    unassignable = entries["K99999"]
    assert is_significant_hit(unassignable, 99999.0, 99999.0) is False


def test_write_mapper_tsv_shape_no_header_tab_separated_sorted(tmp_path):
    # Deliberately unsorted input -- write_mapper_tsv must sort on its own,
    # not rely on the caller having already done so.
    assignments = [
        Assignment("geneB", "K00002", 50.0),
        Assignment("geneA", "K00002", 10.0),
        Assignment("geneA", "K00001", 400.0),
    ]
    out = tmp_path / "out.tsv"
    write_mapper_tsv(assignments, out)
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines == ["geneA\tK00001", "geneA\tK00002", "geneB\tK00002"]


def test_write_mapper_tsv_empty_list_writes_empty_file(tmp_path):
    out = tmp_path / "out.tsv"
    write_mapper_tsv([], out)
    assert out.read_text(encoding="utf-8") == ""


def test_write_mapper_tsv_output_chains_into_load_user_kos_with_zero_adapter_code(tmp_path):
    """The whole point of the mapper-TSV shape: mpph.userdata.load_user_kos
    (unmodified, the real function) must accept it directly."""
    assignments = [
        Assignment("gene1", "K00001", 400.0),
        Assignment("gene2", "K00002", 200.0),
        Assignment("gene3", "K00001", 350.0),
    ]
    out = tmp_path / "sampleA_annotated.tsv"
    write_mapper_tsv(assignments, out)

    kos = load_user_kos(out)
    assert kos == {"sampleA_annotated": {"K00001", "K00002"}}


def test_load_ko_thresholds_reads_plain_ko_list(tmp_path):
    (tmp_path / "ko_list").write_text(_SAMPLE_KO_LIST, encoding="utf-8")
    entries = load_ko_thresholds(tmp_path)
    assert set(entries) == {"K00001", "K00002", "K99999"}


def test_load_ko_thresholds_reads_gzipped_ko_list(tmp_path):
    import gzip
    with gzip.open(tmp_path / "ko_list.gz", "wt", encoding="utf-8") as fh:
        fh.write(_SAMPLE_KO_LIST)
    entries = load_ko_thresholds(tmp_path)
    assert set(entries) == {"K00001", "K00002", "K99999"}


def test_load_ko_thresholds_missing_raises_with_actionable_message(tmp_path):
    with pytest.raises(FileNotFoundError, match="annotate --setup-db"):
        load_ko_thresholds(tmp_path)


def test_discover_profiles_finds_profiles_subdirectory(tmp_path):
    (tmp_path / "profiles").mkdir()
    (tmp_path / "profiles" / "K00001.hmm").write_text("fake hmm")
    (tmp_path / "profiles" / "K00002.hmm").write_text("fake hmm")
    profiles = discover_profiles(tmp_path)
    assert sorted(p.name for p in profiles) == ["K00001.hmm", "K00002.hmm"]


def test_discover_profiles_finds_flat_directory(tmp_path):
    (tmp_path / "K00001.hmm").write_text("fake hmm")
    profiles = discover_profiles(tmp_path)
    assert [p.name for p in profiles] == ["K00001.hmm"]


def test_discover_profiles_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="annotate --setup-db"):
        discover_profiles(tmp_path)


def test_read_ko_subset_accepts_bare_ko_list(tmp_path):
    p = tmp_path / "subset.txt"
    p.write_text("K00001\n# a comment\nK00002\n")
    assert read_ko_subset(p) == {"K00001", "K00002"}


def test_read_ko_subset_accepts_kegg_hal_path_list(tmp_path):
    """KEGG's own prokaryote.hal/eukaryote.hal list profile *paths*, not
    bare ids -- the bare K##### regex scan handles both with zero
    format-specific branching."""
    p = tmp_path / "prokaryote.hal"
    p.write_text("profiles/K00001.hmm\nprofiles/K00002.hmm\n")
    assert read_ko_subset(p) == {"K00001", "K00002"}


def test_read_ko_subset_no_matches_raises(tmp_path):
    p = tmp_path / "empty.txt"
    p.write_text("nothing here\n")
    with pytest.raises(ValueError, match="no KO ids"):
        read_ko_subset(p)
