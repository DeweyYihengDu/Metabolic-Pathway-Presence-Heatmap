"""Tests for the pure logic in mpph.kofam -- no pyhmmer needed (the module
never imports it at module scope, only inside annotate_fasta()), so these
run regardless of whether the optional `annotate` extra is installed.
"""
import pytest

from mpph.kofam import (
    Assignment,
    arbitrate_best_per_gene,
    discover_profiles,
    filter_by_margin,
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


# --- arbitrate_best_per_gene --------------------------------------------

def test_arbitration_keeps_the_largest_margin_not_the_largest_score():
    # The whole point: raw scores are not comparable between KOs because
    # thresholds span ~30 to >2000 bits. A 400-bit hit that only just clears a
    # 395-bit threshold is weaker evidence than a 60-bit hit on a 20-bit
    # threshold, and ranking by score would pick the wrong one.
    assignments = [
        Assignment("geneA", "K_high_score", 400.0, margin=5.0),
        Assignment("geneA", "K_high_margin", 60.0, margin=40.0),
    ]
    kept = arbitrate_best_per_gene(assignments)
    assert [a.ko_id for a in kept] == ["K_high_margin"]


def test_arbitration_leaves_single_ko_genes_untouched():
    assignments = [Assignment("geneA", "K00001", 400.0, margin=300.0),
                   Assignment("geneB", "K00002", 50.0, margin=10.0)]
    assert arbitrate_best_per_gene(assignments) == sorted(
        assignments, key=lambda a: (a.gene_id, a.ko_id))


def test_arbitration_is_independent_of_input_order():
    # Hits arrive in profile order, which is arbitrary; the result must not be.
    a = Assignment("g", "K00001", 100.0, margin=10.0)
    b = Assignment("g", "K00002", 100.0, margin=30.0)
    assert arbitrate_best_per_gene([a, b]) == arbitrate_best_per_gene([b, a])


def test_arbitration_breaks_exact_margin_ties_deterministically():
    a = Assignment("g", "K00002", 100.0, margin=10.0)
    b = Assignment("g", "K00001", 100.0, margin=10.0)
    assert [x.ko_id for x in arbitrate_best_per_gene([a, b])] == ["K00001"]
    assert [x.ko_id for x in arbitrate_best_per_gene([b, a])] == ["K00001"]


def test_min_gap_drops_a_gene_whose_top_two_are_too_close():
    assignments = [Assignment("g", "K00001", 100.0, margin=12.0),
                   Assignment("g", "K00002", 100.0, margin=10.0)]
    assert arbitrate_best_per_gene(assignments, min_gap=10.0) == []
    assert len(arbitrate_best_per_gene(assignments, min_gap=1.0)) == 1


def test_min_gap_never_drops_an_unambiguous_single_ko_gene():
    # A gene with one passing KO has no runner-up; requiring a gap must not
    # silently delete it, which would tank recall for no reason.
    assignments = [Assignment("g", "K00001", 100.0, margin=0.5)]
    assert arbitrate_best_per_gene(assignments, min_gap=1000.0) == assignments


def test_arbitration_output_is_sorted_like_write_mapper_tsv_expects():
    assignments = [Assignment("geneB", "K00002", 50.0, margin=5.0),
                   Assignment("geneA", "K00001", 400.0, margin=300.0)]
    kept = arbitrate_best_per_gene(assignments)
    assert [(a.gene_id, a.ko_id) for a in kept] == [
        ("geneA", "K00001"), ("geneB", "K00002")]


def test_arbitration_can_only_shrink_the_call_set():
    assignments = [Assignment("g1", "K1", 10.0, margin=5.0),
                   Assignment("g1", "K2", 10.0, margin=3.0),
                   Assignment("g2", "K3", 10.0, margin=1.0)]
    kept = arbitrate_best_per_gene(assignments)
    assert len(kept) <= len(assignments)
    assert set(kept).issubset(set(assignments))   # never invents a call


def test_filter_by_margin_drops_only_calls_below_the_cut():
    assignments = [Assignment("g1", "K1", 100.0, margin=25.0),
                   Assignment("g2", "K2", 100.0, margin=5.0),
                   Assignment("g3", "K3", 100.0, margin=20.0)]
    kept = filter_by_margin(assignments, 20.0)
    assert [a.ko_id for a in kept] == ["K1", "K3"]     # >= is inclusive


def test_filter_by_margin_is_a_no_op_at_zero_or_below():
    # Must be the identity, not "drop everything with margin < 0", so the
    # default path is untouched.
    assignments = [Assignment("g", "K1", 100.0, margin=0.0)]
    assert filter_by_margin(assignments, 0.0) == assignments
    assert filter_by_margin(assignments, -5.0) == assignments


def test_filter_by_margin_composes_with_arbitration_without_reordering():
    assignments = [Assignment("gB", "K2", 100.0, margin=30.0),
                   Assignment("gA", "K1", 100.0, margin=50.0),
                   Assignment("gA", "K3", 100.0, margin=10.0)]
    out = filter_by_margin(arbitrate_best_per_gene(assignments), 20.0)
    assert [(a.gene_id, a.ko_id) for a in out] == [("gA", "K1"), ("gB", "K2")]
