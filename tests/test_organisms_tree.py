"""Tests for organism selection and Newick export."""
import numpy as np
import pytest
from scipy.cluster.hierarchy import linkage

from mpph.organisms import select_by_codes, select_by_taxon
from mpph.tree import _safe, linkage_to_newick

GENOMES = [
    ("vch", "Vibrio cholerae O1"),
    ("vibn", "Vibrionimonas magnilacihabitans"),
    ("pmt", "Prochlorococcus marinus MIT 9313"),
]


def test_whole_word_match_excludes_lookalikes():
    hits = {c for c, _ in select_by_taxon(GENOMES, "Vibrio", match="word")}
    assert "vch" in hits
    assert "vibn" not in hits  # "Vibrionimonas" must not match "Vibrio"


def test_prefix_match():
    hits = {c for c, _ in select_by_taxon(GENOMES, "Vibrio cholerae", match="prefix")}
    assert hits == {"vch"}


def test_exact_match():
    assert {c for c, _ in select_by_taxon(GENOMES, "Vibrio cholerae O1",
                                          match="exact")} == {"vch"}
    assert select_by_taxon(GENOMES, "Vibrio", match="exact") == []


def test_select_by_codes_preserves_order():
    picked = select_by_codes(GENOMES, ["pmt", "vch"])
    assert [c for c, _ in picked] == ["pmt", "vch"]


def test_newick_is_wellformed():
    data = np.array([[0.0, 0.0], [0.0, 1.0], [5.0, 5.0]])
    z = linkage(data, method="average")
    nwk = linkage_to_newick(z, ["A", "B", "C"])
    assert nwk.endswith(";")
    assert nwk.count("(") == nwk.count(")")
    for label in ("A", "B", "C"):
        assert label in nwk


def test_safe_quotes_rather_than_strips_special_characters():
    # Stripping brackets/parens instead of quoting can silently collide two
    # genuinely different organism names (a real microbiology convention:
    # "[Eubacterium] rectale" vs. a plain "Eubacterium rectale").
    assert _safe("[Eubacterium] rectale") == "'[Eubacterium] rectale'"
    assert _safe("Eubacterium rectale") == "Eubacterium_rectale"
    assert _safe("[Eubacterium] rectale") != _safe("Eubacterium rectale")
    # A literal single quote inside a label is doubled, not dropped.
    assert _safe("O'Brien strain") == "'O''Brien strain'"


def test_linkage_to_newick_rejects_colliding_labels():
    data = np.array([[0.0, 0.0], [0.0, 1.0], [5.0, 5.0]])
    z = linkage(data, method="average")
    # "A B" and "A_B" both sanitize to "A_B" -- quoting alone can't
    # disambiguate a literal underscore from a space-turned-underscore.
    with pytest.raises(ValueError, match="both sanitize to"):
        linkage_to_newick(z, ["A B", "A_B", "C"])
