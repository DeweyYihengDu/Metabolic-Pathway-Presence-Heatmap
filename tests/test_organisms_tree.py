"""Tests for organism selection and Newick export."""
import numpy as np
from scipy.cluster.hierarchy import linkage

from mpph.organisms import select_by_codes, select_by_taxon
from mpph.tree import linkage_to_newick

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
