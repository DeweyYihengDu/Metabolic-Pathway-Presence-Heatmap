"""Tests for KGML fetching/parsing (KEGG pathway and global-map diagrams)."""
import pytest

from mpph.kgml import fetch_kgml, normalize_map_id, parse_kgml, reactions_for_ortholog

# A small, hand-written KGML snippet covering everything the parser and
# renderer need: compounds, a single-KO ortholog, a multi-KO ortholog (KEGG
# lists alternative isozymes under one entry), a chain of two reactions
# (C00001 -> C00002 -> C00003), and a map-reference box to another pathway.
_SAMPLE_KGML = """<?xml version="1.0"?>
<pathway name="path:ko99999" org="ko" number="99999" title="Test Pathway">
    <entry id="1" name="cpd:C00001" type="compound">
        <graphics name="C00001" fgcolor="#000" bgcolor="#fff"
             type="circle" x="100" y="100" width="8" height="8"/>
    </entry>
    <entry id="2" name="cpd:C00002" type="compound">
        <graphics name="C00002" fgcolor="#000" bgcolor="#fff"
             type="circle" x="200" y="100" width="8" height="8"/>
    </entry>
    <entry id="3" name="cpd:C00003" type="compound">
        <graphics name="C00003" fgcolor="#000" bgcolor="#fff"
             type="circle" x="300" y="100" width="8" height="8"/>
    </entry>
    <entry id="10" name="ko:K00001" type="ortholog" reaction="rn:R00001">
        <graphics name="K00001" fgcolor="#000" bgcolor="#BFBFFF"
             type="rectangle" x="150" y="100" width="46" height="17"/>
    </entry>
    <entry id="11" name="ko:K00002 ko:K00003" type="ortholog" reaction="rn:R00002">
        <graphics name="K00002..." fgcolor="#000" bgcolor="#BFBFFF"
             type="rectangle" x="250" y="100" width="46" height="17"/>
    </entry>
    <entry id="20" name="path:ko00020" type="map">
        <graphics name="Other pathway" fgcolor="#000" bgcolor="#fff"
             type="roundrectangle" x="300" y="300" width="100" height="30"/>
    </entry>
    <reaction id="1" name="rn:R00001" type="irreversible">
        <substrate id="1" name="cpd:C00001"/>
        <product id="2" name="cpd:C00002"/>
    </reaction>
    <reaction id="2" name="rn:R00002" type="reversible">
        <substrate id="2" name="cpd:C00002"/>
        <product id="3" name="cpd:C00003"/>
    </reaction>
</pathway>
"""


@pytest.mark.parametrize("raw,expected", [
    ("01100", "ko01100"), ("ko01100", "ko01100"), ("map01100", "ko01100"),
    ("path:ko01100", "ko01100"), ("00010", "ko00010"),
])
def test_normalize_map_id(raw, expected):
    assert normalize_map_id(raw) == expected


def test_normalize_map_id_rejects_unrecognizable_input():
    with pytest.raises(ValueError, match="not a recognizable"):
        normalize_map_id("glycolysis")


def test_parse_kgml_node_kinds_and_positions():
    pw = parse_kgml(_SAMPLE_KGML)
    assert pw.map_id == "path:ko99999"
    assert pw.title == "Test Pathway"
    kinds = {n.entry_id: n.kind for n in pw.nodes.values()}
    assert kinds == {1: "compound", 2: "compound", 3: "compound",
                     10: "ortholog", 11: "ortholog", 20: "map"}
    assert pw.nodes[10].kos == frozenset({"K00001"})
    assert pw.nodes[11].kos == frozenset({"K00002", "K00003"})  # multi-KO entry
    assert pw.nodes[1].x == 100.0 and pw.nodes[1].y == 100.0
    assert pw.nodes[20].label == "Other pathway"


def test_parse_kgml_reactions_link_compounds():
    pw = parse_kgml(_SAMPLE_KGML)
    assert len(pw.reactions) == 2
    by_name = {next(iter(r.names)): r for r in pw.reactions}
    assert by_name["rn:R00001"].substrate_ids == (1,)
    assert by_name["rn:R00001"].product_ids == (2,)
    assert by_name["rn:R00002"].substrate_ids == (2,)
    assert by_name["rn:R00002"].product_ids == (3,)


def test_reactions_for_ortholog_resolves_by_reaction_name():
    pw = parse_kgml(_SAMPLE_KGML)
    rxns = reactions_for_ortholog(pw, 10)
    assert len(rxns) == 1
    assert rxns[0].names == frozenset({"rn:R00001"})

    rxns_multi_ko = reactions_for_ortholog(pw, 11)
    assert len(rxns_multi_ko) == 1
    assert rxns_multi_ko[0].names == frozenset({"rn:R00002"})


def test_reactions_for_ortholog_empty_for_unknown_entry():
    pw = parse_kgml(_SAMPLE_KGML)
    assert reactions_for_ortholog(pw, 999) == []


class _CountingSession:
    headers: dict = {}

    def __init__(self, text):
        self.calls = 0
        self.text = text

    def get(self, url, timeout=None):
        self.calls += 1
        class _Resp:
            def __init__(self, text):
                self.text = text
            def raise_for_status(self):
                pass
        return _Resp(self.text)


def test_fetch_kgml_normalizes_and_caches(tmp_path):
    session = _CountingSession(_SAMPLE_KGML)
    first = fetch_kgml(session, "map99999", tmp_path)
    second = fetch_kgml(session, "ko99999", tmp_path)  # same normalized id
    assert first == second == _SAMPLE_KGML
    assert session.calls == 1  # second call served from disk
