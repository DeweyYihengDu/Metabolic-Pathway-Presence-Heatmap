"""Fetch and parse KEGG KGML pathway maps: the classic box-and-line pathway
diagram (e.g. ``ko00010`` Glycolysis) and the global metabolic network map
(``ko01100``) that both underlie KEGG's own pathway images.

KGML gives every compound/enzyme node a fixed (x, y) layout position and
lists which reaction(s) each enzyme (KO) entry catalyzes -- this module reads
that structure; :func:`mpph.plot.plot_kgml_map` draws it, colouring enzyme
nodes and their reactions by whether the catalyzing KO is present in one or
two organism/group KO sets you supply (comparative "is this step present"
overlay, not KEGG's own default static colouring).
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import requests

from .kegg import kegg_get

_KO = re.compile(r"K\d{5}")
_MAP_ID = re.compile(r"\d+")


@dataclass(frozen=True)
class KGMLNode:
    """One entry in a KGML map: a compound, an enzyme (ortholog), or a link
    to another pathway map."""
    entry_id: int
    kind: str  # "compound" | "ortholog" | "map"
    kos: frozenset[str]  # non-empty only for kind == "ortholog"
    label: str
    x: float
    y: float
    width: float
    height: float
    shape: str  # graphics "type": "circle" | "rectangle" | "roundrectangle"


@dataclass(frozen=True)
class KGMLReaction:
    """One ``<reaction>`` element: a set of alternative reaction ids (KEGG
    sometimes lists several ``rn:R#####`` under one element) linking
    substrate compound entries to product compound entries."""
    names: frozenset[str]
    substrate_ids: tuple[int, ...]
    product_ids: tuple[int, ...]


@dataclass(frozen=True)
class KGMLPathway:
    map_id: str
    title: str
    nodes: dict[int, KGMLNode]
    reactions: list[KGMLReaction]
    # ortholog entry id -> the reaction id(s) (rn:R#####) it catalyzes, per
    # its own "reaction" attribute -- matched against KGMLReaction.names.
    ortholog_reactions: dict[int, frozenset[str]] = field(default_factory=dict)


def normalize_map_id(map_id: str) -> str:
    """Accept ``"01100"``, ``"ko01100"``, ``"map01100"``, ``"path:ko01100"``.

    KGML is only served for the KO-centric view (``ko#####``), not the bare
    reference map (``map#####``, which 404s on the ``kgml`` endpoint) -- this
    always normalizes to the ``ko#####`` form.
    """
    match = _MAP_ID.search(map_id)
    if not match:
        raise ValueError(f"not a recognizable KEGG map id: {map_id!r}")
    return f"ko{match.group(0)}"


def fetch_kgml(
    session: requests.Session, map_id: str, cache_dir: Path | None,
    *, refresh: bool = False,
) -> str:
    """Fetch the raw KGML XML for a KEGG pathway or the global metabolic map."""
    return kegg_get(session, f"get/{normalize_map_id(map_id)}/kgml",
                    cache_dir, refresh=refresh)


def parse_kgml(xml_text: str) -> KGMLPathway:
    """Parse KGML XML into positioned nodes and the reactions linking them."""
    root = ET.fromstring(xml_text)  # noqa: S314 -- KEGG's own trusted response

    nodes: dict[int, KGMLNode] = {}
    ortholog_reactions: dict[int, frozenset[str]] = {}
    for entry in root.findall("entry"):
        kind = entry.get("type", "")
        if kind not in ("compound", "ortholog", "map"):
            continue
        graphics = entry.find("graphics")
        if graphics is None:
            continue
        entry_id = int(entry.get("id"))
        nodes[entry_id] = KGMLNode(
            entry_id=entry_id, kind=kind,
            kos=frozenset(_KO.findall(entry.get("name", ""))),
            label=graphics.get("name", ""),
            x=float(graphics.get("x", 0.0)), y=float(graphics.get("y", 0.0)),
            width=float(graphics.get("width", 0.0)),
            height=float(graphics.get("height", 0.0)),
            shape=graphics.get("type", "circle"),
        )
        if kind == "ortholog":
            names = frozenset(entry.get("reaction", "").split())
            if names:
                ortholog_reactions[entry_id] = names

    reactions = [
        KGMLReaction(
            names=frozenset(rxn.get("name", "").split()),
            substrate_ids=tuple(int(s.get("id")) for s in rxn.findall("substrate")),
            product_ids=tuple(int(p.get("id")) for p in rxn.findall("product")),
        )
        for rxn in root.findall("reaction")
    ]

    return KGMLPathway(
        map_id=root.get("name", ""), title=root.get("title", ""),
        nodes=nodes, reactions=reactions,
        ortholog_reactions=ortholog_reactions,
    )


def reactions_for_ortholog(pathway: KGMLPathway, entry_id: int) -> list[KGMLReaction]:
    """The ``KGMLReaction`` object(s) a given ortholog entry catalyzes."""
    wanted = pathway.ortholog_reactions.get(entry_id)
    if not wanted:
        return []
    return [r for r in pathway.reactions if r.names & wanted]
