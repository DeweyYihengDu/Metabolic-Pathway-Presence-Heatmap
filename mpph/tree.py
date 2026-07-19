"""Convert a SciPy linkage matrix to a Newick string (for iTOL / FigTree)."""
from __future__ import annotations

import re

import numpy as np

_NEEDS_QUOTING = re.compile(r"[()\[\],:;']")


def _safe(label: str) -> str:
    """Newick-safe form of a label, preserving it losslessly.

    A label containing Newick-structural characters (parens, brackets, comma,
    colon, semicolon, single quote) is wrapped in single quotes per the
    Newick quoting rule (any internal single quote doubled), rather than
    having those characters deleted or replaced -- destructive stripping can
    silently collide two genuinely different labels (e.g. the microbiology
    convention "[Eubacterium] rectale" and a plain "Eubacterium rectale"
    would both become "Eubacterium_rectale"). A plain space with no other
    special character uses the conventional underscore substitution instead
    of quoting, for the widest tool compatibility on the common case.
    """
    if _NEEDS_QUOTING.search(label):
        return "'" + label.replace("'", "''") + "'"
    return label.replace(" ", "_")


def linkage_to_newick(linkage: np.ndarray, labels: list[str]) -> str:
    """Return a Newick tree for a SciPy ``linkage`` matrix over ``labels``.

    Branch lengths are half the merge-height difference between a node and its
    parent, matching the ultrametric heights UPGMA produces.

    Raises ``ValueError`` if two distinct input labels would collide after
    Newick-safe sanitization, rather than silently merging their identities in
    the exported tree.
    """
    n = len(labels)
    safe_labels = [_safe(lab) for lab in labels]
    seen: dict[str, str] = {}
    for original, safe in zip(labels, safe_labels):
        if safe in seen and seen[safe] != original:
            raise ValueError(
                f"Newick export: {original!r} and {seen[safe]!r} both "
                f"sanitize to the same label {safe!r} -- rename one of them "
                "before exporting (the tree would otherwise silently merge "
                "their identities).")
        seen[safe] = original

    if n == 1:
        return f"({safe_labels[0]});"
    heights = {i: 0.0 for i in range(n)}
    nodes: dict[int, str] = {i: safe_labels[i] for i in range(n)}

    for row, (a, b, dist, _count) in enumerate(linkage):
        a, b = int(a), int(b)
        new_id = n + row
        ha = dist / 2.0 - heights.get(a, 0.0)
        hb = dist / 2.0 - heights.get(b, 0.0)
        nodes[new_id] = f"({nodes[a]}:{ha:.5f},{nodes[b]}:{hb:.5f})"
        heights[new_id] = dist / 2.0

    return nodes[n + len(linkage) - 1] + ";"
