"""Convert a SciPy linkage matrix to a Newick string (for iTOL / FigTree)."""
from __future__ import annotations

import numpy as np


def linkage_to_newick(linkage: np.ndarray, labels: list[str]) -> str:
    """Return a Newick tree for a SciPy ``linkage`` matrix over ``labels``.

    Branch lengths are half the merge-height difference between a node and its
    parent, matching the ultrametric heights UPGMA produces.
    """
    n = len(labels)
    if n == 1:
        return f"({_safe(labels[0])});"
    heights = {i: 0.0 for i in range(n)}
    nodes: dict[int, str] = {i: _safe(labels[i]) for i in range(n)}

    for row, (a, b, dist, _count) in enumerate(linkage):
        a, b = int(a), int(b)
        new_id = n + row
        ha = dist / 2.0 - heights.get(a, 0.0)
        hb = dist / 2.0 - heights.get(b, 0.0)
        nodes[new_id] = f"({nodes[a]}:{ha:.5f},{nodes[b]}:{hb:.5f})"
        heights[new_id] = dist / 2.0

    return nodes[n + len(linkage) - 1] + ";"


def _safe(label: str) -> str:
    """Sanitise a label for Newick (strip characters that break parsers)."""
    cleaned = label.replace(" ", "_")
    for ch in "(),:;[]'":
        cleaned = cleaned.replace(ch, "")
    return cleaned
