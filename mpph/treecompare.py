"""Compare a functional dendrogram to a reference tree, and bootstrap support.

Newick parsing here is intentionally small: it extracts the set of non-trivial
bipartitions (clades) so we can compute a Robinson-Foulds distance on the shared
leaf set. Bootstrap support resamples feature columns and re-clusters.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage


# --------------------------------------------------------------------------- #
# Minimal Newick -> clades
# --------------------------------------------------------------------------- #
def _parse_newick(newick: str):
    """Parse a Newick string into nested lists of leaf-name strings."""
    tokens = re.findall(r"[(),]|[^(),;]+", newick.strip().rstrip(";"))
    pos = 0

    def parse():
        nonlocal pos
        node = []
        while pos < len(tokens):
            tok = tokens[pos]
            if tok == "(":
                pos += 1
                node.append(parse())
            elif tok == ")":
                pos += 1
                # optional internal label/branch length after ')'
                if pos < len(tokens) and tokens[pos] not in "(),":
                    pos += 1
                return node
            elif tok == ",":
                pos += 1
            else:
                name = tok.split(":")[0].strip()
                if name:
                    node.append(name)
                pos += 1
        return node

    return parse()


def _leaves(node) -> set[str]:
    if isinstance(node, str):
        return {node}
    out: set[str] = set()
    for child in node:
        out |= _leaves(child)
    return out


def clades(newick: str) -> set[frozenset]:
    """Non-trivial clades (bipartitions) of a Newick tree."""
    tree = _parse_newick(newick)
    all_leaves = _leaves(tree)
    result: set[frozenset] = set()

    def walk(node):
        if isinstance(node, str):
            return
        leafset = _leaves(node)
        if 1 < len(leafset) < len(all_leaves):
            result.add(frozenset(leafset))
        for child in node:
            walk(child)

    walk(tree)
    return result


def robinson_foulds(newick_a: str, newick_b: str) -> dict:
    """Robinson-Foulds distance on the shared leaf set of two Newick trees."""
    ca, cb = clades(newick_a), clades(newick_b)
    shared = _leaves_of(ca) & _leaves_of(cb)

    def restrict(cset):
        out = set()
        for c in cset:
            r = frozenset(c & shared)
            if 1 < len(r) < len(shared):
                out.add(r)
        return out

    ra, rb = restrict(ca), restrict(cb)
    rf = len(ra ^ rb)
    max_rf = len(ra) + len(rb)
    return {"rf_distance": rf, "max_rf": max_rf,
            "normalized_rf": (rf / max_rf) if max_rf else 0.0,
            "n_shared_leaves": len(shared)}


def _leaves_of(clade_set: set[frozenset]) -> set[str]:
    out: set[str] = set()
    for c in clade_set:
        out |= set(c)
    return out


# --------------------------------------------------------------------------- #
# Bootstrap support for the functional dendrogram
# --------------------------------------------------------------------------- #
def _linkage_clades(matrix: np.ndarray, labels: list[str], metric: str):
    link = linkage(matrix, method="average", metric=metric)
    n = len(labels)
    result = set()
    # Each internal node's leaf membership via fcluster at successive heights.
    for k in range(2, n):
        assignments = fcluster(link, k, criterion="maxclust")
        for cid in set(assignments):
            members = frozenset(labels[i] for i in range(n) if assignments[i] == cid)
            if 1 < len(members) < n:
                result.add(members)
    return link, result


def bootstrap_support(
    matrix: pd.DataFrame, *, metric: str = "euclidean",
    n_boot: int = 100, seed: int = 0,
) -> pd.DataFrame:
    """Fraction of feature-resampled replicates recovering each observed clade."""
    labels = list(matrix.index)
    data = matrix.to_numpy(dtype=float)
    _, observed = _linkage_clades(data, labels, metric)
    counts = dict.fromkeys(observed, 0)
    rng = np.random.default_rng(seed)
    n_feat = data.shape[1]
    for _ in range(n_boot):
        cols = rng.integers(0, n_feat, n_feat)
        _, rep = _linkage_clades(data[:, cols], labels, metric)
        for clade in observed:
            if clade in rep:
                counts[clade] += 1
    return pd.DataFrame({
        "clade": ["|".join(sorted(c)) for c in counts],
        "size": [len(c) for c in counts],
        "support": [counts[c] / n_boot for c in counts],
    }).sort_values("support", ascending=False).reset_index(drop=True)
