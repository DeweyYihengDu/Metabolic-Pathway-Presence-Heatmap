"""Comparative and community analyses over the organism x feature matrix.

All functions take a pandas DataFrame (rows = organisms/samples, columns =
features) plus, where relevant, a per-organism grouping. They are pure and
network-free, so they are unit-tested offline.
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

from .modules import module_completeness


# --------------------------------------------------------------------------- #
# Multiple testing
# --------------------------------------------------------------------------- #
def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted q-values."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta effect size in [-1, 1] for two continuous samples."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0 or b.size == 0:
        return float("nan")
    greater = np.sum(a[:, None] > b[None, :])
    less = np.sum(a[:, None] < b[None, :])
    return (greater - less) / (a.size * b.size)


# --------------------------------------------------------------------------- #
# Differential feature analysis between two groups
# --------------------------------------------------------------------------- #
def differential_features(
    matrix: pd.DataFrame,
    groups: pd.Series,
    group_a: str,
    group_b: str,
    *,
    continuous: bool = False,
    present_threshold: float = 1e-9,
) -> pd.DataFrame:
    """Per-feature difference between two groups of organisms.

    Binary (presence / module state): Fisher's exact test, odds ratio and
    prevalence difference. Continuous (completeness): Mann-Whitney U, median
    difference and Cliff's delta. Both add BH-FDR q-values. Not corrected for
    phylogenetic non-independence -- treat as exploratory.
    """
    a_ids = [o for o in matrix.index if groups.get(o) == group_a]
    b_ids = [o for o in matrix.index if groups.get(o) == group_b]
    if not a_ids or not b_ids:
        raise ValueError(f"need members in both groups ({group_a}, {group_b})")

    rows = []
    for feat in matrix.columns:
        a = matrix.loc[a_ids, feat].to_numpy(dtype=float)
        b = matrix.loc[b_ids, feat].to_numpy(dtype=float)
        rec = {"feature_id": feat, "n_a": len(a), "n_b": len(b)}
        if continuous:
            rec["median_a"] = float(np.median(a))
            rec["median_b"] = float(np.median(b))
            rec["median_diff"] = rec["median_a"] - rec["median_b"]
            rec["cliffs_delta"] = cliffs_delta(a, b)
            try:
                rec["p_value"] = float(stats.mannwhitneyu(
                    a, b, alternative="two-sided").pvalue)
            except ValueError:
                rec["p_value"] = 1.0
            rec["test"] = "mann_whitney"
        else:
            ap = int(np.sum(a > present_threshold))
            bp = int(np.sum(b > present_threshold))
            rec["prevalence_a"] = ap / len(a)
            rec["prevalence_b"] = bp / len(b)
            rec["prevalence_diff"] = rec["prevalence_a"] - rec["prevalence_b"]
            table = [[ap, len(a) - ap], [bp, len(b) - bp]]
            odds, p = stats.fisher_exact(table)
            rec["odds_ratio"] = float(odds)
            rec["p_value"] = float(p)
            rec["test"] = "fisher_exact"
        rows.append(rec)

    out = pd.DataFrame(rows)
    out["q_value"] = benjamini_hochberg(out["p_value"].to_numpy())
    return out.sort_values("p_value").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Pan-functional classification
# --------------------------------------------------------------------------- #
def pan_classify(
    matrix: pd.DataFrame,
    *,
    core: float = 1.0,
    soft_core: float = 0.95,
    shell: float = 0.15,
    present_threshold: float = 1e-9,
) -> pd.DataFrame:
    """Classify features as core / soft-core / shell / cloud by prevalence."""
    prevalence = (matrix > present_threshold).mean(axis=0)

    def _cls(p: float) -> str:
        if p >= core:
            return "core"
        if p >= soft_core:
            return "soft-core"
        if p >= shell:
            return "shell"
        return "cloud"

    return pd.DataFrame({
        "feature_id": prevalence.index,
        "prevalence": prevalence.to_numpy(),
        "pan_class": [_cls(p) for p in prevalence],
    }).sort_values("prevalence", ascending=False).reset_index(drop=True)


def accumulation_curve(
    matrix: pd.DataFrame, *, permutations: int = 100, seed: int = 0,
    present_threshold: float = 1e-9,
) -> pd.DataFrame:
    """Bootstrap pan/core accumulation as organisms are added (mean over perms)."""
    present = (matrix.to_numpy() > present_threshold)
    n_org = present.shape[0]
    rng = np.random.default_rng(seed)
    total = np.zeros((permutations, n_org))
    coref = np.zeros((permutations, n_org))
    for p in range(permutations):
        order = rng.permutation(n_org)
        seen = np.zeros(present.shape[1], dtype=bool)
        allrows = np.ones(present.shape[1], dtype=bool)
        for i, idx in enumerate(order):
            seen |= present[idx]
            allrows &= present[idx]
            total[p, i] = seen.sum()
            coref[p, i] = allrows.sum()
    return pd.DataFrame({
        "n_genomes": np.arange(1, n_org + 1),
        "pan_mean": total.mean(axis=0),
        "core_mean": coref.mean(axis=0),
    })


# --------------------------------------------------------------------------- #
# Ordination + PERMANOVA
# --------------------------------------------------------------------------- #
def distance_matrix(matrix: pd.DataFrame, metric: str = "braycurtis") -> np.ndarray:
    """Square distance matrix between organisms for the given metric."""
    return squareform(pdist(matrix.to_numpy(dtype=float), metric=metric))


def pcoa(matrix: pd.DataFrame, metric: str = "braycurtis", n_axes: int = 2):
    """Principal coordinates analysis (classical MDS) of the organisms.

    Returns ``(coords DataFrame, explained-variance array)``.
    """
    d = distance_matrix(matrix, metric)
    n = d.shape[0]
    a = -0.5 * d ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = j @ a @ j
    vals, vecs = np.linalg.eigh(b)
    idx = np.argsort(vals)[::-1]
    vals, vecs = vals[idx], vecs[:, idx]
    pos = vals > 0
    coords = vecs[:, pos] * np.sqrt(vals[pos])
    explained = vals[pos] / vals[pos].sum()
    k = min(n_axes, coords.shape[1])
    cols = [f"PCo{i + 1}" for i in range(k)]
    return (pd.DataFrame(coords[:, :k], index=matrix.index, columns=cols),
            explained[:k])


def permanova(
    matrix: pd.DataFrame, groups: pd.Series, *, metric: str = "braycurtis",
    permutations: int = 999, seed: int = 0,
) -> dict:
    """One-way PERMANOVA (Anderson 2001) pseudo-F and permutation p-value."""
    labels = np.array([groups.get(o) for o in matrix.index])
    d = distance_matrix(matrix, metric)
    n = len(labels)
    uniq = [g for g in pd.unique(labels) if g is not None]
    k = len(uniq)
    if k < 2 or n <= k:
        raise ValueError("PERMANOVA needs >=2 groups and n > number of groups")

    def pseudo_f(lab: np.ndarray) -> float:
        ss_total = (d ** 2).sum() / (2 * n)
        ss_within = 0.0
        for g in uniq:
            members = np.where(lab == g)[0]
            ng = len(members)
            if ng > 1:
                sub = d[np.ix_(members, members)]
                ss_within += (sub ** 2).sum() / (2 * ng)
        ss_among = ss_total - ss_within
        return (ss_among / (k - 1)) / (ss_within / (n - k))

    observed = pseudo_f(labels)
    rng = np.random.default_rng(seed)
    count = 1
    for _ in range(permutations):
        if pseudo_f(rng.permutation(labels)) >= observed:
            count += 1
    return {"pseudo_F": float(observed), "p_value": count / (permutations + 1),
            "n_groups": k, "n_samples": n, "permutations": permutations}


# --------------------------------------------------------------------------- #
# Community metabolic complementarity
# --------------------------------------------------------------------------- #
def pairwise_complementarity(
    org_kos: dict[str, set[str]],
    module_defs: dict[str, tuple[str, str, str]],
    *, max_combination_size: int = 2, only_newly_completed: bool = True,
) -> pd.DataFrame:
    """Find module completions enabled by combining organisms' KO sets.

    For each combination (up to ``max_combination_size``) and each module,
    compare the best single-member completeness with the combined-KO
    completeness. By default reports only modules newly completed by the union.
    Potential complementarity only -- not evidence of metabolite exchange.
    """
    names = list(org_kos)
    rows = []
    for size in range(2, max_combination_size + 1):
        for combo in itertools.combinations(names, size):
            union: set[str] = set().union(*(org_kos[c] for c in combo))
            for mid, (definition, _cat, _t) in module_defs.items():
                indiv = max(module_completeness(definition, org_kos[c], module_defs)
                            for c in combo)
                combined = module_completeness(definition, union, module_defs)
                newly = combined >= 1.0 > indiv
                if only_newly_completed and not newly:
                    continue
                rows.append({
                    "members": "+".join(combo),
                    "module_id": mid,
                    "individual_max": indiv,
                    "combined_score": combined,
                    "newly_completed": newly,
                })
    return pd.DataFrame(rows, columns=["members", "module_id", "individual_max",
                                       "combined_score", "newly_completed"])
