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
    """Benjamini-Hochberg FDR-adjusted q-values.

    Only finite p-values are corrected (and they set the rank denominator);
    non-finite entries stay NaN instead of poisoning every other q-value.
    """
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    n = int(finite.sum())
    if n == 0:
        return q
    pf = p[finite]
    order = np.argsort(pf)
    ranked = pf[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qf = np.empty(n)
    qf[order] = np.clip(ranked, 0, 1)
    q[finite] = qf
    return q


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta effect size in [-1, 1] for two continuous samples.

    Computed from the Mann-Whitney U statistic (``delta = 2U/(n*m) - 1``,
    ties counting as 0.5 in U the same way they count as neither "greater"
    nor "less" here) rather than the full pairwise comparison matrix --
    O(n log n + m log m) via sorting instead of O(n*m) time *and* memory,
    which matters once either group reaches into the thousands (e.g.
    thousands of features/genes rather than the usual handful of organism
    groups). Identical results to the direct pairwise count, ties included.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0 or b.size == 0:
        return float("nan")
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2 * u / (a.size * b.size) - 1


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
    min_known_per_group: int = 1,
) -> pd.DataFrame:
    """Per-feature difference between two groups of organisms.

    Binary (presence / module state): Fisher's exact test, odds ratio and
    prevalence difference. Continuous (completeness): Mann-Whitney U, median
    difference and Cliff's delta. Both add BH-FDR q-values.

    NaN means *unknown*, never absent: it is dropped per feature and per group,
    so prevalences use a known-only denominator and the known/unknown counts are
    reported. A feature with fewer than ``min_known_per_group`` known values in
    either group is not tested (p/q stay NaN) and carries a warning.

    Not corrected for phylogenetic non-independence -- treat as exploratory.
    """
    a_ids = [o for o in matrix.index if groups.get(o) == group_a]
    b_ids = [o for o in matrix.index if groups.get(o) == group_b]
    if not a_ids or not b_ids:
        raise ValueError(f"need members in both groups ({group_a}, {group_b})")

    rows = []
    for feat in matrix.columns:
        a_all = matrix.loc[a_ids, feat].to_numpy(dtype=float)
        b_all = matrix.loc[b_ids, feat].to_numpy(dtype=float)
        a = a_all[np.isfinite(a_all)]
        b = b_all[np.isfinite(b_all)]
        rec = {
            "feature_id": feat,
            "n_a_total": len(a_all), "n_a_known": len(a),
            "n_a_unknown": len(a_all) - len(a),
            "n_b_total": len(b_all), "n_b_known": len(b),
            "n_b_unknown": len(b_all) - len(b),
        }
        testable = len(a) >= min_known_per_group and len(b) >= min_known_per_group
        rec["warning"] = "" if testable else "insufficient_known_values"

        if continuous:
            rec["median_a"] = float(np.median(a)) if len(a) else np.nan
            rec["median_b"] = float(np.median(b)) if len(b) else np.nan
            rec["median_diff"] = rec["median_a"] - rec["median_b"]
            rec["cliffs_delta"] = cliffs_delta(a, b) if testable else np.nan
            rec["test"] = "mann_whitney"
            if testable:
                try:
                    rec["p_value"] = float(stats.mannwhitneyu(
                        a, b, alternative="two-sided").pvalue)
                except ValueError:
                    rec["p_value"] = 1.0
            else:
                rec["p_value"] = np.nan
        else:
            ap = int(np.sum(a > present_threshold))
            bp = int(np.sum(b > present_threshold))
            rec["prevalence_a"] = ap / len(a) if len(a) else np.nan
            rec["prevalence_b"] = bp / len(b) if len(b) else np.nan
            rec["prevalence_diff"] = rec["prevalence_a"] - rec["prevalence_b"]
            rec["test"] = "fisher_exact"
            if testable:
                odds, p = stats.fisher_exact([[ap, len(a) - ap],
                                              [bp, len(b) - bp]])
                rec["odds_ratio"] = float(odds)
                rec["p_value"] = float(p)
            else:
                rec["odds_ratio"] = np.nan
                rec["p_value"] = np.nan
        rows.append(rec)

    out = pd.DataFrame(rows)
    out["q_value"] = benjamini_hochberg(out["p_value"].to_numpy())
    return out.sort_values("p_value", na_position="last").reset_index(drop=True)


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
    min_known_fraction: float = 0.0,
) -> pd.DataFrame:
    """Classify features as core / soft-core / shell / cloud by prevalence.

    Prevalence uses a **known-only denominator**: NaN means unknown, not absent.
    A feature whose known fraction is below ``min_known_fraction`` is labelled
    ``insufficient-data`` rather than assigned a pan class.
    """
    known = matrix.notna()
    n_known = known.sum(axis=0)
    n_present = (matrix > present_threshold).where(known).sum(axis=0)
    n_total = len(matrix.index)
    with np.errstate(invalid="ignore"):
        prevalence = n_present / n_known.replace(0, np.nan)
    known_fraction = n_known / n_total if n_total else n_known * 0.0

    def _cls(p: float, kf: float) -> str:
        if not np.isfinite(p) or kf < min_known_fraction:
            return "insufficient-data"
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
        "n_present": n_present.to_numpy(dtype=int),
        "n_known": n_known.to_numpy(dtype=int),
        "n_unknown": (n_total - n_known).to_numpy(dtype=int),
        "known_fraction": known_fraction.to_numpy(),
        "pan_class": [_cls(p, k) for p, k in zip(prevalence, known_fraction)],
    }).sort_values("prevalence", ascending=False,
                   na_position="last").reset_index(drop=True)


def accumulation_curve(
    matrix: pd.DataFrame, *, permutations: int = 100, seed: int = 0,
    present_threshold: float = 1e-9,
) -> pd.DataFrame:
    """Randomised (permutation / rarefaction-style) pan-core accumulation.

    Organisms are randomly *permuted* (sampled without replacement) and added
    one at a time; this is a rarefaction-style curve, not a with-replacement
    bootstrap. Reports mean, SD and a 95% percentile interval per step. NaN is
    treated as unknown (not present).
    """
    values = matrix.to_numpy(dtype=float)
    present = np.where(np.isfinite(values), values > present_threshold, False)
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
        "pan_sd": total.std(axis=0, ddof=0),
        "pan_ci_lower": np.percentile(total, 2.5, axis=0),
        "pan_ci_upper": np.percentile(total, 97.5, axis=0),
        "core_mean": coref.mean(axis=0),
        "core_sd": coref.std(axis=0, ddof=0),
        "core_ci_lower": np.percentile(coref, 2.5, axis=0),
        "core_ci_upper": np.percentile(coref, 97.5, axis=0),
        "n_replicates": permutations,
        "resampling": "permutation",
    })


# --------------------------------------------------------------------------- #
# Ordination + PERMANOVA
# --------------------------------------------------------------------------- #
def distance_matrix(matrix: pd.DataFrame, metric: str = "braycurtis") -> np.ndarray:
    """Square distance matrix between organisms for the given metric."""
    return squareform(pdist(matrix.to_numpy(dtype=float), metric=metric))


def pcoa(matrix: pd.DataFrame, metric: str = "braycurtis", n_axes: int = 2):
    """Principal coordinates analysis (classical MDS) of the organisms.

    Returns ``(coords DataFrame, explained-variance array, diagnostics dict)``.
    Raises ValueError when no positive axis exists (e.g. every pairwise distance
    is zero) instead of returning an empty result that crashes downstream.
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
    if not pos.any():
        raise ValueError(
            "No positive PCoA axes: all pairwise distances are zero (are all "
            "organisms identical, or is the matrix empty?)")
    coords = vecs[:, pos] * np.sqrt(vals[pos])
    explained = vals[pos] / vals[pos].sum()
    neg = vals < 0
    diagnostics = {
        "n_positive_axes": int(pos.sum()),
        "n_negative_axes": int(neg.sum()),
        "sum_positive_eigenvalues": float(vals[pos].sum()),
        "sum_negative_eigenvalues": float(vals[neg].sum()) if neg.any() else 0.0,
        "negative_fraction": (float(-vals[neg].sum() /
                                    (vals[pos].sum() - vals[neg].sum()))
                              if neg.any() else 0.0),
        "zero_distance_pairs": int(np.isclose(
            d[np.triu_indices(n, k=1)], 0.0, atol=1e-9, rtol=0.0).sum()),
        "metric": metric,
    }
    k = min(n_axes, coords.shape[1])
    cols = [f"PCo{i + 1}" for i in range(k)]
    return (pd.DataFrame(coords[:, :k], index=matrix.index, columns=cols),
            explained[:k], diagnostics)


def permanova(
    matrix: pd.DataFrame, groups: pd.Series, *, metric: str = "braycurtis",
    permutations: int = 999, seed: int = 0,
) -> dict:
    """One-way PERMANOVA (Anderson 2001) pseudo-F and permutation p-value.

    Samples whose group label is missing (None/NaN/empty) are excluded -- a
    missing label is not a group. Requires >=2 groups, each with >=2 samples.
    """
    if permutations < 1:
        raise ValueError("PERMANOVA needs permutations >= 1")
    raw = pd.Series([groups.get(o) for o in matrix.index], index=matrix.index)
    valid = raw.notna() & (raw.astype(str).str.strip() != "")
    if not valid.all():
        matrix = matrix.loc[valid.to_numpy()]
    labels = raw[valid].to_numpy()
    d = distance_matrix(matrix, metric)
    if not np.isfinite(d).all():
        raise ValueError("PERMANOVA requires finite distances (found NaN/Inf)")
    n = len(labels)
    uniq = list(pd.unique(labels))
    k = len(uniq)
    if k < 2 or n <= k:
        raise ValueError("PERMANOVA needs >=2 groups and n > number of groups")
    too_small = [g for g in uniq if int((labels == g).sum()) < 2]
    if too_small:
        raise ValueError(f"PERMANOVA needs >=2 samples per group; too small: "
                         f"{too_small}")
    if np.allclose(d, 0.0, atol=1e-9):
        raise ValueError(
            "PERMANOVA is undefined: every pairwise distance is zero (all "
            "samples are identical under this metric)")

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
            "n_groups": k, "n_samples": n, "permutations": permutations,
            "n_excluded_missing_label": int((~valid).sum()), "metric": metric}


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
