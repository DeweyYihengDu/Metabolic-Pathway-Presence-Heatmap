"""Rank-based gene set enrichment (GSEA / "GSEAPreranked" style).

Unlike `mpph.enrichment` (hypergeometric ORA on a discrete study set), this
scores a **ranked list of every gene** -- typically ranked by a differential
expression statistic -- and tests whether a category's members cluster toward
either end of the ranking, using the weighted running-sum statistic from
Subramanian et al. 2005 (PNAS). More precisely, this is **gene-set permutation
on a fixed ("preranked") ranking**: for each distinct gene-set size, many
random gene sets of that size are drawn from the *same* ranking to build a
null distribution, giving a nominal p-value per category (BH-FDR corrected
across categories tested). This is weaker than the original GSEA's default
*phenotype* permutation (which re-derives the ranking from the raw samples on
every permutation and so also captures gene-gene correlation structure) --
phenotype permutation is not implemented here.

Note the deliberate difference from `mpph.enrichment`: ORA restricts study and
background to genes annotated in the category system before testing. GSEA does
**not** drop unannotated genes from the ranked list -- every gene contributes to
the running sum (as a "miss"), and dropping them would inflate the statistic by
removing genuine background noise.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .analysis import benjamini_hochberg

_EPS = 1e-9


# --------------------------------------------------------------------------- #
# Input: a pre-ranked list, or an expression matrix + two groups
# --------------------------------------------------------------------------- #
def load_ranked_list(path: str | Path) -> pd.Series:
    """Load a ``gene_id<TAB>score`` ranked list, sorted descending by score."""
    rows: list[tuple[str, float]] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        try:
            rows.append((fields[0].strip(), float(fields[1])))
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"no gene/score pairs found in {path}")
    # Check for duplicates on the raw pairs -- dict(rows) would silently keep
    # only the last value for a repeated key, hiding the very thing we're
    # checking for.
    ids = [gene for gene, _score in rows]
    if len(set(ids)) != len(ids):
        seen, dups = set(), []
        for gene in ids:
            if gene in seen:
                dups.append(gene)
            seen.add(gene)
        raise ValueError(f"duplicate gene ids in ranked list: {sorted(set(dups))[:5]}")
    # float() accepts "nan"/"inf"/"-inf" without raising, so these would
    # otherwise sort to one end of the ranking and dominate the running-sum
    # statistic. A common source is a DE tool (e.g. DESeq2) reporting Inf/-Inf
    # log-fold-change or NA for genes with zero counts in one group.
    non_finite = sorted(gene for gene, score in rows if not np.isfinite(score))
    if non_finite:
        raise ValueError(
            f"ranked list has non-finite score(s) (NaN/Inf) for: "
            f"{non_finite[:5]}{' ...' if len(non_finite) > 5 else ''}. Drop or "
            "re-score those genes (e.g. genes with zero counts in one group) "
            "before ranking.")
    series = pd.Series(dict(rows), name="score")
    return series.sort_values(ascending=False)


def load_expression_matrix(path: str | Path) -> pd.DataFrame:
    """Load a ``gene_id<TAB>sample1<TAB>sample2...`` expression table."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    if df.index.duplicated().any():
        dups = df.index[df.index.duplicated()].unique().tolist()
        raise ValueError(f"duplicate gene ids in expression matrix: {dups[:5]}")
    return df.apply(pd.to_numeric, errors="coerce")


def rank_from_expression(
    expr: pd.DataFrame, groups: pd.Series, group_a: str, group_b: str,
    metric: str = "signal2noise",
) -> pd.Series:
    """Rank genes by a differential statistic between two sample groups.

    ``metric='signal2noise'`` is the statistic used in the original GSEA paper:
    ``(mean_a - mean_b) / (std_a + std_b)``. ``metric='log2fc'`` is
    ``log2((mean_a + eps) / (mean_b + eps))``. Both are simple, dependency-free
    and make no distributional assumptions -- this is *not* a substitute for a
    dedicated DE tool (DESeq2/edgeR/limma); rank by that tool's own statistic
    via `load_ranked_list` for a rigorous analysis.
    """
    a_cols = [c for c in expr.columns if groups.get(c) == group_a]
    b_cols = [c for c in expr.columns if groups.get(c) == group_b]
    if len(a_cols) < 2 or len(b_cols) < 2:
        raise ValueError(f"need >=2 samples per group; got {len(a_cols)} in "
                         f"{group_a!r}, {len(b_cols)} in {group_b!r}")
    mean_a, mean_b = expr[a_cols].mean(axis=1), expr[b_cols].mean(axis=1)
    if metric == "log2fc":
        score = np.log2((mean_a + _EPS) / (mean_b + _EPS))
    elif metric == "signal2noise":
        std_a = expr[a_cols].std(axis=1, ddof=1).clip(lower=_EPS)
        std_b = expr[b_cols].std(axis=1, ddof=1).clip(lower=_EPS)
        score = (mean_a - mean_b) / (std_a + std_b)
    else:
        raise ValueError(f"unknown rank metric: {metric}")
    ranked = score.dropna().sort_values(ascending=False)
    # dropna() doesn't catch +-Inf (e.g. from Inf/-Inf already present in the
    # input expression matrix, surviving the epsilon padding above).
    non_finite = sorted(ranked.index[~np.isfinite(ranked.to_numpy())])
    if non_finite:
        raise ValueError(
            f"non-finite score(s) (Inf) computed for: {non_finite[:5]}"
            f"{' ...' if len(non_finite) > 5 else ''} -- check for Inf/-Inf "
            "values in the input expression matrix.")
    return ranked


# --------------------------------------------------------------------------- #
# The weighted running-sum enrichment statistic
# --------------------------------------------------------------------------- #
def _running_sum(hits: np.ndarray, abs_scores: np.ndarray, weight: float) -> np.ndarray:
    """Weighted KS running-sum trajectory for a boolean hit mask."""
    n = hits.size
    n_hit = int(hits.sum())
    n_miss = n - n_hit
    if n_hit == 0 or n_miss == 0:
        return np.zeros(n)
    if weight == 0:
        hit_step = np.where(hits, 1.0 / n_hit, 0.0)
    else:
        w = np.where(hits, abs_scores ** weight, 0.0)
        total = w.sum()
        hit_step = w / total if total > 0 else np.where(hits, 1.0 / n_hit, 0.0)
    miss_step = np.where(hits, 0.0, 1.0 / n_miss)
    return np.cumsum(hit_step - miss_step)


def enrichment_score(
    ranked_genes: list[str], ranked_scores: np.ndarray, gene_set: set[str],
    weight: float = 1.0,
) -> tuple[float, np.ndarray]:
    """Weighted running-sum ES for one gene set against a fixed ranking.

    Returns ``(ES, running_sum)`` where ES is the running sum's maximum-absolute
    deviation from zero (signed) and ``running_sum`` is the full trajectory
    (for plotting / leading-edge extraction).
    """
    hits = np.fromiter((g in gene_set for g in ranked_genes), dtype=bool,
                       count=len(ranked_genes))
    running = _running_sum(hits, np.abs(ranked_scores), weight)
    if not running.any():
        return 0.0, running
    i_max, i_min = int(np.argmax(running)), int(np.argmin(running))
    es_pos, es_neg = running[i_max], running[i_min]
    es = es_pos if abs(es_pos) >= abs(es_neg) else es_neg
    return float(es), running


def leading_edge(ranked_genes: list[str], hits: np.ndarray, running: np.ndarray,
                 es: float) -> list[str]:
    """Genes driving the enrichment: hits up to (ES>0) or from (ES<0) the peak."""
    if es >= 0:
        peak = int(np.argmax(running))
        idx = [i for i in range(peak + 1) if hits[i]]
    else:
        peak = int(np.argmin(running))
        idx = [i for i in range(peak, len(running)) if hits[i]]
    return [ranked_genes[i] for i in idx]


# --------------------------------------------------------------------------- #
# Full analysis: score every category + permutation null + FDR
# --------------------------------------------------------------------------- #
def gsea_analysis(
    ranked: pd.Series,
    category_to_items: dict[str, set[str]],
    category_names: dict[str, str] | None = None,
    *,
    weight: float = 1.0,
    min_size: int = 15,
    max_size: int = 500,
    permutations: int = 1000,
    seed: int = 0,
) -> tuple[pd.DataFrame, dict[str, np.ndarray], list[str]]:
    """Score every category against ``ranked`` and estimate significance.

    Returns ``(results, running_sums, ranked_gene_order)``. ``running_sums``
    maps category id to its running-sum trajectory (for the GSEA-style plot).
    Every gene in ``ranked`` counts toward the statistic (unannotated genes are
    kept as "misses", unlike ORA) -- see the module docstring.
    """
    category_names = category_names or {}
    ranked_genes = list(ranked.index)
    ranked_scores = ranked.to_numpy(dtype=float)
    ranked_set_index = {g: i for i, g in enumerate(ranked_genes)}
    n = len(ranked_genes)
    abs_scores = np.abs(ranked_scores)

    # Restrict each category to genes actually present in the ranking, and
    # drop those outside [min_size, max_size].
    sized: dict[str, set[str]] = {}
    for cat, members in category_to_items.items():
        present = members & ranked_set_index.keys()
        if min_size <= len(present) <= max_size:
            sized[cat] = present

    real: dict[str, tuple[float, np.ndarray, int]] = {}
    for cat, members in sized.items():
        hits = np.fromiter((g in members for g in ranked_genes), dtype=bool, count=n)
        running = _running_sum(hits, abs_scores, weight)
        if not running.any():
            continue
        i_max, i_min = int(np.argmax(running)), int(np.argmin(running))
        es = running[i_max] if abs(running[i_max]) >= abs(running[i_min]) else running[i_min]
        real[cat] = (float(es), running, len(members))

    # Null distribution of ES, batched per distinct size (ES under the null
    # depends only on the *size* of a random gene set, not its identity, for a
    # fixed ranking) -- this reuse is what keeps permutation tractable.
    rng = np.random.default_rng(seed)
    sizes_needed = sorted({size for _es, _run, size in real.values()})
    null_by_size: dict[int, np.ndarray] = {}
    for size in sizes_needed:
        null_by_size[size] = _null_es_batch(abs_scores, n, size, weight, permutations, rng)

    rows = []
    running_sums: dict[str, np.ndarray] = {}
    for cat, (es, running, size) in real.items():
        null = null_by_size[size]
        same_sign = null[null >= 0] if es >= 0 else null[null < 0]
        if same_sign.size == 0:
            nes, p = np.nan, 1.0
        else:
            nes = es / (np.mean(np.abs(same_sign)) + _EPS)
            more_extreme = int(np.sum(np.abs(same_sign) >= abs(es)))
            p = (more_extreme + 1) / (same_sign.size + 1)
        hits = np.fromiter((g in sized[cat] for g in ranked_genes), dtype=bool, count=n)
        edge = leading_edge(ranked_genes, hits, running, es)
        running_sums[cat] = running
        rows.append({
            "category_id": cat,
            "category_name": category_names.get(cat, cat),
            "size": size,
            "ES": es,
            "NES": nes,
            "p_value": float(p),
            "leading_edge_size": len(edge),
            "leading_edge_genes": ",".join(edge),
        })

    out = pd.DataFrame(rows, columns=[
        "category_id", "category_name", "size", "ES", "NES", "p_value",
        "leading_edge_size", "leading_edge_genes",
    ])
    if len(out):
        out["q_value"] = benjamini_hochberg(out["p_value"].to_numpy())
        out = out.sort_values("p_value").reset_index(drop=True)
    else:
        out["q_value"] = pd.Series(dtype=float)
    return out, running_sums, ranked_genes


_NULL_BATCH_SIZE = 200


def _null_es_batch(
    abs_scores: np.ndarray, n: int, size: int, weight: float,
    permutations: int, rng: np.random.Generator,
    *, batch_size: int = _NULL_BATCH_SIZE,
) -> np.ndarray:
    """Vectorised null ES for `permutations` random gene sets of `size`.

    Processed in batches of `batch_size` permutations rather than allocating
    one (permutations, n) array up front: for a large gene universe (tens of
    thousands of genes) and thousands of permutations, the handful of
    same-shaped temporaries this needs (random keys, sort order, hit mask,
    running sum, ...) can reach multiple GB at once. Batching bounds peak
    memory to O(batch_size * n) regardless of the total permutation count,
    with identical results (the same random draws, just requested from the
    generator in smaller chunks -- numpy's `Generator` is stream-based, so
    this doesn't change what's drawn).
    """
    results = np.empty(permutations, dtype=float)
    n_miss = n - size
    for start in range(0, permutations, batch_size):
        b = min(batch_size, permutations - start)
        # One random subset of `size` positions per row.
        rand_keys = rng.random((b, n))
        order = np.argsort(rand_keys, axis=1)
        hit_idx = order[:, :size]  # (b, size) -- indices chosen as "hits"
        hits = np.zeros((b, n), dtype=bool)
        np.put_along_axis(hits, hit_idx, True, axis=1)

        miss_step = np.where(hits, 0.0, 1.0 / n_miss)
        if weight == 0:
            hit_step = np.where(hits, 1.0 / size, 0.0)
        else:
            w = np.where(hits, abs_scores[np.newaxis, :] ** weight, 0.0)
            totals = w.sum(axis=1, keepdims=True)
            totals[totals == 0] = 1.0
            hit_step = w / totals
        running = np.cumsum(hit_step - miss_step, axis=1)
        es_max = running.max(axis=1)
        es_min = running.min(axis=1)
        results[start:start + b] = np.where(
            np.abs(es_max) >= np.abs(es_min), es_max, es_min)
    return results
