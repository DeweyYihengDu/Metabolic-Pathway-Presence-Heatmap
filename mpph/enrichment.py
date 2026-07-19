"""Over-representation (enrichment) analysis for KEGG pathways/modules and GO terms.

Given a *study set* of genes/KOs (e.g. genes unique to a genome, a species'
accessory genome, or an externally supplied DEG list) and a *background*
(the annotated universe it was drawn from), test each category (KEGG pathway,
KEGG module, or a GO term) for over-representation with the hypergeometric
test, one-sided (study set enriched *for* the category, not depleted).

KEGG pathway/module membership comes from the KEGG REST API (public, global
KO<->pathway and KO<->module links -- not organism-specific). GO enrichment
needs a user-supplied gene-to-GO-term mapping: **KEGG does not provide GO
annotations**, so this is bring-your-own-mapping (a GAF-style long table, or
the `GO_terms` column of an eggNOG-mapper ``.annotations`` file). This is
**flat** GO term over-representation -- a gene counts only toward the terms
your mapping lists for it, not their GO-DAG ancestors, unless the mapping
already includes those ancestors explicitly.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import hypergeom

from .analysis import benjamini_hochberg
from .kegg import kegg_get

_GO_ID = re.compile(r"GO:\d{7}")
_KO_PREFIXED = re.compile(r"^ko:(K\d{5})$")


def _normalize_ko_ids(ids: set[str]) -> set[str]:
    """Strip a KEGG ``ko:`` namespace prefix (``ko:K00001`` -> ``K00001``).

    KEGG's own KO-membership endpoints here (and this module's fetchers)
    use bare ``K#####`` ids, but a study/background list copied straight
    from ``link/ko/<org>`` or an eggNOG-mapper ``KEGG_ko`` column keeps the
    ``ko:`` prefix -- which would otherwise silently fail to match any
    category (a plain set difference, not an error) rather than raising.
    Only exact ``ko:K#####`` tokens are touched, so GO-mode gene ids (an
    arbitrary user-supplied namespace) are never affected.
    """
    normalized = set()
    for item in ids:
        match = _KO_PREFIXED.match(item)
        normalized.add(match.group(1) if match else item)
    return normalized


def _odds_ratio_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """2x2 odds ratio + 95% CI (log-scale normal approximation).

    Haldane-Anscombe correction (add 0.5 to every cell) when any cell is
    zero -- standard practice to avoid an undefined or infinite ratio for a
    small/sparse table, at the cost of a slightly conservative estimate.
    """
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    odds_ratio = (a * d) / (b * c)
    log_or = math.log(odds_ratio)
    se_log_or = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return odds_ratio, math.exp(log_or - 1.96 * se_log_or), math.exp(
        log_or + 1.96 * se_log_or)


# --------------------------------------------------------------------------- #
# KEGG category membership (global, not organism-specific)
# --------------------------------------------------------------------------- #
def fetch_ko_pathway_membership(
    session: requests.Session, cache_dir: Path | None, *, refresh: bool = False
) -> dict[str, set[str]]:
    """Return ``{KO: {map_id, ...}}`` from the global ``link/pathway/ko``.

    Each KO links to both a ``map#####`` (reference pathway) and a
    ``ko#####`` (KO-view of the same pathway) entry; both collapse to the same
    5-digit map id.
    """
    text = kegg_get(session, "link/pathway/ko", cache_dir, refresh=refresh)
    membership: dict[str, set[str]] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        ko = parts[0].split(":", 1)[-1]
        map_id = parts[1].split(":", 1)[-1][-5:]
        membership.setdefault(ko, set()).add(map_id)
    return membership


def fetch_ko_module_membership(
    session: requests.Session, cache_dir: Path | None, *, refresh: bool = False
) -> dict[str, set[str]]:
    """Return ``{KO: {module_id, ...}}`` from the global ``link/module/ko``."""
    text = kegg_get(session, "link/module/ko", cache_dir, refresh=refresh)
    membership: dict[str, set[str]] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        ko = parts[0].split(":", 1)[-1]
        mid = parts[1].split(":", 1)[-1]
        membership.setdefault(ko, set()).add(mid)
    return membership


def fetch_pathway_names(
    session: requests.Session, cache_dir: Path | None, *, refresh: bool = False
) -> dict[str, str]:
    """Return ``{map_id: name}`` from the global (non-organism) ``list/pathway``."""
    text = kegg_get(session, "list/pathway", cache_dir, refresh=refresh)
    names: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        map_id = parts[0].split(":", 1)[-1][-5:]
        names[map_id] = parts[1].strip()
    return names


def invert_membership(item_to_categories: dict[str, set[str]]) -> dict[str, set[str]]:
    """Invert ``{item: {category, ...}}`` to ``{category: {item, ...}}``."""
    inverted: dict[str, set[str]] = {}
    for item, cats in item_to_categories.items():
        for cat in cats:
            inverted.setdefault(cat, set()).add(item)
    return inverted


# --------------------------------------------------------------------------- #
# Gene -> GO mapping (bring-your-own; KEGG has no GO annotations)
# --------------------------------------------------------------------------- #
def load_gene_go_map(path: str | Path, fmt: str = "auto") -> dict[str, set[str]]:
    """Load ``{gene_id: {GO:#######, ...}}`` from a file.

    ``fmt``:

    * ``auto`` (default) -- a long table (``gene<TAB>GO:...``, one row per pair
      or per gene with a ``,``/``;``-separated GO list in the second column);
    * ``eggnog`` -- an eggNOG-mapper ``.annotations`` file's ``GO_terms`` column.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"gene-GO map not found: {p}")
    text = p.read_text(encoding="utf-8", errors="ignore")

    if fmt == "eggnog":
        return _parse_eggnog_go(text)

    mapping: dict[str, set[str]] = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        gene = fields[0].strip()
        go_ids = set(_GO_ID.findall(fields[1]))
        if gene and go_ids:
            mapping.setdefault(gene, set()).update(go_ids)
    if not mapping:
        raise ValueError(f"no gene/GO pairs found in {p}")
    return mapping


def _parse_eggnog_go(text: str) -> dict[str, set[str]]:
    header: list[str] | None = None
    idx = None
    mapping: dict[str, set[str]] = {}
    for line in text.splitlines():
        if line.startswith("##") or not line.strip():
            continue
        fields = line.lstrip("#").split("\t")
        if header is None:
            if "GO_terms" in fields:
                header = fields
                idx = header.index("GO_terms")
            continue
        gene = fields[0]
        if idx is not None and idx < len(fields):
            go_ids = set(_GO_ID.findall(fields[idx]))
            if go_ids:
                mapping.setdefault(gene, set()).update(go_ids)
    if not mapping:
        raise ValueError("no GO_terms column with GO ids found (need an "
                         "eggNOG-mapper .annotations file)")
    return mapping


def load_go_names(path: str | Path) -> dict[str, str]:
    """Load an optional ``go_id<TAB>name`` mapping for readable category names."""
    p = Path(path)
    names: dict[str, str] = {}
    for line in p.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) >= 2 and _GO_ID.fullmatch(fields[0].strip()):
            names[fields[0].strip()] = fields[1].strip()
    return names


def read_id_list(path: str | Path) -> set[str]:
    """Read a plain gene/KO id list, one per line (``#`` comments allowed)."""
    ids = set()
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        item = line.split("#", 1)[0].strip()
        if item:
            ids.add(item)
    return ids


# --------------------------------------------------------------------------- #
# Hypergeometric over-representation test
# --------------------------------------------------------------------------- #
def hypergeometric_enrichment(
    study: set[str],
    background: set[str],
    category_to_items: dict[str, set[str]],
    category_names: dict[str, str] | None = None,
    *,
    min_category_size: int = 2,
) -> tuple[pd.DataFrame, dict]:
    """One-sided hypergeometric over-representation test per category.

    Both ``study`` and ``background`` are restricted to genes/KOs that actually
    appear in ``category_to_items`` (the annotated universe) before testing --
    an unannotated gene cannot support or refute enrichment for any category.
    Returns ``(results, stats)`` where ``stats`` records the study/background
    sizes actually used and how many inputs were dropped as unannotated.

    A study gene not present in ``background`` is *also* dropped (a valid
    background must contain the study set) and counted in
    ``n_study_not_in_background``. A ``ko:`` namespace prefix on either input
    is normalized away (see ``_normalize_ko_ids``) before any of this, so it
    never causes a silent all-unannotated result.
    """
    study = _normalize_ko_ids(study)
    background = _normalize_ko_ids(background)
    category_names = category_names or {}
    all_annotated_genes = (set().union(*category_to_items.values())
                          if category_to_items else set())

    not_in_background = study - background
    study_in_bg = study & background
    background_used = background & all_annotated_genes
    study_used = study_in_bg & all_annotated_genes

    stats = {
        "n_study_input": len(study),
        "n_study_not_in_background": len(not_in_background),
        "n_study_used": len(study_used),
        "n_background_input": len(background),
        "n_background_used": len(background_used),
        "n_categories_tested": 0,
    }

    n_bg = len(background_used)
    n_study = len(study_used)
    rows = []
    if n_bg > 0:
        for cat_id, members in category_to_items.items():
            bg_hits = members & background_used
            k_bg = len(bg_hits)
            if k_bg < min_category_size:
                continue
            study_hits = members & study_used
            k_study = len(study_hits)
            fold = ((k_study / n_study) / (k_bg / n_bg)
                    if n_study > 0 and k_bg > 0 else np.nan)
            # P(X >= k_study) with X ~ Hypergeom(N=n_bg, K=k_bg, n=n_study)
            p = hypergeom.sf(k_study - 1, n_bg, k_bg, n_study) if n_study > 0 else 1.0
            # The same 2x2 table as the hypergeometric test itself, with
            # background_used as the universe (study is a subset of it):
            # a=k_study (study & category), b=study \ category,
            # c=(background \ study) & category, d=(background \ study) \ category.
            a, b = k_study, n_study - k_study
            c, d = k_bg - k_study, (n_bg - k_bg) - b
            odds_ratio, or_lo, or_hi = _odds_ratio_ci(a, b, c, d)
            rows.append({
                "category_id": cat_id,
                "category_name": category_names.get(cat_id, cat_id),
                "k_study_hits": k_study,
                "n_study_total": n_study,
                "K_background_hits": k_bg,
                "N_background_total": n_bg,
                "gene_ratio": f"{k_study}/{n_study}" if n_study else "0/0",
                "bg_ratio": f"{k_bg}/{n_bg}",
                "fold_enrichment": fold,
                "odds_ratio": odds_ratio,
                "odds_ratio_ci_low": or_lo,
                "odds_ratio_ci_high": or_hi,
                "p_value": float(p),
                "study_items": ",".join(sorted(study_hits)),
            })
    stats["n_categories_tested"] = len(rows)

    out = pd.DataFrame(rows, columns=[
        "category_id", "category_name", "k_study_hits", "n_study_total",
        "K_background_hits", "N_background_total", "gene_ratio", "bg_ratio",
        "fold_enrichment", "odds_ratio", "odds_ratio_ci_low",
        "odds_ratio_ci_high", "p_value", "study_items",
    ])
    if len(out):
        out["q_value"] = benjamini_hochberg(out["p_value"].to_numpy())
        out = out.sort_values("p_value").reset_index(drop=True)
    else:
        out["q_value"] = pd.Series(dtype=float)
    return out, stats
