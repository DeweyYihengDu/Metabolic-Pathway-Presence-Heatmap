"""Pathway lists per organism and the KEGG functional-category hierarchy."""
from __future__ import annotations

from pathlib import Path

import requests

from .kegg import kegg_get

OVERVIEW_CATEGORY = "Global and overview maps"


def get_pathways(
    session: requests.Session, org_code: str, cache_dir: Path | None,
    *, refresh: bool = False,
) -> dict[str, str]:
    """Return ``{pathway_number: pathway_name}`` for one organism.

    A pathway id such as ``vch01100`` is reduced to its 5-digit KEGG map number
    (``01100``) so the same pathway is comparable across organisms.
    """
    text = kegg_get(session, f"list/pathway/{org_code}", cache_dir, refresh=refresh)
    pathways: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        pathway_number = parts[0][-5:]  # trailing 5 digits are the map id
        pathway_name = parts[1].split(" - ")[0].strip()
        pathways[pathway_number] = pathway_name
    return pathways


def fetch_pathway_categories(
    session: requests.Session, cache_dir: Path | None, *, refresh: bool = False
) -> dict[str, str]:
    """Return ``{map_id: functional_category}`` from KEGG BRITE ``br08901``.

    The category is the second-level ("B") heading of the KEGG pathway
    hierarchy, e.g. ``Carbohydrate metabolism`` or ``Energy metabolism`` -- the
    granularity that is most informative for a microbial comparison.
    """
    text = kegg_get(session, "get/br:br08901", cache_dir, refresh=refresh)
    categories: dict[str, str] = {}
    current_b: str | None = None
    for line in text.splitlines():
        if not line or line[0] in "!+#":
            continue
        level, body = line[0], line[1:].strip()
        if level == "A":
            current_b = None
        elif level == "B":
            current_b = body
        elif level == "C":
            parts = body.split(None, 1)
            if parts and parts[0].isdigit():
                categories[parts[0].zfill(5)] = current_b or "Other"
    return categories
