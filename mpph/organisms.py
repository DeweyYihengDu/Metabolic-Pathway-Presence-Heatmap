"""Selecting KEGG organisms by taxon name or an explicit code list."""
from __future__ import annotations

import re
from pathlib import Path

import requests

from .kegg import kegg_get


def list_genomes(
    session: requests.Session, cache_dir: Path | None, *, refresh: bool = False
) -> list[tuple[str, str]]:
    """Return ``[(org_code, organism_name), ...]`` for every KEGG genome.

    Uses the ``list/genome`` endpoint (the old ``list/organism`` endpoint was
    retired by KEGG and now returns HTTP 400). Each line looks like::

        T00034<TAB>vch; Vibrio cholerae O1 El Tor N16961
    """
    text = kegg_get(session, "list/genome", cache_dir, refresh=refresh)
    genomes: list[tuple[str, str]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2 or "; " not in parts[1]:
            continue
        code, name = parts[1].split("; ", 1)
        genomes.append((code.strip(), name.strip()))
    return genomes


def select_by_taxon(
    genomes: list[tuple[str, str]], taxon: str, exact: bool
) -> list[tuple[str, str]]:
    """Filter genomes whose name matches ``taxon``.

    ``exact=False`` (default) matches ``taxon`` as a whole word, so ``Vibrio``
    no longer accidentally captures ``Vibrionimonas``. ``exact=True`` requires
    the name to *start* with ``taxon`` (useful for a full "Genus species").
    """
    if exact:
        prefix = taxon.lower()
        return [(c, n) for c, n in genomes if n.lower().startswith(prefix)]

    pattern = re.compile(rf"\b{re.escape(taxon)}\b", re.IGNORECASE)
    return [(c, n) for c, n in genomes if pattern.search(n)]


def select_by_codes(
    genomes: list[tuple[str, str]], codes: list[str]
) -> list[tuple[str, str]]:
    """Return genomes whose organism code is in ``codes`` (order preserved)."""
    by_code = {c: n for c, n in genomes}
    wanted = [c.strip() for c in codes if c.strip()]
    return [(c, by_code.get(c, c)) for c in wanted]


def read_code_file(path: str | Path) -> list[str]:
    """Read organism codes from a file (one per line; ``#`` comments allowed)."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    return [ln.split("#", 1)[0].strip() for ln in lines if ln.split("#", 1)[0].strip()]
