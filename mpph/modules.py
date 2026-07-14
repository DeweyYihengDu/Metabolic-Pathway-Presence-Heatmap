"""KEGG modules: definitions, per-organism KO content, completeness scoring.

The KEGG "module completeness" reported here is the fraction of a module's
reaction *steps* that are covered by an organism's KO repertoire. A module
``DEFINITION`` is a boolean expression over KO ids::

    space = AND (sequential steps)   comma = OR (alternatives)
    +     = essential complex        -     = optional component (ignored)

We split the definition into top-level space-separated steps, evaluate each
against the KO set, and report ``complete_steps / total_steps``. This is the
same approximation used by KEGG-Decoder / MicrobeAnnotator.
"""
from __future__ import annotations

import re
from pathlib import Path

import requests

from .kegg import kegg_get

_KO = re.compile(r"K\d{5}")
_SAFE = re.compile(r"^[01()&| ]*$")


def list_modules(
    session: requests.Session, cache_dir: Path | None, *, refresh: bool = False
) -> dict[str, str]:
    """Return ``{module_id: module_name}`` for every KEGG module."""
    text = kegg_get(session, "list/module", cache_dir, refresh=refresh)
    modules: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            modules[parts[0].strip()] = parts[1].strip()
    return modules


def fetch_module_definition(
    session: requests.Session, module_id: str, cache_dir: Path | None,
    *, refresh: bool = False,
) -> tuple[str, str]:
    """Return ``(definition, category)`` for one module from ``get/md:<id>``.

    ``category`` is the second ``;``-separated field of the CLASS line, e.g.
    ``Carbohydrate metabolism``.
    """
    text = kegg_get(session, f"get/md:{module_id}", cache_dir, refresh=refresh)
    definition, category = "", "Other"
    in_def = False
    for line in text.splitlines():
        if line.startswith("DEFINITION"):
            definition = line[len("DEFINITION"):].strip()
            in_def = True
        elif in_def and line.startswith(" "):
            definition += " " + line.strip()  # continuation line
        else:
            in_def = False
        if line.startswith("CLASS"):
            fields = line[len("CLASS"):].strip().split(";")
            if len(fields) >= 2:
                category = fields[1].strip()
    return definition, category


def _parse_entry(entry: str) -> tuple[str, str, str]:
    """Parse one ``get/md:`` entry -> ``(module_id, definition, category)``."""
    module_id, definition, category = "", "", "Other"
    in_def = False
    for line in entry.splitlines():
        if line.startswith("ENTRY"):
            parts = line.split()
            if len(parts) >= 2:
                module_id = parts[1]
        if line.startswith("DEFINITION"):
            definition = line[len("DEFINITION"):].strip()
            in_def = True
        elif in_def and line.startswith(" ") and not line[:12].strip():
            definition += " " + line.strip()
        else:
            in_def = False
        if line.startswith("CLASS"):
            fields = line[len("CLASS"):].strip().split(";")
            if len(fields) >= 2:
                category = fields[1].strip()
    return module_id, definition, category


def fetch_module_definitions(
    session: requests.Session, module_ids: list[str], cache_dir: Path | None,
    *, refresh: bool = False, progress=None,
) -> dict[str, tuple[str, str]]:
    """Fetch ``{module_id: (definition, category)}`` for many modules.

    KEGG's ``get`` accepts up to 10 entries per call, so definitions are fetched
    in batches of 10 (entries are separated by ``///``) -- ~90 requests for the
    full module set instead of ~900.
    """
    defs: dict[str, tuple[str, str]] = {}
    chunks = [module_ids[i:i + 10] for i in range(0, len(module_ids), 10)]
    for n, chunk in enumerate(chunks, 1):
        endpoint = "get/" + "+".join(f"md:{m}" for m in chunk)
        text = kegg_get(session, endpoint, cache_dir, refresh=refresh)
        for entry in text.split("///"):
            if not entry.strip():
                continue
            mid, definition, category = _parse_entry(entry)
            if mid:
                defs[mid] = (definition, category)
        if progress is not None:
            progress(n, len(chunks))
    return defs


def organism_kos(
    session: requests.Session, org_code: str, cache_dir: Path | None,
    *, refresh: bool = False,
) -> set[str]:
    """Return the set of KO ids encoded by an organism (``link/ko/<org>``)."""
    text = kegg_get(session, f"link/ko/{org_code}", cache_dir, refresh=refresh)
    return set(_KO.findall(text))


# --------------------------------------------------------------------------- #
# Completeness evaluation
# --------------------------------------------------------------------------- #
def split_steps(definition: str) -> list[str]:
    """Split a module definition into top-level (paren-depth 0) steps."""
    steps: list[str] = []
    depth = 0
    cur: list[str] = []
    for ch in definition:
        if ch == "(":
            depth += 1
            cur.append(ch)
        elif ch == ")":
            depth -= 1
            cur.append(ch)
        elif ch == " " and depth == 0:
            if cur:
                steps.append("".join(cur))
                cur = []
        else:
            cur.append(ch)
    if cur:
        steps.append("".join(cur))
    return steps


def _strip_optional(step: str) -> str:
    """Remove optional (``-``) components, which do not affect completeness."""
    prev = None
    while prev != step:
        prev = step
        step = re.sub(r"-K\d{5}", "", step)
        step = re.sub(r"-\([^()]*\)", "", step)
        step = re.sub(r"-M\d{5}", "", step)
    return step


def step_complete(step: str, ko_set: set[str]) -> bool:
    """Evaluate whether a single module step is satisfied by ``ko_set``."""
    step = _strip_optional(step).strip()
    if not step:
        return True
    # space (AND) and '+' (complex, AND) -> '&';  ',' (OR) -> '|'
    expr = step.replace(" ", "&").replace("+", "&").replace(",", "|")
    expr = _KO.sub(lambda m: "1" if m.group(0) in ko_set else "0", expr)
    expr = re.sub(r"M\d{5}", "0", expr)  # nested module refs: treat as absent
    if not _SAFE.match(expr) or not expr:
        return False
    try:
        return bool(eval(expr))  # noqa: S307 - expr is validated to [01()&|]
    except SyntaxError:
        return False


def module_completeness(definition: str, ko_set: set[str]) -> float:
    """Fraction of a module's top-level steps satisfied by ``ko_set`` (0..1)."""
    steps = split_steps(definition)
    if not steps:
        return 0.0
    done = sum(1 for s in steps if step_complete(s, ko_set))
    return done / len(steps)
