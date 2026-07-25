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

import math
import re
from dataclasses import dataclass, field
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


def _parse_entry(entry: str) -> tuple[str, str, str, str]:
    """Parse one ``get/md:`` entry -> ``(module_id, definition, category, type)``.

    ``type`` is the ENTRY line's module type (e.g. ``Pathway``, ``Signature``);
    ``category`` is the second ``;``-field of CLASS (e.g. ``Carbohydrate
    metabolism``).
    """
    module_id, definition, category, mtype = "", "", "Other", "Unknown"
    in_def = False
    for line in entry.splitlines():
        if line.startswith("ENTRY"):
            parts = line.split()
            if len(parts) >= 2:
                module_id = parts[1]
            if len(parts) >= 3:
                mtype = parts[2]  # e.g. "Pathway", "Signature", "Reaction"
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
    return module_id, definition, category, mtype


def fetch_module_definitions(
    session: requests.Session, module_ids: list[str], cache_dir: Path | None,
    *, refresh: bool = False, progress=None,
) -> dict[str, tuple[str, str, str]]:
    """Fetch ``{module_id: (definition, category, type)}`` for many modules.

    KEGG's ``get`` accepts up to 10 entries per call, so definitions are fetched
    in batches of 10 (entries are separated by ``///``) -- ~90 requests for the
    full module set instead of ~900.
    """
    defs: dict[str, tuple[str, str, str]] = {}
    chunks = [module_ids[i:i + 10] for i in range(0, len(module_ids), 10)]
    for n, chunk in enumerate(chunks, 1):
        endpoint = "get/" + "+".join(f"md:{m}" for m in chunk)
        text = kegg_get(session, endpoint, cache_dir, refresh=refresh)
        for entry in text.split("///"):
            if not entry.strip():
                continue
            mid, definition, category, mtype = _parse_entry(entry)
            if mid:
                defs[mid] = (definition, category, mtype)
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


def _is_gap(step: str) -> bool:
    """True for a KEGG placeholder step ('--', a reaction with no assigned KO)."""
    s = step.strip()
    return bool(s) and set(s) <= {"-"}


def step_complete(
    step: str, ko_set: set[str],
    module_defs: dict[str, tuple[str, str, str]] | None = None,
    _stack: frozenset[str] = frozenset(),
) -> bool | None:
    """Evaluate whether a single module step is satisfied by ``ko_set``.

    A nested module reference (``M#####``) is resolved recursively: it counts as
    satisfied only if that referenced module is itself fully complete for
    ``ko_set``. ``_stack`` guards against cyclic definitions.

    Returns ``None`` ("unknown") rather than ``False`` when the step's truth
    value genuinely can't be determined -- a nested reference we don't have a
    definition for, a cyclic reference, or syntax this parser doesn't handle --
    instead of silently conflating "can't tell" with "confirmed absent", which
    would understate completeness. This is decided by Kleene/three-valued
    evaluation: every undetermined token is substituted with both 0 and 1: if
    the step's truth value doesn't change either way, it doesn't actually
    depend on that token (e.g. an OR already satisfied by a known KO), so the
    determined result is returned; only a genuine disagreement is "unknown".
    """
    step = _strip_optional(step).strip()
    if not step or _is_gap(step):
        return True

    def _module_value(match) -> str:
        mid = match.group(0)
        if module_defs and mid in module_defs and mid not in _stack:
            sub = module_completeness(module_defs[mid][0], ko_set, module_defs,
                                      _stack | {mid})
            if not math.isnan(sub):  # not NaN: the nested module IS determined
                return "1" if sub >= 1.0 else "0"
        return "?"  # not fetched, a cycle, or itself undetermined

    # space (AND) and '+' (complex, AND) -> '&';  ',' (OR) -> '|'
    expr = step.replace(" ", "&").replace("+", "&").replace(",", "|")
    expr = _KO.sub(lambda m: "1" if m.group(0) in ko_set else "0", expr)
    expr = re.sub(r"M\d{5}", _module_value, expr)

    if "?" not in expr:
        if not _SAFE.match(expr) or not expr:
            return None  # syntax this parser doesn't handle: unknown, not absent
        try:
            return bool(eval(expr))  # expr validated to [01()&|]
        except SyntaxError:
            return None

    pessimistic, optimistic = expr.replace("?", "0"), expr.replace("?", "1")
    if not _SAFE.match(pessimistic) or not pessimistic:
        return None
    try:
        lo, hi = bool(eval(pessimistic)), bool(eval(optimistic))
    except SyntaxError:
        return None
    return lo if lo == hi else None


def module_completeness(
    definition: str, ko_set: set[str],
    module_defs: dict[str, tuple[str, str, str]] | None = None,
    _stack: frozenset[str] = frozenset(),
) -> float:
    """Fraction of a module's real top-level steps satisfied by ``ko_set`` (0..1).

    ``--`` placeholder steps (reactions with no assigned KO) are excluded from
    both numerator and denominator, as are steps whose truth value can't be
    determined (see ``step_complete``) -- consistent with this codebase's NaN
    convention elsewhere, "unknown" is never silently folded into "absent".
    Returns NaN if every real step turns out to be undetermined. Pass
    ``module_defs`` so nested module references resolve; omit it and every
    reference is unknown rather than absent.
    """
    steps = [s for s in split_steps(definition) if not _is_gap(s)]
    if not steps:
        return 0.0
    determined = [step_complete(s, ko_set, module_defs, _stack) for s in steps]
    determined = [r for r in determined if r is not None]
    if not determined:
        return float("nan")
    return sum(determined) / len(determined)


# --------------------------------------------------------------------------- #
# Evidence-aware evaluation (for `mpph explain` and complete/partial states)
# --------------------------------------------------------------------------- #
@dataclass
class StepResult:
    """One module step and whether ``ko_set`` satisfies it.

    ``satisfied`` is ``None`` when the step's truth value is undetermined
    (unresolved nested module reference, cycle, or unsupported syntax) --
    never coerced to ``False``, since that would silently read as "confirmed
    absent" to any caller that doesn't check for ``None`` explicitly.
    """
    expression: str
    satisfied: bool | None
    matched_kos: list[str]
    missing_kos: list[str]


@dataclass
class ModuleEvaluation:
    """Full evidence for one module scored against one KO set."""
    module_id: str
    score: float
    state: str  # complete | partial | absent | unknown
    n_steps: int
    n_satisfied: int
    steps: list[StepResult] = field(default_factory=list)
    matched_kos: set[str] = field(default_factory=set)
    missing_kos: set[str] = field(default_factory=set)
    unresolved_references: set[str] = field(default_factory=set)
    parser_status: str = "valid"


def classify_state(score: float, complete_threshold: float = 1.0,
                   partial_threshold: float = 1e-9) -> str:
    """Map a completeness score to complete / partial / absent."""
    if score >= complete_threshold:
        return "complete"
    if score > partial_threshold:
        return "partial"
    return "absent"


def evaluate_module(
    module_id: str, definition: str, ko_set: set[str],
    module_defs: dict[str, tuple[str, str, str]] | None = None,
    *, complete_threshold: float = 1.0, _stack: frozenset[str] = frozenset(),
) -> ModuleEvaluation:
    """Score ``definition`` against ``ko_set`` with per-step matched/missing KOs."""
    steps_raw = [s for s in split_steps(definition) if not _is_gap(s)]
    if not steps_raw:
        return ModuleEvaluation(module_id, 0.0, "unknown", 0, 0,
                                parser_status="empty_definition")

    results: list[StepResult] = []
    unresolved: set[str] = set()
    matched_all: set[str] = set()
    missing_all: set[str] = set()
    for step in steps_raw:
        stripped = _strip_optional(step)
        kos = set(_KO.findall(stripped))
        matched = sorted(kos & ko_set)
        missing = sorted(kos - ko_set)
        satisfied = step_complete(step, ko_set, module_defs, _stack)
        results.append(StepResult(step, satisfied, matched, missing))
        matched_all.update(matched)
        if satisfied is None:
            # Only the refs this step's own outcome actually hinged on --
            # step_complete may have determined the step regardless (e.g. an
            # OR already satisfied by a known KO), in which case nothing here
            # was truly "unresolved" even though the step mentions a module id.
            for ref in re.findall(r"M\d{5}", stripped):
                if not (module_defs and ref in module_defs and ref not in _stack):
                    unresolved.add(ref)
        elif not satisfied:
            missing_all.update(missing)

    determined = [r.satisfied for r in results if r.satisfied is not None]
    n_sat = sum(1 for v in determined if v)
    score = (n_sat / len(determined)) if determined else float("nan")
    status = "valid" if not unresolved else "unresolved_references"
    state = ("unknown" if math.isnan(score)
             else classify_state(score, complete_threshold))
    return ModuleEvaluation(
        module_id, score, state,
        len(determined), n_sat, results, matched_all, missing_all,
        unresolved, status,
    )
