"""User-defined and built-in metabolic *trait* panels.

A trait is a named set of steps over KO ids, more ecological/biogeochemical than
a raw KEGG module (nitrogen fixation, methanogenesis, ...). Panels are JSON or
YAML (YAML needs PyYAML). Each step is one of:

* ``all_of``  -- every KO must be present (an essential complex/step);
* ``any_of``  -- at least one KO present (alternatives);
* ``optional`` -- ignored in the score (annotation-only).

Trait score = satisfied required steps / total required steps.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

BUILTIN_DIR = Path(__file__).parent / "data" / "traits"


def _resolve_panel_path(path: str | Path) -> Path:
    p = Path(path)
    if not p.exists():
        builtin = BUILTIN_DIR / f"{p.name}.json"
        if builtin.exists():
            return builtin
        raise FileNotFoundError(f"trait panel not found: {path}")
    return p


def _read_panel_data(p: Path) -> dict:
    text = p.read_text(encoding="utf-8")
    if p.suffix in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover
            raise ImportError("YAML trait panels need PyYAML "
                              "(`pip install pyyaml`); or use JSON.") from exc
        return yaml.safe_load(text)
    return json.loads(text)


def load_trait_panel(path: str | Path) -> dict[str, dict]:
    """Load a trait panel from a JSON or YAML file. Returns ``{trait_id: {...}}``."""
    p = _resolve_panel_path(path)
    data = _read_panel_data(p)
    traits = data.get("traits", data)
    if not isinstance(traits, dict):
        raise TypeError("trait panel must map trait ids to definitions")
    return traits


def panel_provenance(path: str | Path) -> dict:
    """Panel-level metadata for a manifest: resolved path, content hash, and
    whatever top-level fields the panel file itself declares (``panel_version``,
    ``schema_version``, ``description``, or any custom field a curator adds).

    The hash lets a later run detect that a panel file changed underneath a
    previously-recorded analysis -- built-in panels are versioned data, not
    frozen constants, and can be revised as marker sets are corrected.
    """
    p = _resolve_panel_path(path)
    data = _read_panel_data(p)
    meta = {k: v for k, v in data.items() if k != "traits"} if isinstance(data, dict) else {}
    return {
        "panel_path": str(p),
        "panel_sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        **meta,
    }


def _step_satisfied(step: dict, ko_set: set[str]) -> bool:
    if "all_of" in step:
        return all(k in ko_set for k in step["all_of"])
    if "any_of" in step:
        return any(k in ko_set for k in step["any_of"])
    return True  # a bare optional step is always "satisfied" but not required


def score_trait(steps: list[dict], ko_set: set[str]) -> float:
    """Fraction of a trait's required steps satisfied by ``ko_set`` (0..1)."""
    required = [s for s in steps if "all_of" in s or "any_of" in s]
    if not required:
        return 0.0
    return sum(_step_satisfied(s, ko_set) for s in required) / len(required)


def score_traits(
    org_kos: dict[str, set[str]], panel: dict[str, dict],
) -> tuple[pd.DataFrame, dict[str, str], dict[str, str]]:
    """Score every organism against every trait.

    Returns ``(matrix, names, categories)`` where ``matrix`` is organisms x
    traits (0..1), ``names`` maps trait id to display name, and ``categories``
    maps trait id to its category (for the figure colour strip).
    """
    traits = sorted(panel)
    names = {t: panel[t].get("name", t) for t in traits}
    categories = {t: panel[t].get("category", "Trait") for t in traits}
    data = {
        org: [score_trait(panel[t].get("steps", []), kos) for t in traits]
        for org, kos in org_kos.items()
    }
    matrix = pd.DataFrame(data, index=traits).T
    return matrix, names, categories
