"""Assemble and filter the organism x feature matrix.

Two matrix flavours share the same downstream code:

* **presence** -- organism x pathway, values in {0, 1}
* **completeness** -- organism x module, values in [0, 1]
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .modules import module_completeness


def build_presence_matrix(
    org_pathways: dict[str, dict[str, str]],
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Organism x pathway presence/absence (0/1) matrix + ``{id: name}`` map."""
    names: dict[str, str] = {}
    for pmap in org_pathways.values():
        names.update(pmap)

    columns = sorted(names)
    rows = list(org_pathways)
    matrix = np.zeros((len(rows), len(columns)), dtype=float)
    col_index = {pid: j for j, pid in enumerate(columns)}
    for i, organism in enumerate(rows):
        for pid in org_pathways[organism]:
            matrix[i, col_index[pid]] = 1.0
    return pd.DataFrame(matrix, index=rows, columns=columns), names


def build_completeness_matrix(
    org_kos: dict[str, set[str]],
    module_defs: dict[str, tuple[str, str, str]],
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Organism x module completeness (0..1) matrix + ``{id: name}`` unused map.

    ``module_defs`` maps ``module_id -> (definition, category, type)``. The full
    map is passed to the scorer so nested module references resolve.
    """
    modules = sorted(module_defs)
    rows = list(org_kos)
    matrix = np.zeros((len(rows), len(modules)), dtype=float)
    for i, organism in enumerate(rows):
        kos = org_kos[organism]
        for j, mid in enumerate(modules):
            definition = module_defs[mid][0]
            matrix[i, j] = module_completeness(definition, kos, module_defs)
    names = {mid: mid for mid in modules}
    return pd.DataFrame(matrix, index=rows, columns=modules), names


def filter_matrix(
    df: pd.DataFrame,
    min_prevalence: float = 0.0,
    max_prevalence: float = 1.0,
    drop_core: bool = False,
    present_threshold: float = 1e-9,
    prevalence_state: str = "any",
    complete_threshold: float = 1.0,
) -> pd.DataFrame:
    """Drop uninformative feature columns by prevalence across organisms.

    With ``prevalence_state='any'`` (default) a feature counts as present when
    its value exceeds ``present_threshold`` (so a partially complete module
    counts). With ``prevalence_state='complete'`` only values >=
    ``complete_threshold`` count -- useful to filter on how often a module is
    *fully* complete rather than merely detectable. ``drop_core`` removes
    features present in every organism.
    """
    if df.shape[0] == 0 or df.shape[1] == 0:
        return df
    if prevalence_state == "complete":
        prevalence = (df >= complete_threshold).mean(axis=0)
    else:
        prevalence = (df > present_threshold).mean(axis=0)
    if drop_core:
        max_prevalence = min(max_prevalence, 1.0 - 1e-9)
    keep = (prevalence >= min_prevalence) & (prevalence <= max_prevalence)
    return df.loc[:, keep]
