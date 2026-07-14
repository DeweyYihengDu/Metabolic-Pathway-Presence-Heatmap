"""MPPH -- Metabolic Pathway Presence Heatmap.

Fetch KEGG metabolic-pathway presence (or module completeness) across genomes
in a taxon -- or across your own MAGs -- and render category-aware heatmaps
with UPGMA dendrograms.

Reference: Y.-H. Du & J.-H. Mu, "Metabolic-Pathway-Presence-Heatmap (MPPH):
Constructing phylogenetic trees based on metabolic pathways", bioRxiv 2023,
doi:10.1101/2023.06.27.546232
"""
from __future__ import annotations

__version__ = "3.0.0"

from .matrix import build_completeness_matrix, build_presence_matrix, filter_matrix
from .modules import module_completeness
from .plot import plot_matrix

__all__ = [
    "__version__",
    "build_presence_matrix",
    "build_completeness_matrix",
    "filter_matrix",
    "module_completeness",
    "plot_matrix",
]
