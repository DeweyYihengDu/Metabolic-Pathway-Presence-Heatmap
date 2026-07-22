"""MPPH -- Metabolic Pathway Presence Heatmap.

Fetch KEGG metabolic-pathway presence (or module completeness) across genomes
in a taxon -- or across your own MAGs -- and render category-aware heatmaps
with UPGMA dendrograms.

The KEGG REST client and all data-access functions are part of the **public
API**: they talk directly to the official KEGG REST endpoint
(``https://rest.kegg.jp``), so anyone can reuse them to connect to KEGG
programmatically. See ``mpph.kegg`` for the low-level connection layer.

Reference: Y.-H. Du & J.-H. Mu, "Metabolic-Pathway-Presence-Heatmap (MPPH):
Constructing phylogenetic trees based on metabolic pathways", bioRxiv 2023,
doi:10.1101/2023.06.27.546232
"""
from __future__ import annotations

__version__ = "3.14.0"

# --- Public KEGG connection interface (open; talks to rest.kegg.jp) ---------
from .kegg import KEGG_API_BASE, kegg_get, kegg_release, make_session

# --- Public data-access interface -------------------------------------------
from .modules import (
    fetch_module_definitions,
    list_modules,
    module_completeness,
    organism_kos,
)
from .organisms import list_genomes, select_by_codes, select_by_taxon
from .pathways import fetch_pathway_categories, get_pathways

# --- Enrichment (KEGG pathway/module ORA; GO needs a user-supplied mapping) -
from .enrichment import (
    fetch_ko_module_membership,
    fetch_ko_pathway_membership,
    fetch_pathway_names,
    hypergeometric_enrichment,
    invert_membership,
    load_gene_go_map,
)

# --- GSEA (rank-based enrichment from expression or a pre-ranked list) ------
from .gsea import (
    enrichment_score,
    gsea_analysis,
    load_expression_matrix,
    load_ranked_list,
    rank_from_expression,
)

# --- KGML pathway/global-map diagrams (comparative two-group overlay) ------
from .kgml import (
    KGMLNode,
    KGMLPathway,
    KGMLReaction,
    fetch_kgml,
    normalize_map_id,
    parse_kgml,
    reactions_for_ortholog,
)

# --- Matrix + figure helpers ------------------------------------------------
from .matrix import build_completeness_matrix, build_presence_matrix, filter_matrix
from .plot import (
    plot_enrichment,
    plot_gsea_running,
    plot_gsea_summary,
    plot_kgml_map,
    plot_matrix,
)

__all__ = [
    "__version__",
    # KEGG connection interface (public)
    "KEGG_API_BASE",
    "make_session",
    "kegg_get",
    "kegg_release",
    # data access
    "list_genomes",
    "select_by_taxon",
    "select_by_codes",
    "get_pathways",
    "fetch_pathway_categories",
    "list_modules",
    "fetch_module_definitions",
    "organism_kos",
    # enrichment
    "fetch_ko_pathway_membership",
    "fetch_ko_module_membership",
    "fetch_pathway_names",
    "invert_membership",
    "hypergeometric_enrichment",
    "load_gene_go_map",
    # gsea
    "load_ranked_list",
    "load_expression_matrix",
    "rank_from_expression",
    "enrichment_score",
    "gsea_analysis",
    # KGML pathway/global-map diagrams
    "fetch_kgml",
    "parse_kgml",
    "normalize_map_id",
    "reactions_for_ortholog",
    "KGMLPathway",
    "KGMLNode",
    "KGMLReaction",
    # analysis + plotting
    "build_presence_matrix",
    "build_completeness_matrix",
    "filter_matrix",
    "module_completeness",
    "plot_matrix",
    "plot_enrichment",
    "plot_gsea_running",
    "plot_gsea_summary",
    "plot_kgml_map",
]
