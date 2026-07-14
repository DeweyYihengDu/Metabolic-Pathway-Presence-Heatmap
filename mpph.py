#!/usr/bin/env python3
"""MPPH -- Metabolic Pathway Presence Heatmap.

Fetch KEGG metabolic-pathway presence/absence across every sequenced genome in
a taxon (genus, family, phylum, ...) and render a heatmap, optionally with
hierarchical clustering that doubles as a pathway-based phylogenetic tree.

Data source : KEGG REST API  (https://www.kegg.jp/kegg/rest/keggapi.html)
Reference   : Y.-H. Du & J.-H. Mu, "Metabolic-Pathway-Presence-Heatmap (MPPH):
              Constructing phylogenetic trees based on metabolic pathways",
              bioRxiv 2023. doi:10.1101/2023.06.27.546232

Example
-------
    python mpph.py Vibrio
    python mpph.py Prochlorococcus --drop-core --cluster --format png pdf
    python mpph.py "Vibrio cholerae" --exact --outdir results
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter

try:  # urllib3 ships with requests; import path differs across versions
    from urllib3.util.retry import Retry
except ImportError:  # pragma: no cover
    from requests.packages.urllib3.util.retry import Retry  # type: ignore

import matplotlib

matplotlib.use("Agg")  # headless-safe; we only save figures, never show them
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, linkage  # noqa: E402

KEGG_API_BASE = "https://rest.kegg.jp"
DEFAULT_CACHE = Path(".mpph_cache")
REQUEST_TIMEOUT = 30  # seconds

# --- Visual design tokens (validated categorical palette; light surface) ----
INK = "#0b0b0b"        # primary text
SECONDARY = "#52514e"  # organism labels
MUTED = "#898781"      # axes / captions
DENDRO = "#a6a49d"     # dendrogram links
ABSENT = "#eceef1"     # "pathway absent" cell
OTHER_COLOR = "#bcbab2"  # overflow / uncategorised pathways
# Fixed hue order, CVD-ordered (see dataviz reference palette):
CATEGORY_PALETTE = [
    "#2a78d6", "#1baf7a", "#eda100", "#008300",
    "#4a3aa7", "#e34948", "#e87ba4", "#eb6834",
]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Helvetica Neue", "Arial", "DejaVu Sans"],
    "svg.fonttype": "none",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.facecolor": "white",
})


# --------------------------------------------------------------------------- #
# KEGG access (with retries + on-disk caching)
# --------------------------------------------------------------------------- #
def make_session() -> requests.Session:
    """A requests session that retries transient failures with backoff."""
    retry = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
    )
    session = requests.Session()
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({"User-Agent": "MPPH/2.0 (+github.com/DeweyYihengDu)"})
    return session


def kegg_get(session: requests.Session, endpoint: str, cache_dir: Path | None) -> str:
    """GET ``{KEGG_API_BASE}/{endpoint}`` with disk caching. Returns raw text."""
    if cache_dir is not None:
        cache_file = cache_dir / (re.sub(r"[^0-9A-Za-z]+", "_", endpoint) + ".tsv")
        if cache_file.exists():
            return cache_file.read_text(encoding="utf-8")

    resp = session.get(f"{KEGG_API_BASE}/{endpoint}", timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    text = resp.text

    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(text, encoding="utf-8")
    return text


# --------------------------------------------------------------------------- #
# Parsing helpers
# --------------------------------------------------------------------------- #
def list_genomes(session: requests.Session, cache_dir: Path | None) -> list[tuple[str, str]]:
    """Return ``[(org_code, organism_name), ...]`` for every KEGG genome.

    Uses the ``list/genome`` endpoint (the old ``list/organism`` endpoint was
    retired by KEGG and now returns HTTP 400). Each line looks like::

        T00034<TAB>vch; Vibrio cholerae O1 El Tor N16961
    """
    text = kegg_get(session, "list/genome", cache_dir)
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


def select_organisms(
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


def get_pathways(
    session: requests.Session, org_code: str, cache_dir: Path | None
) -> dict[str, str]:
    """Return ``{pathway_number: pathway_name}`` for one organism.

    A pathway id such as ``vch01100`` is reduced to its 5-digit KEGG map number
    (``01100``) so the same pathway is comparable across organisms.
    """
    text = kegg_get(session, f"list/pathway/{org_code}", cache_dir)
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
    session: requests.Session, cache_dir: Path | None
) -> dict[str, str]:
    """Return ``{map_id: functional_category}`` from KEGG BRITE ``br08901``.

    The category is the second-level ("B") heading of the KEGG pathway
    hierarchy, e.g. ``Carbohydrate metabolism`` or ``Energy metabolism`` -- the
    granularity that is most informative for a microbial comparison.
    """
    text = kegg_get(session, "get/br:br08901", cache_dir)
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


# --------------------------------------------------------------------------- #
# Matrix construction & filtering
# --------------------------------------------------------------------------- #
def build_matrix(
    org_pathways: dict[str, dict[str, str]],
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Build an organism x pathway presence/absence (0/1) matrix.

    Returns the matrix and a ``{pathway_number: pathway_name}`` lookup.
    """
    pathway_names: dict[str, str] = {}
    for pmap in org_pathways.values():
        pathway_names.update(pmap)

    columns = sorted(pathway_names)
    rows = list(org_pathways)
    matrix = np.zeros((len(rows), len(columns)), dtype=int)
    col_index = {pid: j for j, pid in enumerate(columns)}

    for i, organism in enumerate(rows):
        for pid in org_pathways[organism]:
            matrix[i, col_index[pid]] = 1

    df = pd.DataFrame(matrix, index=rows, columns=columns)
    return df, pathway_names


def filter_matrix(
    df: pd.DataFrame,
    min_prevalence: float = 0.0,
    max_prevalence: float = 1.0,
    drop_core: bool = False,
) -> pd.DataFrame:
    """Drop uninformative pathway columns by prevalence across organisms.

    Prevalence = fraction of organisms possessing the pathway. ``drop_core``
    removes pathways present in *every* organism (equivalent to
    ``max_prevalence < 1``), which sharpens the signal that distinguishes taxa.
    """
    if df.shape[0] == 0:
        return df
    prevalence = df.mean(axis=0)
    if drop_core:
        max_prevalence = min(max_prevalence, 1.0 - 1e-9)
    keep = (prevalence >= min_prevalence) & (prevalence <= max_prevalence)
    return df.loc[:, keep]


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def _assign_category_colors(
    columns: list[str], categories: dict[str, str]
) -> tuple[list[str], dict[str, str], list[tuple[str, str]]]:
    """Map each pathway column to a functional category and a colour.

    The most frequent categories take the fixed palette slots (stable identity);
    any beyond the palette collapse into a muted "Other". Returns the per-column
    category list, a ``category -> hex`` map, and ordered legend entries.
    """
    col_cats = [categories.get(c, "Other") for c in columns]
    counts = Counter(c for c in col_cats if c != "Other")
    top = [cat for cat, _ in counts.most_common(len(CATEGORY_PALETTE))]
    color_of = {cat: CATEGORY_PALETTE[i] for i, cat in enumerate(top)}

    legend = [(cat, color_of[cat]) for cat in top]
    if any(c not in color_of for c in col_cats):
        legend.append(("Other", OTHER_COLOR))
    return col_cats, color_of, legend


def _leaf_order(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray | None]:
    """UPGMA leaf order for the rows of ``matrix`` (and the linkage, or None)."""
    if matrix.shape[0] < 3:
        return np.arange(matrix.shape[0]), None
    link = linkage(matrix, method="average", metric="euclidean")
    order = dendrogram(link, no_plot=True)["leaves"]
    return np.asarray(order), link


def plot_heatmap(
    df: pd.DataFrame,
    categories: dict[str, str],
    outfile: Path,
    title: str,
    subtitle: str,
    cluster: bool,
) -> None:
    """Render a presence/absence heatmap coloured by KEGG functional category.

    Present cells are tinted by the pathway's functional category (with a colour
    strip and legend); absent cells recede to a light neutral. With ``cluster``
    the rows (organisms) and columns (pathways) are UPGMA-ordered and flanked by
    dendrograms -- the pathway-based tree from the reference paper.
    """
    if df.shape[1] == 0:
        raise ValueError("Nothing to plot: the pathway matrix has no columns "
                         "after filtering. Loosen --min-prevalence/--drop-core.")

    binary = df.to_numpy()
    names = list(df.index)
    cols = list(df.columns)
    n_rows, n_cols = binary.shape

    # --- ordering (cluster or grouped-by-category) ---------------------------
    if cluster:
        row_order, row_link = _leaf_order(binary)
        col_order, col_link = _leaf_order(binary.T)
    else:
        row_link = col_link = None
        row_order = np.argsort([n.lower() for n in names])
        col_cats_raw = [categories.get(c, "Other") for c in cols]
        col_order = np.array(sorted(range(n_cols), key=lambda j: (col_cats_raw[j], cols[j])))

    binary = binary[np.ix_(row_order, col_order)]
    names = [names[i] for i in row_order]
    cols = [cols[j] for j in col_order]
    col_cats, color_of, legend = _assign_category_colors(cols, categories)

    def cat_rgb(cat: str) -> tuple[float, float, float]:
        return mcolors.to_rgb(color_of.get(cat, OTHER_COLOR))

    # --- build the RGB image: present -> category colour, absent -> neutral --
    absent_rgb = mcolors.to_rgb(ABSENT)
    rgb = np.empty((n_rows, n_cols, 3))
    for j, cat in enumerate(col_cats):
        present = binary[:, j] == 1
        rgb[present, j] = cat_rgb(cat)
        rgb[~present, j] = absent_rgb
    strip = np.array([[cat_rgb(c) for c in col_cats]])

    # --- figure geometry -----------------------------------------------------
    width = min(30.0, max(7.0, n_cols * 0.14 + 4.5))
    height = min(40.0, max(4.5, n_rows * 0.28 + 3.0))
    show_dendro = cluster and row_link is not None
    left_w = 0.11 if show_dendro else 0.012
    top_h = 0.15 if show_dendro else 0.012

    fig = plt.figure(figsize=(width, height))
    gs = fig.add_gridspec(
        3, 2,
        width_ratios=[left_w, 1.0], height_ratios=[top_h, 0.035, 1.0],
        wspace=0.015, hspace=0.02,
        left=0.02, right=0.995, top=0.88, bottom=0.16,
    )
    ax_top = fig.add_subplot(gs[0, 1])
    ax_strip = fig.add_subplot(gs[1, 1])
    ax_heat = fig.add_subplot(gs[2, 1])
    ax_left = fig.add_subplot(gs[2, 0])

    xmax, ymax = 10 * n_cols, 10 * n_rows

    # heatmap
    ax_heat.imshow(rgb, extent=[0, xmax, 0, ymax], aspect="auto",
                   origin="lower", interpolation="nearest")
    ax_heat.set_xlim(0, xmax)
    ax_heat.set_ylim(0, ymax)
    # white gridlines turn solid colour fields into a legible tile matrix
    if n_cols <= 160:
        ax_heat.set_xticks(np.arange(0, xmax + 1, 10), minor=False)
        for x in range(10, xmax, 10):
            ax_heat.axvline(x, color="white", linewidth=0.6)
    if n_rows <= 90:
        for y in range(10, ymax, 10):
            ax_heat.axhline(y, color="white", linewidth=0.6)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([10 * i + 5 for i in range(n_rows)])
    ax_heat.set_yticklabels(names, fontsize=min(10, 380 / max(n_rows, 1)),
                            color=SECONDARY, style="italic")
    ax_heat.yaxis.tick_right()
    ax_heat.tick_params(length=0)
    for s in ax_heat.spines.values():
        s.set_visible(False)

    # category strip
    ax_strip.imshow(strip, extent=[0, xmax, 0, 1], aspect="auto",
                    origin="lower", interpolation="nearest")
    ax_strip.set_xlim(0, xmax)
    if n_cols <= 160:
        for x in range(10, xmax, 10):
            ax_strip.axvline(x, color="white", linewidth=0.6)
    ax_strip.set_xticks([])
    ax_strip.set_yticks([])
    ax_strip.set_ylabel("category", rotation=0, ha="right", va="center",
                        fontsize=8, color=MUTED, labelpad=8)
    for s in ax_strip.spines.values():
        s.set_visible(False)

    # dendrograms
    if show_dendro:
        dendrogram(col_link, ax=ax_top, orientation="top", no_labels=True,
                   color_threshold=0, above_threshold_color=DENDRO)
        dendrogram(row_link, ax=ax_left, orientation="left", no_labels=True,
                   color_threshold=0, above_threshold_color=DENDRO)
        for coll in list(ax_top.collections) + list(ax_left.collections):
            coll.set_linewidth(0.8)
    for ax in (ax_top, ax_left):
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    if not show_dendro:
        ax_top.set_visible(False)
        ax_left.set_visible(False)

    # titles
    fig.suptitle(title, x=0.5, y=0.99, fontsize=17, fontweight="semibold",
                 color=INK, ha="center")
    fig.text(0.5, 0.925, subtitle, ha="center", fontsize=10.5, color=MUTED)

    # single legend along the bottom: functional categories + "absent"
    handles = [Patch(facecolor=color, edgecolor="none", label=cat)
               for cat, color in legend]
    handles.append(Patch(facecolor=ABSENT, edgecolor=MUTED, linewidth=0.5,
                         label="absent"))
    leg = fig.legend(
        handles=handles,
        title="Present pathway coloured by KEGG functional category  ·  light grey = absent",
        loc="lower center", bbox_to_anchor=(0.5, 0.005),
        ncol=min(len(handles), 5), frameon=False, fontsize=9,
        title_fontsize=10.5, handlelength=1.2, handleheight=1.2,
        columnspacing=1.6, labelspacing=0.7,
    )
    leg.get_title().set_color(INK)
    for txt in leg.get_texts():
        txt.set_color(SECONDARY)

    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run(args: argparse.Namespace) -> int:
    cache_dir = None if args.no_cache else Path(args.cache_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    slug = re.sub(r"[^0-9A-Za-z]+", "_", args.taxon).strip("_") or "taxon"

    session = make_session()

    print(f"[1/4] Fetching KEGG genome list ...", flush=True)
    genomes = list_genomes(session, cache_dir)
    organisms = select_organisms(genomes, args.taxon, args.exact)
    if not organisms:
        print(f"No KEGG genomes matched taxon {args.taxon!r}.", file=sys.stderr)
        print("Try a broader name, drop --exact, or check spelling.", file=sys.stderr)
        return 2
    print(f"      matched {len(organisms)} genome(s) for {args.taxon!r}")

    print(f"[2/4] Fetching pathways for {len(organisms)} genome(s) ...", flush=True)
    org_pathways: dict[str, dict[str, str]] = {}
    for i, (code, name) in enumerate(organisms, 1):
        label = f"{name} ({code})"
        org_pathways[label] = get_pathways(session, code, cache_dir)
        print(f"      [{i}/{len(organisms)}] {label}", flush=True)

    print("[3/4] Building & filtering matrix ...", flush=True)
    categories = fetch_pathway_categories(session, cache_dir)
    df, pathway_names = build_matrix(org_pathways)
    if not args.keep_overview:
        overview = {m for m, c in categories.items()
                    if c == "Global and overview maps"}
        df = df.drop(columns=[c for c in df.columns if c in overview])
    df = filter_matrix(
        df,
        min_prevalence=args.min_prevalence,
        max_prevalence=args.max_prevalence,
        drop_core=args.drop_core,
    )
    print(f"      matrix: {df.shape[0]} organisms x {df.shape[1]} pathways")

    matrix_csv = outdir / f"{slug}_matrix.csv"
    pathways_csv = outdir / f"{slug}_pathways.csv"
    df.to_csv(matrix_csv, index_label="organism")
    kept = sorted(set(df.columns) & set(pathway_names))
    pd.DataFrame(
        {
            "map_id": kept,
            "pathway_name": [pathway_names[m] for m in kept],
            "category": [categories.get(m, "Other") for m in kept],
        }
    ).to_csv(pathways_csv, index=False)

    print("[4/4] Rendering heatmap ...", flush=True)
    title = f"Metabolic pathway presence · {args.taxon}"
    subtitle = (f"{df.shape[0]} genomes × {df.shape[1]} KEGG pathways"
                + ("  ·  UPGMA-clustered" if args.cluster else ""))
    figures = []
    for fmt in args.format:
        fig_path = outdir / f"{slug}_heatmap.{fmt}"
        plot_heatmap(df, categories, fig_path, title, subtitle, cluster=args.cluster)
        figures.append(fig_path.name)

    manifest = {
        "taxon": args.taxon,
        "exact_match": args.exact,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "kegg_api": KEGG_API_BASE,
        "n_organisms": int(df.shape[0]),
        "n_pathways": int(df.shape[1]),
        "filters": {
            "min_prevalence": args.min_prevalence,
            "max_prevalence": args.max_prevalence,
            "drop_core": args.drop_core,
        },
        "clustered": bool(args.cluster),
        "outputs": {
            "matrix": matrix_csv.name,
            "pathways": pathways_csv.name,
            "figures": figures,
        },
    }
    (outdir / f"{slug}_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    print(f"\nDone. Outputs in {outdir}/")
    for p in sorted(outdir.glob(f"{slug}_*")):
        print(f"  - {p.name}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mpph",
        description="Metabolic Pathway Presence Heatmap from the KEGG API.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("taxon", help="Taxon name to query, e.g. 'Vibrio' or 'Prochlorococcus'.")
    p.add_argument("--exact", action="store_true",
                   help="Match organism names that START WITH the taxon "
                        "(e.g. a full 'Genus species'); default is whole-word match.")
    p.add_argument("--outdir", default="output", help="Directory for outputs.")
    p.add_argument("--format", nargs="+", default=["pdf"],
                   choices=["pdf", "png", "svg"], help="Figure format(s).")
    p.add_argument("--cluster", action="store_true",
                   help="Hierarchically cluster rows/columns (UPGMA); adds a "
                        "pathway-based dendrogram alongside the heatmap.")
    p.add_argument("--min-prevalence", type=float, default=0.0,
                   help="Keep pathways present in >= this fraction of organisms.")
    p.add_argument("--max-prevalence", type=float, default=1.0,
                   help="Keep pathways present in <= this fraction of organisms.")
    p.add_argument("--drop-core", action="store_true",
                   help="Drop pathways present in ALL organisms (uninformative).")
    p.add_argument("--keep-overview", action="store_true",
                   help="Keep KEGG 'Global and overview maps' (e.g. 01100 "
                        "Metabolic pathways); dropped by default as aggregate maps.")
    p.add_argument("--cache-dir", default=str(DEFAULT_CACHE),
                   help="Directory for cached KEGG responses.")
    p.add_argument("--no-cache", action="store_true",
                   help="Disable on-disk caching of KEGG responses.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except requests.RequestException as exc:
        print(f"KEGG request failed: {exc}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
