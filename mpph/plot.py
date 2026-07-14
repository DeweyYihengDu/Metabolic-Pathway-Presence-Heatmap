"""The MPPH figure engine: category-aware heatmaps with UPGMA dendrograms."""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe; we only save figures, never show them
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, linkage  # noqa: E402

# --- Visual design tokens (validated categorical palette; light surface) ----
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
DENDRO = "#a6a49d"
ABSENT = "#eceef1"
OTHER_COLOR = "#bcbab2"
CATEGORY_PALETTE = [
    "#2a78d6", "#1baf7a", "#eda100", "#008300",
    "#4a3aa7", "#e34948", "#e87ba4", "#eb6834",
]
COMPLETENESS_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "mpph_seq", ["#eceef1", "#9ec5f4", "#2a78d6", "#184f95"]
)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Helvetica Neue", "Arial", "DejaVu Sans"],
    "svg.fonttype": "none",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.facecolor": "white",
})


def _assign_category_colors(columns, categories):
    """Map each column to a functional category colour + ordered legend."""
    col_cats = [categories.get(c, "Other") for c in columns]
    counts = Counter(c for c in col_cats if c != "Other")
    top = [cat for cat, _ in counts.most_common(len(CATEGORY_PALETTE))]
    color_of = {cat: CATEGORY_PALETTE[i] for i, cat in enumerate(top)}
    legend = [(cat, color_of[cat]) for cat in top]
    if any(c not in color_of for c in col_cats):
        legend.append(("Other", OTHER_COLOR))
    return col_cats, color_of, legend


def _leaf_order(matrix, metric):
    """UPGMA leaf order for the rows of ``matrix`` (and the linkage, or None)."""
    if matrix.shape[0] < 3:
        return np.arange(matrix.shape[0]), None
    link = linkage(matrix, method="average", metric=metric)
    order = dendrogram(link, no_plot=True)["leaves"]
    return np.asarray(order), link


def plot_matrix(
    df: pd.DataFrame,
    categories: dict[str, str],
    outfile: Path,
    title: str,
    subtitle: str,
    *,
    cluster: bool = False,
    mode: str = "presence",
    metric: str = "euclidean",
) -> dict:
    """Render the matrix; return ``{row_link, col_link, row_labels, col_labels}``.

    ``mode='presence'`` colours present cells by functional category;
    ``mode='completeness'`` colours cells by a sequential completeness ramp with
    a colorbar. Rows are organisms, columns are pathways or modules.
    """
    if df.shape[1] == 0:
        raise ValueError("Nothing to plot: the matrix has no columns after "
                         "filtering. Loosen --min-prevalence/--drop-core.")

    values = df.to_numpy(dtype=float)
    names = list(df.index)
    cols = list(df.columns)
    n_rows, n_cols = values.shape
    continuous = mode == "completeness"

    # --- ordering ------------------------------------------------------------
    if cluster:
        row_order, row_link = _leaf_order(values, metric)
        col_order, col_link = _leaf_order(values.T, metric)
    else:
        row_link = col_link = None
        row_order = np.argsort([n.lower() for n in names])
        raw = [categories.get(c, "Other") for c in cols]
        col_order = np.array(sorted(range(n_cols), key=lambda j: (raw[j], cols[j])))

    values = values[np.ix_(row_order, col_order)]
    names = [names[i] for i in row_order]
    cols = [cols[j] for j in col_order]
    col_cats, color_of, legend = _assign_category_colors(cols, categories)

    def cat_rgb(cat):
        return mcolors.to_rgb(color_of.get(cat, OTHER_COLOR))

    # --- RGB image -----------------------------------------------------------
    rgb = np.empty((n_rows, n_cols, 3))
    if continuous:
        rgb = COMPLETENESS_CMAP(values)[:, :, :3]
    else:
        absent_rgb = mcolors.to_rgb(ABSENT)
        for j, cat in enumerate(col_cats):
            present = values[:, j] > 0
            rgb[present, j] = cat_rgb(cat)
            rgb[~present, j] = absent_rgb
    strip = np.array([[cat_rgb(c) for c in col_cats]])

    # --- geometry ------------------------------------------------------------
    width = min(30.0, max(7.0, n_cols * 0.14 + 4.5))
    height = min(40.0, max(4.5, n_rows * 0.28 + 3.0))
    show_dendro = cluster and row_link is not None
    show_xlabels = n_cols <= 50
    left_w = 0.11 if show_dendro else 0.012
    top_h = 0.15 if show_dendro else 0.012

    fig = plt.figure(figsize=(width, height))
    gs = fig.add_gridspec(
        3, 2, width_ratios=[left_w, 1.0], height_ratios=[top_h, 0.035, 1.0],
        wspace=0.015, hspace=0.02, left=0.02, right=0.995, top=0.88,
        bottom=0.22 if show_xlabels else 0.16,
    )
    ax_top = fig.add_subplot(gs[0, 1])
    ax_strip = fig.add_subplot(gs[1, 1])
    ax_heat = fig.add_subplot(gs[2, 1])
    ax_left = fig.add_subplot(gs[2, 0])
    xmax, ymax = 10 * n_cols, 10 * n_rows

    # heatmap + gridlines (tile matrix)
    ax_heat.imshow(rgb, extent=[0, xmax, 0, ymax], aspect="auto",
                   origin="lower", interpolation="nearest")
    ax_heat.set_xlim(0, xmax)
    ax_heat.set_ylim(0, ymax)
    if n_cols <= 160:
        for x in range(10, xmax, 10):
            ax_heat.axvline(x, color="white", linewidth=0.6)
    if n_rows <= 90:
        for y in range(10, ymax, 10):
            ax_heat.axhline(y, color="white", linewidth=0.6)
    # Feature ids on the x-axis when few enough to stay legible; the full
    # ordered list is always in <name>_features.csv for traceability.
    if show_xlabels:
        ax_heat.set_xticks([10 * j + 5 for j in range(n_cols)])
        ax_heat.set_xticklabels(cols, rotation=90, color=MUTED,
                                fontsize=min(8.0, 420 / max(n_cols, 1)))
        ax_heat.xaxis.set_ticks_position("bottom")
    else:
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

    # legend + optional colorbar
    handles = [Patch(facecolor=color, edgecolor="none", label=cat)
               for cat, color in legend]
    if continuous:
        leg_title = "Category strip: KEGG module functional category"
        cax = fig.add_axes([0.055, 0.93, 0.17, 0.018])
        sm = ScalarMappable(norm=mcolors.Normalize(0, 1), cmap=COMPLETENESS_CMAP)
        cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cbar.set_label("module completeness", fontsize=9, color=SECONDARY)
        cbar.ax.tick_params(labelsize=8, color=MUTED, labelcolor=MUTED)
        cbar.outline.set_visible(False)
    else:
        leg_title = ("Present pathway coloured by KEGG functional category"
                     "  ·  light grey = absent")
        handles.append(Patch(facecolor=ABSENT, edgecolor=MUTED, linewidth=0.5,
                             label="absent"))
    leg = fig.legend(
        handles=handles, title=leg_title, loc="lower center",
        bbox_to_anchor=(0.5, 0.005), ncol=min(len(handles), 5), frameon=False,
        fontsize=9, title_fontsize=10.5, handlelength=1.2, handleheight=1.2,
        columnspacing=1.6, labelspacing=0.7,
    )
    leg.get_title().set_color(INK)
    for txt in leg.get_texts():
        txt.set_color(SECONDARY)

    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return {
        "row_link": row_link, "col_link": col_link,
        "row_labels": names, "col_labels": cols,
    }
