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
DENDRO = "#4f4e4a"  # dark enough to stay legible when the figure is scaled down
ABSENT = "#eceef1"
OTHER_COLOR = "#bcbab2"
CATEGORY_PALETTE = [
    "#2a78d6", "#1baf7a", "#eda100", "#008300",
    "#4a3aa7", "#e34948", "#e87ba4", "#eb6834",
]
# Fixed hue slots for the common KEGG categories, so a category keeps the same
# colour across figures (and across prevalence thresholds); anything else takes
# a leftover slot by frequency, then folds into "Other".
FIXED_CATEGORY_COLORS = {
    "Carbohydrate metabolism": "#2a78d6",
    "Energy metabolism": "#1baf7a",
    "Amino acid metabolism": "#eda100",
    "Metabolism of cofactors and vitamins": "#008300",
    "Nucleotide metabolism": "#4a3aa7",
    "Lipid metabolism": "#e34948",
    "Glycan biosynthesis and metabolism": "#e87ba4",
    "Metabolism of terpenoids and polyketides": "#eb6834",
}
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
    """Map each column to a functional category colour + ordered legend.

    Known categories take a fixed hue (stable across figures); remaining
    categories take leftover palette slots by frequency, then fold into "Other".
    """
    col_cats = [categories.get(c, "Other") for c in columns]
    present = [c for c in dict.fromkeys(col_cats) if c != "Other"]

    color_of: dict[str, str] = {}
    for cat in present:
        if cat in FIXED_CATEGORY_COLORS:
            color_of[cat] = FIXED_CATEGORY_COLORS[cat]
    leftover = [c for c in CATEGORY_PALETTE if c not in color_of.values()]
    extras = [cat for cat, _ in Counter(
        c for c in col_cats if c != "Other" and c not in color_of).most_common()]
    for cat, color in zip(extras, leftover):
        color_of[cat] = color

    # Legend ordered as: fixed (in canonical order) then extras, then Other.
    ordered = [c for c in FIXED_CATEGORY_COLORS if c in color_of]
    ordered += [c for c in extras if c in color_of and c not in ordered]
    legend = [(cat, color_of[cat]) for cat in ordered]
    if any(c not in color_of for c in col_cats):
        legend.append(("Other", OTHER_COLOR))
    return col_cats, color_of, legend


def _leaf_order(matrix, metric):
    """UPGMA leaf order for the rows of ``matrix`` (and the linkage, or None)."""
    if matrix.shape[0] < 2:
        return np.arange(matrix.shape[0]), None
    link = linkage(matrix, method="average", metric=metric)
    order = dendrogram(link, no_plot=True)["leaves"]
    return np.asarray(order), link


def _cosmetic_linkage(link):
    """A copy of ``link`` with a minimum *visual* branch height.

    Ties (identical rows/columns, common for small marker panels) merge at
    height 0, which draws as an invisible flat line flush with the leaves --
    it looks broken rather than "all equal". This floor only affects the
    on-screen bracket; leaf order, Newick export and every other consumer keep
    the real (unmodified) linkage returned by ``plot_matrix``.
    """
    disp = link.copy()
    heights = disp[:, 2]
    max_h = heights.max()
    floor = 0.06 * max_h if max_h > 0 else 0.06
    disp[:, 2] = np.maximum.accumulate(np.maximum(heights, floor))
    return disp


def _truncate_label(label: str, max_len: int = 34) -> str:
    """Shorten a long feature id/name for an axis tick; full text stays in CSV."""
    return label if len(label) <= max_len else label[: max_len - 1] + "…"


def plot_ordination(
    coords: pd.DataFrame, explained, outfile: Path, title: str,
    groups: pd.Series | None = None,
) -> None:
    """Scatter of the first two PCoA axes, optionally coloured by group."""
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    # A rank-one solution (e.g. exactly 2 samples, or other degenerate distance
    # structure) has only one positive eigenvalue -- there is no real PCo2.
    rank_one = coords.shape[1] < 2
    x = coords.iloc[:, 0]
    y = pd.Series(0.0, index=coords.index) if rank_one else coords.iloc[:, 1]
    if groups is not None:
        labels = [groups.get(i, "n/a") for i in coords.index]
        uniq = list(dict.fromkeys(labels))
        cmap = {g: CATEGORY_PALETTE[i % len(CATEGORY_PALETTE)]
                for i, g in enumerate(uniq)}
        for g in uniq:
            m = [lab == g for lab in labels]
            ax.scatter(x[m], y[m], s=60, color=cmap[g], label=str(g),
                       edgecolor="white", linewidth=0.6)
        ax.legend(frameon=False, fontsize=9, title_fontsize=10)
    else:
        ax.scatter(x, y, s=60, color=CATEGORY_PALETTE[0],
                   edgecolor="white", linewidth=0.6)
    ev = list(explained) + [0, 0]
    ax.set_xlabel(f"PCo1 ({ev[0] * 100:.1f}%)", color=SECONDARY)
    if rank_one:
        ax.set_ylabel("PCo2 unavailable (rank-one solution)", color=SECONDARY)
        ax.set_yticks([])
    else:
        ax.set_ylabel(f"PCo2 ({ev[1] * 100:.1f}%)", color=SECONDARY)
    ax.set_title(title, fontsize=15, color=INK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_volcano(diff: pd.DataFrame, outfile: Path, title: str,
                 alpha: float = 0.05) -> None:
    """Volcano plot for `compare` output: effect size vs -log10(q-value)."""
    if "prevalence_diff" in diff.columns:
        xcol, xlabel = "prevalence_diff", "prevalence difference (A − B)"
    elif "cliffs_delta" in diff.columns:
        xcol, xlabel = "cliffs_delta", "Cliff's delta (A − B)"
    else:
        xcol, xlabel = diff.columns[0], diff.columns[0]
    q = diff["q_value"].to_numpy(dtype=float)
    qmin = q[q > 0].min() if (q > 0).any() else 1e-300
    y = -np.log10(np.clip(q, qmin, 1.0))
    sig = q < alpha
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.scatter(diff[xcol][~sig], y[~sig], s=28, color="#bcbab2",
               edgecolor="white", linewidth=0.4, label=f"q ≥ {alpha}")
    ax.scatter(diff[xcol][sig], y[sig], s=34, color="#e34948",
               edgecolor="white", linewidth=0.4, label=f"q < {alpha}")
    ax.axhline(-np.log10(alpha), color=MUTED, linewidth=0.8, linestyle="--")
    ax.set_xlabel(xlabel, color=SECONDARY)
    ax.set_ylabel("−log10(q-value)", color=SECONDARY)
    ax.set_title(title, fontsize=15, color=INK)
    ax.legend(frameon=False, fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)


_PAN_COLORS = {"core": "#184f95", "soft-core": "#2a78d6",
               "shell": "#9ec5f4", "cloud": "#eceef1"}


def plot_prevalence(pan: pd.DataFrame, outfile: Path, title: str) -> None:
    """Prevalence histogram coloured by pan-functional class."""
    fig, ax = plt.subplots(figsize=(8, 5))
    order = ["cloud", "shell", "soft-core", "core"]
    bins = np.linspace(0, 1, 21)
    present = [pan.loc[pan["pan_class"] == c, "prevalence"] for c in order]
    ax.hist(present, bins=bins, stacked=True,
            color=[_PAN_COLORS[c] for c in order], label=order,
            edgecolor="white", linewidth=0.4)
    ax.set_xlabel("prevalence across organisms", color=SECONDARY)
    ax.set_ylabel("number of features", color=SECONDARY)
    ax.set_title(title, fontsize=15, color=INK)
    ax.legend(frameon=False, fontsize=9, title="pan class")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_accumulation(acc: pd.DataFrame, outfile: Path, title: str) -> None:
    """Pan / core accumulation curve (mean with 95% interval, if present)."""
    fig, ax = plt.subplots(figsize=(8, 5))
    x = acc["n_genomes"]
    if {"pan_ci_lower", "pan_ci_upper"} <= set(acc.columns):
        ax.fill_between(x, acc["pan_ci_lower"], acc["pan_ci_upper"],
                        color="#2a78d6", alpha=0.15, linewidth=0)
        ax.fill_between(x, acc["core_ci_lower"], acc["core_ci_upper"],
                        color="#e34948", alpha=0.15, linewidth=0)
    ax.plot(x, acc["pan_mean"], "-o", color="#2a78d6",
            markersize=4, label="pan (cumulative distinct)")
    ax.plot(x, acc["core_mean"], "-o", color="#e34948",
            markersize=4, label="core (shared by all)")
    ax.set_xlabel("number of organisms", color=SECONDARY)
    ax.set_ylabel("number of features", color=SECONDARY)
    ax.set_title(title, fontsize=15, color=INK)
    ax.legend(frameon=False, fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_enrichment(
    results: pd.DataFrame, outfile: Path, title: str,
    *, top_n: int = 20, alpha: float = 0.05,
) -> None:
    """Horizontal bar chart of the top enriched categories by q-value.

    Bar length is -log10(q-value); a dashed line marks the significance
    threshold. Each bar is annotated with its gene ratio (study hits / study
    total) so the plot carries both significance and effect size.
    """
    if results.empty:
        raise ValueError("Nothing to plot: no category met --min-category-size.")
    df = results.sort_values("q_value").head(top_n).iloc[::-1]
    q = df["q_value"].to_numpy(dtype=float)
    qmin = q[q > 0].min() if (q > 0).any() else 1e-300
    y = -np.log10(np.clip(q, qmin, 1.0))
    sig = q < alpha

    height = max(3.5, 0.32 * len(df) + 1.2)
    fig, ax = plt.subplots(figsize=(9.5, height))
    colors = [CATEGORY_PALETTE[0] if s else OTHER_COLOR for s in sig]
    labels = [f"{cid}  {name}"[:60] for cid, name in
              zip(df["category_id"], df["category_name"])]
    bars = ax.barh(range(len(df)), y, color=colors, height=0.65)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(labels, fontsize=9, color=SECONDARY)
    ax.axvline(-np.log10(alpha), color=MUTED, linewidth=0.8, linestyle="--")
    ax.set_xlabel("−log10(q-value)", color=SECONDARY)
    ax.set_title(title, fontsize=15, color=INK)
    for bar, ratio in zip(bars, df["gene_ratio"]):
        ax.text(bar.get_width() + max(y) * 0.015, bar.get_y() + bar.get_height() / 2,
                ratio, va="center", fontsize=8, color=MUTED)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    handles = [Patch(facecolor=CATEGORY_PALETTE[0], label=f"q < {alpha}"),
               Patch(facecolor=OTHER_COLOR, label=f"q ≥ {alpha}")]
    ax.legend(handles=handles, frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)


_GSEA_DIVERGING = mcolors.LinearSegmentedColormap.from_list(
    "mpph_gsea_diverging", ["#e34948", "#f4d9d8", "#fcfcfb", "#cfe0f6", "#2a78d6"]
)


def plot_gsea_running(
    ranked_scores: np.ndarray, running: np.ndarray, hits: np.ndarray,
    outfile: Path, title: str, subtitle: str = "",
) -> None:
    """The classic four-panel GSEA running-enrichment-score plot.

    Top to bottom, matching the standard Broad GSEA layout: (1) the running
    enrichment trajectory, filled to zero, peak marked; (2) a diverging
    red-blue strip of the ranking metric at every position, so the eye can
    see at a glance where the ranking crosses zero; (3) a tick-mark rug of
    where category members fall; (4) the ranking metric itself as a filled
    "mountain" chart.
    """
    n = len(ranked_scores)
    peak = int(np.argmax(running)) if abs(running.max()) >= abs(running.min()) \
        else int(np.argmin(running))
    fig, (ax_run, ax_strip, ax_rug, ax_score) = plt.subplots(
        4, 1, figsize=(8.5, 6.3), sharex=True,
        gridspec_kw={"height_ratios": [3, 0.22, 0.42, 1.1], "hspace": 0.06},
    )
    fig.subplots_adjust(top=0.85)

    ax_run.fill_between(range(n), running, 0, color=CATEGORY_PALETTE[0],
                        alpha=0.25, linewidth=0)
    ax_run.plot(range(n), running, color=CATEGORY_PALETTE[0], linewidth=2.2,
               solid_joinstyle="round")
    ax_run.axhline(0, color=MUTED, linewidth=0.6)
    ax_run.axvline(peak, color=SECONDARY, linewidth=0.9, linestyle="--")
    ax_run.plot(peak, running[peak], "o", color=CATEGORY_PALETTE[0],
               markersize=6, markeredgecolor="white", markeredgewidth=1)
    ax_run.set_ylabel("running ES", color=SECONDARY)
    pad = (running.max() - running.min()) * 0.12 or 0.05
    ax_run.set_ylim(min(running.min(), 0) - pad, max(running.max(), 0) + pad)
    fig.suptitle(title, x=0.5, y=0.98, fontsize=15, fontweight="semibold",
                color=INK, ha="center")
    if subtitle:
        fig.text(0.5, 0.905, subtitle, ha="center", fontsize=9.5, color=MUTED)
    for s in ("top", "right"):
        ax_run.spines[s].set_visible(False)

    # Diverging colour strip: the ranking metric's sign/magnitude at every
    # position -- the signature "heatmap" band of a canonical GSEA plot.
    vmax = np.abs(ranked_scores).max() or 1.0
    ax_strip.imshow(ranked_scores[np.newaxis, :], aspect="auto", cmap=_GSEA_DIVERGING,
                    vmin=-vmax, vmax=vmax, extent=[0, n, 0, 1], interpolation="nearest")
    ax_strip.set_yticks([])
    ax_strip.set_xticks([])
    for s in ax_strip.spines.values():
        s.set_visible(False)

    hit_pos = np.nonzero(hits)[0]
    ax_rug.vlines(hit_pos, 0, 1, color=INK, linewidth=0.8)
    ax_rug.set_ylim(0, 1)
    ax_rug.set_yticks([])
    for s in ax_rug.spines.values():
        s.set_visible(False)

    pos = np.clip(ranked_scores, 0, None)
    neg = np.clip(ranked_scores, None, 0)
    ax_score.fill_between(range(n), pos, 0, color="#2a78d6", linewidth=0)
    ax_score.fill_between(range(n), neg, 0, color="#e34948", linewidth=0)
    ax_score.set_xlabel("rank in ordered list", color=SECONDARY)
    ax_score.set_ylabel("ranking\nmetric", color=SECONDARY, fontsize=9)
    ax_score.axhline(0, color=MUTED, linewidth=0.6)
    for s in ("top", "right"):
        ax_score.spines[s].set_visible(False)

    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_gsea_summary(
    results: pd.DataFrame, outfile: Path, title: str,
    *, top_n: int = 20, alpha: float = 0.05,
) -> None:
    """Bar chart of the top GSEA hits by q-value, split by enrichment direction."""
    if results.empty:
        raise ValueError("Nothing to plot: no category met the size filters.")
    df = results.sort_values("q_value").head(top_n).iloc[::-1]
    sig = df["q_value"].to_numpy(dtype=float) < alpha
    pos = df["NES"].to_numpy(dtype=float) >= 0
    colors = [CATEGORY_PALETTE[0] if p else "#e34948" for p in pos]
    colors = [c if s else OTHER_COLOR for c, s in zip(colors, sig)]
    labels = [f"{cid}  {name}"[:60] for cid, name in
              zip(df["category_id"], df["category_name"])]

    height = max(3.5, 0.32 * len(df) + 1.2)
    fig, ax = plt.subplots(figsize=(9.5, height))
    ax.barh(range(len(df)), df["NES"], color=colors, height=0.65)
    ax.axvline(0, color=MUTED, linewidth=0.8)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(labels, fontsize=9, color=SECONDARY)
    ax.set_xlabel("normalized enrichment score (NES)", color=SECONDARY)
    ax.set_title(title, fontsize=15, color=INK)
    handles = [Patch(facecolor=CATEGORY_PALETTE[0], label=f"enriched (q < {alpha})"),
               Patch(facecolor="#e34948", label=f"depleted (q < {alpha})"),
               Patch(facecolor=OTHER_COLOR, label=f"q ≥ {alpha}")]
    ax.legend(handles=handles, frameon=False, fontsize=9, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)


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
    value_label: str = "module completeness",
    strip_label: str = "KEGG functional category",
) -> dict:
    """Render the matrix; return ``{row_link, col_link, row_labels, col_labels,
    row_link_labels, col_link_labels}``.

    ``mode='presence'`` colours present cells by functional category;
    ``mode='completeness'`` colours cells by a sequential completeness ramp with
    a colorbar. Rows are organisms, columns are pathways or modules.
    ``value_label`` / ``strip_label`` name the colorbar and category strip (so a
    trait heatmap reads "trait completeness" / "trait category").

    ``row_labels``/``col_labels`` are in dendrogram (display) order, matching
    the returned/plotted matrix. ``row_link``/``col_link`` index leaves in the
    *pre*-reorder input order instead, so Newick export must pair them with
    ``row_link_labels``/``col_link_labels``, not ``row_labels``/``col_labels``.
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

    # row_link/col_link index leaves in this pre-reorder sequence (the order
    # linkage() was called on) -- keep it for Newick export, which must look
    # up labels by that same indexing, not by the post-dendrogram display order.
    link_row_labels = list(names)
    link_col_labels = list(cols)

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
    show_row_dendro = cluster and row_link is not None
    show_col_dendro = cluster and col_link is not None
    show_xlabels = n_cols <= 50
    tick_labels = [_truncate_label(c) for c in cols] if show_xlabels else cols
    # Give the dendrograms enough of the canvas to actually read their shape.
    left_w = 0.20 if show_row_dendro else 0.012
    top_h = 0.24 if show_col_dendro else 0.012
    # Bottom margin must fit the rotated tick labels (their length varies a lot,
    # from "M00001" to a full trait name) plus the legend below them, or the
    # legend title collides with the labels.
    if show_xlabels:
        label_fontsize = min(8.0, 420 / max(n_cols, 1))
        max_chars = max((len(t) for t in tick_labels), default=0)
        label_frac = (max_chars * label_fontsize * 0.46) / (height * 72.0)
        # 0.06 is the floor proven to clear short labels (e.g. "M00001"); long
        # labels (trait names) need more and scale up from there.
        bottom = min(0.42, 0.16 + max(0.06, label_frac))
    else:
        bottom = 0.16

    fig = plt.figure(figsize=(width, height))
    gs = fig.add_gridspec(
        3, 2, width_ratios=[left_w, 1.0], height_ratios=[top_h, 0.035, 1.0],
        wspace=0.015, hspace=0.02, left=0.02, right=0.995, top=0.88,
        bottom=bottom,
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
        ax_heat.set_xticklabels(tick_labels, rotation=90, color=MUTED,
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

    # dendrograms (each axis drawn only if its linkage exists). The plotted
    # linkage gets a cosmetic minimum branch height (see _cosmetic_linkage);
    # the *returned* row_link/col_link below stay the real, unmodified linkage.
    if show_col_dendro:
        dendrogram(_cosmetic_linkage(col_link), ax=ax_top, orientation="top",
                   no_labels=True, color_threshold=0,
                   above_threshold_color=DENDRO)
        for coll in ax_top.collections:
            coll.set_linewidth(1.5)
    if show_row_dendro:
        dendrogram(_cosmetic_linkage(row_link), ax=ax_left, orientation="left",
                   no_labels=True, color_threshold=0,
                   above_threshold_color=DENDRO)
        for coll in ax_left.collections:
            coll.set_linewidth(1.5)
    for ax in (ax_top, ax_left):
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    if not show_col_dendro:
        ax_top.set_visible(False)
    if not show_row_dendro:
        ax_left.set_visible(False)

    # titles
    fig.suptitle(title, x=0.5, y=0.99, fontsize=17, fontweight="semibold",
                 color=INK, ha="center")
    fig.text(0.5, 0.925, subtitle, ha="center", fontsize=10.5, color=MUTED)

    # legend + optional colorbar
    handles = [Patch(facecolor=color, edgecolor="none", label=cat)
               for cat, color in legend]
    if continuous:
        leg_title = f"Category strip: {strip_label}"
        cax = fig.add_axes([0.055, 0.93, 0.17, 0.018])
        sm = ScalarMappable(norm=mcolors.Normalize(0, 1), cmap=COMPLETENESS_CMAP)
        cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cbar.set_label(value_label, fontsize=9, color=SECONDARY)
        cbar.ax.tick_params(labelsize=8, color=MUTED, labelcolor=MUTED)
        cbar.outline.set_visible(False)
    else:
        leg_title = (f"Present pathway coloured by {strip_label}"
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
        "row_link_labels": link_row_labels, "col_link_labels": link_col_labels,
    }
