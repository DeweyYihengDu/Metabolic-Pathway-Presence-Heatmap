"""Render the benchmark figures from the committed result tables.

Figure 1 (annotation): accuracy against KEGG's own KO assignments, and
pairwise agreement between the three tools, across all three difficulty
tiers.
Figure 2 (enrichment): mpph gsea vs clusterProfiler GSEA, per pathway.

Reads only the CSVs committed alongside this script, so the figures can be
regenerated without re-running the multi-hour benchmark or having the
reference databases on hand.

    python make_figures.py --out-dir figures
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Validated categorical slots 1-3 (light surface): all-pairs CVD dE 9.2,
# normal-vision dE 24.0. Aqua sits below 3:1 contrast on this surface, so
# every bar also carries a printed value -- the documented relief for that.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE = "#fcfcfb"
INK, INK_MUTED = "#0b0b0b", "#52514e"
GRID = "#dcdbd6"

TOOL_COLOR = {"mpph_annotate": BLUE, "kofamscan": ORANGE, "eggnog_mapper": AQUA}
TOOL_LABEL = {"mpph_annotate": "mpph annotate", "kofamscan": "KofamScan",
              "eggnog_mapper": "eggNOG-mapper"}
PAIR_COLOR = {("mpph_annotate", "kofamscan"): BLUE,
              ("mpph_annotate", "eggnog_mapper"): ORANGE,
              ("kofamscan", "eggnog_mapper"): AQUA}
PAIR_LABEL = {("mpph_annotate", "kofamscan"): "mpph ↔ KofamScan",
              ("mpph_annotate", "eggnog_mapper"): "mpph ↔ eggNOG",
              ("kofamscan", "eggnog_mapper"): "KofamScan ↔ eggNOG"}

# Explicit order: the primary validation organism first within its tier,
# then increasing difficulty by tier. Alphabetical would bury E. coli.
ORGANISM_ORDER = ["eco", "bsu", "mja", "vbs", "lbac",
                  "SRR12479784_metabat2_bin17", "SRR13122148_dastool_bin007",
                  "SRR13122166_dastool_bin006"]

# Display names: species in italics, MAGs by bin with their SRA run.
DISPLAY = {
    "eco": "$\\it{E.\\ coli}$ K-12",
    "bsu": "$\\it{B.\\ subtilis}$ 168",
    "mja": "$\\it{M.\\ jannaschii}$",
    "vbs": "$\\it{Verrucomicrobia}$ sp. S94",
    "lbac": "$\\it{Lentisphaerae}$ sp. WC36",
    "SRR12479784_metabat2_bin17": "MAG bin17 (SRR12479784)",
    "SRR13122148_dastool_bin007": "MAG bin007 (SRR13122148)",
    "SRR13122166_dastool_bin006": "MAG bin006 (SRR13122166)",
}
# Deliberately terse: these are rotated into the right margin beside a block
# only 2-3 rows tall, so a longer phrase overruns into the neighbouring tier.
TIER_LABEL = {"1_model_organism": "Model", "2_nonmodel_genome": "Non-model",
              "3_mag": "MAGs"}


def _style_axis(ax) -> None:
    """Recessive grid and axes; the data carries the ink."""
    ax.set_facecolor(SURFACE)
    ax.xaxis.grid(True, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.yaxis.grid(False)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, length=0)


def _grouped_barh(ax, groups, series, values, colors, labels, bar_h=0.26):
    """Horizontal grouped bars: `groups` on y, one bar per series."""
    n = len(series)
    for i, key in enumerate(series):
        offsets = [g - (n - 1) / 2 * bar_h + i * bar_h for g in range(len(groups))]
        vals = [values[(grp, key)] for grp in groups]
        ax.barh(offsets, vals, height=bar_h * 0.86, color=colors[key],
                label=labels[key], zorder=3)
        for y, v in zip(offsets, vals):
            ax.text(v + 0.012, y, f"{v:.3f}", va="center", ha="left",
                    fontsize=7, color=INK_MUTED, zorder=4)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([DISPLAY.get(g, g) for g in groups], fontsize=8.5, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.18)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])


def _tier_dividers(ax, groups, tiers) -> None:
    """A hairline between tiers, with the tier named on the right."""
    seen, start = [], 0
    for i, t in enumerate(tiers):
        if i and t != tiers[i - 1]:
            ax.axhline(i - 0.5, color=GRID, linewidth=0.9, zorder=2)
            seen.append((tiers[i - 1], start, i - 1))
            start = i
    seen.append((tiers[-1], start, len(groups) - 1))
    for tier, lo, hi in seen:
        ax.text(1.155, (lo + hi) / 2, TIER_LABEL[tier], rotation=270,
                va="center", ha="center", fontsize=7.5, color=INK_MUTED)


def figure_annotation(accuracy: pd.DataFrame, agreement: pd.DataFrame,
                      out_path: Path) -> None:
    def _ordered(df):
        present = [o for o in ORGANISM_ORDER if o in set(df["organism"])]
        tiers = [df.loc[df["organism"] == o, "tier"].iloc[0] for o in present]
        return present, tiers

    acc = accuracy.dropna(subset=["f1"]).copy()
    acc_groups, acc_tiers = _ordered(acc)
    acc_vals = {(r.organism, r.tool): r.f1 for r in acc.itertuples()}

    agr = agreement.copy()
    agr_groups, agr_tiers = _ordered(agr)
    agr_vals = {(r.organism, (r.tool_a, r.tool_b)): r.jaccard for r in agr.itertuples()}

    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, figsize=(8.4, 9.2), facecolor=SURFACE,
        gridspec_kw={"height_ratios": [len(acc_groups), len(agr_groups)],
                     "hspace": 0.30})

    tools = ["mpph_annotate", "kofamscan", "eggnog_mapper"]
    _style_axis(ax_a)
    _grouped_barh(ax_a, acc_groups, tools, acc_vals, TOOL_COLOR, TOOL_LABEL)
    _tier_dividers(ax_a, acc_groups, acc_tiers)
    ax_a.set_xlabel("F1 vs KEGG's own KO assignments", fontsize=8.5, color=INK_MUTED)
    ax_a.set_title("a   Accuracy against an independently-derived reference",
                   loc="left", fontsize=10.5, color=INK, pad=26)
    # Legend above the axes: every bar reaches ~1.0, so there is no interior
    # space it could occupy without covering data.
    ax_a.legend(frameon=False, fontsize=8, ncol=3, loc="lower left",
                bbox_to_anchor=(0, 1.005), handlelength=1.4, columnspacing=1.6)

    pairs = list(PAIR_COLOR)
    _style_axis(ax_b)
    _grouped_barh(ax_b, agr_groups, pairs, agr_vals, PAIR_COLOR, PAIR_LABEL)
    _tier_dividers(ax_b, agr_groups, agr_tiers)
    ax_b.set_xlabel("Jaccard overlap of (gene, KO) calls", fontsize=8.5, color=INK_MUTED)
    ax_b.set_title("b   Agreement between tools, including where no reference exists",
                   loc="left", fontsize=10.5, color=INK, pad=26)
    ax_b.legend(frameon=False, fontsize=8, ncol=3, loc="lower left",
                bbox_to_anchor=(0, 1.005), handlelength=1.4, columnspacing=1.6)

    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {out_path} and {out_path.with_suffix('.pdf')}")


def figure_enrichment(mpph: pd.DataFrame, cp: pd.DataFrame, out_path: Path) -> None:
    merged = mpph.merge(cp, on="category_id", how="inner")
    rho, _ = stats.spearmanr(merged["NES_mpph"], merged["NES_cp"])

    fig, ax = plt.subplots(figsize=(5.0, 5.0), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, length=0)

    lim = [min(merged["NES_mpph"].min(), merged["NES_cp"].min()) - 0.25,
           max(merged["NES_mpph"].max(), merged["NES_cp"].max()) + 0.25]
    ax.plot(lim, lim, color=INK_MUTED, linewidth=1.0, linestyle=(0, (4, 3)),
            zorder=2, label="y = x")
    ax.scatter(merged["NES_mpph"], merged["NES_cp"], s=34, color=BLUE,
               edgecolor=SURFACE, linewidth=1.2, zorder=3, alpha=0.9)

    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("mpph gsea  NES", fontsize=9, color=INK_MUTED)
    ax.set_ylabel("clusterProfiler GSEA  NES", fontsize=9, color=INK_MUTED)
    ax.set_title(f"Same ranked list, same gene sets\n{len(merged)} KEGG metabolic "
                 f"pathways   ·   Spearman ρ = {rho:.3f}",
                 loc="left", fontsize=10, color=INK, pad=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {out_path} and {out_path.with_suffix('.pdf')} (rho={rho:.4f})")


def main() -> int:
    here = Path(__file__).parent
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", default=str(here / "figures"))
    args = p.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    accuracy = pd.read_csv(here / "annotation" / "results" / "summary_accuracy.csv")
    agreement = pd.read_csv(here / "annotation" / "results" / "summary_agreement.csv")
    figure_annotation(accuracy, agreement, out_dir / "fig1_annotation_benchmark.png")

    mpph = pd.read_csv(here / "enrichment" / "results_mpph_gsea.csv",
                       dtype={"category_id": str})[["category_id", "NES"]].rename(
        columns={"NES": "NES_mpph"})
    cp = pd.read_csv(here / "enrichment" / "results_clusterprofiler_gsea.csv",
                     dtype={"ID": str})[["ID", "NES"]].rename(
        columns={"ID": "category_id", "NES": "NES_cp"})
    figure_enrichment(mpph, cp, out_dir / "fig2_enrichment_benchmark.png")

    calib_path = here / "phylo" / "calibration.csv"
    if calib_path.exists():
        figure_phylo_calibration(pd.read_csv(calib_path),
                                 out_dir / "fig3_phylo_calibration.png")
    return 0



def figure_phylo_calibration(calib: pd.DataFrame, out_path: Path) -> None:
    """Type I error, uncorrected vs corrected, one point per simulated cell.

    Every point is a scenario with no group effect present, so every
    rejection is a false positive: a calibrated test sits at or below the
    nominal 0.05 line.
    """
    df = calib.dropna(subset=["phylo_type1"]).copy()
    fig, ax = plt.subplots(figsize=(5.4, 5.0), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, length=0)

    ax.axhline(0.05, color=INK_MUTED, linewidth=1.0, linestyle=(0, (4, 3)),
               zorder=2, label="nominal 0.05")
    ax.axvline(0.05, color=INK_MUTED, linewidth=1.0, linestyle=(0, (4, 3)),
               zorder=2)
    markers = {"balanced": "o", "unbalanced": "^"}
    for shape, marker in markers.items():
        sub = df[df["shape"] == shape]
        ax.scatter(sub["fisher_type1"], sub["phylo_type1"], s=42, marker=marker,
                   color=BLUE, edgecolor=SURFACE, linewidth=1.2, alpha=0.85,
                   zorder=3, label=f"{shape} tree")

    ax.set_xlim(-0.04, 1.02)
    # Tight to the data: nothing exceeds 0.072, and a 0.35 ceiling
    # would spend most of the panel on empty space above the story.
    ax.set_ylim(-0.005, 0.09)
    ax.set_xlabel("Type I error, uncorrected (Fisher)", fontsize=9, color=INK_MUTED)
    ax.set_ylabel("Type I error, corrected (--tree)", fontsize=9, color=INK_MUTED)
    ax.set_title(f"No group effect present: every rejection is a false positive\n"
                 f"{len(df)} simulated scenarios   ·   worst "
                 f"{df['fisher_type1'].max():.2f} → {df['phylo_type1'].max():.2f}",
                 loc="left", fontsize=10, color=INK, pad=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {out_path} and {out_path.with_suffix('.pdf')}")

if __name__ == "__main__":
    raise SystemExit(main())
