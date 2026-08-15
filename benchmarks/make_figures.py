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
import numpy as np
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
# Within the model tier, ordered by domain (Bacteria, Archaea, Eukaryota) so
# the domain labels down the right margin come out as contiguous blocks.
ORGANISM_ORDER = ["eco", "bsu", "mja", "sce", "ath", "vbs", "lbac",
                  "SRR12479784_metabat2_bin17", "SRR13122148_dastool_bin007",
                  "SRR13122166_dastool_bin006"]

# Display names: species in italics, MAGs by bin with their SRA run. Proteome
# size is shown for the eukaryotes only, where it is the point -- Arabidopsis
# is 11x E. coli and is what makes this a scale test, not just a scope one.
DISPLAY = {
    "eco": "$\\it{E.\\ coli}$ K-12",
    "bsu": "$\\it{B.\\ subtilis}$ 168",
    "mja": "$\\it{M.\\ jannaschii}$",
    "sce": "$\\it{S.\\ cerevisiae}$ (6,021 prot.)",
    "ath": "$\\it{A.\\ thaliana}$ (48,265 prot.)",
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
# Abbreviated for the same reason -- see _domain_labels. An em dash for the
# MAGs: their domain is genuinely unknown, not omitted for space.
DOMAIN_LABEL = {"Bacteria": "Bact.", "Archaea": "Arch.",
                "Eukaryota": "Euk.", "unassigned": "—"}

# Right margin holds two annotation columns outside the data range. The gap
# has to clear the value labels, which sit at v + 0.012 and reach x ~ 1.06 on
# panel b where every mpph<->KofamScan bar is exactly 1.000 -- at a narrower
# X_MAX those labels collide with the domain column.
X_MAX = 1.30
X_DOMAIN, X_TIER = 1.135, 1.25


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
    ax.set_xlim(0, X_MAX)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])


def _runs(labels):
    """Contiguous runs as (label, first_index, last_index)."""
    out, start = [], 0
    for i, t in enumerate(labels):
        if i and t != labels[i - 1]:
            out.append((labels[i - 1], start, i - 1))
            start = i
    out.append((labels[-1], start, len(labels) - 1))
    return out


def _tier_dividers(ax, groups, tiers) -> None:
    """A hairline between tiers, with the tier named on the right."""
    for tier, lo, hi in _runs(tiers):
        if lo:
            ax.axhline(lo - 0.5, color=GRID, linewidth=0.9, zorder=2)
        ax.text(X_TIER, (lo + hi) / 2, TIER_LABEL[tier], rotation=270,
                va="center", ha="center", fontsize=7.5, color=INK_MUTED)


def _domain_labels(ax, domains) -> None:
    """Domain of life in a second, inner margin column.

    Deliberately unrotated and abbreviated: `_tier_dividers` can rotate
    because its blocks are 2-3 rows tall, but Archaea is a single row here
    and a rotated word would overrun into its neighbours. No dividers --
    the tier hairlines already segment the axis, and a second set of rules
    would read as a grid rather than as an annotation.
    """
    for domain, lo, hi in _runs(domains):
        ax.text(X_DOMAIN, (lo + hi) / 2, DOMAIN_LABEL[domain], va="center",
                ha="center", fontsize=6.8, color=INK_MUTED, zorder=4)


def figure_annotation(accuracy: pd.DataFrame, agreement: pd.DataFrame,
                      out_path: Path) -> None:
    def _ordered(df):
        known = set(ORGANISM_ORDER)
        missing = sorted(set(df["organism"]) - known)
        if missing:
            # Silently dropping a benchmarked genome from the figure would be
            # the worst possible failure mode here -- it would look like a
            # clean result rather than a missing one.
            raise KeyError(f"organism(s) not in ORGANISM_ORDER: {missing}")
        present = [o for o in ORGANISM_ORDER if o in set(df["organism"])]

        def first(organism, column):
            # Boolean mask, not .set_index().loc[] -- the latter returns a
            # scalar rather than a Series when an organism has exactly one
            # row, and .iloc[0] would then raise.
            return df.loc[df["organism"] == organism, column].iloc[0]

        return (present,
                [first(o, "tier") for o in present],
                [first(o, "domain") for o in present])

    acc = accuracy.dropna(subset=["f1"]).copy()
    acc_groups, acc_tiers, acc_domains = _ordered(acc)
    acc_vals = {(r.organism, r.tool): r.f1 for r in acc.itertuples()}

    agr = agreement.copy()
    agr_groups, agr_tiers, agr_domains = _ordered(agr)
    agr_vals = {(r.organism, (r.tool_a, r.tool_b)): r.jaccard for r in agr.itertuples()}

    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, figsize=(8.4, 9.2), facecolor=SURFACE,
        gridspec_kw={"height_ratios": [len(acc_groups), len(agr_groups)],
                     "hspace": 0.30})

    tools = ["mpph_annotate", "kofamscan", "eggnog_mapper"]
    _style_axis(ax_a)
    _grouped_barh(ax_a, acc_groups, tools, acc_vals, TOOL_COLOR, TOOL_LABEL)
    _tier_dividers(ax_a, acc_groups, acc_tiers)
    _domain_labels(ax_a, acc_domains)
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
    _domain_labels(ax_b, agr_domains)
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

    curve_path = here / "threshold_audit" / "results" / "curve.csv"
    summary_path = here / "threshold_audit" / "results" / "summary.csv"
    if curve_path.exists() and summary_path.exists():
        figure_margin_precision(pd.read_csv(curve_path),
                                pd.read_csv(summary_path),
                                out_dir / "fig4_margin_precision.png")
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

DISPLAY_SHORT = {
    "eco": r"$\it{E.\ coli}$", "bsu": r"$\it{B.\ subtilis}$",
    "mja": r"$\it{M.\ jannaschii}$", "vbs": r"$\it{Verrucomicrobia}$ S94",
    "lbac": r"$\it{Lentisphaerae}$ WC36", "sce": r"$\it{S.\ cerevisiae}$",
    "ath": r"$\it{A.\ thaliana}$",
}


def figure_margin_precision(curve: pd.DataFrame, summary: pd.DataFrame,
                            out_path: Path) -> None:
    """The reliability of a KO call as a function of its margin above threshold.

    This is a property of the field's method, not of this implementation:
    every tool built on KOfam emits all of these as the same `1`. Panel (a) is
    the curve per genome; panel (b) is the consequence -- the small, knowable
    set of low-margin calls carries a wildly disproportionate share of errors.
    """
    order = ["eco", "bsu", "mja", "vbs", "lbac", "sce", "ath"]
    present = [o for o in order if o in set(curve["organism"])]

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(11.0, 4.6), facecolor=SURFACE,
        gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.28})

    for ax in (ax_a, ax_b):
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=INK_MUTED, length=0)

    # (a) precision vs margin. x is the band index so the open-ended top bin
    # is drawn at a finite position; a log axis would misrepresent it.
    bands = (curve[curve["organism"] == present[0]]
             .sort_values("margin_lo")[["margin_lo", "margin_hi"]])
    labels = [f"{int(lo)}–{int(hi)}" if np.isfinite(hi) else f"≥{int(lo)}"
              for lo, hi in bands.itertuples(index=False)]
    xs = range(len(labels))
    palette = [BLUE, ORANGE, AQUA, "#7b5cd6", "#c2185b", "#00796b", "#8d6e63"]
    for org, colour in zip(present, palette):
        sub = curve[curve["organism"] == org].sort_values("margin_lo")
        ax_a.plot(xs, sub["precision"], marker="o", markersize=4.5,
                  linewidth=1.7, color=colour, zorder=3,
                  label=DISPLAY_SHORT.get(org, org))
    ax_a.set_xticks(list(xs))
    ax_a.set_xticklabels(labels, fontsize=7.5, rotation=30, ha="right")
    ax_a.set_ylim(0, 1.02)
    ax_a.set_xlabel("margin above that KO's own threshold (bits)",
                    fontsize=8.5, color=INK_MUTED)
    ax_a.set_ylabel("precision vs KEGG's own assignments",
                    fontsize=8.5, color=INK_MUTED)
    ax_a.set_title("a   Every one of these is written as the same `1`",
                   loc="left", fontsize=10.5, color=INK, pad=26)
    ax_a.legend(frameon=False, fontsize=7.5, loc="lower right", ncol=2,
                handlelength=1.4)

    # (b) share of calls vs share of errors, low-margin band.
    s = summary.set_index("organism").loc[present]
    y = np.arange(len(present))
    h = 0.36
    ax_b.barh(y - h / 2, 100 * s["frac_calls_low_margin"], height=h,
              color=BLUE, zorder=3, label="share of all calls")
    ax_b.barh(y + h / 2, 100 * s["frac_fp_from_low_margin"], height=h,
              color=ORANGE, zorder=3, label="share of all false positives")
    for i, org in enumerate(present):
        ax_b.text(100 * s.loc[org, "frac_fp_from_low_margin"] + 1.0, i + h / 2,
                  f"{s.loc[org, 'fp_enrichment_low_margin']:.1f}×",
                  va="center", fontsize=7.5, color=INK_MUTED, zorder=4)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([DISPLAY_SHORT.get(o, o) for o in present], fontsize=8.5)
    ax_b.invert_yaxis()
    # Headroom for the "N.Nx" enrichment labels, which are drawn just past the
    # end of the longest bar and would otherwise be clipped at the axis edge.
    ax_b.set_xlim(0, 100 * s["frac_fp_from_low_margin"].max() + 9)
    ax_b.set_xlabel("percent", fontsize=8.5, color=INK_MUTED)
    ax_b.set_title("b   Calls within 20 bits of threshold",
                   loc="left", fontsize=10.5, color=INK, pad=26)
    # Above the axes, not inside it: the y-axis is inverted so the last rows
    # sit at the bottom right, exactly where a "lower right" legend lands.
    ax_b.legend(frameon=False, fontsize=7.5, ncol=2, loc="lower left",
                bbox_to_anchor=(0, 1.005), handlelength=1.4)

    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight",
                facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {out_path} and {out_path.with_suffix('.pdf')}")

if __name__ == "__main__":
    raise SystemExit(main())
