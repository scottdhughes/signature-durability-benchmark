"""
Generate publication-quality paired forest plot for signature-durability-benchmark paper.

Figure 1: IFN-γ evaluated two ways — cross-context (all 30 cohorts) vs within-program
(6 interferon cohorts) — producing opposite conclusions.

Outputs:
  paper/figure1.pdf  (vector, for submission)
  paper/figure1.png  (300 DPI raster, for preview)
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from signature_durability_benchmark.meta_analysis import dl_random_effects_meta

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(BASE, "..")
EFFECTS_CSV = os.path.join(ROOT, "outputs", "canonical_v7", "per_cohort_effects.csv")
MANIFEST_TSV = os.path.join(ROOT, "config", "cohort_manifest.tsv")
OUT_PDF = os.path.join(ROOT, "paper", "figure1.pdf")
OUT_PNG = os.path.join(ROOT, "paper", "figure1.png")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
effects = pd.read_csv(EFFECTS_CSV)
manifest = pd.read_csv(MANIFEST_TSV, sep="\t")

cohort_program = dict(zip(manifest.cohort_id, manifest.biological_program))

SIG = "hallmark_ifng_response"

ifng = effects[effects.signature_id == SIG].copy()
ifng["biological_program"] = ifng.cohort_id.map(cohort_program)
ifng["ci_lo"] = ifng.cohens_d - 1.96 * np.sqrt(ifng.cohens_d_var)
ifng["ci_hi"] = ifng.cohens_d + 1.96 * np.sqrt(ifng.cohens_d_var)

# ---------------------------------------------------------------------------
# Color palette — monochrome + cool blue accents, program colors
# ---------------------------------------------------------------------------
PROGRAM_COLORS = {
    "interferon":    "#2171b5",   # strong blue
    "inflammation":  "#d62728",   # muted red
    "proliferation": "#2ca02c",   # green
    "hypoxia":       "#e6820e",   # orange
    "emt":           "#7b2d8b",   # purple
}
POOLED_COLOR_ALL = "#555555"       # neutral gray for cross-context pooled
POOLED_COLOR_WITHIN = "#2171b5"    # blue for within-program pooled

# ---------------------------------------------------------------------------
# Short cohort name helper
# ---------------------------------------------------------------------------
LABEL_MAP = {
    "sepsis_mortality_gse65682":   "Sepsis mortality (GSE65682)",
    "breast_prognosis_gse2034":    "Breast prognosis (GSE2034)",
    "breast_er_gse7390":           "Breast ER status (GSE7390)",
    "hypoxia_cellline_gse53012":   "Hypoxia cell line (GSE53012)",
    "sepsis_blood_gse28750":       "Sepsis blood (GSE28750)",
    "tb_blood_gse19491":           "TB blood (GSE19491)",
    "ipf_lung_gse47460":           "IPF lung (GSE47460)",
    "copd_lung_gse47460":          "COPD lung (GSE47460)",
    "influenza_pbmc_gse101702":    "Influenza PBMC (GSE101702)",
    "trauma_blood_gse36809":       "Trauma blood (GSE36809)",
    "hcc_liver_gse6764":           "HCC liver (GSE6764)",
    "crohn_intestine_gse112366":   "Crohn intestine (GSE112366)",
    "rsv_blood_gse34205":          "RSV blood (GSE34205)",
    "viral_challenge_gse73072":    "Viral challenge (GSE73072)",
    "breast_relapse_gse1456":      "Breast relapse (GSE1456)",
    "hypoxia_mcf7_gse3188":        "Hypoxia MCF7 (GSE3188)",
    "hypoxia_timecourse_gse47533": "Hypoxia timecourse (GSE47533)",
    "hypoxia_multicell_gse18494":  "Hypoxia multicell (GSE18494)",
    "emt_tgfb_gse17708":           "EMT TGF-β (GSE17708)",
    "emt_hmle_gse24202":           "EMT HMLE (GSE24202)",
    "emt_mammary_gse43495":        "EMT mammary (GSE43495)",
    "sepsis_shock_gse95233":       "Sepsis shock (GSE95233)",
    "melioidosis_blood_gse69528":  "Melioidosis blood (GSE69528)",
    "lung_tumor_gse19188":         "Lung tumor (GSE19188)",
    "influenza_challenge_gse68310":"Influenza challenge (GSE68310)",
    "ccrcc_kidney_gse36895":       "ccRCC kidney (GSE36895)",
    "gbm_brain_gse4290":           "GBM brain (GSE4290)",
    "emt_arpe19_gse12548":         "EMT ARPE-19 (GSE12548)",
    "ipf_lung_gse53845":           "IPF lung (GSE53845)",
    "influenza_severe_gse111368":  "Influenza severe (GSE111368)",
}

# ---------------------------------------------------------------------------
# Compute pooled effects via DL random-effects
# ---------------------------------------------------------------------------
def dl_pooled(df):
    ds = df.cohens_d.values.tolist()
    vs = df.cohens_d_var.values.tolist()
    res = dl_random_effects_meta(ds, vs)
    se = res["pooled_se"]
    mu = res["pooled_effect"]
    i2 = res["i_squared"]
    return mu, mu - 1.96 * se, mu + 1.96 * se, i2

# All 30 cohorts
mu_all, lo_all, hi_all, i2_all = dl_pooled(ifng)

# Within interferon program
ifng_ifn = ifng[ifng.biological_program == "interferon"].copy()
mu_ifn, lo_ifn, hi_ifn, i2_ifn = dl_pooled(ifng_ifn)

# ---------------------------------------------------------------------------
# Sort cohorts for display
# Order: by biological_program (so colors cluster), then by effect size within program
# For the ALL panel we want a visually compelling order — mixed programs, sorted by d
# For the WITHIN panel just sort by d ascending
# ---------------------------------------------------------------------------
PROG_ORDER = ["interferon", "inflammation", "proliferation", "hypoxia", "emt"]

ifng_sorted_all = (
    ifng
    .assign(prog_rank=ifng.biological_program.map({p: i for i, p in enumerate(PROG_ORDER)}).fillna(99))
    .sort_values(["prog_rank", "cohens_d"])
    .reset_index(drop=True)
)

ifng_sorted_within = ifng_ifn.sort_values("cohens_d").reset_index(drop=True)

# ---------------------------------------------------------------------------
# Diamond helper: draws a filled diamond at the pooled estimate row
# ---------------------------------------------------------------------------
def draw_diamond(ax, x_center, x_lo, x_hi, y_center, height=0.45, color="#444444", zorder=5):
    """Draw a horizontal diamond (rhombus) centered at (x_center, y_center)."""
    diamond_x = [x_lo, x_center, x_hi, x_center, x_lo]
    diamond_y = [y_center, y_center + height / 2, y_center, y_center - height / 2, y_center]
    ax.fill(diamond_x, diamond_y, color=color, zorder=zorder, clip_on=True)
    ax.plot(diamond_x, diamond_y, color=color, linewidth=0.8, zorder=zorder + 1, clip_on=True)

# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "lines.linewidth": 1.0,
    "pdf.fonttype": 42,   # embed fonts as TrueType for journal submission
    "ps.fonttype": 42,
})

N_ALL = len(ifng_sorted_all)    # 30
N_IFN = len(ifng_sorted_within)  # 6

# Heights: each cohort row = 1 unit, plus 2 rows gap+diamond at bottom
ROW_H_ALL = 0.24   # inches per row
ROW_H_IFN = 0.50   # taller rows for the 6-cohort panel

# Extra space at bottom of each panel: separator + diamond row + xlabel
BOTTOM_EXTRA = 1.5  # row units

HEIGHT_ALL = (N_ALL + BOTTOM_EXTRA) * ROW_H_ALL + 0.6   # +0.6 for title
HEIGHT_IFN = (N_IFN + BOTTOM_EXTRA) * ROW_H_IFN + 0.6

FIG_HEIGHT = max(HEIGHT_ALL, HEIGHT_IFN) + 0.8  # common height, title space
FIG_WIDTH = 16.0   # inches (two panels)

fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT), facecolor="white")

# Leave room at top for suptitle, and bottom margin
TOP_MARGIN = 0.88
BOTTOM_MARGIN = 0.09
LEFT_MARGIN = 0.02
MID_GAP = 0.04
RIGHT_MARGIN = 0.01

panel_width = (1.0 - LEFT_MARGIN - MID_GAP - RIGHT_MARGIN) / 2

ax_left = fig.add_axes([LEFT_MARGIN, BOTTOM_MARGIN,
                         panel_width, TOP_MARGIN - BOTTOM_MARGIN])
ax_right = fig.add_axes([LEFT_MARGIN + panel_width + MID_GAP, BOTTOM_MARGIN,
                          panel_width, TOP_MARGIN - BOTTOM_MARGIN])

# ---------------------------------------------------------------------------
# Shared helper to draw one forest panel
# ---------------------------------------------------------------------------
def draw_forest_panel(
    ax, df_sorted, mu, lo, hi, i2,
    panel_label, pooled_color,
    all_blue=False,    # if True, override all point colors to blue
    x_range=(-4.5, 4.5),
):
    n = len(df_sorted)

    # y coords: top cohort at n-1, bottom cohort at 0, then gap, then diamond
    y_cohorts = np.arange(n - 1, -1, -1, dtype=float)  # n-1 down to 0
    y_diamond = -1.5  # below the separator line

    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Vertical reference line at 0
    ax.axvline(0, color="#cccccc", linewidth=0.9, zorder=1, linestyle="--", dashes=(4, 3))

    # Plot each cohort
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        y = y_cohorts[i]
        prog = row.biological_program
        color = PROGRAM_COLORS.get(prog, "#888888") if not all_blue else PROGRAM_COLORS["interferon"]

        # CI whisker
        ax.plot([row.ci_lo, row.ci_hi], [y, y],
                color=color, linewidth=1.1, zorder=3, solid_capstyle="round")

        # Study point — size proportional to inverse variance (precision)
        weight = 1.0 / (row.cohens_d_var + 0.01)
        # Clamp size range for visual clarity
        pt_size = np.clip(weight * 0.8, 15, 120)
        ax.scatter(row.cohens_d, y, s=pt_size, color=color, zorder=4,
                   edgecolors="white", linewidths=0.5)

    # Separator line above diamond
    ax.axhline(-0.75, color="#aaaaaa", linewidth=0.7, zorder=2, xmin=0.05, xmax=0.95)

    # Diamond for pooled estimate
    draw_diamond(ax, mu, lo, hi, y_diamond, height=0.55, color=pooled_color)

    # I² annotation on diamond row — prominently to the right of diamond
    i2_pct = i2 * 100
    xspan = x_range[1] - x_range[0]
    ax.text(x_range[0] + xspan * 0.88, y_diamond,
            f"I² = {i2_pct:.0f}%",
            ha="right", va="center", fontsize=9, color="#222222",
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="#f0f0f0",
                      edgecolor="#bbbbbb", linewidth=0.7))

    # Pooled g annotation — left of center axis
    ax.text(x_range[0] + xspan * 0.02, y_diamond,
            f"g = {mu:.3f}  [{lo:.3f}, {hi:.3f}]",
            ha="left", va="center", fontsize=8, color="#222222",
            style="italic")

    # y-axis labels (cohort names)
    ytick_positions = list(y_cohorts) + [y_diamond]
    ytick_labels = [LABEL_MAP.get(row.cohort_id, row.cohort_id)
                    for _, row in df_sorted.iterrows()] + ["Pooled (RE)"]

    ax.set_yticks(ytick_positions)
    ax.set_yticklabels(ytick_labels, fontsize=7.8)
    ax.tick_params(axis="y", length=0, pad=4)

    # x-axis
    ax.set_xlim(x_range)
    ax.set_xlabel("Hedges' g (effect size)", fontsize=9, labelpad=6)
    xticks = np.arange(int(x_range[0]), int(x_range[1]) + 1, 1)
    ax.set_xticks(xticks)
    ax.xaxis.set_tick_params(length=3)

    # y limits
    ax.set_ylim(y_diamond - 1.0, n - 0.3)

    # Panel title
    ax.set_title(panel_label, fontsize=10, fontweight="bold", pad=8, loc="left")

    return ax


# ---------------------------------------------------------------------------
# Draw left panel: all 30 cohorts
# ---------------------------------------------------------------------------
draw_forest_panel(
    ax=ax_left,
    df_sorted=ifng_sorted_all,
    mu=mu_all,
    lo=lo_all,
    hi=hi_all,
    i2=i2_all,
    panel_label="(A)  Cross-context evaluation  (all 30 cohorts)",
    pooled_color=POOLED_COLOR_ALL,
    all_blue=False,
    x_range=(-5.0, 5.0),
)

# Program legend on left panel
legend_handles = [
    mpatches.Patch(color=PROGRAM_COLORS["interferon"],    label="Interferon"),
    mpatches.Patch(color=PROGRAM_COLORS["inflammation"],  label="Inflammation"),
    mpatches.Patch(color=PROGRAM_COLORS["proliferation"], label="Proliferation"),
    mpatches.Patch(color=PROGRAM_COLORS["hypoxia"],       label="Hypoxia"),
    mpatches.Patch(color=PROGRAM_COLORS["emt"],           label="EMT"),
]
ax_left.legend(
    handles=legend_handles,
    loc="lower right",
    fontsize=7.5,
    frameon=True,
    framealpha=0.92,
    edgecolor="#cccccc",
    ncol=1,
    title="Biological program",
    title_fontsize=7.5,
)

# ---------------------------------------------------------------------------
# Draw right panel: 6 interferon cohorts
# ---------------------------------------------------------------------------
draw_forest_panel(
    ax=ax_right,
    df_sorted=ifng_sorted_within,
    mu=mu_ifn,
    lo=lo_ifn,
    hi=hi_ifn,
    i2=i2_ifn,
    panel_label="(B)  Within-program evaluation  (6 interferon cohorts)",
    pooled_color=POOLED_COLOR_WITHIN,
    all_blue=True,
    x_range=(-0.3, 2.3),
)

# ---------------------------------------------------------------------------
# Suptitle
# ---------------------------------------------------------------------------
fig.text(
    0.5, 0.955,
    "IFN-γ: same signature, different evaluation framework",
    ha="center", va="top",
    fontsize=12, fontweight="bold", color="#111111",
)
fig.text(
    0.5, 0.935,
    "hallmark_ifng_response  |  DerSimonian-Laird random-effects meta-analysis  |  95% CI",
    ha="center", va="top",
    fontsize=8, color="#555555",
)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
os.makedirs(os.path.join(ROOT, "paper"), exist_ok=True)

plt.savefig(OUT_PDF, bbox_inches="tight", dpi=300, facecolor="white")
plt.savefig(OUT_PNG, bbox_inches="tight", dpi=300, facecolor="white")
print(f"Saved:\n  {OUT_PDF}\n  {OUT_PNG}")

# Report computed pooled values
print(f"\nCross-context (all {len(ifng)} cohorts): g = {mu_all:.4f}  [{lo_all:.4f}, {hi_all:.4f}]  I² = {i2_all:.4f}")
print(f"Within-program ({len(ifng_ifn)} interferon cohorts): g = {mu_ifn:.4f}  [{lo_ifn:.4f}, {hi_ifn:.4f}]  I² = {i2_ifn:.4f}")
