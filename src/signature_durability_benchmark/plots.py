"""Visualization functions."""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def forest_plot(per_cohort_effects: pd.DataFrame, out_path: str | Path) -> None:
    """Forest plot of per-cohort effect sizes for each signature."""
    out_path = Path(out_path)
    sigs = per_cohort_effects["signature_id"].unique()
    n_sigs = len(sigs)
    fig, ax = plt.subplots(figsize=(12, max(6, n_sigs * 0.4)))

    y_positions = []
    y_labels = []
    for i, sig in enumerate(sorted(sigs)):
        subset = per_cohort_effects[per_cohort_effects["signature_id"] == sig]
        effects = subset["cohens_d"].values
        y_pos = n_sigs - i
        ax.scatter(effects, [y_pos] * len(effects), alpha=0.5, s=20, color="steelblue")
        mean_eff = np.mean(effects)
        ax.plot(mean_eff, y_pos, "D", color="darkred", markersize=8, zorder=5)
        y_positions.append(y_pos)
        y_labels.append(sig[:35])

    ax.axvline(0, color="gray", linestyle="--", alpha=0.5)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=7)
    ax.set_xlabel("Cohen's d")
    ax.set_title("Per-Cohort Effect Sizes by Signature")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def null_separation_plot(null_summary: pd.DataFrame, out_path: str | Path) -> None:
    """Observed vs null effect distribution."""
    out_path = Path(out_path)
    fig, ax = plt.subplots(figsize=(10, 6))

    if "null_mean" in null_summary.columns and "observed_effect" in null_summary.columns:
        null_means = null_summary["null_mean"].values
        null_sds = null_summary["null_sd"].values
        observed = null_summary["observed_effect"].values
        sigs = null_summary["signature_id"].values

        x = np.arange(len(sigs))
        ax.bar(x, null_means, yerr=null_sds, alpha=0.3, color="gray", label="Null mean ± SD")
        ax.scatter(x, observed, color="darkred", s=40, zorder=5, label="Observed")
        ax.set_xticks(x)
        ax.set_xticklabels([s[:20] for s in sigs], rotation=90, fontsize=6)
        ax.axhline(0, color="black", linestyle="-", alpha=0.3)
        ax.set_ylabel("Pooled Effect Size")
        ax.set_title("Observed vs Null Distribution")
        ax.legend()

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def stability_heatmap(loo_df: pd.DataFrame, out_path: str | Path) -> None:
    """LOO stability heatmap: signatures × dropped cohorts."""
    out_path = Path(out_path)

    if loo_df.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No LOO data", ha="center", va="center")
        ax.set_axis_off()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return

    pivot = loo_df.pivot_table(index="signature_id", columns="dropped_cohort", values="pooled_effect", aggfunc="first")
    fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns) * 0.8), max(6, len(pivot) * 0.4)))

    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([c[:15] for c in pivot.columns], rotation=90, fontsize=6)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([s[:30] for s in pivot.index], fontsize=6)
    ax.set_title("Leave-One-Cohort-Out Pooled Effects")
    plt.colorbar(im, ax=ax, label="Pooled Effect")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def platform_transfer_panel(platform_holdout_df: pd.DataFrame, out_path: str | Path) -> None:
    """Effect sizes split by held-out platform."""
    out_path = Path(out_path)

    if platform_holdout_df.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No platform holdout data", ha="center", va="center")
        ax.set_axis_off()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return

    platforms = platform_holdout_df["held_out_platform"].unique()
    sigs = platform_holdout_df["signature_id"].unique()

    fig, axes = plt.subplots(1, len(platforms), figsize=(5 * len(platforms), max(6, len(sigs) * 0.3)), sharey=True)
    if len(platforms) == 1:
        axes = [axes]

    for ax_idx, plat in enumerate(sorted(platforms)):
        subset = platform_holdout_df[platform_holdout_df["held_out_platform"] == plat]
        y_pos = range(len(subset))
        ax = axes[ax_idx]
        ax.barh(list(y_pos), subset["pooled_effect"].values, color="steelblue", alpha=0.7)
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels([s[:25] for s in subset["signature_id"].values], fontsize=6)
        ax.axvline(0, color="gray", linestyle="--", alpha=0.5)
        ax.set_title(f"Held out: {plat}", fontsize=9)
        ax.set_xlabel("Pooled Effect")

    plt.suptitle("Platform Holdout Analysis", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
