"""Data-backed diagnostic workflow figure generation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from .utils import read_json


def build_workflow_figure_data(root: Path) -> dict[str, Any]:
    triage = read_json(root / "outputs" / "triage_ifng" / "diagnostic.json")
    failure = read_json(root / "outputs" / "canonical_v8" / "failure_mode_analysis.json")
    prospective_v1 = read_json(root / "outputs" / "canonical_v8" / "prospective_holdout_validation.json")
    prospective_v2 = read_json(
        root
        / "data"
        / "prospective_holdout"
        / "external_timestamps"
        / "prospective_holdout_v2"
        / "declaration_receipt.json"
    )
    external_transfer = read_json(root / "outputs" / "canonical_v8" / "external_rnaseq_validation.json")

    ifn_focus = triage["program_diagnostics"]["interferon"]
    inflam_fail = failure["ambiguous_program_examples"]["hallmark_inflammatory_response"]
    external_ifn = external_transfer["primary_bulk_rnaseq_panel"]["ifn_focus_summary"]

    return {
        "input_signature": "30-gene IFN-gamma core",
        "frozen_panel": "35 frozen cohorts\n5 programs / 14 platforms",
        "within_outside": (
            f"within IFN g={ifn_focus['within']['pooled_effect']:.3f}\n"
            f"outside IFN g={ifn_focus['outside']['pooled_effect']:.3f}"
        ),
        "permutation": (
            f"Q_B/Q_tot={triage['permutation_program_structure']['observed_Q_B_fraction']:.3f}\n"
            f"perm p={triage['permutation_program_structure']['p_value']:.3g}"
        ),
        "durable_branch": (
            f"best program: {triage['inferred_program']['program']}\n"
            f"within-program: {triage['classification']['within_program']}\n"
            f"full model: {triage['classification']['full_model']}"
        ),
        "failure_branch": (
            "coarse-label / failure-mode\n"
            f"inflammatory largest LOPO drop when hiding {inflam_fail['largest_drop_when_hiding']}\n"
            f"ΔQ_B/Q_tot={inflam_fail['largest_drop_value']:.3f}"
        ),
        "external_transfer": (
            "held-out external transfer\n"
            f"4 IFN signatures positive in all {external_transfer['primary_bulk_rnaseq_panel']['n_cohorts']} bulk RNA-seq cohorts"
        ),
        "prospective_v1": (
            "prospective v1 scored challenge\n"
            f"overall_success={prospective_v1['prediction_summary']['overall_success']}"
        ),
        "prospective_v2": (
            "prospective v2 declared future challenge\n"
            f"RFC3161 TSA reply: {prospective_v2['tsa_reply_time']}"
        ),
    }


def _box(ax: plt.Axes, xy: tuple[float, float], text: str, *, width: float = 0.2, height: float = 0.12, color: str = "#edf2f7") -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.015,rounding_size=0.02",
        linewidth=1.2,
        facecolor=color,
        edgecolor="#1f2937",
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=10)


def _arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=1.3,
            color="#374151",
        )
    )


def render_workflow_figure(data: dict[str, Any], pdf_path: Path, png_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    _box(ax, (0.04, 0.78), f"1. Input Signature\n{data['input_signature']}", color="#dbeafe")
    _box(ax, (0.28, 0.78), f"2. Score Across Frozen Cohorts\n{data['frozen_panel']}", color="#e0f2fe")
    _box(ax, (0.52, 0.78), f"3. Within / Outside Split\n{data['within_outside']}", color="#dcfce7")
    _box(ax, (0.76, 0.78), f"4. Permutation Check\n{data['permutation']}", color="#fef3c7")

    _box(ax, (0.20, 0.48), f"5a. Durable-Within-Context\n{data['durable_branch']}", width=0.26, color="#d1fae5")
    _box(ax, (0.56, 0.48), f"5b. Coarse-Label / Failure Mode\n{data['failure_branch']}", width=0.30, color="#fee2e2")

    _box(ax, (0.08, 0.17), f"6a. Audit Branch\n{data['external_transfer']}", width=0.24, color="#ede9fe")
    _box(ax, (0.38, 0.17), f"6b. Audit Branch\n{data['prospective_v1']}", width=0.24, color="#ede9fe")
    _box(ax, (0.68, 0.17), f"6c. Audit Branch\n{data['prospective_v2']}", width=0.24, color="#ede9fe")

    _arrow(ax, (0.24, 0.84), (0.28, 0.84))
    _arrow(ax, (0.48, 0.84), (0.52, 0.84))
    _arrow(ax, (0.72, 0.84), (0.76, 0.84))
    _arrow(ax, (0.80, 0.78), (0.73, 0.60))
    _arrow(ax, (0.80, 0.78), (0.40, 0.60))
    _arrow(ax, (0.33, 0.48), (0.20, 0.29))
    _arrow(ax, (0.33, 0.48), (0.50, 0.29))
    _arrow(ax, (0.71, 0.48), (0.80, 0.29))

    ax.text(0.5, 0.95, "Program-Conditioned Diagnostic Workflow", ha="center", va="center", fontsize=16, fontweight="bold")
    ax.text(
        0.5,
        0.91,
        "Input signature -> frozen benchmark structure -> durable-vs-ambiguous interpretation -> external and prospective audit layers",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#374151",
    )

    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=180)
    plt.close(fig)
