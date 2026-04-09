"""Certificate builders."""
from __future__ import annotations
from typing import Any
import pandas as pd
from .constants import DURABLE_CLASSES

def _verdict(score: float) -> str:
    if score >= 0.80:
        return "passed"
    if score >= 0.50:
        return "mixed"
    return "failed"

def generate_durability_certificate(scores_df: pd.DataFrame, aggregate: dict[str, Any]) -> dict[str, Any]:
    """Direction consistency, aggregate effect, LOO stability on durable subset."""
    fm_primary = scores_df[(scores_df["model"] == "full_model") & (scores_df["split"] == "primary")]
    durable_rows = fm_primary[fm_primary["expected_class"].isin(DURABLE_CLASSES)]
    if durable_rows.empty:
        return {"certificate_type": "durability", "verdict": "failed", "summary": "No durable signatures", "metrics": {}}

    mean_dir = float(durable_rows["direction_consistency"].mean())
    mean_loo = float(durable_rows["loo_stability"].mean())
    mean_effect = float(durable_rows["aggregate_effect"].abs().mean())
    score = (mean_dir + mean_loo) / 2.0

    return {
        "certificate_type": "durability",
        "verdict": _verdict(score),
        "summary": f"Mean direction consistency {mean_dir:.3f}, LOO stability {mean_loo:.3f}, |effect| {mean_effect:.3f}",
        "metrics": {"mean_direction_consistency": mean_dir, "mean_loo_stability": mean_loo, "mean_abs_effect": mean_effect, "composite_score": score},
    }

def generate_platform_transfer_certificate(scores_df: pd.DataFrame, aggregate: dict[str, Any]) -> dict[str, Any]:
    """Platform holdout consistency across models."""
    fm_primary = scores_df[(scores_df["model"] == "full_model") & (scores_df["split"] == "primary")]
    fm_primary = fm_primary[fm_primary["expected_class"] != "insufficient_coverage"]
    if fm_primary.empty or "platform_holdout_consistency" not in fm_primary.columns:
        return {"certificate_type": "platform_transfer", "verdict": "failed", "summary": "No data", "metrics": {}}

    mean_plat = float(fm_primary["platform_holdout_consistency"].mean())
    return {
        "certificate_type": "platform_transfer",
        "verdict": _verdict(mean_plat),
        "summary": f"Mean platform holdout consistency {mean_plat:.3f}",
        "metrics": {"mean_platform_holdout_consistency": mean_plat},
    }

def generate_confounder_rejection_certificate(scores_df: pd.DataFrame, aggregate: dict[str, Any]) -> dict[str, Any]:
    """Confounded subset correctly identified."""
    fm_primary = scores_df[(scores_df["model"] == "full_model") & (scores_df["split"] == "primary")]
    confounded_rows = fm_primary[fm_primary["expected_class"] == "confounded"]
    if confounded_rows.empty:
        return {"certificate_type": "confounder_rejection", "verdict": "failed", "summary": "No confounded signatures", "metrics": {}}

    correct = float((confounded_rows["predicted_class"] == "confounded").mean())
    return {
        "certificate_type": "confounder_rejection",
        "verdict": _verdict(correct),
        "summary": f"Confounded rejection accuracy {correct:.3f}",
        "metrics": {"confounded_rejection_accuracy": correct},
    }

def generate_coverage_certificate(scores_df: pd.DataFrame, aggregate: dict[str, Any]) -> dict[str, Any]:
    """Coverage sentinels correctly identified."""
    fm_primary = scores_df[(scores_df["model"] == "full_model") & (scores_df["split"] == "primary")]
    lowcov_rows = fm_primary[fm_primary["expected_class"] == "insufficient_coverage"]
    if lowcov_rows.empty:
        return {"certificate_type": "coverage", "verdict": "passed", "summary": "No coverage sentinels", "metrics": {}}

    correct = float((lowcov_rows["predicted_class"] == "insufficient_coverage").mean())
    mean_cov = float(fm_primary["mean_coverage"].mean())
    return {
        "certificate_type": "coverage",
        "verdict": _verdict(correct),
        "summary": f"Coverage sentinel accuracy {correct:.3f}, mean coverage {mean_cov:.3f}",
        "metrics": {"coverage_sentinel_accuracy": correct, "mean_coverage": mean_cov},
    }
