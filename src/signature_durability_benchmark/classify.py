"""Classification rules for 6 models."""
from __future__ import annotations
from typing import Any

THRESHOLDS = {
    "coverage_min": 0.60,
    "direction_consistency_min": 0.50,
    "aggregate_p_max": 0.10,
    "i_squared_max": 0.75,
    "loo_stability_min": 0.60,
    "null_separation_p_max": 0.10,
}

def classify_signature(profile: dict[str, Any], model_name: str) -> dict[str, Any]:
    t = THRESHOLDS
    coverage = profile.get("mean_coverage", 0.0)
    aggregate_effect = abs(profile.get("aggregate_effect", 0.0))
    aggregate_p = profile.get("aggregate_p", 1.0)
    dir_consistency = profile.get("direction_consistency", 0.0)
    i_sq = profile.get("i_squared", 0.0)
    loo = profile.get("loo_stability", 0.0)
    max_conf = abs(profile.get("max_confounder_effect", 0.0))
    null_p = profile.get("null_separation_p", 1.0)
    within_p = profile.get("within_p", None)
    within_d = profile.get("within_d", None)

    if coverage < t["coverage_min"]:
        predicted = "insufficient_coverage"
    elif model_name == "overlap_only":
        if dir_consistency < t["direction_consistency_min"]:
            predicted = "brittle"
        else:
            predicted = "durable"
    elif model_name == "effect_only":
        if aggregate_p > t["aggregate_p_max"] or dir_consistency < t["direction_consistency_min"]:
            predicted = "brittle"
        else:
            predicted = "durable"
    elif model_name == "null_aware":
        if null_p > t["null_separation_p_max"]:
            predicted = "brittle"
        elif dir_consistency < t["direction_consistency_min"]:
            predicted = "brittle"
        elif i_sq > t["i_squared_max"] or loo < t["loo_stability_min"]:
            predicted = "mixed"
        else:
            predicted = "durable"
    elif model_name == "no_confounder":
        if aggregate_p > t["aggregate_p_max"] or dir_consistency < t["direction_consistency_min"]:
            predicted = "brittle"
        elif i_sq > t["i_squared_max"] or loo < t["loo_stability_min"]:
            predicted = "mixed"
        else:
            predicted = "durable"
    elif model_name == "within_program":
        # Confounder check first (same as full model)
        if max_conf >= aggregate_effect:
            predicted = "confounded"
        # Within-program check: if signature has within-program data and is significant
        elif within_p is not None and within_p < 0.05 and abs(within_d) > 0.2:
            predicted = "durable"
        # Fall back to cross-context checks
        elif aggregate_p > t["aggregate_p_max"] or dir_consistency < t["direction_consistency_min"]:
            predicted = "brittle"
        elif i_sq > t["i_squared_max"] or loo < t["loo_stability_min"]:
            predicted = "mixed"
        else:
            predicted = "durable"
    else:  # full_model
        if max_conf >= aggregate_effect:
            predicted = "confounded"
        elif aggregate_p > t["aggregate_p_max"] or dir_consistency < t["direction_consistency_min"]:
            predicted = "brittle"
        elif i_sq > t["i_squared_max"] or loo < t["loo_stability_min"]:
            predicted = "mixed"
        else:
            predicted = "durable"

    return {
        "model": model_name,
        "predicted_class": predicted,
        "aggregate_effect": float(aggregate_effect),
        "max_confounder_effect": float(max_conf),
        "direction_consistency": float(dir_consistency),
        "i_squared": float(i_sq),
        "loo_stability": float(loo),
        "coverage": float(coverage),
    }
