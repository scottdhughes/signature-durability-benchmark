"""Cross-cohort meta-analysis: fixed-effect pooling and heterogeneity."""
from __future__ import annotations
from typing import Any
import numpy as np
from scipy import stats

def fixed_effect_meta(effects: list[float], variances: list[float]) -> dict[str, Any]:
    effects_arr = np.array(effects, dtype=float)
    variances_arr = np.array(variances, dtype=float)
    variances_arr = np.where(variances_arr < 1e-10, 1e-10, variances_arr)
    weights = 1.0 / variances_arr
    pooled = float(np.sum(weights * effects_arr) / np.sum(weights))
    pooled_se = float(1.0 / np.sqrt(np.sum(weights)))
    z = pooled / pooled_se if pooled_se > 0 else 0.0
    p = float(2 * (1 - stats.norm.cdf(abs(z))))
    return {"pooled_effect": pooled, "pooled_se": pooled_se, "z_score": float(z), "pooled_p": p}

def i_squared(effects: list[float], variances: list[float]) -> float:
    if len(effects) < 2:
        return 0.0
    meta = fixed_effect_meta(effects, variances)
    effects_arr = np.array(effects)
    variances_arr = np.where(np.array(variances) < 1e-10, 1e-10, np.array(variances))
    weights = 1.0 / variances_arr
    q = float(np.sum(weights * (effects_arr - meta["pooled_effect"]) ** 2))
    df = len(effects) - 1
    return max(0.0, float((q - df) / q)) if q > 0 else 0.0

def leave_one_out(effects: list[float], variances: list[float]) -> list[dict[str, Any]]:
    results = []
    for i in range(len(effects)):
        sub_effects = effects[:i] + effects[i + 1:]
        sub_variances = variances[:i] + variances[i + 1:]
        if len(sub_effects) < 1:
            continue
        meta = fixed_effect_meta(sub_effects, sub_variances)
        results.append({"dropped_index": i, **meta})
    return results

def loo_stability(effects: list[float], variances: list[float], direction_threshold: float = 0.0) -> float:
    """Fraction of LOO iterations where pooled effect stays positive."""
    loo = leave_one_out(effects, variances)
    if not loo:
        return 0.0
    consistent = sum(1 for r in loo if r["pooled_effect"] > direction_threshold)
    return float(consistent / len(loo))

def platform_holdout(effects: list[float], variances: list[float], platforms: list[str]) -> dict[str, Any]:
    unique = sorted(set(platforms))
    results = {}
    for platform in unique:
        idx = [i for i, p in enumerate(platforms) if p != platform]
        if len(idx) < 1:
            continue
        sub_effects = [effects[i] for i in idx]
        sub_variances = [variances[i] for i in idx]
        meta = fixed_effect_meta(sub_effects, sub_variances)
        results[platform] = meta
    return results

def platform_holdout_consistency(effects: list[float], variances: list[float], platforms: list[str]) -> float:
    """Fraction of platform-holdout iterations where pooled effect stays positive."""
    holdout = platform_holdout(effects, variances, platforms)
    if not holdout:
        return 0.0
    consistent = sum(1 for r in holdout.values() if r["pooled_effect"] > 0)
    return float(consistent / len(holdout))

def direction_consistency(effects: list[float]) -> float:
    """Fraction of cohort effects that are positive."""
    if not effects:
        return 0.0
    return float(sum(1 for e in effects if e > 0) / len(effects))
