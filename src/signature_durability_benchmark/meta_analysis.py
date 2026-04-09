"""Cross-cohort meta-analysis: fixed-effect pooling, DerSimonian-Laird random-effects, and heterogeneity."""
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
    """Fraction of LOO iterations where pooled effect retains the same sign as the full pooled estimate."""
    if len(effects) < 2:
        return 0.0
    full_meta = fixed_effect_meta(effects, variances)
    full_sign = 1 if full_meta["pooled_effect"] > 0 else -1
    loo = leave_one_out(effects, variances)
    if not loo:
        return 0.0
    consistent = sum(1 for r in loo if (r["pooled_effect"] > 0) == (full_sign > 0))
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
    """Fraction of platform-holdout iterations where pooled effect retains full-data sign."""
    full_meta = fixed_effect_meta(effects, variances)
    full_sign = 1 if full_meta["pooled_effect"] > 0 else -1
    holdout = platform_holdout(effects, variances, platforms)
    if not holdout:
        return 0.0
    consistent = sum(1 for r in holdout.values() if (r["pooled_effect"] > 0) == (full_sign > 0))
    return float(consistent / len(holdout))

def direction_consistency(effects: list[float]) -> float:
    """Fraction of cohort effects matching the majority direction."""
    if not effects:
        return 0.0
    n_pos = sum(1 for e in effects if e > 0)
    n_neg = len(effects) - n_pos
    return float(max(n_pos, n_neg) / len(effects))


def dl_random_effects_meta(effects: list[float], variances: list[float]) -> dict[str, Any]:
    """DerSimonian-Laird random-effects meta-analysis."""
    effects_arr = np.array(effects, dtype=float)
    variances_arr = np.where(np.array(variances) < 1e-10, 1e-10, np.array(variances))

    # Fixed-effect weights for Q statistic
    w_fe = 1.0 / variances_arr
    pooled_fe = np.sum(w_fe * effects_arr) / np.sum(w_fe)

    # Q statistic and tau-squared
    Q = float(np.sum(w_fe * (effects_arr - pooled_fe) ** 2))
    k = len(effects)
    df = k - 1
    C = float(np.sum(w_fe) - np.sum(w_fe**2) / np.sum(w_fe))
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0

    # Random-effects weights
    w_re = 1.0 / (variances_arr + tau2)
    pooled_re = float(np.sum(w_re * effects_arr) / np.sum(w_re))
    se_re = float(1.0 / np.sqrt(np.sum(w_re)))

    z = pooled_re / se_re if se_re > 0 else 0.0
    p = float(2 * (1 - stats.norm.cdf(abs(z))))

    # I-squared
    i2 = max(0.0, (Q - df) / Q) if Q > 0 and df > 0 else 0.0

    return {
        "pooled_effect": pooled_re,
        "pooled_se": se_re,
        "pooled_p": p,
        "tau2": tau2,
        "i_squared": i2,
        "Q": Q,
        "k": k,
        "method": "DL_random_effects",
    }


def within_program_meta(
    effects: list[float],
    variances: list[float],
    cohort_programs: list[str],
    signature_program: str,
) -> dict[str, Any]:
    """Compute within-program and outside-program meta-analysis separately.

    Uses DerSimonian-Laird random-effects for both subsets, which is the
    appropriate method given observed heterogeneity within programs.
    """
    within_idx = [i for i, p in enumerate(cohort_programs) if p == signature_program]
    outside_idx = [i for i, p in enumerate(cohort_programs) if p != signature_program]

    within_effects = [effects[i] for i in within_idx]
    within_vars = [variances[i] for i in within_idx]
    outside_effects = [effects[i] for i in outside_idx]
    outside_vars = [variances[i] for i in outside_idx]

    within_result = dl_random_effects_meta(within_effects, within_vars) if len(within_effects) >= 2 else None
    outside_result = dl_random_effects_meta(outside_effects, outside_vars) if len(outside_effects) >= 2 else None

    return {
        "within": within_result,
        "outside": outside_result,
        "within_k": len(within_effects),
        "outside_k": len(outside_effects),
    }
