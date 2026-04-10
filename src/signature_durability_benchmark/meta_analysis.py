"""Cross-cohort meta-analysis: pooling, heterogeneity, and program-structure diagnostics."""
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
    p = float(2 * stats.norm.sf(abs(z)))
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
    p = float(2 * stats.norm.sf(abs(z)))

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


def guarded_hksj_random_effects_meta(effects: list[float], variances: list[float]) -> dict[str, Any]:
    """DerSimonian-Laird random-effects meta-analysis with guarded HKSJ inference.

    The pooled effect estimate is the same DL random-effects estimate used in the
    canonical v8 analysis scripts. The primary p-value is computed from the
    Hartung-Knapp-Sidik-Jonkman standard error with the guarded variant
    `max(SE_DL, SE_HKSJ)`.
    """
    dl = dl_random_effects_meta(effects, variances)
    k = dl["k"]
    if k < 2:
        return dl

    effects_arr = np.array(effects, dtype=float)
    variances_arr = np.where(np.array(variances) < 1e-10, 1e-10, np.array(variances))
    tau2 = float(dl["tau2"])
    pooled = float(dl["pooled_effect"])

    w_re = 1.0 / (variances_arr + tau2)
    q_star = float(np.sum(w_re * (effects_arr - pooled) ** 2))
    df = k - 1
    se_hksj = float(np.sqrt(q_star / (df * np.sum(w_re)))) if df > 0 and np.sum(w_re) > 0 else float(dl["pooled_se"])
    se_guarded = max(float(dl["pooled_se"]), se_hksj)
    t_hksj = pooled / se_hksj if se_hksj > 0 else 0.0
    p_hksj = float(2 * stats.t.sf(abs(t_hksj), df=df)) if df > 0 else float(dl["pooled_p"])
    t_guarded = pooled / se_guarded if se_guarded > 0 else 0.0
    p_guarded = float(2 * stats.t.sf(abs(t_guarded), df=df)) if df > 0 else float(dl["pooled_p"])

    return {
        **dl,
        "hksj_se": se_hksj,
        "hksj_p": p_hksj,
        "guarded_se": se_guarded,
        "guarded_p": p_guarded,
        "pooled_se": se_guarded,
        "pooled_p": p_guarded,
        "method": "DL_HKSJ_guarded_random_effects",
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

    within_result = guarded_hksj_random_effects_meta(within_effects, within_vars) if len(within_effects) >= 2 else None
    outside_result = guarded_hksj_random_effects_meta(outside_effects, outside_vars) if len(outside_effects) >= 2 else None

    return {
        "within": within_result,
        "outside": outside_result,
        "within_k": len(within_effects),
        "outside_k": len(outside_effects),
    }


def q_decomposition(
    effects: list[float],
    variances: list[float],
    groups: list[str],
) -> dict[str, Any]:
    """Decompose weighted heterogeneity into within-group and between-group Q.

    This is the standard subgroup Q partition:
      Q_total = Q_within + Q_between
    """
    if len(effects) != len(variances) or len(effects) != len(groups):
        raise ValueError("effects, variances, and groups must have the same length")
    if len(effects) < 2:
        return {
            "Q_total": 0.0,
            "Q_within": 0.0,
            "Q_between": 0.0,
            "Q_B_fraction": 0.0,
            "I2_total": 0.0,
            "I2_within": 0.0,
            "p_between": 1.0,
            "k_groups": len(set(groups)),
            "k_effects": len(effects),
        }

    effects_arr = np.array(effects, dtype=float)
    variances_arr = np.where(np.array(variances, dtype=float) < 1e-10, 1e-10, np.array(variances, dtype=float))
    groups_arr = np.array(groups, dtype=str)
    weights = 1.0 / variances_arr

    pooled = float(np.sum(weights * effects_arr) / np.sum(weights))
    q_total = float(np.sum(weights * (effects_arr - pooled) ** 2))

    q_within = 0.0
    k_groups = 0
    for group in np.unique(groups_arr):
        mask = groups_arr == group
        if mask.sum() == 0:
            continue
        k_groups += 1
        if mask.sum() == 1:
            continue
        g_weights = weights[mask]
        g_effects = effects_arr[mask]
        g_pooled = float(np.sum(g_weights * g_effects) / np.sum(g_weights))
        q_within += float(np.sum(g_weights * (g_effects - g_pooled) ** 2))

    q_between = float(q_total - q_within)
    df_total = len(effects) - 1
    df_within = len(effects) - k_groups
    df_between = k_groups - 1
    i2_total = max(0.0, (q_total - df_total) / q_total) if q_total > 0 else 0.0
    i2_within = max(0.0, (q_within - df_within) / q_within) if q_within > 0 and df_within > 0 else 0.0
    p_between = float(stats.chi2.sf(q_between, df_between)) if df_between > 0 else 1.0

    return {
        "Q_total": q_total,
        "Q_within": q_within,
        "Q_between": q_between,
        "Q_B_fraction": (q_between / q_total) if q_total > 0 else 0.0,
        "I2_total": i2_total,
        "I2_within": i2_within,
        "p_between": p_between,
        "k_groups": k_groups,
        "k_effects": len(effects),
    }


def permutation_q_decomposition(
    effects: list[float],
    variances: list[float],
    groups: list[str],
    n_permutations: int,
    seed: int,
) -> dict[str, Any]:
    """Empirical test for whether observed Q_between exceeds random labelings."""
    observed = q_decomposition(effects, variances, groups)
    observed_fraction = observed["Q_B_fraction"]
    groups_arr = np.array(groups, dtype=str)
    rng = np.random.default_rng(seed)
    null_fractions = np.empty(n_permutations, dtype=float)

    for i in range(n_permutations):
        permuted = rng.permutation(groups_arr)
        null_fractions[i] = q_decomposition(effects, variances, permuted.tolist())["Q_B_fraction"]

    p_value = float((np.sum(null_fractions >= observed_fraction) + 1) / (n_permutations + 1))
    return {
        "observed_Q_B_fraction": observed_fraction,
        "null_mean": float(np.mean(null_fractions)),
        "null_95th": float(np.percentile(null_fractions, 95)),
        "null_99th": float(np.percentile(null_fractions, 99)),
        "p_value": p_value,
        "n_permutations": n_permutations,
        "seed": seed,
    }
