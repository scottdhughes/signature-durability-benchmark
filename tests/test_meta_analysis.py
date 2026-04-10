import numpy as np
import pandas as pd
import pytest
from signature_durability_benchmark.meta_analysis import (
    dl_random_effects_meta, fixed_effect_meta, i_squared, leave_one_out,
    loo_stability, direction_consistency, platform_holdout_consistency,
    permutation_q_decomposition, q_decomposition, within_program_meta,
)
from signature_durability_benchmark.null_model import generate_null_signatures, null_separation_p

def test_fixed_effect_consistent_studies():
    effects = [0.5, 0.6, 0.4, 0.55]
    variances = [0.1, 0.1, 0.1, 0.1]
    pooled = fixed_effect_meta(effects, variances)
    assert 0.4 < pooled["pooled_effect"] < 0.7
    assert pooled["pooled_p"] < 0.05

def test_i_squared_homogeneous():
    effects = [0.5, 0.5, 0.5]
    variances = [0.1, 0.1, 0.1]
    assert i_squared(effects, variances) < 0.1

def test_i_squared_heterogeneous():
    effects = [2.0, -1.0, 3.0, -2.0]
    variances = [0.01, 0.01, 0.01, 0.01]
    assert i_squared(effects, variances) > 0.5

def test_loo_stability_all_positive():
    effects = [0.5, 0.6, 0.4, 0.55]
    variances = [0.1, 0.1, 0.1, 0.1]
    assert loo_stability(effects, variances) == 1.0

def test_direction_consistency():
    assert direction_consistency([0.5, 0.3, -0.1, 0.2]) == 0.75

def test_null_signatures_correct_count():
    universe = [f"GENE_{i}" for i in range(1000)]
    sigs = generate_null_signatures(30, universe, 10, seed=42)
    assert len(sigs) == 10
    assert all(len(s) == 30 for s in sigs)

def test_null_separation_extreme():
    assert null_separation_p(10.0, [0.1, 0.2, 0.3, 0.4]) == 0.0

def test_null_separation_within():
    assert null_separation_p(0.2, [0.1, 0.2, 0.3, 0.4]) > 0.5

def test_platform_holdout_consistency_uniform():
    effects = [0.5, 0.6, 0.4, 0.55]
    variances = [0.1, 0.1, 0.1, 0.1]
    platforms = ["affy", "affy", "illumina", "illumina"]
    assert platform_holdout_consistency(effects, variances, platforms) == 1.0


# --- DerSimonian-Laird random-effects ---

def test_dl_random_effects_homogeneous():
    """With homogeneous effects, DL should give tau2 ~ 0 and match fixed-effect."""
    effects = [0.5, 0.5, 0.5]
    variances = [0.1, 0.1, 0.1]
    result = dl_random_effects_meta(effects, variances)
    assert abs(result["pooled_effect"] - 0.5) < 0.01
    assert result["tau2"] < 0.01
    assert result["i_squared"] < 0.1
    assert result["k"] == 3
    assert result["method"] == "DL_random_effects"

def test_dl_random_effects_heterogeneous():
    """With heterogeneous effects, DL should estimate positive tau2."""
    effects = [2.0, -1.0, 3.0, -2.0]
    variances = [0.01, 0.01, 0.01, 0.01]
    result = dl_random_effects_meta(effects, variances)
    assert result["tau2"] > 0
    assert result["i_squared"] > 0.5
    assert result["k"] == 4

def test_dl_random_effects_two_studies():
    """DL should work with just 2 studies."""
    effects = [0.5, 0.8]
    variances = [0.05, 0.05]
    result = dl_random_effects_meta(effects, variances)
    assert 0.4 < result["pooled_effect"] < 0.9
    assert result["k"] == 2


# --- Within-program meta-analysis ---

def test_within_program_meta_basic():
    """Within-program should partition effects correctly."""
    effects = [0.5, 0.6, 0.1, 0.05, -0.1, 0.02]
    variances = [0.05, 0.05, 0.1, 0.1, 0.1, 0.1]
    programs = ["ifn", "ifn", "prolif", "prolif", "emt", "emt"]
    result = within_program_meta(effects, variances, programs, "ifn")
    assert result["within_k"] == 2
    assert result["outside_k"] == 4
    assert result["within"] is not None
    assert result["outside"] is not None
    # Within-program effect should be higher
    assert result["within"]["pooled_effect"] > result["outside"]["pooled_effect"]

def test_within_program_meta_no_match():
    """If signature program matches no cohort, within should be empty."""
    effects = [0.5, 0.6]
    variances = [0.05, 0.05]
    programs = ["ifn", "ifn"]
    result = within_program_meta(effects, variances, programs, "emt")
    assert result["within_k"] == 0
    assert result["within"] is None
    assert result["outside_k"] == 2
    assert result["outside"] is not None

def test_within_program_meta_single_cohort():
    """If only 1 within-program cohort, within should be None (need >=2)."""
    effects = [0.5, 0.6, 0.1]
    variances = [0.05, 0.05, 0.1]
    programs = ["ifn", "prolif", "prolif"]
    result = within_program_meta(effects, variances, programs, "ifn")
    assert result["within_k"] == 1
    assert result["within"] is None
    assert result["outside_k"] == 2
    assert result["outside"] is not None


def test_q_decomposition_all_between():
    effects = [1.0, 1.1, -0.9, -1.0]
    variances = [0.1, 0.1, 0.1, 0.1]
    groups = ["a", "a", "b", "b"]
    result = q_decomposition(effects, variances, groups)
    assert result["Q_between"] > result["Q_within"]
    assert result["Q_B_fraction"] > 0.5


def test_permutation_q_decomposition_detects_structure():
    effects = [1.1, 1.0, 1.2, -1.0, -1.1, -1.2]
    variances = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    groups = ["a", "a", "a", "b", "b", "b"]
    result = permutation_q_decomposition(effects, variances, groups, n_permutations=200, seed=42)
    assert result["observed_Q_B_fraction"] > result["null_mean"]
    assert result["p_value"] < 0.1
