import numpy as np
import pandas as pd
import pytest
from signature_durability_benchmark.meta_analysis import (
    fixed_effect_meta, i_squared, leave_one_out, loo_stability,
    direction_consistency, platform_holdout_consistency,
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
