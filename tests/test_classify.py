import pytest
from signature_durability_benchmark.classify import classify_signature

def _base_profile(**overrides):
    profile = {
        "mean_coverage": 0.85,
        "aggregate_effect": 0.8,
        "aggregate_p": 0.001,
        "direction_consistency": 0.90,
        "i_squared": 0.30,
        "loo_stability": 0.85,
        "max_confounder_effect": 0.2,
        "null_separation_p": 0.01,
    }
    profile.update(overrides)
    return profile

def test_durable():
    assert classify_signature(_base_profile(), "full_model")["predicted_class"] == "durable"

def test_confounded():
    assert classify_signature(_base_profile(max_confounder_effect=0.9), "full_model")["predicted_class"] == "confounded"

def test_brittle_low_direction():
    assert classify_signature(_base_profile(direction_consistency=0.3), "full_model")["predicted_class"] == "brittle"

def test_brittle_high_p():
    assert classify_signature(_base_profile(aggregate_p=0.5), "full_model")["predicted_class"] == "brittle"

def test_mixed_high_heterogeneity():
    assert classify_signature(_base_profile(i_squared=0.9), "full_model")["predicted_class"] == "mixed"

def test_mixed_low_loo():
    assert classify_signature(_base_profile(loo_stability=0.3), "full_model")["predicted_class"] == "mixed"

def test_insufficient_coverage():
    assert classify_signature(_base_profile(mean_coverage=0.2), "full_model")["predicted_class"] == "insufficient_coverage"

def test_overlap_only_never_confounded():
    result = classify_signature(_base_profile(max_confounder_effect=0.9), "overlap_only")
    assert result["predicted_class"] != "confounded"

def test_overlap_only_never_mixed():
    result = classify_signature(_base_profile(i_squared=0.9, loo_stability=0.3), "overlap_only")
    assert result["predicted_class"] in ("durable", "brittle", "insufficient_coverage")

def test_effect_only_never_confounded():
    result = classify_signature(_base_profile(max_confounder_effect=0.9), "effect_only")
    assert result["predicted_class"] != "confounded"

def test_null_aware_never_confounded():
    result = classify_signature(_base_profile(max_confounder_effect=0.9), "null_aware")
    assert result["predicted_class"] != "confounded"

def test_no_confounder_never_confounded():
    result = classify_signature(_base_profile(max_confounder_effect=0.9), "no_confounder")
    assert result["predicted_class"] != "confounded"
