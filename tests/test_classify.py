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


# --- within_program model tests ---

def test_within_program_durable_with_significant_within_data():
    """Signature with significant within-program data should be classified durable."""
    result = classify_signature(
        _base_profile(within_p=0.001, within_d=0.85),
        "within_program",
    )
    assert result["predicted_class"] == "durable"


def test_within_program_confounded_overrides_within_data():
    """Confounder check takes priority even when within-program data is significant."""
    result = classify_signature(
        _base_profile(max_confounder_effect=0.9, within_p=0.001, within_d=0.85),
        "within_program",
    )
    assert result["predicted_class"] == "confounded"


def test_within_program_no_within_data_falls_back():
    """Without within-program data, should fall back to cross-context checks."""
    result = classify_signature(
        _base_profile(),  # no within_p or within_d
        "within_program",
    )
    # With the base profile (good metrics, low confounder), falls through to durable
    assert result["predicted_class"] == "durable"


def test_within_program_nonsignificant_within_falls_back_to_mixed():
    """Non-significant within-program data with high heterogeneity -> mixed."""
    result = classify_signature(
        _base_profile(within_p=0.3, within_d=0.5, i_squared=0.9),
        "within_program",
    )
    assert result["predicted_class"] == "mixed"


def test_within_program_small_effect_falls_back():
    """Within-program d below 0.2 threshold should not trigger durable."""
    result = classify_signature(
        _base_profile(within_p=0.01, within_d=0.1, i_squared=0.9),
        "within_program",
    )
    assert result["predicted_class"] == "mixed"


def test_within_program_insufficient_coverage():
    """Coverage check applies before within-program logic."""
    result = classify_signature(
        _base_profile(mean_coverage=0.2, within_p=0.001, within_d=0.85),
        "within_program",
    )
    assert result["predicted_class"] == "insufficient_coverage"
