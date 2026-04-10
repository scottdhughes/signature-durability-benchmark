from signature_durability_benchmark.diagnostic import _infer_best_program


def test_infer_best_program_prefers_within_outside_gap():
    diagnostics = {
        "interferon": {
            "within_k": 11,
            "outside_k": 24,
            "within": {"pooled_effect": 1.2, "pooled_p": 0.001},
            "outside": {"pooled_effect": 0.2, "pooled_p": 0.4},
        },
        "hypoxia": {
            "within_k": 6,
            "outside_k": 29,
            "within": {"pooled_effect": 0.7, "pooled_p": 0.02},
            "outside": {"pooled_effect": 0.5, "pooled_p": 0.04},
        },
    }
    result = _infer_best_program(diagnostics)
    assert result is not None
    assert result["best"]["program"] == "interferon"


def test_infer_best_program_ignores_negative_within_signal():
    diagnostics = {
        "interferon": {
            "within_k": 11,
            "outside_k": 24,
            "within": {"pooled_effect": -0.8, "pooled_p": 0.001},
            "outside": {"pooled_effect": 0.1, "pooled_p": 0.7},
        },
        "emt": {
            "within_k": 6,
            "outside_k": 29,
            "within": {"pooled_effect": 0.6, "pooled_p": 0.01},
            "outside": {"pooled_effect": 0.1, "pooled_p": 0.8},
        },
    }
    result = _infer_best_program(diagnostics)
    assert result is not None
    assert result["best"]["program"] == "emt"
