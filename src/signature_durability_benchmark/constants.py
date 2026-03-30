"""Pipeline constants."""

DURABLE_CLASSES = frozenset({"durable"})
NONDURABLE_CLASSES = frozenset({"brittle", "mixed", "confounded"})
ALL_CLASSES = frozenset({"durable", "brittle", "mixed", "confounded", "insufficient_coverage"})

MODEL_NAMES = [
    "full_model",
    "overlap_only",
    "effect_only",
    "null_aware",
    "no_confounder",
]

FREEZE_AUDIT_NAME = "freeze_audit.json"

REQUIRED_OUTPUTS = [
    "manifest.json",
    "normalization_audit.json",
    "cohort_overlap_summary.csv",
    "per_cohort_effects.csv",
    "aggregate_durability_scores.csv",
    "matched_null_summary.csv",
    "leave_one_cohort_out.csv",
    "platform_holdout_summary.csv",
    "durability_certificate.json",
    "platform_transfer_certificate.json",
    "confounder_rejection_certificate.json",
    "coverage_certificate.json",
    "benchmark_protocol.json",
    "verification.json",
    "public_summary.md",
    "forest_plot.png",
    "null_separation_plot.png",
    "stability_heatmap.png",
    "platform_transfer_panel.png",
]
