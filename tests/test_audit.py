from pathlib import Path

from signature_durability_benchmark.audit import build_provenance_audit


ROOT = Path(__file__).resolve().parents[1]


def test_provenance_audit_core_checks_pass() -> None:
    report = build_provenance_audit(ROOT)
    active = report["active_cohort_panel"]
    signatures = report["signature_panel"]
    outputs = report["canonical_outputs"]

    assert active["active_cohorts"] == 35
    assert active["active_sample_sum"] == 5922
    assert active["active_manifest_matches_freeze_manifest"] is True
    assert active["phenotype_counts_match_manifest"] is True
    assert active["matrix_columns_match_manifest"] is True

    assert signatures["freeze_signature_id_count"] == 30
    assert signatures["config_manifest_matches_freeze_manifest"] is True
    assert signatures["missing_from_config_manifest"] == []
    assert signatures["missing_from_freeze_manifest"] == []
    assert signatures["paper_target_all_non_synthetic"] is True

    assert outputs["per_cohort_rows_match_expected"] is True
    assert outputs["per_cohort_effect_rows"] == 1050


def test_provenance_audit_has_no_unexpected_keyword_hits() -> None:
    report = build_provenance_audit(ROOT)
    assert report["keyword_audit"]["unexpected_hits"] == []
