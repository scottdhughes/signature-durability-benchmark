import json
from pathlib import Path

import yaml

from signature_durability_benchmark.generalization import EXPECTED_GENERALIZATION_WINNER, rank_generalization_candidates
from signature_durability_benchmark.workflow_figure import build_workflow_figure_data, render_workflow_figure


ROOT = Path(__file__).resolve().parents[1]


def test_rank_generalization_candidates_prefers_expected_hypoxia_case():
    ranked = rank_generalization_candidates(
        [
            {
                "signature_id": "hallmark_emt",
                "passes_p_bonf": False,
                "loo_perfect": True,
                "split_perfect": True,
                "mean_abs_foreign_effect": 0.50,
                "within_outside_gap": 2.06,
            },
            {
                "signature_id": "hallmark_hypoxia",
                "passes_p_bonf": True,
                "loo_perfect": True,
                "split_perfect": True,
                "mean_abs_foreign_effect": 1.10,
                "within_outside_gap": 2.73,
            },
        ]
    )
    assert ranked[0]["signature_id"] == EXPECTED_GENERALIZATION_WINNER


def test_workflow_figure_data_contains_required_branches(tmp_path: Path):
    data = build_workflow_figure_data(ROOT)
    assert "input_signature" in data
    assert "durable_branch" in data
    assert "failure_branch" in data
    assert "external_transfer" in data
    assert "prospective_v1" in data
    assert "prospective_v2" in data

    pdf = tmp_path / "workflow.pdf"
    png = tmp_path / "workflow.png"
    render_workflow_figure(data, pdf, png)
    assert pdf.exists() and pdf.stat().st_size > 0
    assert png.exists() and png.stat().st_size > 0


def test_archive_metadata_files_parse():
    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    assert zenodo["license"] == "MIT"
    assert citation["license"] == "MIT"
    assert citation["title"] == "signature-durability-benchmark"


def test_declaration_bundle_contains_required_files():
    bundle = ROOT / "submission" / "archive_bundles" / "prospective_holdout_v2_declaration"
    assert (bundle / "prediction_registry_v2.tsv").exists()
    assert (bundle / "PREDICTION_PROTOCOL_v2.md").exists()
    assert (bundle / "prospective_holdout_v2" / "declaration_receipt.json").exists()
    assert (bundle / "CHECKSUMS.sha256").exists()
    assert (bundle / "RELEASE_NOTES.md").exists()


def test_public_json_artifacts_do_not_embed_local_absolute_paths():
    targets = [
        ROOT / "outputs" / "canonical_v8" / "prospective_holdout_validation.json",
        ROOT / "outputs" / "canonical_v8" / "external_hypoxia_validation.json",
        ROOT / "outputs" / "canonical_v8" / "generalization_case_study.json",
        ROOT / "outputs" / "canonical_v8" / "rescued_signature_case_study.json",
        ROOT / "outputs" / "triage_ifng" / "diagnostic.json",
    ]
    for path in targets:
        text = path.read_text(encoding="utf-8")
        assert "/Users/scott/" not in text, path
