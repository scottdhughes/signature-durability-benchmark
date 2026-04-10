import json
from pathlib import Path

import pandas as pd
import pytest

from signature_durability_benchmark.prospective import _build_pheno, verify_round_receipt
from signature_durability_benchmark.utils import sha256_file


def _write_registry(path: Path, round_id: str = "prospective_holdout_v2") -> None:
    path.write_text(
        "\n".join(
            [
                "round_id\tdeclared_at_utc\texternal_cohort_id\tgeo_accession\tprediction_group\ttissue\tcontrast_description\tphenotype_field\tcase_label\tcontrol_label\texpected_case_n\texpected_control_n\tacquisition_mode\tseries_matrix_url\tmatrix_url",
                f"{round_id}\t2026-04-10T05:32:41Z\tcohort_a\tGSE000001\tprimary\tblood\tdesc\ttreatment\tcase\tcontrol\t1\t1\tdirect_normalized_matrix\thttps://example.org/series.txt.gz\thttps://example.org/matrix.txt.gz",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_protocol(path: Path) -> None:
    path.write_text("# Protocol\n", encoding="utf-8")


def _write_receipt(path: Path, registry_path: Path, protocol_path: Path, round_id: str = "prospective_holdout_v2") -> None:
    payload = {
        "round_id": round_id,
        "verification_status": "OK",
        "registry_sha256": sha256_file(registry_path),
        "protocol_sha256": sha256_file(protocol_path),
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def test_verify_round_receipt_passes(tmp_path: Path):
    registry = tmp_path / "registry.tsv"
    protocol = tmp_path / "protocol.md"
    receipt = tmp_path / "receipt.json"
    _write_registry(registry)
    _write_protocol(protocol)
    _write_receipt(receipt, registry, protocol)

    verified = verify_round_receipt(registry, protocol, receipt)
    assert verified["round_id"] == "prospective_holdout_v2"


def test_verify_round_receipt_fails_for_tampered_registry(tmp_path: Path):
    registry = tmp_path / "registry.tsv"
    protocol = tmp_path / "protocol.md"
    receipt = tmp_path / "receipt.json"
    _write_registry(registry)
    _write_protocol(protocol)
    _write_receipt(receipt, registry, protocol)
    registry.write_text(registry.read_text(encoding="utf-8").replace("cohort_a", "cohort_b"), encoding="utf-8")

    with pytest.raises(ValueError, match="Registry SHA256"):
        verify_round_receipt(registry, protocol, receipt)


def test_verify_round_receipt_fails_for_tampered_protocol(tmp_path: Path):
    registry = tmp_path / "registry.tsv"
    protocol = tmp_path / "protocol.md"
    receipt = tmp_path / "receipt.json"
    _write_registry(registry)
    _write_protocol(protocol)
    _write_receipt(receipt, registry, protocol)
    protocol.write_text("# Changed protocol\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Protocol SHA256"):
        verify_round_receipt(registry, protocol, receipt)


def test_verify_round_receipt_fails_for_wrong_round_id(tmp_path: Path):
    registry = tmp_path / "registry.tsv"
    protocol = tmp_path / "protocol.md"
    receipt = tmp_path / "receipt.json"
    _write_registry(registry, round_id="prospective_holdout_v2")
    _write_protocol(protocol)
    _write_receipt(receipt, registry, protocol, round_id="prospective_holdout_v1")

    with pytest.raises(ValueError, match="Receipt round_id"):
        verify_round_receipt(registry, protocol, receipt)


def test_build_pheno_supports_subset_and_multilabel_case_values():
    meta = pd.DataFrame(
        {
            "title": ["S1", "S2", "S3", "S4"],
            "treatment": ["Acute_Survivor", "Acute_Non-survivor", "Control", "Acute_Survivor"],
            "timepoint": ["D3", "D3", "D3", "D0"],
        }
    )
    spec = pd.Series(
        {
            "external_cohort_id": "cohort_a",
            "phenotype_field": "treatment",
            "case_label": "Acute_Non-survivor|Acute_Survivor",
            "control_label": "Control",
            "subset_field": "timepoint",
            "subset_value": "D3",
            "expected_case_n": 2,
            "expected_control_n": 1,
        }
    )

    pheno = _build_pheno(spec, meta, ["S1", "S2", "S3", "S4"])
    assert pheno["sample"].tolist() == ["S1", "S2", "S3"]
    assert pheno["phenotype"].tolist() == ["case", "case", "control"]
