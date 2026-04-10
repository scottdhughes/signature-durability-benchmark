#!/usr/bin/env python3
"""Backward-compatible v1 wrapper around the generic prospective evaluator."""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from signature_durability_benchmark.prospective import evaluate_round, write_legacy_v1_outputs


ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "prospective_holdout"
REGISTRY = DATA_DIR / "prediction_registry_v1.tsv"
PROTOCOL = DATA_DIR / "PREDICTION_PROTOCOL.md"
OUT_DIR = ROOT / "outputs" / "canonical_v8"
ROUND_DIR = OUT_DIR / "prospective_rounds" / "prospective_holdout_v1"
LEGACY_JSON = OUT_DIR / "prospective_holdout_validation.json"
LEGACY_MD = OUT_DIR / "prospective_holdout_validation.md"
LEGACY_CSV = OUT_DIR / "prospective_holdout_per_cohort_effects.csv"


def main() -> None:
    evaluation = evaluate_round(
        REGISTRY,
        PROTOCOL,
        out_dir=ROUND_DIR,
    )
    shutil.copy2(ROUND_DIR / "per_cohort_effects.csv", LEGACY_CSV)
    write_legacy_v1_outputs(
        evaluation,
        registry_path=REGISTRY,
        protocol_path=PROTOCOL,
        out_json=LEGACY_JSON,
        out_md=LEGACY_MD,
    )
    print(json.dumps({"overall_success": evaluation["prediction_summary"]["overall_success"], "output_json": str(LEGACY_JSON)}, indent=2))


if __name__ == "__main__":
    main()
