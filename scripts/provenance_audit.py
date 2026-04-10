#!/usr/bin/env python3
"""Run the reproducibility and provenance audit."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from signature_durability_benchmark.audit import write_provenance_audit


def main() -> None:
    report = write_provenance_audit(ROOT)
    summary = {
        "status": report["status"],
        "active_cohorts": report["active_cohort_panel"]["active_cohorts"],
        "active_samples": report["active_cohort_panel"]["active_sample_sum"],
        "paper_target_all_non_synthetic": report["signature_panel"]["paper_target_all_non_synthetic"],
        "unexpected_keyword_hits": len(report["keyword_audit"]["unexpected_hits"]),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
