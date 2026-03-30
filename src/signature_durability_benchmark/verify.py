"""Verification checks."""
from __future__ import annotations
from pathlib import Path
from typing import Any
from .config import SkillConfig
from .constants import REQUIRED_OUTPUTS
from .utils import read_json, write_json

def run_verification(config: SkillConfig, run_dir: str | Path) -> dict[str, Any]:
    run_path = Path(run_dir)
    checks = []

    # Check 1: manifest exists
    manifest_path = run_path / "manifest.json"
    checks.append({"name": "manifest_exists", "passed": manifest_path.exists()})

    # Check 2: all required outputs exist
    missing = [f for f in REQUIRED_OUTPUTS if not (run_path / f).exists()]
    checks.append({"name": "artifacts_exist", "passed": len(missing) == 0, "missing": missing})

    # Check 3: all required outputs are nonempty
    empty = [f for f in REQUIRED_OUTPUTS if (run_path / f).exists() and (run_path / f).stat().st_size == 0]
    checks.append({"name": "artifacts_nonempty", "passed": len(empty) == 0, "empty": empty})

    if not manifest_path.exists():
        result = {"status": "failed", "checks": checks}
        write_json(run_path / "verification.json", result)
        return result

    manifest = read_json(manifest_path)

    # Check 4: success rule
    rule = manifest.get("aggregate_metrics", {}).get("success_rule", {})
    success = rule.get("success", False)
    checks.append({"name": "freeze_ready_success_rule", "passed": success, "details": rule})

    # Check 5: AUPRC margin positive
    margin = rule.get("auprc_margin", 0.0)
    checks.append({"name": "auprc_margin_positive", "passed": margin > 0})

    # Check 6: secondary wins >= 2
    sec_wins = rule.get("secondary_wins", 0)
    checks.append({"name": "secondary_wins_sufficient", "passed": sec_wins >= 2})

    # Check 7: signature count matches config
    cfg_panel = config.path("signature_panel")
    if cfg_panel.exists():
        import pandas as pd
        panel = pd.read_csv(cfg_panel, sep="\t")
        expected = len(panel)
        actual = manifest.get("config", {}).get("signature_count", 0)
        checks.append({"name": "signature_count_matches", "passed": expected == actual})

    # Check 8: cohort count matches
    cfg_cohort = config.path("cohort_manifest")
    if cfg_cohort.exists():
        import pandas as pd
        cohorts = pd.read_csv(cfg_cohort, sep="\t")
        expected = len(cohorts)
        actual = manifest.get("config", {}).get("cohort_count", 0)
        checks.append({"name": "cohort_count_matches", "passed": expected == actual})

    # Check 9: certificates exist and have verdicts
    for cert_name in ["durability_certificate.json", "platform_transfer_certificate.json",
                       "confounder_rejection_certificate.json", "coverage_certificate.json"]:
        cert_path = run_path / cert_name
        if cert_path.exists():
            cert = read_json(cert_path)
            checks.append({"name": f"{cert_name}_valid", "passed": "verdict" in cert})
        else:
            checks.append({"name": f"{cert_name}_valid", "passed": False, "reason": "file missing"})

    # Check 10: runtime environment
    rt = manifest.get("runtime", {})
    checks.append({"name": "runtime_environment", "passed": "python" in str(rt)})

    all_passed = all(c["passed"] for c in checks)
    status = "passed" if all_passed else "failed"

    result = {"status": status, "checks": checks}
    write_json(run_path / "verification.json", result)
    return result
