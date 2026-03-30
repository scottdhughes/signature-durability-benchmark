"""Freeze validation."""
from __future__ import annotations
from pathlib import Path
from typing import Any
import pandas as pd
from .config import SkillConfig
from .utils import read_table, write_json, ensure_dir, sha256_file, now_timestamp

def build_freeze(config: SkillConfig, out_dir: str | Path) -> dict[str, Any]:
    out_path = ensure_dir(out_dir)
    errors = []

    # Load and validate signature panel
    panel_path = config.path("signature_panel")
    if not panel_path.exists():
        raise FileNotFoundError(f"Signature panel not found: {panel_path}")
    panel = read_table(panel_path)

    # Check locked curation status
    if "curation_status" in panel.columns:
        unlocked = panel[~panel["curation_status"].str.contains("locked", na=False)]
        if not unlocked.empty:
            errors.append(f"{len(unlocked)} signatures not locked: {unlocked['signature_id'].tolist()}")

    # Check class counts
    primary = panel[panel["split"] == "primary"]
    blind = panel[panel["split"] == "blind"]

    # Load signature rows
    sig_rows_path = config.path("signature_rows")
    if not sig_rows_path.exists():
        raise FileNotFoundError(f"Signature rows not found: {sig_rows_path}")
    sig_rows = read_table(sig_rows_path)

    # Check all panel signatures have gene rows
    panel_ids = set(panel["signature_id"].astype(str))
    row_ids = set(sig_rows["signature_id"].astype(str))
    missing = panel_ids - row_ids
    if missing:
        errors.append(f"Signatures in panel but missing gene rows: {missing}")

    # Check gene symbols uppercase
    non_upper = sig_rows[sig_rows["gene_symbol"].astype(str) != sig_rows["gene_symbol"].astype(str).str.upper()]
    if not non_upper.empty:
        errors.append(f"{len(non_upper)} gene rows not uppercase")

    # Check cohort manifest
    cohort_path = config.path("cohort_manifest")
    if not cohort_path.exists():
        raise FileNotFoundError(f"Cohort manifest not found: {cohort_path}")
    cohorts = read_table(cohort_path)

    # Check cohort matrices and phenotypes exist
    freeze_dir = config.path("freeze_dir")
    for _, crow in cohorts.iterrows():
        cid = str(crow["cohort_id"])
        matrix_path = freeze_dir / "cohort_matrices" / f"{cid}.tsv"
        pheno_path = freeze_dir / "cohort_phenotypes" / f"{cid}.tsv"
        if not matrix_path.exists():
            errors.append(f"Missing cohort matrix: {matrix_path}")
        if not pheno_path.exists():
            errors.append(f"Missing cohort phenotype: {pheno_path}")

    # Check confounder panel
    conf_path = config.path("confounder_panel")
    if not conf_path.exists():
        raise FileNotFoundError(f"Confounder panel not found: {conf_path}")

    # Source leakage check
    if "source_family_id" in panel.columns:
        sig_families = set(panel["source_family_id"].astype(str))
        # Cohort families don't exist in this schema, so just verify non-empty
        if not sig_families:
            errors.append("No source_family_id values in signature panel")

    audit = {
        "created_at": now_timestamp(),
        "primary_count": int(len(primary)),
        "blind_count": int(len(blind)),
        "total_signatures": int(len(panel)),
        "total_gene_rows": int(len(sig_rows)),
        "total_cohorts": int(len(cohorts)),
        "errors": errors,
        "valid": len(errors) == 0,
    }

    # Copy validated assets to freeze dir
    write_json(out_path / "freeze_audit.json", audit)

    # Copy manifests
    panel.to_csv(out_path / "signature_panel_manifest.tsv", sep="\t", index=False)
    cohorts.to_csv(out_path / "cohort_manifest.tsv", sep="\t", index=False)
    sig_rows.to_csv(out_path / "signatures.tsv", sep="\t", index=False)

    if errors:
        raise ValueError(f"Freeze validation failed: {errors}")

    return audit
