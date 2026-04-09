"""Gene normalization and coverage auditing."""
from __future__ import annotations
from typing import Any
import pandas as pd

def normalize_signature(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    normalized = frame.copy()
    normalized["gene_symbol"] = normalized["gene_symbol"].astype(str).str.upper().str.strip()
    normalized["direction"] = normalized["direction"].astype(str).str.lower().str.strip()
    normalized["weight"] = pd.to_numeric(normalized["weight"], errors="coerce").fillna(1.0).abs()
    normalized = normalized.loc[normalized["gene_symbol"] != ""].copy()
    deduplicated = (
        normalized.sort_values(["weight", "gene_symbol"], ascending=[False, True])
        .drop_duplicates(subset=["gene_symbol"], keep="first")
        .sort_values("gene_symbol").reset_index(drop=True)
    )
    audit = {
        "raw_row_count": int(frame.shape[0]),
        "normalized_gene_count": int(deduplicated.shape[0]),
        "dropped_blank_rows": int(frame.shape[0] - normalized.shape[0]),
        "duplicate_rows_removed": int(normalized.shape[0] - deduplicated.shape[0]),
    }
    return deduplicated, audit

def compute_coverage(signature: pd.DataFrame, cohort_genes: set[str]) -> float:
    if signature.empty:
        return 0.0
    matched = signature["gene_symbol"].isin(cohort_genes).sum()
    return float(matched / len(signature))
