"""Per-cohort signature scoring."""
from __future__ import annotations
from typing import Any
import numpy as np
import pandas as pd
from .normalize import normalize_signature, compute_coverage

def _sample_scores(signature: pd.DataFrame, expr: pd.DataFrame) -> pd.Series:
    """Compute weighted signed mean of signature genes per sample."""
    shared = signature.loc[signature["gene_symbol"].isin(expr.index)].copy()
    if shared.empty:
        return pd.Series(dtype=float)
    values = expr.loc[shared["gene_symbol"].values].copy()
    weights = shared.set_index("gene_symbol")["weight"]
    signs = shared.set_index("gene_symbol")["direction"].map({"up": 1.0, "down": -1.0}).fillna(0.0)
    weighted = values.multiply(weights * signs, axis=0)
    return weighted.sum(axis=0) / weights.sum()

def _cohens_d(case_scores: np.ndarray, control_scores: np.ndarray) -> float:
    n1, n2 = len(case_scores), len(control_scores)
    if n1 < 2 or n2 < 2:
        return 0.0
    mean_diff = float(np.mean(case_scores) - np.mean(control_scores))
    pooled_var = ((n1 - 1) * np.var(case_scores, ddof=1) + (n2 - 1) * np.var(control_scores, ddof=1)) / (n1 + n2 - 2)
    pooled_sd = float(np.sqrt(pooled_var)) if pooled_var > 0 else 1e-10
    return mean_diff / pooled_sd

def _cohens_d_variance(d: float, n1: int, n2: int) -> float:
    """Approximate variance of Cohen's d for meta-analysis weighting."""
    if n1 < 2 or n2 < 2:
        return 1.0
    return (n1 + n2) / (n1 * n2) + (d * d) / (2 * (n1 + n2))

def score_signature_in_cohort(
    signature: pd.DataFrame,
    expr: pd.DataFrame,
    pheno: pd.DataFrame,
    phenotype_column: str,
    case_label: str,
    control_label: str,
) -> dict[str, Any]:
    normalized, audit = normalize_signature(signature)
    cohort_genes = set(expr.index.astype(str))
    coverage = compute_coverage(normalized, cohort_genes)
    if coverage < 0.01:
        return {"cohens_d": 0.0, "cohens_d_var": 1.0, "direction_consistent": False,
                "coverage_fraction": coverage, "case_n": 0, "control_n": 0,
                "mean_case": 0.0, "mean_control": 0.0, **audit}
    scores = _sample_scores(normalized, expr)
    sample_col = "sample" if "sample" in pheno.columns else pheno.columns[0]
    case_samples = pheno.loc[pheno[phenotype_column] == case_label, sample_col].values
    control_samples = pheno.loc[pheno[phenotype_column] == control_label, sample_col].values
    case_scores = scores.reindex(case_samples).dropna().values
    control_scores = scores.reindex(control_samples).dropna().values
    d = _cohens_d(case_scores, control_scores)
    d_var = _cohens_d_variance(d, len(case_scores), len(control_scores))
    return {
        "cohens_d": float(d),
        "cohens_d_var": float(d_var),
        "direction_consistent": d > 0,
        "coverage_fraction": float(coverage),
        "case_n": int(len(case_scores)),
        "control_n": int(len(control_scores)),
        "mean_case": float(np.mean(case_scores)) if len(case_scores) else 0.0,
        "mean_control": float(np.mean(control_scores)) if len(control_scores) else 0.0,
        **audit,
    }
