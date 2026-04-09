"""Confounder scoring: score signature against nuisance gene sets."""
from __future__ import annotations
from typing import Any
from pathlib import Path
import yaml
import pandas as pd
from .scoring import score_signature_in_cohort

def load_confounder_sets(path: str | Path) -> dict[str, pd.DataFrame]:
    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    tables = {}
    for name, value in payload["sets"].items():
        rows = [{"gene_symbol": str(g).upper(), "direction": "up", "weight": 1.0} for g in value.get("genes_up", [])]
        tables[name] = pd.DataFrame(rows)
    return tables

def score_confounders_in_cohort(
    confounder_sets: dict[str, pd.DataFrame],
    expr: pd.DataFrame,
    pheno: pd.DataFrame,
    phenotype_column: str,
    case_label: str,
    control_label: str,
) -> dict[str, float]:
    scores = {}
    for name, gene_set in confounder_sets.items():
        result = score_signature_in_cohort(gene_set, expr, pheno, phenotype_column, case_label, control_label)
        scores[name] = abs(float(result["cohens_d"]))
    return scores
