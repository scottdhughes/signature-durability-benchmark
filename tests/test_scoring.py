import numpy as np
import pandas as pd
import pytest
from signature_durability_benchmark.scoring import score_signature_in_cohort
from signature_durability_benchmark.normalize import normalize_signature, compute_coverage

def test_perfect_separation_gives_large_effect():
    genes = ["GENE_A", "GENE_B", "GENE_C"]
    sig = pd.DataFrame({"gene_symbol": genes, "direction": ["up", "up", "up"], "weight": [1.0, 1.0, 1.0]})
    expr = pd.DataFrame({
        "gene_symbol": genes,
        "S1": [10.0, 10.0, 10.0], "S2": [10.0, 10.0, 10.0],
        "S3": [0.0, 0.0, 0.0], "S4": [0.0, 0.0, 0.0],
    }).set_index("gene_symbol")
    pheno = pd.DataFrame({"sample": ["S1", "S2", "S3", "S4"], "phenotype": ["case", "case", "control", "control"]})
    result = score_signature_in_cohort(sig, expr, pheno, "phenotype", "case", "control")
    assert result["cohens_d"] > 2.0
    assert result["direction_consistent"]
    assert result["coverage_fraction"] == 1.0

def test_no_overlap_gives_zero_coverage():
    sig = pd.DataFrame({"gene_symbol": ["MISSING_A", "MISSING_B"], "direction": ["up", "up"], "weight": [1.0, 1.0]})
    expr = pd.DataFrame({"gene_symbol": ["GENE_X"], "S1": [5.0], "S2": [5.0]}).set_index("gene_symbol")
    pheno = pd.DataFrame({"sample": ["S1", "S2"], "phenotype": ["case", "control"]})
    result = score_signature_in_cohort(sig, expr, pheno, "phenotype", "case", "control")
    assert result["coverage_fraction"] == 0.0

def test_down_direction_inverts_score():
    sig = pd.DataFrame({"gene_symbol": ["G1"], "direction": ["down"], "weight": [1.0]})
    expr = pd.DataFrame({"gene_symbol": ["G1"], "S1": [10.0], "S2": [9.0], "S3": [0.0], "S4": [1.0]}).set_index("gene_symbol")
    pheno = pd.DataFrame({"sample": ["S1", "S2", "S3", "S4"], "phenotype": ["case", "case", "control", "control"]})
    result = score_signature_in_cohort(sig, expr, pheno, "phenotype", "case", "control")
    assert result["cohens_d"] < 0  # inverted because direction=down

def test_normalize_deduplicates():
    sig = pd.DataFrame({"gene_symbol": ["A", "A", "B"], "direction": ["up", "up", "up"], "weight": [1.0, 0.5, 1.0]})
    normalized, audit = normalize_signature(sig)
    assert len(normalized) == 2
    assert audit["duplicate_rows_removed"] == 1

def test_coverage_partial():
    sig = pd.DataFrame({"gene_symbol": ["A", "B", "C", "D"], "direction": ["up"]*4, "weight": [1.0]*4})
    coverage = compute_coverage(sig, {"A", "B"})
    assert coverage == 0.5
