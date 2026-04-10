#!/usr/bin/env python3
"""Compute per-cohort effect sizes for:
1. All 35 cohorts (30 original + 5 new IFN expansion)
2. All 30 frozen benchmark signatures, including the orthogonal Schoggins IRG panel

Outputs: outputs/canonical_v8/per_cohort_effects.csv

The Schoggins IRG signature is built from Schoggins et al. 2011 (Nature
472:481-485, doi:10.1038/nature09907) — genes validated as antiviral via
overexpression screens. Uses a 76-gene core with 25-36% overlap with the
frozen MSigDB-anchored Hallmark IFN cores, making it a genuinely orthogonal
validation panel.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from signature_durability_benchmark.scoring import _cohens_d, _cohens_d_variance

ROOT = Path(__file__).resolve().parent.parent
FREEZE = ROOT / "data" / "freeze"
MATRICES = FREEZE / "cohort_matrices"
PHENO_DIR = FREEZE / "cohort_phenotypes"
OUT_DIR = ROOT / "outputs" / "canonical_v8"
OUT_DIR.mkdir(exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════════
# Schoggins 2011 IRG Core Panel (orthogonal to MSigDB)
# ═══════════════════════════════════════════════════════════════════════════════

SCHOGGINS_2011_CORE = [
    "IRF1", "IRF7", "C6orf150", "CMPK2", "DDX60", "HERC6", "HERC5",
    "HES4", "IFI6", "IFIH1", "IFIT1", "IFIT3", "IFITM1", "IFITM2", "IFITM3",
    "IRGM", "ISG15", "ISG20", "LAMP3", "LGALS9", "LY6E", "MOV10",
    "MX1", "MX2", "NAMPT", "NT5C3", "OAS1", "OAS2", "OAS3", "OASL",
    "PHF11", "PML", "PNPT1", "RNF19B", "RSAD2", "RTP4", "SAMD9L",
    "SERPING1", "TDRD7", "TREX1", "TRIM5", "TRIM21", "TRIM22",
    "TRIM25", "TRIM34", "TRIM56", "UBD", "USP18", "ZBP1",
    "CXCL10", "CXCL11", "DHX58", "EIF2AK2", "GBP1", "GBP2", "GBP3", "GBP4", "GBP5",
    "IFI27", "IFI35", "IFI44", "IFI44L", "IFIT2", "IFIT5", "PARP9", "PARP10",
    "PARP12", "PARP14", "PLSCR1", "SAMD9", "SP100", "STAT1", "STAT2", "TAP1",
    "XAF1", "ZC3HAV1",
]

# ═══════════════════════════════════════════════════════════════════════════════
# Load signatures without mutating the frozen release asset
# ═══════════════════════════════════════════════════════════════════════════════

print("Loading signatures...")
sigs_path = FREEZE / "signatures.tsv"
sigs = pd.read_csv(sigs_path, sep="\t")
print(f"  Existing signatures: {sigs['signature_id'].nunique()}")

# Legacy fallback: keep the scored path runnable even if an older freeze is
# missing Schoggins, but never rewrite the frozen asset on disk.
if "schoggins_2011_irg" not in sigs["signature_id"].unique():
    print("  WARNING: schoggins_2011_irg missing from frozen signatures.tsv; adding it in memory only")
    new_rows = [
        {"signature_id": "schoggins_2011_irg", "gene_symbol": g, "direction": "up", "weight": 1.0}
        for g in SCHOGGINS_2011_CORE
    ]
    sigs = pd.concat([sigs, pd.DataFrame(new_rows)], ignore_index=True)
    print(f"  Added schoggins_2011_irg signature with {len(new_rows)} genes in memory")

# ═══════════════════════════════════════════════════════════════════════════════
# Score a cohort for a signature (z-score → Cohen's d)
# ═══════════════════════════════════════════════════════════════════════════════

def score_cohort_simple(expr: pd.DataFrame, pheno: pd.DataFrame, genes: list, case_label: str, control_label: str) -> dict:
    """Compute per-sample signature score = mean z-score across signature genes,
    then Cohen's d (Hedges' g) between case and control."""

    # Map gene symbols (uppercase)
    expr_upper = expr.copy()
    expr_upper.index = expr_upper.index.str.upper()
    genes_upper = [g.upper() for g in genes]
    genes_in_matrix = [g for g in genes_upper if g in expr_upper.index]
    coverage = len(genes_in_matrix) / len(genes) if genes else 0

    if len(genes_in_matrix) < 5:
        return {"cohens_d": 0.0, "cohens_d_var": 1.0, "direction_consistent": False,
                "coverage_fraction": coverage, "case_n": 0, "control_n": 0}

    # Subset expression matrix to signature genes
    sig_expr = expr_upper.loc[genes_in_matrix]

    # Z-score per gene across samples
    gene_mean = sig_expr.mean(axis=1)
    gene_std = sig_expr.std(axis=1).replace(0, 1e-9)
    sig_expr_z = sig_expr.sub(gene_mean, axis=0).div(gene_std, axis=0)

    # Per-sample score = mean z across signature genes
    per_sample = sig_expr_z.mean(axis=0)

    # Get case and control sample IDs
    case_samples = pheno[pheno["phenotype"] == case_label]["sample_id"].tolist()
    ctrl_samples = pheno[pheno["phenotype"] == control_label]["sample_id"].tolist()

    case_scores = per_sample[per_sample.index.isin(case_samples)].values
    ctrl_scores = per_sample[per_sample.index.isin(ctrl_samples)].values

    if len(case_scores) < 3 or len(ctrl_scores) < 3:
        return {"cohens_d": 0.0, "cohens_d_var": 1.0, "direction_consistent": False,
                "coverage_fraction": coverage, "case_n": len(case_scores), "control_n": len(ctrl_scores)}

    d = _cohens_d(case_scores, ctrl_scores)
    var = _cohens_d_variance(d, len(case_scores), len(ctrl_scores))

    return {
        "cohens_d": round(float(d), 4),
        "cohens_d_var": round(float(var), 6),
        "direction_consistent": bool(d > 0),
        "coverage_fraction": round(coverage, 4),
        "case_n": len(case_scores),
        "control_n": len(ctrl_scores),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Load cohort manifest
# ═══════════════════════════════════════════════════════════════════════════════

cohort_manifest = pd.read_csv(ROOT / "config" / "cohort_manifest.tsv", sep="\t")
print(f"\nCohorts: {len(cohort_manifest)}")
print(f"Programs: {dict(cohort_manifest['biological_program'].value_counts())}")

# ═══════════════════════════════════════════════════════════════════════════════
# Score every (signature, cohort) pair
# ═══════════════════════════════════════════════════════════════════════════════

all_results = []
sig_ids = sigs["signature_id"].unique()
print(f"\nSignatures: {len(sig_ids)}")
print(f"  Will score {len(sig_ids)} × {len(cohort_manifest)} = {len(sig_ids) * len(cohort_manifest)} (sig, cohort) pairs")

for _, cohort in cohort_manifest.iterrows():
    cohort_id = cohort["cohort_id"]
    expr_path = MATRICES / f"{cohort_id}.tsv"
    pheno_path = PHENO_DIR / f"{cohort_id}.tsv"

    if not expr_path.exists() or not pheno_path.exists():
        print(f"  SKIP {cohort_id}: files not found")
        continue

    try:
        expr = pd.read_csv(expr_path, sep="\t", index_col=0)
        pheno = pd.read_csv(pheno_path, sep="\t")
    except Exception as e:
        print(f"  SKIP {cohort_id}: {e}")
        continue

    for sig_id in sig_ids:
        sig_genes = sigs[sigs["signature_id"] == sig_id]["gene_symbol"].tolist()
        result = score_cohort_simple(
            expr, pheno, sig_genes,
            case_label=cohort["case_label"],
            control_label=cohort["control_label"],
        )
        result["signature_id"] = sig_id
        result["cohort_id"] = cohort_id
        all_results.append(result)

    print(f"  Scored {cohort_id}: {len(sig_ids)} signatures")

# ═══════════════════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════════════════

df = pd.DataFrame(all_results)
# Reorder columns to match original schema
df = df[["signature_id", "cohort_id", "cohens_d", "cohens_d_var", "direction_consistent", "coverage_fraction", "case_n", "control_n"]]
df.to_csv(OUT_DIR / "per_cohort_effects.csv", index=False)
print(f"\n✓ Saved {len(df)} rows to {OUT_DIR}/per_cohort_effects.csv")

# Quick summary: Schoggins in IFN cohorts
print("\n=== Schoggins IRG in IFN cohorts ===")
ifn_cohorts = cohort_manifest[cohort_manifest["biological_program"] == "interferon"]["cohort_id"].tolist()
schoggins_df = df[(df["signature_id"] == "schoggins_2011_irg") & (df["cohort_id"].isin(ifn_cohorts))]
print(schoggins_df[["cohort_id", "cohens_d", "coverage_fraction"]].to_string(index=False))

print("\n=== IFN-γ in IFN cohorts (comparison) ===")
ifng_df = df[(df["signature_id"] == "hallmark_ifng_response") & (df["cohort_id"].isin(ifn_cohorts))]
print(ifng_df[["cohort_id", "cohens_d"]].to_string(index=False))
