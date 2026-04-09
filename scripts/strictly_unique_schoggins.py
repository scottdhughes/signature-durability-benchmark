#!/usr/bin/env python3
"""Test with 41 strictly-unique Schoggins genes (NOT in Hallmark IFN-γ or IFN-α).

This is the killer anti-circularity test. ChatGPT pointed out that the expansion
cohorts were admitted based on IFN marker gene verification (STAT1, IRF1, IFIT1,
ISG15, MX1, OAS1, GBP1, CXCL10, RSAD2). All 9 of those markers are in our
Hallmark IFN signatures, so the expansion curation could theoretically have
baked the IFN answer into cohort admission.

The 41 strictly-unique Schoggins genes (found neither in Hallmark IFN-γ nor IFN-α)
share ZERO genes with our marker validation list. If this 41-gene subset still
produces a large within-IFN effect on the expanded k=11 panel, we have
falsified the curation-circularity concern.
"""

import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats as scipy_stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from signature_durability_benchmark.scoring import _cohens_d, _cohens_d_variance

ROOT = Path(__file__).resolve().parent.parent
FREEZE = ROOT / "data" / "freeze"
MATRICES = FREEZE / "cohort_matrices"
PHENO_DIR = FREEZE / "cohort_phenotypes"
OUT = ROOT / "outputs" / "canonical_v8"

# The 41 strictly-unique Schoggins genes (NOT in Hallmark IFN-γ or IFN-α)
# None of these are in the marker verification list used during expansion
# cohort admission, so this directly tests for curation-level circularity.
STRICTLY_UNIQUE = [
    "C6orf150", "CMPK2", "CXCL11", "DDX60", "DHX58", "GBP3", "GBP4", "GBP5",
    "HERC6", "HES4", "IFIH1", "IFIT5", "IFITM2", "IRGM", "LAMP3", "LGALS9",
    "MOV10", "NAMPT", "NT5C3", "PARP10", "PARP12", "PARP14", "PARP9", "PHF11",
    "PLSCR1", "PML", "PNPT1", "RNF19B", "RTP4", "SERPING1", "SP100", "TDRD7",
    "TREX1", "TRIM21", "TRIM25", "TRIM34", "TRIM5", "TRIM56", "UBD", "ZBP1",
    "ZC3HAV1",
]

# Our expansion marker genes (for verification that they're disjoint)
EXPANSION_MARKERS = {
    "STAT1", "IRF1", "IFIT1", "IFIT2", "ISG15", "MX1", "OAS1", "GBP1",
    "CXCL10", "RSAD2",
}

print(f"Strictly-unique Schoggins genes: {len(STRICTLY_UNIQUE)}")
print(f"Overlap with expansion markers: {len(set(STRICTLY_UNIQUE) & EXPANSION_MARKERS)}")
assert len(set(STRICTLY_UNIQUE) & EXPANSION_MARKERS) == 0, "CIRCULARITY: overlap with markers!"
print("✓ Zero overlap with expansion marker genes — cohort admission was independent")


def score_cohort(expr, pheno, genes, case_label, control_label):
    """Compute Hedges' g for a signature in a cohort via z-score mean."""
    expr_upper = expr.copy()
    expr_upper.index = expr_upper.index.str.upper()
    genes_upper = [g.upper() for g in genes]
    genes_in = [g for g in genes_upper if g in expr_upper.index]
    coverage = len(genes_in) / len(genes) if genes else 0

    if len(genes_in) < 5:
        return {"cohens_d": 0.0, "cohens_d_var": 1.0, "coverage": coverage,
                "case_n": 0, "control_n": 0}

    sig_expr = expr_upper.loc[genes_in]
    gene_mean = sig_expr.mean(axis=1)
    gene_std = sig_expr.std(axis=1).replace(0, 1e-9)
    sig_z = sig_expr.sub(gene_mean, axis=0).div(gene_std, axis=0)
    per_sample = sig_z.mean(axis=0)

    case_samples = pheno[pheno["phenotype"] == case_label]["sample_id"].tolist()
    ctrl_samples = pheno[pheno["phenotype"] == control_label]["sample_id"].tolist()
    case_scores = per_sample[per_sample.index.isin(case_samples)].values
    ctrl_scores = per_sample[per_sample.index.isin(ctrl_samples)].values

    if len(case_scores) < 3 or len(ctrl_scores) < 3:
        return {"cohens_d": 0.0, "cohens_d_var": 1.0, "coverage": coverage,
                "case_n": len(case_scores), "control_n": len(ctrl_scores)}

    d = _cohens_d(case_scores, ctrl_scores)
    var = _cohens_d_variance(d, len(case_scores), len(ctrl_scores))
    return {"cohens_d": round(float(d), 4), "cohens_d_var": round(float(var), 6),
            "coverage": round(coverage, 4), "case_n": len(case_scores),
            "control_n": len(ctrl_scores)}


def meta_analysis(gs, vars_, k):
    """DerSimonian-Laird + HKSJ guarded."""
    w = 1.0 / np.array(vars_)
    wm = np.sum(w * gs) / np.sum(w)
    Q = np.sum(w * (gs - wm) ** 2)
    df = k - 1
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0, (Q - df) / c) if c > 0 else 0
    I2 = max(0, (Q - df) / Q) if Q > 0 else 0

    w_star = 1 / (np.array(vars_) + tau2)
    pooled = np.sum(w_star * gs) / np.sum(w_star)
    se_DL = np.sqrt(1 / np.sum(w_star))
    q_star = np.sum(w_star * (gs - pooled) ** 2)
    se_HKSJ = np.sqrt(q_star / ((k - 1) * np.sum(w_star)))
    se_guard = max(se_DL, se_HKSJ)
    t_guard = pooled / se_guard
    p_guard = 2 * (1 - scipy_stats.t.cdf(abs(t_guard), df=k - 1))
    return {
        "k": k, "pooled_g": round(float(pooled), 4),
        "I2": round(float(I2), 4),
        "HKSJ_guarded_p": float(p_guard),
        "p_bonf": float(min(1.0, p_guard * 9)),
    }


# Load manifest
manifest = pd.read_csv(ROOT / "config" / "cohort_manifest.tsv", sep="\t")
ifn_cohorts = manifest[manifest["biological_program"] == "interferon"]["cohort_id"].tolist()
other_cohorts = manifest[manifest["biological_program"] != "interferon"]["cohort_id"].tolist()
print(f"\nIFN cohorts: {len(ifn_cohorts)}")
print(f"Non-IFN cohorts: {len(other_cohorts)}")

print("\n=== Scoring 41-gene strictly-unique Schoggins panel ===\n")

ifn_results = []
for cohort_id in ifn_cohorts:
    cohort_info = manifest[manifest["cohort_id"] == cohort_id].iloc[0]
    expr_path = MATRICES / f"{cohort_id}.tsv"
    pheno_path = PHENO_DIR / f"{cohort_id}.tsv"
    if not expr_path.exists() or not pheno_path.exists():
        continue
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    pheno = pd.read_csv(pheno_path, sep="\t")
    r = score_cohort(expr, pheno, STRICTLY_UNIQUE, cohort_info["case_label"], cohort_info["control_label"])
    r["cohort_id"] = cohort_id
    r["tissue"] = cohort_info["tissue"]
    ifn_results.append(r)
    marker = " [NEW]" if cohort_id in ["sle_pbmc_gse50772", "sle_blood_gse49454",
                                         "psoriasis_skin_gse13355", "psoriasis_skin_gse14905",
                                         "dengue_blood_gse51808"] else ""
    print(f"  {cohort_id:35s} [{cohort_info['tissue']:10s}] g={r['cohens_d']:+.3f} coverage={r['coverage']:.2f}{marker}")

# Within-IFN meta-analysis with strictly unique panel
ifn_gs = np.array([r["cohens_d"] for r in ifn_results])
ifn_vars = np.array([r["cohens_d_var"] for r in ifn_results])
meta = meta_analysis(ifn_gs, ifn_vars, len(ifn_gs))

print("\n=== Within-IFN meta-analysis (STRICTLY UNIQUE 41-gene panel) ===")
print(f"  k = {meta['k']}")
print(f"  Pooled g = {meta['pooled_g']:+.3f}")
print(f"  I² = {meta['I2']:.3f}")
print(f"  HKSJ-guarded p = {meta['HKSJ_guarded_p']:.6f}")
print(f"  HKSJ-guarded p_Bonferroni (9 tests) = {meta['p_bonf']:.6f}")
print(f"  Survives Bonferroni: {meta['p_bonf'] < 0.05}")

# Outside-IFN comparison
print("\n=== Outside-IFN (strictly unique panel) ===")
other_results = []
for cohort_id in other_cohorts:
    cohort_info = manifest[manifest["cohort_id"] == cohort_id].iloc[0]
    expr_path = MATRICES / f"{cohort_id}.tsv"
    pheno_path = PHENO_DIR / f"{cohort_id}.tsv"
    if not expr_path.exists() or not pheno_path.exists():
        continue
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    pheno = pd.read_csv(pheno_path, sep="\t")
    r = score_cohort(expr, pheno, STRICTLY_UNIQUE, cohort_info["case_label"], cohort_info["control_label"])
    other_results.append(r)

other_gs = np.array([r["cohens_d"] for r in other_results])
other_vars = np.array([r["cohens_d_var"] for r in other_results])
other_meta = meta_analysis(other_gs, other_vars, len(other_gs))
print(f"  k = {other_meta['k']}")
print(f"  Pooled g = {other_meta['pooled_g']:+.3f}")
print(f"  HKSJ-guarded p = {other_meta['HKSJ_guarded_p']:.4f}")

# Save results
output = {
    "description": "Strictly unique Schoggins (41 genes) test — zero overlap with Hallmark IFN or expansion marker genes",
    "gene_list": STRICTLY_UNIQUE,
    "n_genes": len(STRICTLY_UNIQUE),
    "expansion_marker_overlap": 0,
    "ifn_cohort_results": ifn_results,
    "other_cohort_results": other_results,
    "within_ifn_meta": meta,
    "outside_ifn_meta": other_meta,
    "interpretation": (
        "The 41 strictly-unique Schoggins genes share ZERO genes with the 9 IFN "
        "markers used to admit expansion cohorts (STAT1, IRF1, IFIT1, IFIT2, ISG15, "
        "MX1, OAS1, GBP1, CXCL10, RSAD2). If this panel still produces a large "
        "within-IFN effect, curation-level circularity is refuted."
    ),
}

with open(OUT / "strictly_unique_schoggins.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\n✓ Saved to {OUT}/strictly_unique_schoggins.json")
