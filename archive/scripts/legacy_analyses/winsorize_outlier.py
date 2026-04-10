"""Winsorize the GSE47533 outlier in hypoxia within-program analysis.

Computes three versions of the within-program DerSimonian-Laird meta-analysis
for hallmark_hypoxia within hypoxia cohorts (k=6):
  1. Original: all 6 cohorts as-is
  2. Leave-one-out: drop GSE47533 (k=5)
  3. Winsorized: replace GSE47533's d with the next-largest d in the program

This directly addresses the reviewer concern that the large d=16.3 (now
Hedges' g corrected) from GSE47533 may be driving the hypoxia result.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json
import sys

sys.path.insert(0, "src")
from signature_durability_benchmark.meta_analysis import dl_random_effects_meta

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
effects = pd.read_csv("outputs/canonical_v7/per_cohort_effects.csv")
manifest = pd.read_csv("config/cohort_manifest.tsv", sep="\t")
cohort_program = dict(zip(manifest.cohort_id, manifest.biological_program))

# Filter: hallmark_hypoxia within hypoxia cohorts
mask = (effects.signature_id == "hallmark_hypoxia") & (
    effects.cohort_id.map(cohort_program) == "hypoxia"
)
hyp = effects[mask].copy().reset_index(drop=True)

print(f"Hypoxia within hypoxia: k={len(hyp)} cohorts")
print(hyp[["cohort_id", "cohens_d", "cohens_d_var", "case_n", "control_n"]].to_string(index=False))

# Bonferroni correction: 9 within-program tests (from within_program_durability.csv)
BONF_N = 9
BONF_THRESHOLD = 0.05 / BONF_N

# ---------------------------------------------------------------------------
# Helper: DL meta with Bonferroni
# ---------------------------------------------------------------------------
def dl_with_bonf(ds, vs, label):
    result = dl_random_effects_meta(list(ds), list(vs))
    p_bonf = min(1.0, result["pooled_p"] * BONF_N)
    print(f"\n{label} (k={len(ds)}):")
    print(f"  Pooled g = {result['pooled_effect']:.4f}")
    print(f"  SE = {result['pooled_se']:.4f}")
    print(f"  tau2 = {result['tau2']:.4f}")
    print(f"  I2 = {result['i_squared']:.4f}")
    print(f"  p (DL) = {result['pooled_p']:.2e}")
    print(f"  p (Bonf) = {p_bonf:.2e}")
    print(f"  Survives Bonferroni (p < {BONF_THRESHOLD:.4f})? {'YES' if p_bonf < 0.05 else 'NO'}")
    return {
        "k": len(ds),
        "pooled_g": round(result["pooled_effect"], 4),
        "pooled_se": round(result["pooled_se"], 4),
        "tau2": round(result["tau2"], 4),
        "i_squared": round(result["i_squared"], 4),
        "p_DL": result["pooled_p"],
        "p_Bonferroni": p_bonf,
        "survives_bonferroni": p_bonf < 0.05,
        "cohort_gs": [round(float(d), 4) for d in ds],
    }


# ---------------------------------------------------------------------------
# Version 1: Original (all 6)
# ---------------------------------------------------------------------------
ds_all = hyp.cohens_d.values
vs_all = hyp.cohens_d_var.values
original = dl_with_bonf(ds_all, vs_all, "ORIGINAL (all 6 cohorts)")

# ---------------------------------------------------------------------------
# Version 2: Leave-one-out (drop GSE47533)
# ---------------------------------------------------------------------------
outlier_mask = hyp.cohort_id == "hypoxia_timecourse_gse47533"
ds_loo = hyp.loc[~outlier_mask, "cohens_d"].values
vs_loo = hyp.loc[~outlier_mask, "cohens_d_var"].values
loo = dl_with_bonf(ds_loo, vs_loo, "LEAVE-ONE-OUT (k=5, drop GSE47533)")

# ---------------------------------------------------------------------------
# Version 3: Winsorized (cap GSE47533 at next-largest d)
# ---------------------------------------------------------------------------
outlier_idx = hyp.index[outlier_mask][0]
outlier_d = float(hyp.loc[outlier_idx, "cohens_d"])

# Find next-largest d (by absolute value, excluding the outlier)
other_ds = hyp.loc[~outlier_mask, "cohens_d"].values
next_largest = float(other_ds[np.argmax(np.abs(other_ds))])

print(f"\nOutlier GSE47533: g = {outlier_d:.4f}")
print(f"Next-largest |g| in program: {next_largest:.4f} (|g|={abs(next_largest):.4f})")
print(f"Winsorizing: replace {outlier_d:.4f} -> {next_largest:.4f}")

ds_winsor = hyp.cohens_d.values.copy()
ds_winsor[outlier_idx] = next_largest
# Also recompute variance for the Winsorized value using the same n
outlier_n1 = int(hyp.loc[outlier_idx, "case_n"])
outlier_n2 = int(hyp.loc[outlier_idx, "control_n"])
# Var(g) ~ (n1+n2)/(n1*n2) + g^2/(2*(n1+n2))
# Use a simple approximation: scale variance by (g_new/g_old)^2
# More precisely: recompute from formula
var_winsor = (outlier_n1 + outlier_n2) / (outlier_n1 * outlier_n2) + (next_largest ** 2) / (2 * (outlier_n1 + outlier_n2))
# Apply Hedges J correction to variance
J = 1.0 - 3.0 / (4 * (outlier_n1 + outlier_n2) - 9)
var_winsor_g = var_winsor * J * J

vs_winsor = hyp.cohens_d_var.values.copy()
vs_winsor[outlier_idx] = var_winsor_g

winsorized = dl_with_bonf(ds_winsor, vs_winsor, "WINSORIZED (GSE47533 capped at next-largest)")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY: Hypoxia within-program robustness")
print("-" * 70)
print(f"{'Version':<30} {'g':>8} {'p_DL':>12} {'p_Bonf':>12} {'Survives':>10}")
print("-" * 70)
for label, res in [("Original (k=6)", original), ("LOO w/o GSE47533 (k=5)", loo), ("Winsorized (k=6)", winsorized)]:
    surv = "YES" if res["survives_bonferroni"] else "NO"
    print(f"{label:<30} {res['pooled_g']:8.4f} {res['p_DL']:12.2e} {res['p_Bonferroni']:12.2e} {surv:>10}")

print(f"\nKey finding: Outlier effect on pooled g = {original['pooled_g']:.4f} -> {loo['pooled_g']:.4f} (LOO) / {winsorized['pooled_g']:.4f} (Winsorized)")
print(f"Direction: {'preserved' if np.sign(original['pooled_g']) == np.sign(loo['pooled_g']) else 'REVERSED'}")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
output = {
    "analysis": "winsorize_outlier",
    "signature": "hallmark_hypoxia",
    "program": "hypoxia",
    "outlier_cohort": "hypoxia_timecourse_gse47533",
    "outlier_g": round(outlier_d, 4),
    "winsorize_cap_g": round(next_largest, 4),
    "bonferroni_n_tests": BONF_N,
    "bonferroni_threshold": BONF_THRESHOLD,
    "original": original,
    "leave_one_out": loo,
    "winsorized": winsorized,
}

with open("outputs/canonical_v7/winsorize_outlier.json", "w") as f:
    json.dump(output, f, indent=2)
print("\nSaved to outputs/canonical_v7/winsorize_outlier.json")
