#!/usr/bin/env python3
"""Descriptive within-IFN heterogeneity decomposition.

Reviewer concern: within-IFN I² remains 0.91-0.93 even after the expansion,
suggesting substantial residual heterogeneity. Decompose this into:
- Platform (Affymetrix vs Illumina)
- Tissue (blood/PBMC vs skin)
- Etiology (viral vs autoimmune)

This script uses inverse-variance weighted subgroup partitions of Q for the
IFN-gamma core. With only 11 cohorts and several sparse combined cells, the
results are descriptive rather than a definitive mixed-effects meta-regression.
"""

import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats as scipy_stats

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs" / "canonical_v8"

# Load data
per_cohort = pd.read_csv(OUT / "per_cohort_effects.csv")
manifest = pd.read_csv(ROOT / "config" / "cohort_manifest.tsv", sep="\t")

# Filter to IFN cohorts
ifn_cohorts = manifest[manifest["biological_program"] == "interferon"].copy()
print(f"IFN cohorts: {len(ifn_cohorts)}")

# Categorize tissue and etiology
def categorize_tissue(row):
    t = row["tissue"].lower()
    if "skin" in t:
        return "skin"
    else:
        return "blood_pbmc"

def categorize_etiology(row):
    cid = row["cohort_id"]
    if "sle" in cid or "psoriasis" in cid:
        return "autoimmune"
    else:
        return "viral"

def categorize_platform_family(row):
    p = row["platform_id"]
    if "GPL570" in p or "GPL571" in p or "GPL13158" in p:
        return "Affymetrix"
    elif "GPL10558" in p or "GPL6947" in p:
        return "Illumina"
    else:
        return "Other"

ifn_cohorts["tissue_cat"] = ifn_cohorts.apply(categorize_tissue, axis=1)
ifn_cohorts["etiology"] = ifn_cohorts.apply(categorize_etiology, axis=1)
ifn_cohorts["platform_family"] = ifn_cohorts.apply(categorize_platform_family, axis=1)

print(f"\nTissue: {ifn_cohorts['tissue_cat'].value_counts().to_dict()}")
print(f"Etiology: {ifn_cohorts['etiology'].value_counts().to_dict()}")
print(f"Platform: {ifn_cohorts['platform_family'].value_counts().to_dict()}")

# Get IFN-γ core effects for each IFN cohort
ifng_data = per_cohort[per_cohort["signature_id"] == "hallmark_ifng_response"].merge(
    ifn_cohorts[["cohort_id", "tissue_cat", "etiology", "platform_family"]], on="cohort_id"
)


def partition_summary(gs, vars_):
    """Inverse-variance summary used for weighted-Q partitions, not inferential pooling."""
    gs = np.array(gs)
    vars_ = np.array(vars_)
    k = len(gs)
    w = 1.0 / vars_
    wm = np.sum(w * gs) / np.sum(w)
    Q = np.sum(w * (gs - wm) ** 2)
    df = k - 1
    I2 = max(0, (Q - df) / Q) if Q > 0 else 0
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0, (Q - df) / c) if c > 0 else 0
    return {"k": k, "pooled_g": round(wm, 4), "Q": round(Q, 3), "I2": round(I2, 4), "tau2": round(tau2, 4)}


# Full IFN meta
print("\n=== Full IFN (k=11) ===")
full = partition_summary(ifng_data["cohens_d"], ifng_data["cohens_d_var"])
print(f"  inverse-variance mean g = {full['pooled_g']}, Q = {full['Q']}, I² = {full['I2']}, tau² = {full['tau2']}")

# Subgroup by tissue
print("\n=== Subgroup: Blood/PBMC vs Skin ===")
subgroup_results = {}
for cat_col, cat_name in [("tissue_cat", "Tissue"), ("etiology", "Etiology"), ("platform_family", "Platform")]:
    subgroup_results[cat_col] = {}
    print(f"\n--- {cat_name} ---")
    for level in ifng_data[cat_col].unique():
        sub = ifng_data[ifng_data[cat_col] == level]
        if len(sub) >= 2:
            r = partition_summary(sub["cohens_d"], sub["cohens_d_var"])
            subgroup_results[cat_col][level] = r
            print(f"  {level:15s}: k={r['k']}, g={r['pooled_g']:+.3f}, Q={r['Q']:.2f}, I²={r['I2']:.3f}")


# Descriptive weighted Q partition: how much of Q is aligned with each moderator?
def meta_regression_r2(gs, vars_, moderator_dummy):
    """Simple inverse-variance weighted Q partition R² for a moderator."""
    w = 1.0 / np.array(vars_)
    gs = np.array(gs)
    mod = np.array(moderator_dummy)
    k = len(gs)

    # Weighted mean
    mu = np.sum(w * gs) / np.sum(w)
    # Total SS
    ss_total = np.sum(w * (gs - mu) ** 2)
    # Weighted mean within each level
    levels = np.unique(mod)
    ss_within = 0
    for lv in levels:
        idx = mod == lv
        if idx.sum() >= 1:
            sub_mu = np.sum(w[idx] * gs[idx]) / np.sum(w[idx])
            ss_within += np.sum(w[idx] * (gs[idx] - sub_mu) ** 2)
    ss_between = ss_total - ss_within
    r2 = ss_between / ss_total if ss_total > 0 else 0
    return {"Q_total": round(ss_total, 3), "Q_within": round(ss_within, 3),
            "Q_between": round(ss_between, 3), "R2": round(r2, 4)}


print("\n\n=== Descriptive Decomposition: Variance Explained by Moderators ===")
regression_results = {}
for col, name in [("tissue_cat", "Tissue (blood/PBMC vs skin)"),
                   ("etiology", "Etiology (viral vs autoimmune)"),
                   ("platform_family", "Platform (Affymetrix vs Illumina)")]:
    r = meta_regression_r2(ifng_data["cohens_d"], ifng_data["cohens_d_var"], ifng_data[col])
    regression_results[col] = r
    print(f"\n{name}:")
    print(f"  Q_total = {r['Q_total']}, Q_between = {r['Q_between']} ({r['R2']*100:.1f}% of variance)")


# Combined model: all three moderators
# Simple approach: iterate levels
print("\n=== Within vs Between: Combined tissue × etiology × platform ===")
ifng_data["combined"] = ifng_data["tissue_cat"] + "_" + ifng_data["etiology"] + "_" + ifng_data["platform_family"]
print(ifng_data.groupby("combined")["cohens_d"].agg(["mean", "std", "count"]).to_string())
r_combined = meta_regression_r2(ifng_data["cohens_d"], ifng_data["cohens_d_var"], ifng_data["combined"])
print(f"\nCombined Q_between = {r_combined['Q_between']} ({r_combined['R2']*100:.1f}% of variance)")

# Residual after all moderators
residual_I2 = 1 - r_combined["R2"]
print(f"Residual variance (not explained by tissue/etiology/platform) = {residual_I2*100:.1f}%")

output = {
    "description": "Descriptive weighted heterogeneity decomposition for within-IFN IFN-gamma effects",
    "full_ifn_partition_summary": full,
    "full_ifn_partition_summary_note": (
        "This inverse-variance weighted mean is used only for descriptive weighted-Q "
        "partitions in this file. The inferential guarded DL/HKSJ pooled IFN-gamma "
        "effect is reported separately in outputs/canonical_v8/hartung_knapp_expanded.json."
    ),
    "subgroup_analyses": subgroup_results,
    "variance_explained_by_moderator": regression_results,
    "combined_r2": r_combined,
    "interpretation": (
        f"Within-IFN I² of {full['I2']} for the IFN-gamma core is descriptively partitioned into: "
        f"tissue explains {regression_results['tissue_cat']['R2']*100:.1f}%, "
        f"etiology explains {regression_results['etiology']['R2']*100:.1f}%, "
        f"platform explains {regression_results['platform_family']['R2']*100:.1f}%. "
        f"Combined moderators explain {r_combined['R2']*100:.1f}% of total within-IFN Q. "
        "Because k=11 and the combined partition contains sparse cells, this is "
        "a descriptive heterogeneity decomposition rather than a definitive "
        "inferential meta-regression."
    ),
}

with open(OUT / "within_ifn_metaregression.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\n✓ Saved to {OUT}/within_ifn_metaregression.json")
