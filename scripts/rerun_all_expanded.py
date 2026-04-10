#!/usr/bin/env python3
"""Re-run all meta-analyses on the expanded panel:
- 35 cohorts (30 original + 5 IFN expansion)
- 30 signatures (29 + 1 new Schoggins 2011 IRG)
- 11 interferon cohorts (6 original + 5 new)

Outputs to outputs/canonical_v8/:
- i2_decomposition_expanded.json
- within_program_expanded.csv
- hartung_knapp_expanded.json
- permutation_validation_expanded.json
- lopo_cross_validation_expanded.json
"""

import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats as scipy_stats
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import LeaveOneOut

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs" / "canonical_v8"
OUT.mkdir(exist_ok=True)

SEED = 42
N_PERMS = 10000

# ═══════════════════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════════════════

per_cohort = pd.read_csv(OUT / "per_cohort_effects.csv")
manifest = pd.read_csv(ROOT / "config" / "cohort_manifest.tsv", sep="\t")
cohort_to_program = dict(zip(manifest["cohort_id"], manifest["biological_program"]))

print(f"Loaded {len(per_cohort)} effect rows for {per_cohort['signature_id'].nunique()} signatures × {per_cohort['cohort_id'].nunique()} cohorts")
print(f"Programs: {dict(pd.Series(list(cohort_to_program.values())).value_counts())}")

# Target signatures for the paper
TARGET_SIGS = pd.read_csv(ROOT / "config" / "paper_target_signatures.tsv", sep="\t")["signature_id"].tolist()

PROGRAMS = ["inflammation", "interferon", "proliferation", "hypoxia", "emt"]

SIG_PROGRAM = {
    "hallmark_ifng_response": "interferon",
    "hallmark_ifna_response": "interferon",
    "schoggins_2011_irg": "interferon",
    "blind_durable_ifn_composite": "interferon",
    "hallmark_inflammatory_response": "inflammation",
    "hallmark_tnfa_nfkb": "inflammation",
    "hallmark_hypoxia": "hypoxia",
    "hallmark_e2f_targets": "proliferation",
    "hallmark_emt": "emt",
}


# ═══════════════════════════════════════════════════════════════════════════════
# I² Decomposition
# ═══════════════════════════════════════════════════════════════════════════════

def i2_decompose(df, sig):
    """Return (Q_total, Q_within, Q_between, Q_B_frac, I2_total, I2_within, p_between)."""
    d = df[df["signature_id"] == sig].copy()
    d = d.merge(manifest[["cohort_id", "biological_program"]], on="cohort_id")

    if len(d) < 4:
        return None

    d["w"] = 1.0 / d["cohens_d_var"].replace(0, 1e-9)
    d["wd"] = d["w"] * d["cohens_d"]

    grand_mean = d["wd"].sum() / d["w"].sum()
    Q_total = (d["w"] * (d["cohens_d"] - grand_mean) ** 2).sum()

    Q_within = 0.0
    for prog, group in d.groupby("biological_program"):
        if len(group) < 2:
            continue
        pm = (group["w"] * group["cohens_d"]).sum() / group["w"].sum()
        Q_within += (group["w"] * (group["cohens_d"] - pm) ** 2).sum()

    Q_between = Q_total - Q_within
    K = d["biological_program"].nunique()
    df_Q = len(d) - 1
    df_W = len(d) - K
    df_B = K - 1

    I2_total = max(0, (Q_total - df_Q) / Q_total) if Q_total > 0 else 0
    I2_within = max(0, (Q_within - df_W) / Q_within) if Q_within > 0 else 0
    p_between = 1 - scipy_stats.chi2.cdf(Q_between, df_B) if df_B > 0 else np.nan

    return {
        "Q_total": round(Q_total, 3),
        "Q_within": round(Q_within, 3),
        "Q_between": round(Q_between, 3),
        "Q_B_fraction": round(Q_between / Q_total, 4) if Q_total > 0 else 0,
        "I2_total": round(I2_total, 4),
        "I2_within": round(I2_within, 4),
        "p_between": float(p_between),
        "k_cohorts": len(d),
        "k_programs": K,
    }


print("\n=== I² Decomposition (Expanded Panel) ===")
i2_results = {}
for sig in TARGET_SIGS:
    r = i2_decompose(per_cohort, sig)
    if r:
        i2_results[sig] = r
        print(f"  {sig:35s}: Q_B/Q_tot = {r['Q_B_fraction']:.3f}  (k={r['k_cohorts']}, programs={r['k_programs']})")

with open(OUT / "i2_decomposition_expanded.json", "w") as f:
    json.dump(i2_results, f, indent=2)


# ═══════════════════════════════════════════════════════════════════════════════
# Within-program meta-analysis (DL + HKSJ)
# ═══════════════════════════════════════════════════════════════════════════════

def meta_analysis(gs, vars_, k):
    """DerSimonian-Laird + HKSJ (guarded) meta-analysis."""
    w = 1.0 / np.array(vars_)
    weighted_mean = np.sum(w * gs) / np.sum(w)
    Q = np.sum(w * (gs - weighted_mean) ** 2)
    df = k - 1
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0, (Q - df) / c) if c > 0 else 0
    I2 = max(0, (Q - df) / Q) if Q > 0 else 0

    # DL
    w_star = 1 / (np.array(vars_) + tau2)
    pooled_g = np.sum(w_star * gs) / np.sum(w_star)
    se_DL = np.sqrt(1 / np.sum(w_star))
    z_DL = pooled_g / se_DL
    p_DL = 2 * (1 - scipy_stats.norm.cdf(abs(z_DL)))

    # HKSJ
    q_star = np.sum(w_star * (gs - pooled_g) ** 2)
    se_HKSJ = np.sqrt(q_star / ((k - 1) * np.sum(w_star)))
    se_HKSJ_guarded = max(se_DL, se_HKSJ)

    t_HKSJ = pooled_g / se_HKSJ
    p_HKSJ = 2 * (1 - scipy_stats.t.cdf(abs(t_HKSJ), df=k - 1))

    t_guarded = pooled_g / se_HKSJ_guarded
    p_guarded = 2 * (1 - scipy_stats.t.cdf(abs(t_guarded), df=k - 1))

    return {
        "k": k,
        "pooled_g": round(float(pooled_g), 4),
        "tau2": round(float(tau2), 4),
        "I2": round(float(I2), 4),
        "DL_se": round(float(se_DL), 4),
        "DL_p": float(p_DL),
        "HKSJ_se": round(float(se_HKSJ), 4),
        "HKSJ_p": float(p_HKSJ),
        "HKSJ_guarded_se": round(float(se_HKSJ_guarded), 4),
        "HKSJ_guarded_p": float(p_guarded),
        "p_bonf_HKSJ_guarded": float(min(1.0, p_guarded * 9)),  # 9 tests
    }


def within_program_meta(df, sig, program):
    """Pool within a specific program only."""
    d = df[df["signature_id"] == sig].merge(manifest[["cohort_id", "biological_program"]], on="cohort_id")
    d = d[d["biological_program"] == program]
    if len(d) < 3:
        return None
    return meta_analysis(d["cohens_d"].values, d["cohens_d_var"].values, len(d))


print("\n=== Within-Interferon Meta-Analysis (Expanded k=11) ===")
hksj_results = {}
for sig in TARGET_SIGS:
    r = within_program_meta(per_cohort, sig, "interferon")
    if r:
        hksj_results[sig] = r
        bonf = r["p_bonf_HKSJ_guarded"]
        sig_marker = "***" if bonf < 0.001 else "**" if bonf < 0.01 else "*" if bonf < 0.05 else ""
        print(f"  {sig:35s}: g={r['pooled_g']:+.3f}  HKSJ-guarded p_bonf={bonf:.4f} {sig_marker}  (k={r['k']}, I²={r['I2']:.2f})")

# Also compute outside-program
print("\n=== Outside-Interferon (for comparison) ===")
for sig in TARGET_SIGS:
    d = per_cohort[per_cohort["signature_id"] == sig].merge(manifest[["cohort_id", "biological_program"]], on="cohort_id")
    d_out = d[d["biological_program"] != "interferon"]
    if len(d_out) >= 3:
        r = meta_analysis(d_out["cohens_d"].values, d_out["cohens_d_var"].values, len(d_out))
        hksj_results[sig + "__outside_interferon"] = r
        print(f"  {sig:35s}: g={r['pooled_g']:+.3f}  p={r['HKSJ_p']:.4f}  (k={r['k']})")

with open(OUT / "hartung_knapp_expanded.json", "w") as f:
    json.dump(hksj_results, f, indent=2)


# ═══════════════════════════════════════════════════════════════════════════════
# Permutation test (10K iterations)
# ═══════════════════════════════════════════════════════════════════════════════

def permutation_test(df, sig, n_perms=N_PERMS, seed=SEED):
    d = df[df["signature_id"] == sig].merge(manifest[["cohort_id", "biological_program"]], on="cohort_id")
    if len(d) < 4:
        return None

    # Observed
    obs = i2_decompose(df, sig)
    obs_frac = obs["Q_B_fraction"]

    # Permute program labels across cohorts
    rng = np.random.default_rng(seed)
    programs = d["biological_program"].values.copy()
    null_fracs = []
    for _ in range(n_perms):
        perm_programs = rng.permutation(programs)
        d["biological_program_perm"] = perm_programs

        d["w"] = 1.0 / d["cohens_d_var"].replace(0, 1e-9)
        grand = (d["w"] * d["cohens_d"]).sum() / d["w"].sum()
        Q_total = (d["w"] * (d["cohens_d"] - grand) ** 2).sum()

        Q_w = 0.0
        for _, grp in d.groupby("biological_program_perm"):
            if len(grp) < 2:
                continue
            pm = (grp["w"] * grp["cohens_d"]).sum() / grp["w"].sum()
            Q_w += (grp["w"] * (grp["cohens_d"] - pm) ** 2).sum()
        Q_b = Q_total - Q_w
        null_fracs.append(Q_b / Q_total if Q_total > 0 else 0)

    null_arr = np.array(null_fracs)
    p = (np.sum(null_arr >= obs_frac) + 1) / (n_perms + 1)
    return {
        "observed_Q_B_fraction": obs_frac,
        "null_mean": round(float(np.mean(null_arr)), 4),
        "null_95th": round(float(np.percentile(null_arr, 95)), 4),
        "null_99th": round(float(np.percentile(null_arr, 99)), 4),
        "p_value": float(p),
    }


print("\n=== Permutation Test (10K iterations) ===")
perm_results = {}
for sig in TARGET_SIGS:
    r = permutation_test(per_cohort, sig)
    if r:
        perm_results[sig] = r
        sig_marker = "***" if r["p_value"] < 0.001 else "**" if r["p_value"] < 0.01 else "*" if r["p_value"] < 0.05 else ""
        print(f"  {sig:35s}: obs={r['observed_Q_B_fraction']:.3f}  null_mean={r['null_mean']:.3f}  p={r['p_value']:.4f} {sig_marker}")

with open(OUT / "permutation_validation_expanded.json", "w") as f:
    json.dump(perm_results, f, indent=2)


# ═══════════════════════════════════════════════════════════════════════════════
# LOPO-CV
# ═══════════════════════════════════════════════════════════════════════════════

print("\n=== Leave-One-Program-Out Q_B Stability ===")
lopo_results = {}
for sig in TARGET_SIGS:
    if sig not in i2_results:
        continue
    full_qb = i2_results[sig]["Q_B_fraction"]
    lopo_results[sig] = {"full": full_qb, "loo": {}}
    for held_out in PROGRAMS:
        cohorts_kept = [c for c, p in cohort_to_program.items() if p != held_out]
        d_loo = per_cohort[(per_cohort["signature_id"] == sig) & (per_cohort["cohort_id"].isin(cohorts_kept))]
        if len(d_loo) < 4:
            continue
        # compute Q_B/Q_tot on subset
        d_loo = d_loo.merge(manifest[["cohort_id", "biological_program"]], on="cohort_id")
        d_loo["w"] = 1.0 / d_loo["cohens_d_var"].replace(0, 1e-9)
        grand = (d_loo["w"] * d_loo["cohens_d"]).sum() / d_loo["w"].sum()
        Q_total = (d_loo["w"] * (d_loo["cohens_d"] - grand) ** 2).sum()
        Q_w = 0.0
        for _, grp in d_loo.groupby("biological_program"):
            if len(grp) < 2:
                continue
            pm = (grp["w"] * grp["cohens_d"]).sum() / grp["w"].sum()
            Q_w += (grp["w"] * (grp["cohens_d"] - pm) ** 2).sum()
        qb_loo = (Q_total - Q_w) / Q_total if Q_total > 0 else 0
        lopo_results[sig]["loo"][held_out] = round(qb_loo, 4)

    sig_prog = SIG_PROGRAM.get(sig, "?")
    print(f"\n  {sig} (on-program = {sig_prog}):")
    print(f"    Full: {full_qb:.3f}")
    for held_out, v in lopo_results[sig]["loo"].items():
        delta = v - full_qb
        marker = " <-- ON-PROGRAM" if held_out == sig_prog else ""
        print(f"    Hide {held_out:15s}: {v:.3f} (Δ={delta:+.3f}){marker}")

with open(OUT / "lopo_cross_validation_expanded.json", "w") as f:
    json.dump(lopo_results, f, indent=2)

print("\n\n=== DONE ===")
print(f"Results in: {OUT}")
