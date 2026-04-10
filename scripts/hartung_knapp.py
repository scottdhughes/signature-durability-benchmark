"""Hartung-Knapp-Sidik-Jonkman (HKSJ) correction for small-k meta-analyses.

The standard DerSimonian-Laird random-effects meta-analysis uses a normal
distribution for the pooled effect CI. With small k (number of studies),
this substantially underestimates the CI width and inflates Type I errors.

The HKSJ correction replaces the DL standard error with a t-distribution-
based estimate that accounts for the uncertainty in tau-squared estimation.
This produces wider but more honest confidence intervals.

We apply this to IFN-gamma and IFN-alpha within the interferon program (k=6),
which are the key results. The question: do they still survive Bonferroni
under the more conservative HKSJ inference?
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

BONF_N = 9
BONF_THRESHOLD = 0.05 / BONF_N

# ---------------------------------------------------------------------------
# HKSJ meta-analysis
# ---------------------------------------------------------------------------
def hksj_meta(ds, vs, label):
    """Compute DL random-effects pooled effect, then apply HKSJ correction.

    Steps:
    1. Compute DL tau-squared and random-effects weights w_RE = 1/(v_i + tau2)
    2. Compute RE pooled effect: d_RE = sum(w_RE * d_i) / sum(w_RE)
    3. HKSJ SE: q = sum(w_RE * (d_i - d_RE)^2); SE_HK = sqrt(q / ((k-1) * sum(w_RE)))
    4. Use t(k-1) distribution for CI and p-value
    """
    ds = np.array(ds, dtype=float)
    vs = np.where(np.array(vs, dtype=float) < 1e-10, 1e-10, np.array(vs, dtype=float))
    k = len(ds)

    # Step 1: DL tau-squared estimation
    dl_result = dl_random_effects_meta(list(ds), list(vs))
    tau2 = dl_result["tau2"]

    # Step 2: Random-effects weights and pooled effect
    w_re = 1.0 / (vs + tau2)
    d_re = np.sum(w_re * ds) / np.sum(w_re)
    se_dl = 1.0 / np.sqrt(np.sum(w_re))

    # Step 3: HKSJ SE
    q_hk = np.sum(w_re * (ds - d_re) ** 2)
    se_hk = np.sqrt(q_hk / ((k - 1) * np.sum(w_re)))

    # Ensure HKSJ SE is not smaller than DL SE (conservative guard)
    # Some implementations enforce this; we report both
    se_hk_guarded = max(se_hk, se_dl)

    # Step 4: t-distribution inference
    df = k - 1
    t_stat = d_re / se_hk if se_hk > 0 else 0.0
    t_stat_guarded = d_re / se_hk_guarded if se_hk_guarded > 0 else 0.0

    p_hk = float(2 * (1 - stats.t.cdf(abs(t_stat), df)))
    p_hk_guarded = float(2 * (1 - stats.t.cdf(abs(t_stat_guarded), df)))

    # CIs
    t_crit = stats.t.ppf(0.975, df)
    ci_hk = (d_re - t_crit * se_hk, d_re + t_crit * se_hk)
    ci_hk_guarded = (d_re - t_crit * se_hk_guarded, d_re + t_crit * se_hk_guarded)
    ci_dl = (d_re - 1.96 * se_dl, d_re + 1.96 * se_dl)

    # Bonferroni
    p_dl = dl_result["pooled_p"]
    p_bonf_dl = min(1.0, p_dl * BONF_N)
    p_bonf_hk = min(1.0, p_hk * BONF_N)
    p_bonf_hk_g = min(1.0, p_hk_guarded * BONF_N)

    print(f"\n{label} (k={k}):")
    print(f"  Pooled g = {d_re:.4f}")
    print(f"  tau2 = {tau2:.4f}, I2 = {dl_result['i_squared']:.4f}")
    print(f"  DL:   SE={se_dl:.4f},  z={d_re/se_dl:.3f},  p={p_dl:.2e},  95% CI=({ci_dl[0]:.4f}, {ci_dl[1]:.4f})")
    print(f"  HKSJ: SE={se_hk:.4f},  t={t_stat:.3f},  p={p_hk:.2e},  95% CI=({ci_hk[0]:.4f}, {ci_hk[1]:.4f})")
    if se_hk_guarded > se_hk:
        print(f"  HKSJ (guarded): SE={se_hk_guarded:.4f},  t={t_stat_guarded:.3f},  p={p_hk_guarded:.2e},  95% CI=({ci_hk_guarded[0]:.4f}, {ci_hk_guarded[1]:.4f})")
    print(f"  Bonferroni (n={BONF_N}): DL p={p_bonf_dl:.2e}, HKSJ p={p_bonf_hk:.2e}")
    print(f"  Survives Bonferroni? DL={'YES' if p_bonf_dl < 0.05 else 'NO'}, HKSJ={'YES' if p_bonf_hk < 0.05 else 'NO'}")

    return {
        "k": k,
        "pooled_g": round(float(d_re), 4),
        "tau2": round(float(tau2), 4),
        "i_squared": round(float(dl_result["i_squared"]), 4),
        "DL": {
            "se": round(float(se_dl), 4),
            "z": round(float(d_re / se_dl), 4),
            "p": p_dl,
            "ci_95": [round(ci_dl[0], 4), round(ci_dl[1], 4)],
            "p_Bonferroni": p_bonf_dl,
            "survives_bonferroni": p_bonf_dl < 0.05,
        },
        "HKSJ": {
            "se": round(float(se_hk), 4),
            "t_stat": round(float(t_stat), 4),
            "df": df,
            "p": p_hk,
            "ci_95": [round(ci_hk[0], 4), round(ci_hk[1], 4)],
            "p_Bonferroni": p_bonf_hk,
            "survives_bonferroni": p_bonf_hk < 0.05,
        },
        "HKSJ_guarded": {
            "se": round(float(se_hk_guarded), 4),
            "t_stat": round(float(t_stat_guarded), 4),
            "df": df,
            "p": p_hk_guarded,
            "ci_95": [round(ci_hk_guarded[0], 4), round(ci_hk_guarded[1], 4)],
            "p_Bonferroni": p_bonf_hk_g,
            "survives_bonferroni": p_bonf_hk_g < 0.05,
        },
        "cohort_gs": None,  # filled below
    }


# ---------------------------------------------------------------------------
# IFN-gamma within interferon
# ---------------------------------------------------------------------------
sigs = {
    "hallmark_ifng_response": "IFN-gamma",
    "hallmark_ifna_response": "IFN-alpha",
}

results = {}

for sig_id, short_name in sigs.items():
    mask = (effects.signature_id == sig_id) & (
        effects.cohort_id.map(cohort_program) == "interferon"
    )
    sub = effects[mask].copy().reset_index(drop=True)

    print(f"\n{'='*70}")
    print(f"{short_name} within interferon (k={len(sub)})")
    print(sub[["cohort_id", "cohens_d", "cohens_d_var", "case_n", "control_n"]].to_string(index=False))

    ds = sub.cohens_d.values
    vs = sub.cohens_d_var.values

    res = hksj_meta(ds, vs, f"{short_name} within interferon")
    res["cohort_gs"] = [
        {"cohort": row.cohort_id, "g": round(float(row.cohens_d), 4), "var": round(float(row.cohens_d_var), 4)}
        for _, row in sub.iterrows()
    ]
    results[sig_id] = res


# ---------------------------------------------------------------------------
# Also run for ALL 7 hallmarks within their home program
# ---------------------------------------------------------------------------
hallmarks = {
    "hallmark_ifng_response": "interferon",
    "hallmark_ifna_response": "interferon",
    "hallmark_inflammatory_response": "inflammation",
    "hallmark_hypoxia": "hypoxia",
    "hallmark_e2f_targets": "proliferation",
    "hallmark_emt": "emt",
    "hallmark_tnfa_nfkb": "inflammation",
}

print("\n\n" + "=" * 70)
print("ALL HALLMARKS: DL vs HKSJ comparison")
print("=" * 70)
print(f"{'Signature':<30} {'k':>3} {'g':>7} {'p_DL':>12} {'p_HKSJ':>12} {'p_Bonf_DL':>12} {'p_Bonf_HK':>12} {'DL?':>4} {'HK?':>4}")
print("-" * 105)

for sig_id, prog in hallmarks.items():
    mask = (effects.signature_id == sig_id) & (
        effects.cohort_id.map(cohort_program) == prog
    )
    sub = effects[mask]
    if len(sub) < 2:
        continue

    ds = sub.cohens_d.values
    vs = sub.cohens_d_var.values

    # Quick DL + HKSJ
    dl = dl_random_effects_meta(list(ds), list(vs))
    tau2 = dl["tau2"]
    w_re = 1.0 / (np.where(vs < 1e-10, 1e-10, vs) + tau2)
    d_re = np.sum(w_re * ds) / np.sum(w_re)
    q_hk = np.sum(w_re * (ds - d_re) ** 2)
    k = len(ds)
    se_hk = np.sqrt(q_hk / ((k - 1) * np.sum(w_re)))
    se_dl = 1.0 / np.sqrt(np.sum(w_re))
    se_hk = max(se_hk, se_dl)  # guarded
    t_stat = d_re / se_hk if se_hk > 0 else 0.0
    p_hk = float(2 * (1 - stats.t.cdf(abs(t_stat), k - 1)))

    p_dl = dl["pooled_p"]
    p_bonf_dl = min(1.0, p_dl * BONF_N)
    p_bonf_hk = min(1.0, p_hk * BONF_N)

    short = sig_id.replace("hallmark_", "")
    dl_surv = "YES" if p_bonf_dl < 0.05 else "no"
    hk_surv = "YES" if p_bonf_hk < 0.05 else "no"
    print(f"{short:<30} {k:3d} {d_re:7.3f} {p_dl:12.2e} {p_hk:12.2e} {p_bonf_dl:12.2e} {p_bonf_hk:12.2e} {dl_surv:>4} {hk_surv:>4}")

    if sig_id not in results:
        results[sig_id] = {
            "k": k,
            "pooled_g": round(float(d_re), 4),
            "program": prog,
            "DL_p": p_dl,
            "DL_p_Bonferroni": p_bonf_dl,
            "HKSJ_guarded_p": p_hk,
            "HKSJ_guarded_p_Bonferroni": p_bonf_hk,
            "survives_DL_Bonferroni": p_bonf_dl < 0.05,
            "survives_HKSJ_Bonferroni": p_bonf_hk < 0.05,
        }


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
output = {
    "analysis": "hartung_knapp_sidik_jonkman",
    "bonferroni_n_tests": BONF_N,
    "bonferroni_threshold": BONF_THRESHOLD,
    "note": "All effect sizes are Hedges' g (small-sample-corrected Cohen's d)",
    "results": results,
}

with open("outputs/canonical_v7/hartung_knapp.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to outputs/canonical_v7/hartung_knapp.json")
