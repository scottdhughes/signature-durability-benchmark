"""I-squared decomposition: between-program vs within-program heterogeneity.

Uses the standard Cochrane subgroup Q decomposition (Borenstein et al. 2009):
  Q_total = Q_within + Q_between

This decomposes the widely-reported I-squared into components attributable to
within-program variation (genuine measurement noise, cohort-level differences)
versus between-program variation (signatures activating differently in different
biological contexts).

Finding: if Q_between dominates Q_total, then published I-squared values for
gene signatures are primarily a *context artifact* rather than evidence of
poor signature quality.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json
import sys

sys.path.insert(0, "src")

effects = pd.read_csv("outputs/canonical_v6/per_cohort_effects.csv")
manifest = pd.read_csv("config/cohort_manifest.tsv", sep="\t")

# Map cohorts to programs
cohort_program = dict(zip(manifest.cohort_id, manifest.biological_program))

# The 7 Hallmark signatures and their "home" programs
hallmarks = {
    "hallmark_ifng_response": "interferon",
    "hallmark_ifna_response": "interferon",
    "hallmark_inflammatory_response": "inflammation",
    "hallmark_hypoxia": "hypoxia",
    "hallmark_e2f_targets": "proliferation",
    "hallmark_emt": "emt",
    "hallmark_tnfa_nfkb": "inflammation",
}

results = {}
print(f"{'Signature':<25} {'Q_total':>8} {'Q_within':>9} {'Q_between':>10} {'frac_btw':>9} {'I2_total':>9} {'I2_within':>10} {'p_between':>10}")
print("-" * 100)

for sig_id, sig_program in hallmarks.items():
    sig_effects = effects[effects.signature_id == sig_id].copy()

    ds = sig_effects.cohens_d.values
    vs = np.where(sig_effects.cohens_d_var.values < 1e-10, 1e-10, sig_effects.cohens_d_var.values)
    ws = 1.0 / vs
    k = len(ds)

    # Total pooled effect and Q
    d_pooled = np.sum(ws * ds) / np.sum(ws)
    Q_total = np.sum(ws * (ds - d_pooled) ** 2)
    df_total = k - 1
    I2_total = max(0, (Q_total - df_total) / Q_total) if Q_total > 0 else 0

    # Assign programs
    programs = np.array([cohort_program.get(c, "unknown") for c in sig_effects.cohort_id])
    unique_programs = sorted(set(programs))

    # Within-program Q: sum of subgroup Q values
    Q_within = 0.0
    df_within = 0
    program_stats = {}

    for prog in unique_programs:
        mask = programs == prog
        d_prog = ds[mask]
        w_prog = ws[mask]
        k_prog = len(d_prog)
        if k_prog < 2:
            program_stats[prog] = {
                "pooled": float(d_prog[0]) if k_prog == 1 else 0.0,
                "k": k_prog,
                "Q": 0.0,
                "weight_sum": float(np.sum(w_prog)),
            }
            continue

        d_prog_pooled = np.sum(w_prog * d_prog) / np.sum(w_prog)
        Q_prog = np.sum(w_prog * (d_prog - d_prog_pooled) ** 2)
        Q_within += Q_prog
        df_within += k_prog - 1
        program_stats[prog] = {
            "pooled": float(d_prog_pooled),
            "k": k_prog,
            "Q": float(Q_prog),
            "weight_sum": float(np.sum(w_prog)),
        }

    # Between-program Q
    Q_between = Q_total - Q_within
    n_subgroups = len([p for p in unique_programs if program_stats[p]["k"] >= 1])
    df_between = n_subgroups - 1

    # I-squared within (pooled across subgroups)
    I2_within = max(0, (Q_within - df_within) / Q_within) if Q_within > 0 and df_within > 0 else 0

    # Fraction of total Q explained by between-program
    frac_between = Q_between / Q_total if Q_total > 0 else 0

    # P-value for between-program heterogeneity (chi-squared test)
    p_between = 1.0 - stats.chi2.cdf(Q_between, df_between) if df_between > 0 else 1.0

    results[sig_id] = {
        "Q_total": round(float(Q_total), 2),
        "Q_within": round(float(Q_within), 2),
        "Q_between": round(float(Q_between), 2),
        "df_total": int(df_total),
        "df_within": int(df_within),
        "df_between": int(df_between),
        "I2_total": round(float(I2_total), 4),
        "I2_within": round(float(I2_within), 4),
        "fraction_Q_between": round(float(frac_between), 4),
        "p_between": float(p_between),
        "program_means": {
            k: {kk: round(vv, 4) for kk, vv in v.items()}
            for k, v in program_stats.items()
        },
    }

    short = sig_id.replace("hallmark_", "")
    print(
        f"{short:<25} {Q_total:8.1f} {Q_within:9.1f} {Q_between:10.1f} "
        f"{frac_between:9.1%} {I2_total:9.4f} {I2_within:10.4f} {p_between:10.2e}"
    )

# Summary statistics
fracs = [r["fraction_Q_between"] for r in results.values()]
print("\n" + "=" * 100)
print(f"Mean fraction of Q_total explained by between-program: {np.mean(fracs):.1%}")
print(f"Median: {np.median(fracs):.1%}")
print(f"Range: {np.min(fracs):.1%} - {np.max(fracs):.1%}")

i2_totals = [r["I2_total"] for r in results.values()]
i2_withins = [r["I2_within"] for r in results.values()]
print(f"\nMean I2_total (reported): {np.mean(i2_totals):.2f}")
print(f"Mean I2_within (context-adjusted): {np.mean(i2_withins):.2f}")
print(f"Reduction: {np.mean(i2_totals) - np.mean(i2_withins):.2f} ({(np.mean(i2_totals) - np.mean(i2_withins))/np.mean(i2_totals):.0%} of total)")

# Check p-values
sig_between = sum(1 for r in results.values() if r["p_between"] < 0.05)
print(f"\nSignatures with significant between-program heterogeneity (p<0.05): {sig_between}/{len(results)}")

with open("outputs/canonical_v6/i2_decomposition.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nSaved to outputs/canonical_v6/i2_decomposition.json")
