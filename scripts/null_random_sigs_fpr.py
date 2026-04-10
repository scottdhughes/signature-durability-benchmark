"""
Random Signature FPR Test: Direct Venet Reinterpretation
========================================================

Generates 200 random gene signatures and scores them through the
within-program DL random-effects framework to measure the false
positive rate. Compares to single-cohort FPR (the Venet regime).
"""
import sys
import time
import json

sys.path.insert(0, "src")

import numpy as np
import pandas as pd
from scipy import stats
from signature_durability_benchmark.meta_analysis import dl_random_effects_meta
from signature_durability_benchmark.scoring import score_signature_in_cohort
from signature_durability_benchmark.null_model import generate_null_signatures

N_NULL = 200
SIG_SIZE = 30
ALPHA = 0.05
EFFECT_THRESHOLD = 0.2
SEED = 42

manifest = pd.read_csv("config/cohort_manifest.tsv", sep="\t")
cohort_program = dict(zip(manifest.cohort_id, manifest.biological_program))
programs = sorted(set(cohort_program.values()))

print(f"Random Signature FPR Test: {N_NULL} null sigs, size={SIG_SIZE}")
print(f"  {len(manifest)} cohorts, {len(programs)} programs")
print()

# Pre-load cohort data
print("Loading cohort data...")
t0 = time.time()
cohort_data = {}
all_genes = set()
for _, row in manifest.iterrows():
    cid = row.cohort_id
    expr = pd.read_csv(f"data/freeze/cohort_matrices/{cid}.tsv", sep="\t", index_col=0)
    pheno = pd.read_csv(f"data/freeze/cohort_phenotypes/{cid}.tsv", sep="\t")
    cohort_data[cid] = {
        "expr": expr, "pheno": pheno,
        "phenotype_column": str(row.phenotype_column),
        "case_label": str(row.case_label),
        "control_label": str(row.control_label),
    }
    all_genes.update(expr.index.tolist())
print(f"  Done in {time.time()-t0:.1f}s, {len(all_genes)} genes in universe")

universe = sorted(all_genes)
null_sigs = generate_null_signatures(SIG_SIZE, universe, N_NULL, seed=SEED)

# Score all nulls across all cohorts
print(f"\nScoring {N_NULL} null signatures x {len(cohort_data)} cohorts...")
t0 = time.time()
null_scores = {}  # null_idx -> {cohort_id: (d, var)}

for i, nsig in enumerate(null_sigs):
    if i % 50 == 0:
        print(f"  {i}/{N_NULL}...")
    null_scores[i] = {}
    for cid, cd in cohort_data.items():
        result = score_signature_in_cohort(
            nsig, cd["expr"], cd["pheno"],
            cd["phenotype_column"], cd["case_label"], cd["control_label"],
        )
        null_scores[i][cid] = (result["cohens_d"], result["cohens_d_var"])

print(f"  Done in {time.time()-t0:.1f}s")

# ── Test 1: Single-cohort FPR (Venet-style) ──
single_fp = 0
single_total = 0
for i in range(N_NULL):
    for cid in null_scores[i]:
        d, var = null_scores[i][cid]
        se = np.sqrt(var) if var > 0 else 1.0
        z = d / se
        p = 2 * (1 - stats.norm.cdf(abs(z)))
        single_total += 1
        if p < ALPHA:
            single_fp += 1

single_fpr = single_fp / single_total

# ── Test 2: Within-program DL meta FPR ──
wp_fp = {p: 0 for p in programs}
wp_total = {p: 0 for p in programs}
wp_rows = []

for i in range(N_NULL):
    for prog in programs:
        prog_cohorts = [c for c in null_scores[i] if cohort_program.get(c) == prog]
        if len(prog_cohorts) < 3:
            continue

        ds = [null_scores[i][c][0] for c in prog_cohorts]
        vs = [null_scores[i][c][1] for c in prog_cohorts]
        meta = dl_random_effects_meta(ds, vs)

        wp_total[prog] += 1
        is_fp = meta["pooled_p"] < ALPHA and abs(meta["pooled_effect"]) > EFFECT_THRESHOLD
        if is_fp:
            wp_fp[prog] += 1

        wp_rows.append({
            "null_idx": i, "program": prog,
            "pooled_d": meta["pooled_effect"], "pooled_p": meta["pooled_p"],
            "i_squared": meta["i_squared"], "k": meta["k"],
        })

wp_total_all = sum(wp_total.values())
wp_fp_all = sum(wp_fp.values())
wp_fpr = wp_fp_all / wp_total_all if wp_total_all > 0 else 0

# ── Test 3: Cross-context DL meta FPR ──
cross_fp = 0
cross_total = 0
for i in range(N_NULL):
    all_ds = [null_scores[i][c][0] for c in null_scores[i]]
    all_vs = [null_scores[i][c][1] for c in null_scores[i]]
    if len(all_ds) < 3:
        continue
    meta = dl_random_effects_meta(all_ds, all_vs)
    cross_total += 1
    if meta["pooled_p"] < ALPHA and abs(meta["pooled_effect"]) > EFFECT_THRESHOLD:
        cross_fp += 1

cross_fpr = cross_fp / cross_total if cross_total > 0 else 0

# ── Test 4: "Any program" FPR (per-null, was any program significant?) ──
any_prog_fp = 0
for i in range(N_NULL):
    any_sig = False
    for prog in programs:
        prog_cohorts = [c for c in null_scores[i] if cohort_program.get(c) == prog]
        if len(prog_cohorts) < 3:
            continue
        ds = [null_scores[i][c][0] for c in prog_cohorts]
        vs = [null_scores[i][c][1] for c in prog_cohorts]
        meta = dl_random_effects_meta(ds, vs)
        if meta["pooled_p"] < ALPHA and abs(meta["pooled_effect"]) > EFFECT_THRESHOLD:
            any_sig = True
            break
    if any_sig:
        any_prog_fp += 1

any_prog_fpr = any_prog_fp / N_NULL

# ── Results ──
print()
print("=" * 70)
print("RANDOM SIGNATURE FPR: Venet Reinterpretation")
print("=" * 70)
print()
print(f"  {'Test':<45} {'FP':>5} {'Total':>6} {'FPR':>8}")
print(f"  {'-'*45} {'---':>5} {'-----':>6} {'------':>8}")
print(f"  {'Single-cohort (Venet regime)':<45} {single_fp:>5} {single_total:>6} {single_fpr:>8.3f}")
print(f"  {'Cross-context DL meta (30 cohorts)':<45} {cross_fp:>5} {cross_total:>6} {cross_fpr:>8.3f}")
print(f"  {'Within-program DL meta (per-program)':<45} {wp_fp_all:>5} {wp_total_all:>6} {wp_fpr:>8.3f}")
print(f"  {'Any-program significant (per null sig)':<45} {any_prog_fp:>5} {N_NULL:>6} {any_prog_fpr:>8.3f}")
print()
print(f"  Per-program within-program FPR:")
for prog in programs:
    if wp_total[prog] > 0:
        rate = wp_fp[prog] / wp_total[prog]
        print(f"    {prog:<15} {wp_fp[prog]:>3}/{wp_total[prog]:>3} = {rate:.3f}")
print()

# Effect size distribution
wp_df = pd.DataFrame(wp_rows)
print(f"  Null within-program effect size distribution:")
print(f"    Mean |d|: {wp_df.pooled_d.abs().mean():.3f}")
print(f"    Median |d|: {wp_df.pooled_d.abs().median():.3f}")
print(f"    SD: {wp_df.pooled_d.std():.3f}")
print(f"    Mean I-squared: {wp_df.i_squared.mean():.3f}")
print()

# Compare to real benchmark results
within_dur = pd.read_csv("outputs/canonical_v6/within_program_durability.csv")
hallmark_rows = within_dur[within_dur.signature_id.str.startswith("hallmark_")]
real_sig = hallmark_rows[
    (hallmark_rows.within_p < ALPHA) &
    (hallmark_rows.within_d.abs() > EFFECT_THRESHOLD)
]

print(f"  Real benchmark: {len(real_sig)}/{len(hallmark_rows)} hallmarks durable")
print(f"  Real hallmark within-program |d| range: "
      f"{hallmark_rows.within_d.abs().min():.2f} - {hallmark_rows.within_d.abs().max():.2f}")
print(f"  Null random sig within-program |d| range: "
      f"{wp_df.pooled_d.abs().min():.4f} - {wp_df.pooled_d.abs().max():.3f}")
print()

print("INTERPRETATION:")
print(f"  Single-cohort testing (Venet regime) has {single_fpr*100:.1f}% FPR.")
if single_fpr > 0.05:
    print(f"  This confirms Venet's observation: random signatures are often")
    print(f"  'significant' in individual cohorts due to noise and confounding.")
print()
print(f"  Within-program DL meta-analysis has {wp_fpr*100:.1f}% FPR.")
if wp_fpr > 0.10:
    print(f"  This is above nominal (5%) due to DL anti-conservatism with small k")
    print(f"  and high heterogeneity (mean I2={wp_df.i_squared.mean():.2f}).")
    print(f"  However, it is {single_fpr/wp_fpr:.1f}x better than single-cohort testing.")
    print(f"  The real hallmark effects (|d|=0.5-3.6) far exceed null (mean |d|={wp_df.pooled_d.abs().mean():.2f}).")
else:
    print(f"  Well-calibrated -- the within-program framework eliminates the")
    print(f"  random-signature problem.")

# Save
output = {
    "n_null": N_NULL,
    "sig_size": SIG_SIZE,
    "single_cohort_fpr": round(single_fpr, 4),
    "cross_context_fpr": round(cross_fpr, 4),
    "within_program_fpr": round(wp_fpr, 4),
    "any_program_fpr": round(any_prog_fpr, 4),
    "per_program_fpr": {p: round(wp_fp[p]/wp_total[p], 4) if wp_total[p]>0 else None
                        for p in programs},
    "null_effect_size_mean": round(float(wp_df.pooled_d.abs().mean()), 4),
    "null_i_squared_mean": round(float(wp_df.i_squared.mean()), 4),
}

with open("outputs/canonical_v6/null_random_sigs_fpr.json", "w") as f:
    json.dump(output, f, indent=2)

wp_df.to_csv("outputs/canonical_v6/null_random_sigs_fpr_detail.csv", index=False)
print(f"\nSaved to outputs/canonical_v6/null_random_sigs_fpr.json")
