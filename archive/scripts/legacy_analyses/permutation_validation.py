import sys
sys.path.insert(0, 'src')
import pandas as pd
import numpy as np
import json

effects = pd.read_csv("outputs/canonical_v7/per_cohort_effects.csv")
manifest = pd.read_csv("config/cohort_manifest.tsv", sep="\t")
cohort_program = dict(zip(manifest.cohort_id, manifest.biological_program))

hallmarks = [
    "hallmark_ifng_response", "hallmark_ifna_response",
    "hallmark_inflammatory_response", "hallmark_hypoxia",
    "hallmark_e2f_targets", "hallmark_emt", "hallmark_tnfa_nfkb"
]

def compute_qb_fraction(ds, vs, programs):
    """Compute Q_B/Q_total for a given program assignment."""
    ds, vs = np.array(ds), np.array(vs)
    vs = np.where(vs < 1e-10, 1e-10, vs)
    ws = 1.0 / vs

    # Total pooled and Q_total
    d_pooled = np.sum(ws * ds) / np.sum(ws)
    Q_total = np.sum(ws * (ds - d_pooled)**2)

    # Q_within
    Q_within = 0
    unique_progs = sorted(set(programs))
    for prog in unique_progs:
        mask = np.array([p == prog for p in programs])
        if mask.sum() < 2:
            continue
        d_prog = ds[mask]
        w_prog = ws[mask]
        d_prog_pooled = np.sum(w_prog * d_prog) / np.sum(w_prog)
        Q_within += np.sum(w_prog * (d_prog - d_prog_pooled)**2)

    Q_between = Q_total - Q_within
    return Q_between / Q_total if Q_total > 0 else 0

np.random.seed(42)
N_PERMS = 10000

results = {}
for sig_id in hallmarks:
    sig_eff = effects[effects.signature_id == sig_id].sort_values("cohort_id")
    cohorts = sig_eff.cohort_id.tolist()
    ds = sig_eff.cohens_d.values
    vs = sig_eff.cohens_d_var.values
    programs = [cohort_program[c] for c in cohorts]

    # Observed Q_B/Q_total
    observed = compute_qb_fraction(ds, vs, programs)

    # Permutation null: shuffle program labels preserving group sizes
    null_dist = np.empty(N_PERMS)
    for i in range(N_PERMS):
        perm_programs = list(np.random.permutation(programs))
        null_dist[i] = compute_qb_fraction(ds, vs, perm_programs)

    p_value = np.mean(null_dist >= observed)

    results[sig_id] = {
        "observed_qb_fraction": round(float(observed), 4),
        "null_mean": round(float(np.mean(null_dist)), 4),
        "null_median": round(float(np.median(null_dist)), 4),
        "null_95th": round(float(np.percentile(null_dist, 95)), 4),
        "null_99th": round(float(np.percentile(null_dist, 99)), 4),
        "p_value": round(float(p_value), 4),
        "null_distribution": null_dist.tolist(),  # for histograms
    }

    short = sig_id.replace("hallmark_", "")
    print(f"{short:<25} observed={observed:.3f}  null_95th={np.percentile(null_dist, 95):.3f}  null_99th={np.percentile(null_dist, 99):.3f}  p={p_value:.4f}")

# Summary
sig_below_001 = sum(1 for r in results.values() if r["p_value"] < 0.01)
sig_below_005 = sum(1 for r in results.values() if r["p_value"] < 0.05)
print(f"\nSignificant at p<0.01: {sig_below_001}/7")
print(f"Significant at p<0.05: {sig_below_005}/7")

with open("outputs/canonical_v7/permutation_validation.json", "w") as f:
    json.dump({"n_permutations": N_PERMS, "seed": 42, "results": results}, f, indent=2)
print("Saved to outputs/canonical_v7/permutation_validation.json")
