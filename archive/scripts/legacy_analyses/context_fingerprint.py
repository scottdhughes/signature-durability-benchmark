"""Context fingerprinting: predict biological program from per-cohort effect profiles.

For each signature, compute a 5-element "context fingerprint" = mean Cohen's d
within each biological program's cohorts. If Hallmark signatures cluster by their
expected program, this validates that cross-context heterogeneity carries
*information* about biological context rather than being noise.

Practical implication: a researcher with an unknown signature could test it across
a multi-program panel and read off which biology it captures from the fingerprint
shape -- a diagnostic tool for signature interpretation.
"""
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list
import json

effects = pd.read_csv("outputs/canonical_v6/per_cohort_effects.csv")
manifest = pd.read_csv("config/cohort_manifest.tsv", sep="\t")
cohort_program = dict(zip(manifest.cohort_id, manifest.biological_program))

programs = sorted(set(cohort_program.values()))

# All signatures including confounded/stealth for comparison
all_sigs = sorted(effects.signature_id.unique())

# Focus on interpretable subset
focus_sigs = [
    "hallmark_ifng_response",
    "hallmark_ifna_response",
    "hallmark_inflammatory_response",
    "hallmark_hypoxia",
    "hallmark_e2f_targets",
    "hallmark_emt",
    "hallmark_tnfa_nfkb",
    "stealth_confounded_inflam_immune",
    "stealth_confounded_hypoxia_prolif",
    "confounded_proliferation",
    "confounded_immune_infiltrate",
    "mixed_emt_inflammation",
    "mixed_hypoxia_stress",
    "brittle_small_study_inflammation",
    "brittle_overfit_noise",
]

# Expected program for hallmarks
expected = {
    "hallmark_ifng_response": "interferon",
    "hallmark_ifna_response": "interferon",
    "hallmark_inflammatory_response": "inflammation",
    "hallmark_hypoxia": "hypoxia",
    "hallmark_e2f_targets": "proliferation",
    "hallmark_emt": "emt",
    "hallmark_tnfa_nfkb": "inflammation",
}

# Compute fingerprints
fingerprints = {}
for sig_id in focus_sigs:
    sig_eff = effects[effects.signature_id == sig_id]
    fp = {}
    for prog in programs:
        prog_cohorts = [c for c, p in cohort_program.items() if p == prog]
        prog_effects = sig_eff[sig_eff.cohort_id.isin(prog_cohorts)].cohens_d.values
        fp[prog] = round(float(np.mean(prog_effects)), 3) if len(prog_effects) > 0 else 0
    fingerprints[sig_id] = fp

# Print fingerprint table
print(f"{'Signature':<40} " + " | ".join(f"{p[:5]:>7}" for p in programs))
print("-" * (40 + 11 * len(programs)))

for sig_id in focus_sigs:
    fp = fingerprints[sig_id]
    short = sig_id.replace("hallmark_", "H:").replace("stealth_confounded_", "S:").replace("confounded_", "C:").replace("mixed_", "M:").replace("brittle_", "B:")
    vals = " | ".join(f"{fp[p]:+7.3f}" for p in programs)
    print(f"{short:<40} {vals}")

# Prediction accuracy for hallmarks
print("\n" + "=" * 80)
print("Program prediction from fingerprint (argmax of mean |d| per program):")
print("-" * 80)

correct = 0
total = 0
for sig_id in expected:
    fp = fingerprints[sig_id]
    # For prediction, use absolute value of mean effect (some programs have
    # negative effects that are still strong signals)
    predicted = max(fp, key=lambda p: abs(fp[p]))
    actual = expected[sig_id]
    match = predicted == actual
    correct += match
    total += 1
    mark = "CORRECT" if match else "WRONG"
    short = sig_id.replace("hallmark_", "")
    print(f"  {mark:<8} {short:<25} predicted={predicted:<15} actual={actual:<15} (d={fp[predicted]:+.3f} vs d={fp[actual]:+.3f})")

print(f"\nAccuracy: {correct}/{total} = {correct/total:.0%}")

# Compute pairwise distances between hallmark fingerprints
print("\n" + "=" * 80)
print("Euclidean distance matrix between Hallmark fingerprints:")
print("(Signatures from same program should cluster together)")
print("-" * 80)

hallmark_ids = list(expected.keys())
fp_matrix = np.array([[fingerprints[s][p] for p in programs] for s in hallmark_ids])
dist_matrix = squareform(pdist(fp_matrix, metric="euclidean"))

# Print distance matrix
short_names = [s.replace("hallmark_", "") for s in hallmark_ids]
print(f"{'':>15}", end="")
for sn in short_names:
    print(f"  {sn[:10]:>10}", end="")
print()
for i, sn1 in enumerate(short_names):
    print(f"{sn1:>15}", end="")
    for j, sn2 in enumerate(short_names):
        print(f"  {dist_matrix[i,j]:10.3f}", end="")
    print()

# Check: do same-program pairs have smaller distances than different-program pairs?
same_program_dists = []
diff_program_dists = []
for i in range(len(hallmark_ids)):
    for j in range(i + 1, len(hallmark_ids)):
        d = dist_matrix[i, j]
        if expected[hallmark_ids[i]] == expected[hallmark_ids[j]]:
            same_program_dists.append(d)
        else:
            diff_program_dists.append(d)

print(f"\nMean distance, same-program pairs:      {np.mean(same_program_dists):.3f} (n={len(same_program_dists)})")
print(f"Mean distance, different-program pairs:  {np.mean(diff_program_dists):.3f} (n={len(diff_program_dists)})")
print(f"Ratio (diff/same): {np.mean(diff_program_dists)/np.mean(same_program_dists):.2f}x")

# Confounded signature analysis
print("\n" + "=" * 80)
print("Confounded signature fingerprints reveal their confound:")
print("-" * 80)
confounded = {
    "stealth_confounded_hypoxia_prolif": "hypoxia + proliferation",
    "stealth_confounded_inflam_immune": "inflammation + immune",
    "confounded_proliferation": "pure proliferation",
    "confounded_immune_infiltrate": "immune infiltrate",
}
for sig_id, desc in confounded.items():
    fp = fingerprints[sig_id]
    top2 = sorted(fp, key=lambda p: abs(fp[p]), reverse=True)[:2]
    short = sig_id.replace("stealth_confounded_", "S:").replace("confounded_", "C:")
    print(f"  {short:<35} top programs: {top2[0]} ({fp[top2[0]]:+.3f}), {top2[1]} ({fp[top2[1]]:+.3f})  [{desc}]")

# Save
output = {
    "fingerprints": fingerprints,
    "programs": programs,
    "prediction_accuracy": f"{correct}/{total}",
    "same_program_mean_dist": round(float(np.mean(same_program_dists)), 4),
    "diff_program_mean_dist": round(float(np.mean(diff_program_dists)), 4),
}
with open("outputs/canonical_v6/context_fingerprints.json", "w") as f:
    json.dump(output, f, indent=2)
print("\nSaved to outputs/canonical_v6/context_fingerprints.json")
