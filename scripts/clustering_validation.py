import sys
sys.path.insert(0, 'src')
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import adjusted_rand_score
import json

effects = pd.read_csv("outputs/canonical_v7/per_cohort_effects.csv")
manifest = pd.read_csv("config/cohort_manifest.tsv", sep="\t")
cohort_program = dict(zip(manifest.cohort_id, manifest.biological_program))

hallmarks = [
    "hallmark_ifng_response", "hallmark_ifna_response",
    "hallmark_inflammatory_response", "hallmark_hypoxia",
    "hallmark_e2f_targets", "hallmark_emt", "hallmark_tnfa_nfkb"
]

# Build cohort x hallmark effect matrix
cohort_ids = sorted(manifest.cohort_id.tolist())
matrix = np.zeros((len(cohort_ids), len(hallmarks)))

for j, sig_id in enumerate(hallmarks):
    sig_eff = effects[effects.signature_id == sig_id]
    for i, cid in enumerate(cohort_ids):
        row = sig_eff[sig_eff.cohort_id == cid]
        if len(row) > 0:
            matrix[i, j] = row.cohens_d.values[0]

true_labels = [cohort_program[c] for c in cohort_ids]
unique_programs = sorted(set(true_labels))
true_numeric = [unique_programs.index(l) for l in true_labels]

# K-means (K=5)
np.random.seed(42)
kmeans = KMeans(n_clusters=5, random_state=42, n_init=10)
km_labels = kmeans.fit_predict(matrix)
ari_kmeans = adjusted_rand_score(true_numeric, km_labels)

# Ward's hierarchical
Z = linkage(matrix, method='ward')
ward_labels = fcluster(Z, t=5, criterion='maxclust') - 1
ari_ward = adjusted_rand_score(true_numeric, ward_labels)

print(f"K-means ARI: {ari_kmeans:.3f}")
print(f"Ward's ARI: {ari_ward:.3f}")

# Show assignments
print(f"\n{'Cohort':<40} {'True':<15} {'K-means':<10} {'Ward':<10}")
print("-" * 75)
for i, cid in enumerate(cohort_ids):
    true_prog = true_labels[i]
    km_cluster = km_labels[i]
    ward_cluster = ward_labels[i]
    # Map cluster numbers to most common true label in that cluster
    match_km = "?"
    match_ward = "?"
    print(f"{cid:<40} {true_prog:<15} {km_cluster:<10} {ward_cluster:<10}")

# Check IPF lung clustering
print("\n=== IPF LUNG COHORTS ===")
for cid in ["ipf_lung_gse47460", "ipf_lung_gse53845"]:
    if cid in cohort_ids:
        idx = cohort_ids.index(cid)
        print(f"{cid}: true={true_labels[idx]}, km_cluster={km_labels[idx]}, ward_cluster={ward_labels[idx]}")
        # What other cohorts share this cluster?
        km_same = [cohort_ids[j] for j in range(len(cohort_ids)) if km_labels[j] == km_labels[idx] and j != idx]
        print(f"  K-means cluster mates: {km_same[:5]}")

results = {
    "ari_kmeans": round(float(ari_kmeans), 4),
    "ari_ward": round(float(ari_ward), 4),
    "cohort_assignments": [
        {"cohort": cohort_ids[i], "true_program": true_labels[i],
         "kmeans_cluster": int(km_labels[i]), "ward_cluster": int(ward_labels[i])}
        for i in range(len(cohort_ids))
    ],
}

with open("outputs/canonical_v7/clustering_validation.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to outputs/canonical_v7/clustering_validation.json")
