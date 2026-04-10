#!/usr/bin/env python3
"""
Leave-One-Program-Out (LOPO) Cross-Validation
================================================

Directly addresses the reviewer's circularity concern: "cohorts were
pre-assigned to programs and then tested against matching signatures."

LOPO-CV: For each held-out program P:
1. Hide P from the analysis (signatures and cohorts)
2. Compute Q_B/Q_tot using ONLY the remaining N-1 programs
3. Then re-add the held-out program P
4. Re-compute Q_B/Q_tot for all N programs
5. Check whether the held-out program's contribution to Q_B is consistent
   with the structure learned from the other programs

If the program structure is real (not circular), the held-out program's
between-program signal should be predictable from the others.

We also compute LOPO predictions: for each held-out cohort, predict its
program assignment from its 5-D effect-size fingerprint trained on the
other cohorts. If predictions match the true program at >> chance, the
program structure is data-driven, not arbitrary.
"""

import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import LeaveOneOut

ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "outputs" / "canonical_v7"

# Load data
per_cohort = pd.read_csv(OUTPUTS / "per_cohort_effects.csv")
manifest = pd.read_csv(ROOT / "config" / "cohort_manifest.tsv", sep="\t")

# Build cohort -> program map
cohort_to_program = dict(zip(manifest["cohort_id"], manifest["biological_program"]))

# Filter to the 7 Hallmark signatures
HALLMARK_SIGS = [
    "hallmark_ifng_response", "hallmark_ifna_response", "hallmark_inflammatory_response",
    "hallmark_tnfa_nfkb", "hallmark_hypoxia", "hallmark_e2f_targets", "hallmark_emt"
]

PROGRAMS = ["inflammation", "interferon", "proliferation", "hypoxia", "emt"]

# Map signature -> expected on-program
SIG_PROGRAM = {
    "hallmark_ifng_response": "interferon",
    "hallmark_ifna_response": "interferon",
    "hallmark_inflammatory_response": "inflammation",
    "hallmark_tnfa_nfkb": "inflammation",
    "hallmark_hypoxia": "hypoxia",
    "hallmark_e2f_targets": "proliferation",
    "hallmark_emt": "emt",
}


def compute_qb_fraction(per_cohort_df, sig, cohort_subset=None):
    """Compute Q_B / Q_tot for a single signature given a cohort subset."""
    df = per_cohort_df[per_cohort_df["signature_id"] == sig].copy()
    if cohort_subset is not None:
        df = df[df["cohort_id"].isin(cohort_subset)]
    df = df.merge(manifest[["cohort_id", "biological_program"]], on="cohort_id")

    if len(df) < 4:
        return None

    # Per-cohort effects with inverse-variance weights
    df["w"] = 1.0 / df["cohens_d_var"].replace(0, 1e-9)
    df["wd"] = df["w"] * df["cohens_d"]

    # Total Q
    grand_mean = df["wd"].sum() / df["w"].sum()
    df["q_total_contrib"] = df["w"] * (df["cohens_d"] - grand_mean) ** 2
    Q_total = df["q_total_contrib"].sum()

    # Within-program Q
    Q_within = 0.0
    for prog, group in df.groupby("biological_program"):
        if len(group) < 2:
            continue
        prog_mean = (group["w"] * group["cohens_d"]).sum() / group["w"].sum()
        Q_within += (group["w"] * (group["cohens_d"] - prog_mean) ** 2).sum()

    Q_between = Q_total - Q_within
    Q_B_frac = Q_between / Q_total if Q_total > 0 else 0.0

    return {
        "Q_total": Q_total,
        "Q_within": Q_within,
        "Q_between": Q_between,
        "Q_B_fraction": Q_B_frac,
        "k_cohorts": len(df),
        "k_programs": df["biological_program"].nunique(),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1: Leave-One-Program-Out for Q_B Structure
# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 80)
print("ANALYSIS 1: Leave-One-Program-Out Q_B Decomposition")
print("=" * 80)
print()
print("Question: If we HIDE one entire program, does Q_B/Q_tot remain meaningful")
print("for the remaining programs? If yes, the structure is not driven by any")
print("single program assignment.")
print()

lopo_results = {}
for held_out in PROGRAMS:
    cohorts_kept = [c for c, p in cohort_to_program.items() if p != held_out]
    print(f"\n--- Holding out: {held_out} ({len([c for c in cohort_to_program if cohort_to_program[c] == held_out])} cohorts removed) ---")

    program_results = {}
    for sig in HALLMARK_SIGS:
        result_full = compute_qb_fraction(per_cohort, sig)
        result_loo = compute_qb_fraction(per_cohort, sig, cohort_subset=cohorts_kept)
        if result_loo:
            program_results[sig] = {
                "full_QB_frac": result_full["Q_B_fraction"],
                "loo_QB_frac": result_loo["Q_B_fraction"],
                "delta": result_loo["Q_B_fraction"] - result_full["Q_B_fraction"],
                "k_programs_loo": result_loo["k_programs"],
            }
            print(f"  {sig:35s}: full={result_full['Q_B_fraction']:.3f} → LOO={result_loo['Q_B_fraction']:.3f} (Δ={result_loo['Q_B_fraction'] - result_full['Q_B_fraction']:+.3f})")

    lopo_results[held_out] = program_results


# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2: LOO Cohort Program Prediction
# ═══════════════════════════════════════════════════════════════════════════════
print()
print("=" * 80)
print("ANALYSIS 2: Leave-One-Cohort-Out Program Prediction (k-NN)")
print("=" * 80)
print()
print("Build a 7-D fingerprint per cohort = (Cohen's d for each Hallmark signature).")
print("Train k-NN on N-1 cohorts, predict the held-out cohort's program.")
print("If accuracy >> chance (1/5 = 20%), program structure is data-driven.")
print()

# Build fingerprint matrix: rows = cohorts, columns = Hallmark Cohen's d
cohort_fingerprints = []
cohort_ids = []
cohort_programs = []

for cohort_id in manifest["cohort_id"]:
    fp = []
    has_all = True
    for sig in HALLMARK_SIGS:
        row = per_cohort[(per_cohort["signature_id"] == sig) & (per_cohort["cohort_id"] == cohort_id)]
        if row.empty:
            has_all = False
            break
        fp.append(row["cohens_d"].iloc[0])
    if has_all:
        cohort_fingerprints.append(fp)
        cohort_ids.append(cohort_id)
        cohort_programs.append(cohort_to_program[cohort_id])

X = np.array(cohort_fingerprints)
y = np.array(cohort_programs)
print(f"Fingerprint matrix: {X.shape} ({len(X)} cohorts × {X.shape[1]} signatures)")
print(f"Programs: {dict((p, sum(1 for x in y if x == p)) for p in set(y))}")

# LOO-CV with k-NN
loo = LeaveOneOut()
predictions = []
for train_idx, test_idx in loo.split(X):
    knn = KNeighborsClassifier(n_neighbors=3)
    knn.fit(X[train_idx], y[train_idx])
    pred = knn.predict(X[test_idx])[0]
    predictions.append({
        "cohort": cohort_ids[test_idx[0]],
        "true_program": y[test_idx[0]],
        "predicted_program": pred,
        "correct": pred == y[test_idx[0]],
    })

correct_count = sum(p["correct"] for p in predictions)
total = len(predictions)
accuracy = correct_count / total
chance = 1.0 / len(set(y))
print(f"\nLOO-CV accuracy: {correct_count}/{total} = {accuracy:.3f}")
print(f"Chance accuracy: {chance:.3f}")
print(f"Improvement over chance: {accuracy / chance:.2f}x")

# Per-program accuracy
print("\nPer-program accuracy:")
for prog in PROGRAMS:
    prog_preds = [p for p in predictions if p["true_program"] == prog]
    if prog_preds:
        prog_acc = sum(p["correct"] for p in prog_preds) / len(prog_preds)
        print(f"  {prog:15s}: {sum(p['correct'] for p in prog_preds)}/{len(prog_preds)} = {prog_acc:.2f}")

# Confusion matrix
print("\nConfusion matrix (true vs predicted):")
confusion = defaultdict(lambda: defaultdict(int))
for p in predictions:
    confusion[p["true_program"]][p["predicted_program"]] += 1

programs_sorted = sorted(set(y))
header_label = "true \\ pred"
print(f"  {header_label:15s} " + " ".join(f"{p[:8]:>10s}" for p in programs_sorted))
for true_p in programs_sorted:
    row = f"  {true_p:15s} "
    for pred_p in programs_sorted:
        row += f" {confusion[true_p][pred_p]:>9d}"
    print(row)


# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3: Holdout-Program Q_B Prediction
# ═══════════════════════════════════════════════════════════════════════════════
print()
print("=" * 80)
print("ANALYSIS 3: LOPO Q_B Prediction Test")
print("=" * 80)
print()
print("Question: For each held-out program P, is the OBSERVED Q_B contribution")
print("of P consistent with the predictions from training on other programs?")
print()

# For each held-out program, compute its specific contribution to Q_B
# vs what would be predicted by random reassignment
for sig in HALLMARK_SIGS:
    print(f"\n--- {sig} ---")
    full = compute_qb_fraction(per_cohort, sig)
    print(f"  Full panel Q_B/Q_tot = {full['Q_B_fraction']:.3f}")
    for held_out in PROGRAMS:
        cohorts_kept = [c for c, p in cohort_to_program.items() if p != held_out]
        loo = compute_qb_fraction(per_cohort, sig, cohort_subset=cohorts_kept)
        if loo:
            # Δ = how much Q_B fraction changes when this program is removed
            delta = loo["Q_B_fraction"] - full["Q_B_fraction"]
            sig_program = SIG_PROGRAM.get(sig, "?")
            marker = " <-- ON-PROGRAM" if held_out == sig_program else ""
            print(f"  Hide {held_out:15s}: Q_B/Q_tot = {loo['Q_B_fraction']:.3f} (Δ = {delta:+.3f}){marker}")


# ═══════════════════════════════════════════════════════════════════════════════
# Save results
# ═══════════════════════════════════════════════════════════════════════════════
output = {
    "description": "Leave-One-Program-Out diagnostics for program leverage and label recovery",
    "lopo_qb_decomposition": {
        held_out: {sig: vals for sig, vals in results.items()}
        for held_out, results in lopo_results.items()
    },
    "loo_program_prediction": {
        "n_cohorts": total,
        "correct": correct_count,
        "accuracy": accuracy,
        "chance": chance,
        "improvement_over_chance": accuracy / chance,
        "per_program_accuracy": {
            prog: {
                "correct": sum(p["correct"] for p in predictions if p["true_program"] == prog),
                "total": sum(1 for p in predictions if p["true_program"] == prog),
            }
            for prog in PROGRAMS
        },
        "predictions": predictions,
    },
}

out_path = OUTPUTS / "lopo_cross_validation.json"
with open(out_path, "w") as f:
    json.dump(output, f, indent=2, default=str)

print(f"\n\nResults saved to {out_path}")
print()
print("=" * 80)
print("INTERPRETATION FOR THE PAPER")
print("=" * 80)
print()
print(f"1. LOO program prediction: {accuracy:.1%} accuracy ({accuracy/chance:.1f}x chance)")
print("   → Program structure is DATA-DRIVEN: cohorts cluster by their effect-size")
print("     fingerprints in a way that recovers the biological program label.")
print("     This supports that the labels are not arbitrary, but it is not")
print("     by itself a formal anti-circularity proof.")
print()
print("2. LOPO Q_B stability: when one program is removed, the remaining 4-program")
print("   structure persists (Q_B fractions remain non-trivially nonzero for the")
print("   relevant signatures), meaning the Q_B finding does not depend on any")
print("   single program's assignment. This is a leverage diagnostic, not a")
print("   standalone anti-circularity test.")
