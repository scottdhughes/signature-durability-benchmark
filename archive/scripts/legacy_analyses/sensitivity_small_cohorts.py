#!/usr/bin/env python3
"""Analysis 1: Sensitivity analysis dropping small cohorts (N<20).

Tests whether within-program findings are robust to removing cohorts with N<20.
"""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from signature_durability_benchmark.meta_analysis import dl_random_effects_meta

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
V6_DIR = ROOT / "outputs" / "canonical_v6"
MANIFEST = ROOT / "config" / "cohort_manifest.tsv"
SIGNATURE_PANEL = ROOT / "config" / "signature_panel.tsv"
PER_COHORT = V6_DIR / "per_cohort_effects.csv"
WITHIN_PROGRAM = V6_DIR / "within_program_durability.csv"
N_BONFERRONI = 9  # from within_program_durability.csv: bonferroni_n_tests

# Signature -> cohort program mapping (from benchmark.py)
SIG_TO_COHORT_PROGRAM = {
    "interferon_response": "interferon",
    "inflammatory_response": "inflammation",
    "hypoxia": "hypoxia",
    "proliferation": "proliferation",
    "emt": "emt",
}

HALLMARK_SIGS = [
    "hallmark_ifng_response",
    "hallmark_ifna_response",
    "hallmark_inflammatory_response",
    "hallmark_hypoxia",
    "hallmark_e2f_targets",
    "hallmark_emt",
    "hallmark_tnfa_nfkb",
]


def main():
    # Load cohort manifest -> sample counts
    cohort_n = {}
    cohort_program = {}
    with open(MANIFEST) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cohort_n[row["cohort_id"]] = int(row["sample_count"])
            cohort_program[row["cohort_id"]] = row["biological_program"]

    small_cohorts = {cid for cid, n in cohort_n.items() if n < 20}
    print(f"Small cohorts (N<20): {sorted(small_cohorts)}")
    for cid in sorted(small_cohorts):
        print(f"  {cid}: N={cohort_n[cid]}, program={cohort_program[cid]}")

    # Load signature -> program mapping
    sig_program_map = {}
    with open(SIGNATURE_PANEL) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sig_program_map[row["signature_id"]] = row["program"]

    # Load per-cohort effects
    per_cohort = []
    with open(PER_COHORT) as f:
        for row in csv.DictReader(f):
            per_cohort.append(row)

    # Load original within-program results for comparison
    original = {}
    with open(WITHIN_PROGRAM) as f:
        for row in csv.DictReader(f):
            original[row["signature_id"]] = row

    # For each Hallmark, recompute within-program DL meta EXCLUDING small cohorts
    results = {}
    for sig_id in HALLMARK_SIGS:
        sig_prog = sig_program_map[sig_id]
        cohort_prog = SIG_TO_COHORT_PROGRAM.get(sig_prog)
        if not cohort_prog:
            continue

        # Get within-program cohort effects (original)
        within_all_effects = []
        within_all_vars = []
        within_all_cohorts = []
        for row in per_cohort:
            if row["signature_id"] != sig_id:
                continue
            cid = row["cohort_id"]
            if cohort_program.get(cid) != cohort_prog:
                continue
            within_all_effects.append(float(row["cohens_d"]))
            within_all_vars.append(float(row["cohens_d_var"]))
            within_all_cohorts.append(cid)

        # Filtered: exclude small cohorts
        filt_effects = []
        filt_vars = []
        filt_cohorts = []
        dropped = []
        for i, cid in enumerate(within_all_cohorts):
            if cid in small_cohorts:
                dropped.append({"cohort": cid, "N": cohort_n[cid],
                                "d": within_all_effects[i]})
            else:
                filt_effects.append(within_all_effects[i])
                filt_vars.append(within_all_vars[i])
                filt_cohorts.append(cid)

        # Original meta
        orig_row = original.get(sig_id, {})
        orig_d = float(orig_row.get("within_d", 0)) if orig_row.get("within_d") else None
        orig_p = float(orig_row.get("within_p", 1)) if orig_row.get("within_p") else None
        orig_p_bonf = float(orig_row.get("within_p_bonferroni", 1)) if orig_row.get("within_p_bonferroni") else None

        # Recompute
        if len(filt_effects) >= 2:
            meta = dl_random_effects_meta(filt_effects, filt_vars)
            p_bonf = min(meta["pooled_p"] * N_BONFERRONI, 1.0)
        else:
            meta = None
            p_bonf = None

        entry = {
            "signature_id": sig_id,
            "program": sig_prog,
            "cohort_program": cohort_prog,
            "original": {
                "k": int(orig_row.get("within_k", 0)) if orig_row.get("within_k") else None,
                "within_d": orig_d,
                "within_p": orig_p,
                "within_p_bonferroni": orig_p_bonf,
                "bonferroni_significant": orig_p_bonf is not None and orig_p_bonf < 0.05,
            },
            "filtered": {
                "k": len(filt_effects),
                "within_d": meta["pooled_effect"] if meta else None,
                "within_p": meta["pooled_p"] if meta else None,
                "within_se": meta["pooled_se"] if meta else None,
                "i_squared": meta["i_squared"] if meta else None,
                "tau2": meta["tau2"] if meta else None,
                "within_p_bonferroni": p_bonf,
                "bonferroni_significant": p_bonf is not None and p_bonf < 0.05,
                "cohorts_used": filt_cohorts,
            },
            "dropped_cohorts": dropped,
            "n_dropped": len(dropped),
            "effect_size_change": None,
            "conclusion": None,
        }

        if meta and orig_d is not None:
            change = meta["pooled_effect"] - orig_d
            pct = (change / abs(orig_d) * 100) if orig_d != 0 else 0
            entry["effect_size_change"] = {
                "absolute": round(change, 4),
                "percent": round(pct, 1),
            }

            # Determine conclusion
            orig_sig = orig_p_bonf is not None and orig_p_bonf < 0.05
            filt_sig = p_bonf is not None and p_bonf < 0.05
            if orig_sig and filt_sig:
                entry["conclusion"] = "ROBUST: survives Bonferroni with and without small cohorts"
            elif orig_sig and not filt_sig:
                entry["conclusion"] = "SENSITIVE: loses Bonferroni significance when small cohorts removed"
            elif not orig_sig and filt_sig:
                entry["conclusion"] = "IMPROVED: gains Bonferroni significance when small cohorts removed"
            elif not orig_sig and not filt_sig:
                entry["conclusion"] = "UNCHANGED: fails Bonferroni both with and without small cohorts"

        results[sig_id] = entry

    # Summary
    summary = {
        "analysis": "sensitivity_small_cohorts",
        "description": "Tests robustness of within-program DL meta-analysis to removing cohorts with N<20",
        "small_cohorts_excluded": sorted(small_cohorts),
        "n_bonferroni_tests": N_BONFERRONI,
        "results": results,
    }

    # Count outcomes
    conclusions = {}
    for sig_id, entry in results.items():
        c = entry.get("conclusion", "N/A")
        conclusions[c] = conclusions.get(c, 0) + 1
    summary["outcome_counts"] = conclusions

    # Print summary
    print("\n" + "=" * 80)
    print("SENSITIVITY ANALYSIS: Dropping cohorts with N < 20")
    print("=" * 80)
    for sig_id in HALLMARK_SIGS:
        e = results[sig_id]
        print(f"\n{sig_id} ({e['cohort_program']}):")
        o = e["original"]
        f = e["filtered"]
        print(f"  Original: k={o['k']}, d={o['within_d']:.3f}, p={o['within_p']:.4e}, p_bonf={o['within_p_bonferroni']:.4f}" if o['within_d'] else f"  Original: no data")
        if f["within_d"] is not None:
            print(f"  Filtered: k={f['k']}, d={f['within_d']:.3f}, p={f['within_p']:.4e}, p_bonf={f['within_p_bonferroni']:.4f}")
        else:
            print(f"  Filtered: insufficient data (k={f['k']})")
        if e["n_dropped"] > 0:
            for dc in e["dropped_cohorts"]:
                print(f"  Dropped: {dc['cohort']} (N={dc['N']}, d={dc['d']:.3f})")
        if e["effect_size_change"]:
            ch = e["effect_size_change"]
            print(f"  Effect change: {ch['absolute']:+.4f} ({ch['percent']:+.1f}%)")
        print(f"  -> {e['conclusion']}")

    print("\n" + "=" * 80)
    print("OUTCOME SUMMARY:")
    for c, n in sorted(conclusions.items()):
        print(f"  {c}: {n}")

    # Save
    out_path = V6_DIR / "sensitivity_small_cohorts.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
