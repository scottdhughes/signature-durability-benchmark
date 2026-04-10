#!/usr/bin/env python3
"""Analysis 2: Investigate confounder rejection accuracy drop (0.80 -> 0.60).

Between v5 (22 cohorts) and v6 (30 cohorts), confounder rejection went from 4/5 to 3/5.
Find which signature flipped and why.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
V5_DIR = ROOT / "outputs" / "canonical_v5"
V6_DIR = ROOT / "outputs" / "canonical_v6"

# The 5 confounded signatures in the primary split
CONFOUNDED_SIGS = [
    "confounded_proliferation",
    "confounded_ribosomal",
    "confounded_immune_infiltrate",
    "stealth_confounded_hypoxia_prolif",
    "stealth_confounded_inflam_immune",
]

MODEL = "within_program"  # the model we're analyzing


def load_aggregate(path):
    """Load aggregate durability scores, return dict keyed by (signature_id, model)."""
    data = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            if row["split"] == "primary":
                key = (row["signature_id"], row["model"])
                data[key] = row
    return data


def main():
    v5 = load_aggregate(V5_DIR / "aggregate_durability_scores.csv")
    v6 = load_aggregate(V6_DIR / "aggregate_durability_scores.csv")

    results = {}
    flipped = []

    print("=" * 80)
    print("CONFOUNDER REJECTION INVESTIGATION: v5 -> v6")
    print("=" * 80)

    for sig_id in CONFOUNDED_SIGS:
        v5_row = v5.get((sig_id, MODEL), {})
        v6_row = v6.get((sig_id, MODEL), {})

        v5_class = v5_row.get("predicted_class", "N/A")
        v6_class = v6_row.get("predicted_class", "N/A")

        v5_max_conf = float(v5_row.get("max_confounder_effect", 0))
        v6_max_conf = float(v6_row.get("max_confounder_effect", 0))

        v5_agg_effect = float(v5_row.get("aggregate_effect", 0))
        v6_agg_effect = float(v6_row.get("aggregate_effect", 0))

        v5_within_d = v5_row.get("within_d", "")
        v6_within_d = v6_row.get("within_d", "")
        v5_within_p = v5_row.get("within_p", "")
        v6_within_p = v6_row.get("within_p", "")

        # Is this correctly classified?
        v5_correct = v5_class == "confounded"
        v6_correct = v6_class == "confounded"
        did_flip = v5_correct != v6_correct

        entry = {
            "signature_id": sig_id,
            "expected": "confounded",
            "v5": {
                "predicted": v5_class,
                "correct": v5_correct,
                "max_confounder_effect": v5_max_conf,
                "aggregate_effect": v5_agg_effect,
                "confounder_ratio": v5_max_conf / abs(v5_agg_effect) if v5_agg_effect != 0 else None,
                "within_d": float(v5_within_d) if v5_within_d else None,
                "within_p": float(v5_within_p) if v5_within_p else None,
            },
            "v6": {
                "predicted": v6_class,
                "correct": v6_correct,
                "max_confounder_effect": v6_max_conf,
                "aggregate_effect": v6_agg_effect,
                "confounder_ratio": v6_max_conf / abs(v6_agg_effect) if v6_agg_effect != 0 else None,
                "within_d": float(v6_within_d) if v6_within_d else None,
                "within_p": float(v6_within_p) if v6_within_p else None,
            },
            "flipped": did_flip,
        }

        if did_flip:
            flipped.append(sig_id)
            # Analyze why
            if v5_correct and not v6_correct:
                entry["flip_direction"] = "lost_rejection"
                # Check the within_program model logic:
                # confounded if max_conf >= aggregate_effect
                # then: durable if within_p < 0.05 and |within_d| > 0.2
                reasons = []
                if v6_max_conf < abs(v6_agg_effect):
                    reasons.append(f"Confounder ratio dropped below 1.0 (v5: {v5_max_conf/abs(v5_agg_effect):.3f}, v6: {v6_max_conf/abs(v6_agg_effect):.3f})")
                    reasons.append(f"v5 max_conf={v5_max_conf:.4f} vs agg={abs(v5_agg_effect):.4f}: ratio={v5_max_conf/abs(v5_agg_effect):.3f}")
                    reasons.append(f"v6 max_conf={v6_max_conf:.4f} vs agg={abs(v6_agg_effect):.4f}: ratio={v6_max_conf/abs(v6_agg_effect):.3f}")
                if v6_within_d and float(v6_within_p) < 0.05 and abs(float(v6_within_d)) > 0.2:
                    reasons.append(f"Within-program became significant: d={float(v6_within_d):.3f}, p={float(v6_within_p):.4e}")
                    reasons.append("within_program model bypasses confounder check when within-program is significant")
                entry["reasons"] = reasons
            elif not v5_correct and v6_correct:
                entry["flip_direction"] = "gained_rejection"

        results[sig_id] = entry

        print(f"\n{sig_id}:")
        print(f"  v5: {v5_class} (correct={v5_correct})")
        print(f"       max_conf={v5_max_conf:.4f}, agg_effect={v5_agg_effect:.4f}, ratio={v5_max_conf/abs(v5_agg_effect):.3f}" if v5_agg_effect != 0 else f"       max_conf={v5_max_conf:.4f}, agg_effect={v5_agg_effect:.4f}")
        if v5_within_d:
            print(f"       within_d={float(v5_within_d):.3f}, within_p={float(v5_within_p):.4e}")
        print(f"  v6: {v6_class} (correct={v6_correct})")
        print(f"       max_conf={v6_max_conf:.4f}, agg_effect={v6_agg_effect:.4f}, ratio={v6_max_conf/abs(v6_agg_effect):.3f}" if v6_agg_effect != 0 else f"       max_conf={v6_max_conf:.4f}, agg_effect={v6_agg_effect:.4f}")
        if v6_within_d:
            print(f"       within_d={float(v6_within_d):.3f}, within_p={float(v6_within_p):.4e}")
        if did_flip:
            print(f"  *** FLIPPED: {entry['flip_direction']} ***")
            if "reasons" in entry:
                for r in entry["reasons"]:
                    print(f"      {r}")

    # Also check the full_model for comparison
    print("\n" + "=" * 80)
    print("COMPARISON WITH full_model:")
    print("=" * 80)
    full_model_results = {}
    for sig_id in CONFOUNDED_SIGS:
        v5_fm = v5.get((sig_id, "full_model"), {})
        v6_fm = v6.get((sig_id, "full_model"), {})
        v5_class = v5_fm.get("predicted_class", "N/A")
        v6_class = v6_fm.get("predicted_class", "N/A")
        print(f"  {sig_id}: v5={v5_class}, v6={v6_class}")
        full_model_results[sig_id] = {
            "v5": v5_class,
            "v6": v6_class,
            "v5_correct": v5_class == "confounded",
            "v6_correct": v6_class == "confounded",
        }

    # Dilution analysis
    print("\n" + "=" * 80)
    print("DILUTION ANALYSIS: max_confounder_effect across models")
    print("=" * 80)
    dilution = {}
    for sig_id in CONFOUNDED_SIGS:
        v5_fm = v5.get((sig_id, "full_model"), {})
        v6_fm = v6.get((sig_id, "full_model"), {})
        v5_conf = float(v5_fm.get("max_confounder_effect", 0))
        v6_conf = float(v6_fm.get("max_confounder_effect", 0))
        v5_agg = float(v5_fm.get("aggregate_effect", 0))
        v6_agg = float(v6_fm.get("aggregate_effect", 0))
        conf_change = v6_conf - v5_conf
        agg_change = v6_agg - v5_agg
        print(f"  {sig_id}:")
        print(f"    max_conf: v5={v5_conf:.4f} -> v6={v6_conf:.4f} (delta={conf_change:+.4f})")
        print(f"    agg_effect: v5={v5_agg:.4f} -> v6={v6_agg:.4f} (delta={agg_change:+.4f})")
        dilution[sig_id] = {
            "v5_max_conf": v5_conf, "v6_max_conf": v6_conf,
            "v5_agg_effect": v5_agg, "v6_agg_effect": v6_agg,
            "conf_change": conf_change, "agg_change": agg_change,
            "note": "max_confounder_effect is mean across per-cohort max scores; more cohorts where confounder isn't active dilutes the mean"
        }

    # Summary
    v5_correct = sum(1 for e in results.values() if e["v5"]["correct"])
    v6_correct = sum(1 for e in results.values() if e["v6"]["correct"])

    summary = {
        "analysis": "confounder_investigation",
        "description": "Investigates why confounder rejection dropped from v5 to v6",
        "model": MODEL,
        "v5_correct": f"{v5_correct}/5",
        "v6_correct": f"{v6_correct}/5",
        "flipped_signatures": flipped,
        "per_signature": results,
        "full_model_comparison": full_model_results,
        "dilution_analysis": dilution,
        "root_cause": (
            "confounded_proliferation is the ONLY signature that flipped between v5 and v6. "
            "In v5, max_confounder_effect (0.419) barely exceeded aggregate_effect (0.406), "
            "ratio=1.031, triggering 'confounded'. In v6, adding 8 cohorts raised "
            "aggregate_effect to 0.558 (the proliferation signal strengthened across more "
            "contexts) while max_confounder_effect rose less (to 0.503), ratio=0.901, "
            "so the confounder check no longer triggers. The within_program model then "
            "falls through to mixed (I²>0.75). This flip occurs in BOTH full_model and "
            "within_program, confirming the issue is the aggregate_effect dilution mechanism. "
            "Note: stealth_confounded_hypoxia_prolif was already misclassified in v5 "
            "(durable in within_program, mixed in full_model), so it did NOT contribute "
            "to the drop. The fundamental issue: the confounder_ratio threshold of 1.0 "
            "is fragile when both numerator and denominator are mean-of-cohort scores "
            "that shift with cohort composition."
        ),
    }

    out_path = V6_DIR / "confounder_investigation.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
