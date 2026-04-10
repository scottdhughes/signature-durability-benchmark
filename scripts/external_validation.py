#!/usr/bin/env python3
"""Held-out validation of within-program findings on the v8 expanded panel.

For each target signature, this script evaluates:
  1. Leave-one-cohort-out prediction within the signature's home program.
  2. Split-half replication within the same program.

The key question is whether a program-conditioned finding is robust to withholding
entire cohorts, not just significant on the full pooled analysis.
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from signature_durability_benchmark.meta_analysis import guarded_hksj_random_effects_meta


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs" / "canonical_v8"
PER_COHORT = OUT / "per_cohort_effects.csv"
MANIFEST = ROOT / "config" / "cohort_manifest.tsv"
OUT_JSON = OUT / "external_validation.json"

TARGETS = {
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


def loo_within_program(sig_df: pd.DataFrame) -> dict:
    rows = []
    for i in range(len(sig_df)):
        held_out = sig_df.iloc[i]
        remaining = sig_df.drop(sig_df.index[i])
        meta = guarded_hksj_random_effects_meta(remaining["cohens_d"].tolist(), remaining["cohens_d_var"].tolist())
        held_effect = float(held_out["cohens_d"])
        pooled_effect = float(meta["pooled_effect"])
        sign_match = (held_effect > 0) == (pooled_effect > 0)
        rows.append(
            {
                "held_out_cohort": held_out["cohort_id"],
                "held_out_effect": round(held_effect, 4),
                "pooled_effect_remaining": round(pooled_effect, 4),
                "pooled_p_remaining": float(meta["pooled_p"]),
                "i_squared_remaining": round(float(meta["i_squared"]), 4),
                "sign_match": bool(sign_match),
            }
        )
    matches = sum(row["sign_match"] for row in rows)
    return {
        "n_holdouts": len(rows),
        "sign_matches": matches,
        "prediction_accuracy": round(matches / len(rows), 4) if rows else 0.0,
        "details": rows,
    }


def split_half_within_program(sig_df: pd.DataFrame) -> dict:
    k = len(sig_df)
    half = k // 2
    rows = []
    for combo in itertools.combinations(range(k), half):
        other = tuple(i for i in range(k) if i not in combo)
        if k % 2 == 0 and combo[0] > other[0]:
            continue
        if k % 2 == 1 and combo[0] > other[0]:
            continue

        group_a = sig_df.iloc[list(combo)]
        group_b = sig_df.iloc[list(other)]
        if len(group_a) < 2 or len(group_b) < 2:
            continue

        meta_a = guarded_hksj_random_effects_meta(group_a["cohens_d"].tolist(), group_a["cohens_d_var"].tolist())
        meta_b = guarded_hksj_random_effects_meta(group_b["cohens_d"].tolist(), group_b["cohens_d_var"].tolist())
        sign_agree = (meta_a["pooled_effect"] > 0) == (meta_b["pooled_effect"] > 0)
        both_significant = meta_a["pooled_p"] < 0.05 and meta_b["pooled_p"] < 0.05

        rows.append(
            {
                "group_a": group_a["cohort_id"].tolist(),
                "group_b": group_b["cohort_id"].tolist(),
                "pooled_effect_a": round(float(meta_a["pooled_effect"]), 4),
                "pooled_p_a": float(meta_a["pooled_p"]),
                "pooled_effect_b": round(float(meta_b["pooled_effect"]), 4),
                "pooled_p_b": float(meta_b["pooled_p"]),
                "sign_agree": bool(sign_agree),
                "both_significant": bool(both_significant),
            }
        )

    agree = sum(row["sign_agree"] for row in rows)
    both = sum(row["both_significant"] for row in rows)
    return {
        "n_splits": len(rows),
        "sign_agreements": agree,
        "sign_agreement_rate": round(agree / len(rows), 4) if rows else 0.0,
        "both_significant": both,
        "both_significant_rate": round(both / len(rows), 4) if rows else 0.0,
        "details": rows,
    }


def main() -> None:
    effects = pd.read_csv(PER_COHORT)
    manifest = pd.read_csv(MANIFEST, sep="\t")
    cohort_program = dict(zip(manifest["cohort_id"], manifest["biological_program"]))

    results = {
        "description": (
            "Held-out validation of within-program findings on the expanded v8 panel "
            "using leave-one-cohort-out prediction and split-half replication."
        ),
        "signatures": {},
    }

    print("=" * 80)
    print("HELD-OUT VALIDATION ON THE EXPANDED V8 PANEL")
    print("=" * 80)

    for sig, program in TARGETS.items():
        sig_df = effects[(effects["signature_id"] == sig) & (effects["cohort_id"].map(cohort_program) == program)].copy()
        if len(sig_df) < 3:
            continue

        loo = loo_within_program(sig_df)
        split = split_half_within_program(sig_df) if len(sig_df) >= 4 else None
        results["signatures"][sig] = {
            "program": program,
            "within_k": int(len(sig_df)),
            "leave_one_out": loo,
            "split_half": split,
        }

        print(
            f"{sig:35s} program={program:13s} "
            f"LOO={loo['prediction_accuracy']:.1%} "
            f"split_agree={split['sign_agreement_rate']:.1%} "
            f"split_both_sig={split['both_significant_rate']:.1%}"
        )

    ifn_sigs = [
        "hallmark_ifng_response",
        "hallmark_ifna_response",
        "schoggins_2011_irg",
        "blind_durable_ifn_composite",
    ]
    results["ifn_focus_summary"] = {
        "all_four_loo_prediction_accuracy": {
            sig: results["signatures"][sig]["leave_one_out"]["prediction_accuracy"] for sig in ifn_sigs
        },
        "all_four_split_half_sign_agreement": {
            sig: results["signatures"][sig]["split_half"]["sign_agreement_rate"] for sig in ifn_sigs
        },
        "all_four_split_half_both_significant": {
            sig: results["signatures"][sig]["split_half"]["both_significant_rate"] for sig in ifn_sigs
        },
    }

    with open(OUT_JSON, "w") as handle:
        json.dump(results, handle, indent=2)
    print(f"\nSaved to {OUT_JSON}")


if __name__ == "__main__":
    main()
