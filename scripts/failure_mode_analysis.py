#!/usr/bin/env python3
"""Program-boundary and cross-talk analysis on the expanded v8 panel.

This script formalizes the benchmark's second methodological claim:
the same diagnostics that cleanly validate IFN also reveal when a program
partition is too coarse or biologically mixed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs" / "canonical_v8"
PER_COHORT = OUT / "per_cohort_effects.csv"
MANIFEST = ROOT / "config" / "cohort_manifest.tsv"
LOPO = OUT / "lopo_cross_validation_expanded.json"
I2 = OUT / "i2_decomposition_expanded.json"
HK = OUT / "hartung_knapp_expanded.json"
OUT_JSON = OUT / "failure_mode_analysis.json"

PROGRAM_MAP = {
    "hallmark_inflammatory_response": "inflammation",
    "hallmark_tnfa_nfkb": "inflammation",
    "hallmark_hypoxia": "hypoxia",
    "hallmark_e2f_targets": "proliferation",
    "hallmark_emt": "emt",
}


def top_off_program_effects(df: pd.DataFrame, sig: str, home_program: str, top_n: int = 5) -> list[dict]:
    subset = df[(df["signature_id"] == sig) & (df["biological_program"] != home_program)].copy()
    subset["abs_d"] = subset["cohens_d"].abs()
    subset = subset.sort_values("abs_d", ascending=False).head(top_n)
    return [
        {
            "cohort_id": row["cohort_id"],
            "foreign_program": row["biological_program"],
            "cohens_d": round(float(row["cohens_d"]), 4),
            "case_label": row["case_label"],
            "control_label": row["control_label"],
        }
        for _, row in subset.iterrows()
    ]


def own_vs_foreign_summary(df: pd.DataFrame, sig: str, home_program: str) -> dict:
    own = df[(df["signature_id"] == sig) & (df["biological_program"] == home_program)]
    foreign = df[(df["signature_id"] == sig) & (df["biological_program"] != home_program)]
    return {
        "home_program": home_program,
        "mean_abs_effect_home": round(float(own["cohens_d"].abs().mean()), 4),
        "mean_abs_effect_foreign": round(float(foreign["cohens_d"].abs().mean()), 4),
        "top_off_program_effects": top_off_program_effects(df, sig, home_program),
    }


def e2f_boundary_analysis(df: pd.DataFrame) -> dict:
    subset = df[(df["signature_id"] == "hallmark_e2f_targets") & (df["biological_program"] == "proliferation")].copy()
    subset = subset.sort_values("cohens_d", ascending=False)
    rows = [
        {
            "cohort_id": row["cohort_id"],
            "cohens_d": round(float(row["cohens_d"]), 4),
            "case_label": row["case_label"],
            "control_label": row["control_label"],
            "tissue": row["tissue"],
        }
        for _, row in subset.iterrows()
    ]
    return {
        "within_proliferation_effect_span": {
            "min": round(float(subset["cohens_d"].min()), 4),
            "max": round(float(subset["cohens_d"].max()), 4),
            "span": round(float(subset["cohens_d"].max() - subset["cohens_d"].min()), 4),
        },
        "within_proliferation_details": rows,
        "interpretation": (
            "E2F targets have high within-program heterogeneity not because the signature is random, "
            "but because the coarse 'proliferation' bucket mixes tumor-vs-normal, prognosis, and subtype contrasts."
        ),
    }


def ambiguous_lopo_signature(sig: str, lopo_data: dict) -> dict:
    full = float(lopo_data[sig]["full"])
    drops = {program: round(full - float(value), 4) for program, value in lopo_data[sig]["loo"].items()}
    biggest = max(drops, key=drops.get)
    return {
        "full_Q_B_fraction": round(full, 4),
        "drop_by_hidden_program": drops,
        "largest_drop_when_hiding": biggest,
        "largest_drop_value": drops[biggest],
    }


def main() -> None:
    df = pd.read_csv(PER_COHORT).merge(
        pd.read_csv(MANIFEST, sep="\t")[["cohort_id", "biological_program", "case_label", "control_label", "tissue"]],
        on="cohort_id",
        how="left",
    )
    lopo = json.loads(LOPO.read_text())
    i2 = json.loads(I2.read_text())
    hk = json.loads(HK.read_text())

    results = {
        "description": (
            "Failure-mode analysis showing that ambiguous diagnostics correspond to "
            "coarse or mixed biological program labels rather than broken signatures."
        ),
        "cross_talk_summary": {
            sig: own_vs_foreign_summary(df, sig, home_program)
            for sig, home_program in PROGRAM_MAP.items()
        },
        "ambiguous_program_examples": {
            "hallmark_inflammatory_response": ambiguous_lopo_signature("hallmark_inflammatory_response", lopo),
            "hallmark_tnfa_nfkb": ambiguous_lopo_signature("hallmark_tnfa_nfkb", lopo),
            "hallmark_e2f_targets": {
                **ambiguous_lopo_signature("hallmark_e2f_targets", lopo),
                "q_between_fraction": i2["hallmark_e2f_targets"]["Q_B_fraction"],
                "within_program_p_bonf": hk["hallmark_e2f_targets"]["p_bonf_HKSJ_guarded"],
                "within_program_i2": hk["hallmark_e2f_targets"]["I2"],
                **e2f_boundary_analysis(df),
            },
        },
    }

    with open(OUT_JSON, "w") as handle:
        json.dump(results, handle, indent=2)

    print("=" * 80)
    print("FAILURE-MODE ANALYSIS (V8)")
    print("=" * 80)
    for sig in ["hallmark_inflammatory_response", "hallmark_tnfa_nfkb"]:
        item = results["ambiguous_program_examples"][sig]
        print(
            f"{sig:35s} full={item['full_Q_B_fraction']:.3f} "
            f"largest_drop={item['largest_drop_when_hiding']} "
            f"(Δ={item['largest_drop_value']:+.3f})"
        )
    e2f = results["ambiguous_program_examples"]["hallmark_e2f_targets"]
    print(
        f"hallmark_e2f_targets{'':15s} full={e2f['full_Q_B_fraction']:.3f} "
        f"span={e2f['within_proliferation_effect_span']['span']:.3f} "
        f"within_I2={e2f['within_program_i2']:.3f}"
    )
    print(f"\nSaved to {OUT_JSON}")


if __name__ == "__main__":
    main()
