#!/usr/bin/env python3
"""Deterministic second positive generalization case study."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

import sys

sys.path.insert(0, str(ROOT / "src"))

from signature_durability_benchmark.config import load_config
from signature_durability_benchmark.diagnostic import run_triage
from signature_durability_benchmark.generalization import (
    EXPECTED_GENERALIZATION_WINNER,
    GENERALIZATION_CANDIDATES,
    rank_generalization_candidates,
    render_generalization_markdown,
    summarize_foreign_effects,
)
from signature_durability_benchmark.meta_analysis import guarded_hksj_random_effects_meta
from signature_durability_benchmark.utils import ensure_dir, read_json, read_table, write_json, write_text


CONFIG = ROOT / "config" / "benchmark_config.yaml"
OUT_DIR = ROOT / "outputs" / "canonical_v8"
PER_COHORT = OUT_DIR / "per_cohort_effects.csv"
MANIFEST = ROOT / "config" / "cohort_manifest.tsv"
HARTUNG_KNAPP = OUT_DIR / "hartung_knapp_expanded.json"
HELD_OUT = OUT_DIR / "external_validation.json"
EXTERNAL_HYPOXIA = OUT_DIR / "external_hypoxia_validation.json"
SIGNATURES = ROOT / "data" / "freeze" / "signatures.tsv"
OUT_JSON = OUT_DIR / "generalization_case_study.json"
OUT_MD = OUT_DIR / "generalization_case_study.md"
ASSET_DIR = OUT_DIR / "generalization_case_study_assets"


def _home_program_cohorts(per_cohort: pd.DataFrame, *, signature_id: str, home_program: str, limit: int = 6) -> list[dict[str, object]]:
    subset = per_cohort.loc[
        (per_cohort["signature_id"] == signature_id)
        & (per_cohort["biological_program"] == home_program)
    ].copy()
    subset = subset.sort_values("cohens_d", ascending=False).head(limit)
    return [
        {
            "cohort_id": str(row["cohort_id"]),
            "cohens_d": float(row["cohens_d"]),
            "case_label": str(row["case_label"]),
            "control_label": str(row["control_label"]),
        }
        for _, row in subset.iterrows()
    ]


def main() -> None:
    per_cohort = read_table(PER_COHORT)
    manifest = read_table(MANIFEST)[["cohort_id", "biological_program", "case_label", "control_label"]]
    per_cohort = per_cohort.merge(manifest, on="cohort_id", how="left")
    hartung_knapp = read_json(HARTUNG_KNAPP)
    held_out = read_json(HELD_OUT)["signatures"]
    if not EXTERNAL_HYPOXIA.exists():
        subprocess.run([sys.executable, str(ROOT / "scripts" / "external_hypoxia_validation.py")], check=True)
    external_hypoxia = read_json(EXTERNAL_HYPOXIA) if EXTERNAL_HYPOXIA.exists() else None

    candidate_rows: list[dict[str, object]] = []
    for signature_id in GENERALIZATION_CANDIDATES:
        home_program = str(held_out[signature_id]["program"])
        within = per_cohort.loc[
            (per_cohort["signature_id"] == signature_id)
            & (per_cohort["biological_program"] == home_program)
        ].copy()
        outside = per_cohort.loc[
            (per_cohort["signature_id"] == signature_id)
            & (per_cohort["biological_program"] != home_program)
        ].copy()
        within_meta = guarded_hksj_random_effects_meta(within["cohens_d"].tolist(), within["cohens_d_var"].tolist())
        outside_meta = guarded_hksj_random_effects_meta(outside["cohens_d"].tolist(), outside["cohens_d_var"].tolist())
        p_bonf = float(hartung_knapp[signature_id]["p_bonf_HKSJ_guarded"])
        candidate_rows.append(
            {
                "signature_id": signature_id,
                "home_program": home_program,
                "frozen_p_bonf_hksj_guarded": p_bonf,
                "passes_p_bonf": p_bonf < 0.05,
                "loo_accuracy": float(held_out[signature_id]["leave_one_out"]["prediction_accuracy"]),
                "loo_perfect": float(held_out[signature_id]["leave_one_out"]["prediction_accuracy"]) == 1.0,
                "split_sign_agreement": float(held_out[signature_id]["split_half"]["sign_agreement_rate"]),
                "split_perfect": float(held_out[signature_id]["split_half"]["sign_agreement_rate"]) == 1.0,
                "mean_abs_foreign_effect": float(outside["cohens_d"].abs().mean()),
                "within_outside_gap": float(within_meta["pooled_effect"] - abs(outside_meta["pooled_effect"])),
                "within_effect": float(within_meta["pooled_effect"]),
                "outside_effect": float(outside_meta["pooled_effect"]),
                "within_p": float(within_meta["pooled_p"]),
                "outside_p": float(outside_meta["pooled_p"]),
            }
        )

    ranked = rank_generalization_candidates(candidate_rows)
    for index, row in enumerate(ranked, start=1):
        row["rank"] = index
    chosen = dict(ranked[0])

    signatures = read_table(SIGNATURES)
    chosen_signature = signatures.loc[
        signatures["signature_id"] == chosen["signature_id"],
        ["gene_symbol", "direction", "weight"],
    ].copy()
    triage_dir = ensure_dir(ASSET_DIR / f"triage_{chosen['signature_id']}")
    input_tsv = triage_dir / f"{chosen['signature_id']}.tsv"
    chosen_signature.to_csv(input_tsv, sep="\t", index=False)
    config = load_config(CONFIG)
    run_triage(config, input_tsv, triage_dir, declared_program=str(chosen["home_program"]))
    triage_json = read_json(triage_dir / "diagnostic.json")
    triage_summary = (triage_dir / "diagnostic_summary.md").read_text(encoding="utf-8")

    chosen["triage_best_program_matches_home"] = bool(
        triage_json.get("inferred_program", {}).get("program") == chosen["home_program"]
    )
    chosen["triage_full_model_class"] = str(triage_json["classification"]["full_model"])
    chosen["triage_within_program_class"] = str(triage_json["classification"]["within_program"])

    payload = {
        "description": (
            "Deterministic selection of the second non-IFN positive generalization case "
            "from hallmark_hypoxia versus hallmark_emt using frozen benchmark artifacts."
        ),
        "selection_rule": {
            "candidates": list(GENERALIZATION_CANDIDATES),
            "ranking": [
                "passes p_bonf_HKSJ_guarded < 0.05",
                "perfect leave-one-out sign prediction",
                "perfect split-half sign agreement",
                "lower mean absolute foreign-program effect",
                "higher within-minus-abs(outside) separation",
            ],
            "expected_winner": EXPECTED_GENERALIZATION_WINNER,
        },
        "selection_rule_satisfied": ranked[0]["signature_id"] == EXPECTED_GENERALIZATION_WINNER,
        "ranked_candidates": ranked,
        "selected_case": chosen,
        "home_program_cohorts": _home_program_cohorts(
            per_cohort,
            signature_id=str(chosen["signature_id"]),
            home_program=str(chosen["home_program"]),
        ),
        "largest_foreign_effects": summarize_foreign_effects(
            per_cohort,
            signature_id=str(chosen["signature_id"]),
            home_program=str(chosen["home_program"]),
        ),
        "triage_summary": triage_summary,
        "triage_output_dir": str(triage_dir),
        "external_validation": external_hypoxia if chosen["signature_id"] == "hallmark_hypoxia" else None,
    }
    write_json(OUT_JSON, payload)
    write_text(OUT_MD, render_generalization_markdown(payload))
    print(json.dumps({"selected_case": chosen["signature_id"], "selection_rule_satisfied": payload["selection_rule_satisfied"]}, indent=2))


if __name__ == "__main__":
    main()
