"""Helpers for selecting and rendering the second positive generalization case."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


GENERALIZATION_CANDIDATES = ("hallmark_hypoxia", "hallmark_emt")
EXPECTED_GENERALIZATION_WINNER = "hallmark_hypoxia"


def rank_generalization_candidates(
    candidate_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    ranked = sorted(
        candidate_rows,
        key=lambda row: (
            bool(row["passes_p_bonf"]),
            bool(row["loo_perfect"]),
            bool(row["split_perfect"]),
            -float(row["mean_abs_foreign_effect"]),
            float(row["within_outside_gap"]),
            row["signature_id"],
        ),
        reverse=True,
    )
    return ranked


def summarize_foreign_effects(
    per_cohort: pd.DataFrame,
    *,
    signature_id: str,
    home_program: str,
    limit: int = 6,
) -> list[dict[str, Any]]:
    foreign = per_cohort.loc[
        (per_cohort["signature_id"] == signature_id)
        & (per_cohort["biological_program"] != home_program)
    ].copy()
    foreign["abs_d"] = foreign["cohens_d"].abs()
    foreign = foreign.sort_values(["abs_d", "cohens_d"], ascending=[False, False]).head(limit)
    return [
        {
            "cohort_id": str(row["cohort_id"]),
            "foreign_program": str(row["biological_program"]),
            "cohens_d": float(row["cohens_d"]),
            "case_label": str(row["case_label"]),
            "control_label": str(row["control_label"]),
        }
        for _, row in foreign.iterrows()
    ]


def render_generalization_markdown(payload: dict[str, Any]) -> str:
    chosen = payload["selected_case"]
    lines = [
        "# Second Positive Generalization Case",
        "",
        f"- Selected signature: `{chosen['signature_id']}`",
        f"- Home program: `{chosen['home_program']}`",
        f"- Ranking expectation: `{EXPECTED_GENERALIZATION_WINNER}`",
        f"- Selection rule satisfied: `{payload['selection_rule_satisfied']}`",
        "",
        "## Candidate Ranking",
        "",
        "| Rank | Signature | p_Bonf<0.05 | LOO=1.0 | Split=1.0 | mean abs foreign d | within-minus-abs(outside) gap |",
        "|---:|---|---:|---:|---:|---:|---:|",
    ]
    for row in payload["ranked_candidates"]:
        lines.append(
            f"| {row['rank']} | {row['signature_id']} | {row['passes_p_bonf']} | "
            f"{row['loo_perfect']} | {row['split_perfect']} | "
            f"{row['mean_abs_foreign_effect']:.3f} | {row['within_outside_gap']:.3f} |"
        )

    lines.extend(
        [
            "",
            "## Selected Case Checks",
            "",
            f"- Inferred program matches home program: `{chosen['triage_best_program_matches_home']}`",
            f"- Within-program class: `{chosen['triage_within_program_class']}`",
            f"- Full-model class: `{chosen['triage_full_model_class']}`",
            f"- Positive statistically supported pooled effect (frozen artifact): `{chosen['passes_p_bonf']}`",
            f"- Leave-one-out sign prediction: `{chosen['loo_accuracy']:.3f}`",
            f"- Split-half sign agreement: `{chosen['split_sign_agreement']:.3f}`",
            "",
            "## Home-Program Cohorts",
            "",
            "| Cohort | Hedges' g | Contrast |",
            "|---|---:|---|",
        ]
    )
    for row in payload["home_program_cohorts"]:
        lines.append(
            f"| {row['cohort_id']} | {row['cohens_d']:+.3f} | {row['case_label']} vs {row['control_label']} |"
        )

    lines.extend(
        [
            "",
            "## Largest Foreign-Program Effects",
            "",
            "| Cohort | Foreign program | Hedges' g | Contrast |",
            "|---|---|---:|---|",
        ]
    )
    for row in payload["largest_foreign_effects"]:
        lines.append(
            f"| {row['cohort_id']} | {row['foreign_program']} | {row['cohens_d']:+.3f} | {row['case_label']} vs {row['control_label']} |"
        )

    external = payload.get("external_validation")
    if external:
        lines.extend(
            [
                "",
                "## External Support",
                "",
                f"- Cohort: `{external['cohort']['external_cohort_id']}` ({external['cohort']['geo_accession']})",
                f"- Chosen-signature external effect: `{external['chosen_signature_support']['effect_size_g']:+.3f}`",
                f"- Rank among externally scored signatures: `{external['chosen_signature_support']['rank_among_scored_signatures']}`",
                f"- Positive direction: `{external['chosen_signature_support']['positive_direction']}`",
                "",
            ]
        )

    lines.extend(
        [
            "## Triage Interpretation",
            "",
            payload["triage_summary"].strip(),
            "",
        ]
    )
    return "\n".join(lines)
