"""Arbitrary-signature diagnostic / triage workflow."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .benchmark import _build_profile, _compute_null_effects, _score_signature_across_cohorts
from .classify import classify_signature
from .confounders import load_confounder_sets
from .config import SkillConfig
from .meta_analysis import permutation_q_decomposition, within_program_meta
from .normalize import normalize_signature
from .utils import (
    ensure_columns,
    ensure_dir,
    repo_relpath,
    read_table,
    stable_name_seed,
    write_json,
    write_table,
    write_text,
)

PROGRAM_DISPLAY = {
    "interferon": "interferon",
    "inflammation": "inflammation",
    "proliferation": "proliferation",
    "hypoxia": "hypoxia",
    "emt": "EMT",
}


def _load_cohort_data(
    cohort_manifest: pd.DataFrame,
    freeze_dir: Path,
) -> dict[str, dict[str, Any]]:
    cohort_data: dict[str, dict[str, Any]] = {}
    for _, crow in cohort_manifest.iterrows():
        cid = str(crow["cohort_id"])
        expr = read_table(freeze_dir / "cohort_matrices" / f"{cid}.tsv").set_index("gene_symbol")
        pheno = read_table(freeze_dir / "cohort_phenotypes" / f"{cid}.tsv")
        cohort_data[cid] = {"expr": expr, "pheno": pheno, "manifest": crow}
    return cohort_data


def _build_gene_universe(cohort_data: dict[str, dict[str, Any]]) -> list[str]:
    all_genes: set[str] = set()
    for cdata in cohort_data.values():
        all_genes.update(cdata["expr"].index.astype(str).tolist())
    return sorted(all_genes)


def _read_signature_input(path: str | Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    frame = read_table(path)
    if frame.shape[1] == 1:
        frame = pd.DataFrame({"gene_symbol": frame.iloc[:, 0], "direction": "up", "weight": 1.0})
    else:
        if "gene_symbol" not in frame.columns:
            frame = frame.rename(columns={frame.columns[0]: "gene_symbol"})
        if "direction" not in frame.columns:
            frame["direction"] = "up"
        if "weight" not in frame.columns:
            frame["weight"] = 1.0
    ensure_columns(frame, ["gene_symbol", "direction", "weight"], "triage signature")
    normalized, audit = normalize_signature(frame[["gene_symbol", "direction", "weight"]])
    return normalized, audit


def _infer_best_program(program_diagnostics: dict[str, dict[str, Any]]) -> dict[str, Any] | None:
    ranked: list[dict[str, Any]] = []
    for program, diag in program_diagnostics.items():
        within = diag.get("within")
        outside = diag.get("outside")
        if not within or diag.get("within_k", 0) < 3:
            continue
        within_effect = float(within["pooled_effect"])
        outside_effect = float(outside["pooled_effect"]) if outside else 0.0
        separation = within_effect - abs(outside_effect)
        if within_effect <= 0:
            separation -= 1.0
        ranked.append(
            {
                "program": program,
                "separation_score": separation,
                "within_effect": within_effect,
                "outside_effect": outside_effect,
                "within_p": float(within["pooled_p"]),
                "outside_p": float(outside["pooled_p"]) if outside else 1.0,
            }
        )
    if not ranked:
        return None
    ranked.sort(key=lambda row: (row["separation_score"], -row["within_p"]), reverse=True)
    return {"best": ranked[0], "ranked": ranked}


def _interpret_diagnostic(
    signature_name: str,
    focus_program: str,
    profile: dict[str, Any],
    focus_diag: dict[str, Any] | None,
    permutation: dict[str, Any],
    dominant_confounder: dict[str, Any] | None,
) -> str:
    coverage = float(profile.get("mean_coverage", 0.0))
    predicted = profile.get("predicted_class", "unclassified")
    if coverage < 0.60:
        return (
            f"{signature_name} has insufficient coverage across the frozen cohorts "
            f"(mean coverage {coverage:.2f}), so the benchmark cannot make a strong durability call."
        )

    if dominant_confounder and dominant_confounder["weighted_effect"] >= abs(float(profile.get("aggregate_effect", 0.0))):
        return (
            f"{signature_name} looks dominated by the {dominant_confounder['name']} confounder set "
            f"(weighted confounder effect {dominant_confounder['weighted_effect']:.3f} versus aggregate "
            f"effect {abs(float(profile.get('aggregate_effect', 0.0))):.3f}), so the signature is better "
            f"described as confounded than durable."
        )

    if focus_diag and focus_diag.get("within") and focus_diag.get("outside"):
        within = focus_diag["within"]
        outside = focus_diag["outside"]
        if (
            within["pooled_effect"] > 0.5
            and within["pooled_p"] < 0.05
            and abs(outside["pooled_effect"]) < 0.35
            and outside["pooled_p"] >= 0.05
            and permutation["p_value"] < 0.05
        ):
            return (
                f"{signature_name} behaves like a context-matched {PROGRAM_DISPLAY.get(focus_program, focus_program)} "
                f"signature: the within-program pooled effect is {within['pooled_effect']:.3f} "
                f"(p={within['pooled_p']:.4g}) while the outside-program effect shrinks to "
                f"{outside['pooled_effect']:.3f} (p={outside['pooled_p']:.4g}), and the observed "
                f"program structure exceeds random labelings (permutation p={permutation['p_value']:.4g})."
            )

    if predicted == "mixed":
        return (
            f"{signature_name} shows real signal but does not reduce cleanly to one program. "
            f"The aggregate effect is {profile['aggregate_effect']:.3f} with I²={profile['i_squared']:.3f}, "
            f"suggesting either cross-context activation or overly broad program boundaries."
        )

    if predicted == "brittle":
        return (
            f"{signature_name} does not show a stable, interpretable cross-cohort pattern in this benchmark. "
            f"The aggregate effect is {profile['aggregate_effect']:.3f} with aggregate p={profile['aggregate_p']:.4g} "
            f"and direction consistency {profile['direction_consistency']:.2f}."
        )

    return (
        f"{signature_name} receives a {predicted} label under the within-program benchmark rules, "
        f"but the main value is the program-conditioned breakdown rather than the top-line class alone."
    )


def run_triage(
    config: SkillConfig,
    input_path: str | Path,
    out_dir: str | Path,
    declared_program: str | None = None,
) -> Path:
    signature_path = Path(input_path).resolve()
    output_dir = ensure_dir(out_dir)
    signature_name = signature_path.stem

    signature_df, audit = _read_signature_input(signature_path)
    manifest = read_table(config.path("cohort_manifest"))
    freeze_dir = config.path("freeze_dir")
    cohort_data = _load_cohort_data(manifest, freeze_dir)
    confounders = load_confounder_sets(config.path("confounder_panel"))
    universe = _build_gene_universe(cohort_data)

    scored = _score_signature_across_cohorts(
        signature_name,
        signature_df,
        cohort_data,
        confounders,
    )
    per_cohort_df = pd.DataFrame(scored["per_cohort_records"]).merge(
        manifest[["cohort_id", "biological_program", "tissue", "platform", "platform_id"]],
        on="cohort_id",
        how="left",
    )

    null_effects = _compute_null_effects(
        signature_df,
        universe,
        int(config.scoring.get("null_draws", 200)),
        int(config.runtime.get("random_seed", 42)) + stable_name_seed(signature_name),
        cohort_data,
    )
    profile = _build_profile(
        scored["cohort_effects"],
        scored["cohort_variances"],
        scored["cohort_platforms"],
        scored["cohort_confounder_maxes"],
        scored["coverages"],
        null_effects,
    )

    groups = [str(cdata["manifest"]["biological_program"]) for cdata in cohort_data.values()]
    permutation = permutation_q_decomposition(
        scored["cohort_effects"],
        scored["cohort_variances"],
        groups,
        n_permutations=1000,
        seed=int(config.runtime.get("random_seed", 42)),
    )

    program_diagnostics: dict[str, dict[str, Any]] = {}
    for program in sorted(set(groups)):
        meta = within_program_meta(
            scored["cohort_effects"],
            scored["cohort_variances"],
            groups,
            program,
        )
        within = meta["within"]
        outside = meta["outside"]
        gap = None
        if within and outside:
            gap = float(within["pooled_effect"] - abs(outside["pooled_effect"]))
        program_diagnostics[program] = {
            "within_k": meta["within_k"],
            "outside_k": meta["outside_k"],
            "within": within,
            "outside": outside,
            "within_outside_gap": gap,
        }

    inferred = _infer_best_program(program_diagnostics)
    focus_program = declared_program or (inferred["best"]["program"] if inferred else None)
    focus_diag = program_diagnostics.get(focus_program) if focus_program else None

    profile["predicted_class"] = classify_signature(profile, "full_model")["predicted_class"]
    if focus_diag and focus_diag.get("within"):
        profile["within_p"] = focus_diag["within"]["pooled_p"]
        profile["within_d"] = focus_diag["within"]["pooled_effect"]
    within_program_class = classify_signature(profile, "within_program")["predicted_class"]

    overlap_by_confounder = {}
    signature_genes = set(signature_df["gene_symbol"].astype(str).str.upper())
    for name, cdf in confounders.items():
        conf_genes = set(cdf["gene_symbol"].astype(str).str.upper())
        overlap_by_confounder[name] = len(signature_genes & conf_genes) / len(signature_genes) if signature_genes else 0.0
    dominant_name = max(overlap_by_confounder, key=overlap_by_confounder.get) if overlap_by_confounder else None
    dominant_confounder = None
    if dominant_name:
        dominant_confounder = {
            "name": dominant_name,
            "overlap_fraction": overlap_by_confounder[dominant_name],
            "weighted_effect": float(profile.get("max_confounder_effect", 0.0)),
        }

    interpretation = _interpret_diagnostic(
        signature_name,
        focus_program or "unknown",
        profile,
        focus_diag,
        permutation,
        dominant_confounder,
    )

    result = {
        "signature_name": signature_name,
        "input_path": repo_relpath(signature_path),
        "declared_program": declared_program,
        "inferred_program": inferred["best"] if inferred else None,
        "program_ranking": inferred["ranked"] if inferred else [],
        "normalization_audit": audit,
        "aggregate_profile": profile,
        "permutation_program_structure": permutation,
        "program_diagnostics": program_diagnostics,
        "classification": {
            "full_model": profile["predicted_class"],
            "within_program": within_program_class,
            "focus_program": focus_program,
        },
        "confounder_overlap": overlap_by_confounder,
        "dominant_confounder": dominant_confounder,
        "interpretation": interpretation,
    }

    summary_lines = [
        f"# Diagnostic Report: {signature_name}",
        "",
        f"- Input: `{signature_path.name}`",
        f"- Inferred best-supported program: `{focus_program}`" if focus_program else "- Inferred best-supported program: none",
        f"- Full-model classification: `{profile['predicted_class']}`",
        f"- Within-program classification: `{within_program_class}`",
        f"- Mean coverage: `{profile['mean_coverage']:.3f}`",
        f"- Aggregate effect: `{profile['aggregate_effect']:.3f}` (p=`{profile['aggregate_p']:.4g}`)",
        f"- Program-structure permutation p: `{permutation['p_value']:.4g}`",
        "",
        interpretation,
        "",
        "## Program Ranking",
        "",
    ]
    if inferred and inferred["ranked"]:
        for row in inferred["ranked"]:
            summary_lines.append(
                f"- `{row['program']}`: separation={row['separation_score']:.3f}, "
                f"within={row['within_effect']:.3f} (p={row['within_p']:.4g}), "
                f"outside={row['outside_effect']:.3f} (p={row['outside_p']:.4g})"
            )
    else:
        summary_lines.append("- No program had enough within-program cohorts for ranking.")

    write_json(output_dir / "diagnostic.json", result)
    write_table(per_cohort_df, output_dir / "per_cohort_effects.csv")
    write_text(output_dir / "diagnostic_summary.md", "\n".join(summary_lines) + "\n")
    return output_dir
