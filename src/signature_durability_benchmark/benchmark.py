"""Full benchmark pipeline runner.

Ties together scoring, meta-analysis, null model, confounders, and
classification to produce a complete durability assessment of every
signature in the panel across all cohorts.
"""
from __future__ import annotations

import importlib
import logging
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

from .classify import classify_signature
from .confounders import load_confounder_sets, score_confounders_in_cohort
from .constants import ALL_CLASSES, DURABLE_CLASSES, MODEL_NAMES, REQUIRED_OUTPUTS
from .config import SkillConfig
from .meta_analysis import (
    direction_consistency,
    fixed_effect_meta,
    i_squared,
    leave_one_out,
    loo_stability,
    platform_holdout,
    platform_holdout_consistency,
)
from .normalize import normalize_signature
from .null_model import generate_null_signatures, null_separation_p
from .scoring import score_signature_in_cohort
from .utils import (
    ensure_dir,
    now_timestamp,
    read_json,
    read_table,
    runtime_summary,
    set_runtime_environment,
    sha256_file,
    write_json,
    write_table,
    write_text,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_cohort_data(
    cohort_manifest: pd.DataFrame,
    freeze_dir: Path,
) -> dict[str, dict[str, Any]]:
    """Load expression matrices and phenotype tables for every cohort."""
    cohort_data: dict[str, dict[str, Any]] = {}
    for _, crow in cohort_manifest.iterrows():
        cid = str(crow["cohort_id"])
        expr = read_table(freeze_dir / "cohort_matrices" / f"{cid}.tsv")
        expr = expr.set_index(expr.columns[0])  # gene_symbol as index
        pheno = read_table(freeze_dir / "cohort_phenotypes" / f"{cid}.tsv")
        cohort_data[cid] = {"expr": expr, "pheno": pheno, "manifest": crow}
    return cohort_data


def _build_gene_universe(cohort_data: dict[str, dict[str, Any]]) -> list[str]:
    """Union of all gene symbols across cohorts, sorted for determinism."""
    all_genes: set[str] = set()
    for cdata in cohort_data.values():
        all_genes.update(cdata["expr"].index.astype(str).tolist())
    return sorted(all_genes)


def _score_signature_across_cohorts(
    sig_id: str,
    raw_sig: pd.DataFrame,
    cohort_data: dict[str, dict[str, Any]],
    confounder_sets: dict[str, pd.DataFrame],
) -> dict[str, Any]:
    """Score one signature in every cohort and return per-cohort records + aggregates."""
    per_cohort_records: list[dict[str, Any]] = []
    cohort_effects: list[float] = []
    cohort_variances: list[float] = []
    cohort_platforms: list[str] = []
    cohort_confounder_maxes: list[float] = []
    coverages: list[float] = []

    for cid, cdata in cohort_data.items():
        crow = cdata["manifest"]
        result = score_signature_in_cohort(
            raw_sig,
            cdata["expr"],
            cdata["pheno"],
            str(crow["phenotype_column"]),
            str(crow["case_label"]),
            str(crow["control_label"]),
        )

        per_cohort_records.append({
            "signature_id": sig_id,
            "cohort_id": cid,
            "cohens_d": result["cohens_d"],
            "cohens_d_var": result["cohens_d_var"],
            "direction_consistent": result["direction_consistent"],
            "coverage_fraction": result["coverage_fraction"],
            "case_n": result["case_n"],
            "control_n": result["control_n"],
        })

        cohort_effects.append(result["cohens_d"])
        cohort_variances.append(result["cohens_d_var"])
        cohort_platforms.append(str(crow["platform"]))
        coverages.append(result["coverage_fraction"])

        # Score confounders in this cohort
        conf_scores = score_confounders_in_cohort(
            confounder_sets,
            cdata["expr"],
            cdata["pheno"],
            str(crow["phenotype_column"]),
            str(crow["case_label"]),
            str(crow["control_label"]),
        )
        cohort_confounder_maxes.append(
            max(conf_scores.values()) if conf_scores else 0.0
        )

    return {
        "per_cohort_records": per_cohort_records,
        "cohort_effects": cohort_effects,
        "cohort_variances": cohort_variances,
        "cohort_platforms": cohort_platforms,
        "cohort_confounder_maxes": cohort_confounder_maxes,
        "coverages": coverages,
    }


def _compute_null_effects(
    raw_sig: pd.DataFrame,
    universe: list[str],
    null_draws: int,
    seed: int,
    cohort_data: dict[str, dict[str, Any]],
) -> list[float]:
    """Generate null signatures and compute their pooled meta-analytic effects."""
    null_sigs = generate_null_signatures(len(raw_sig), universe, null_draws, seed)
    null_effects: list[float] = []
    for nsig in null_sigs:
        null_cohort_effects: list[float] = []
        for cid, cdata in cohort_data.items():
            crow = cdata["manifest"]
            nr = score_signature_in_cohort(
                nsig,
                cdata["expr"],
                cdata["pheno"],
                str(crow["phenotype_column"]),
                str(crow["case_label"]),
                str(crow["control_label"]),
            )
            null_cohort_effects.append(nr["cohens_d"])
        # Use uniform variance weighting for null (no real variance estimate)
        null_meta = fixed_effect_meta(
            null_cohort_effects, [1.0] * len(null_cohort_effects)
        )
        null_effects.append(null_meta["pooled_effect"])
    return null_effects


def _build_profile(
    cohort_effects: list[float],
    cohort_variances: list[float],
    cohort_platforms: list[str],
    cohort_confounder_maxes: list[float],
    coverages: list[float],
    null_effects: list[float],
) -> dict[str, float]:
    """Compute the aggregate profile dict used for classification."""
    mean_coverage = float(np.mean(coverages)) if coverages else 0.0
    meta = fixed_effect_meta(cohort_effects, cohort_variances)
    i_sq = i_squared(cohort_effects, cohort_variances)
    loo = loo_stability(cohort_effects, cohort_variances)
    dir_cons = direction_consistency(cohort_effects)
    plat_cons = platform_holdout_consistency(
        cohort_effects, cohort_variances, cohort_platforms
    )
    max_conf = (
        float(np.mean(cohort_confounder_maxes))
        if cohort_confounder_maxes
        else 0.0
    )
    null_p = null_separation_p(meta["pooled_effect"], null_effects)

    return {
        "mean_coverage": mean_coverage,
        "aggregate_effect": meta["pooled_effect"],
        "aggregate_p": meta["pooled_p"],
        "direction_consistency": dir_cons,
        "i_squared": i_sq,
        "loo_stability": loo,
        "max_confounder_effect": max_conf,
        "null_separation_p": null_p,
        "platform_holdout_consistency": plat_cons,
    }


def _classify_all_models(
    sig_id: str,
    srow: pd.Series,
    profile: dict[str, float],
) -> list[dict[str, Any]]:
    """Classify a signature under all 5 models and return records."""
    records: list[dict[str, Any]] = []
    for model_name in MODEL_NAMES:
        result = classify_signature(profile, model_name)
        records.append({
            "signature_id": sig_id,
            "split": str(srow["split"]),
            "expected_class": str(srow["expected_class"]),
            "predicted_class": result["predicted_class"],
            "model": model_name,
            **profile,
        })
    return records


def _build_loo_table(
    sig_id: str,
    cohort_effects: list[float],
    cohort_variances: list[float],
    cohort_ids: list[str],
) -> list[dict[str, Any]]:
    """Build leave-one-cohort-out records for one signature."""
    loo_results = leave_one_out(cohort_effects, cohort_variances)
    rows: list[dict[str, Any]] = []
    for rec in loo_results:
        dropped_idx = rec["dropped_index"]
        rows.append({
            "signature_id": sig_id,
            "dropped_cohort": cohort_ids[dropped_idx] if dropped_idx < len(cohort_ids) else "unknown",
            "pooled_effect": rec["pooled_effect"],
            "pooled_se": rec["pooled_se"],
            "pooled_p": rec["pooled_p"],
        })
    return rows


def _build_platform_holdout_table(
    sig_id: str,
    cohort_effects: list[float],
    cohort_variances: list[float],
    cohort_platforms: list[str],
) -> list[dict[str, Any]]:
    """Build platform-holdout records for one signature."""
    holdout_results = platform_holdout(cohort_effects, cohort_variances, cohort_platforms)
    rows: list[dict[str, Any]] = []
    for plat, meta_result in holdout_results.items():
        rows.append({
            "signature_id": sig_id,
            "held_out_platform": plat,
            "pooled_effect": meta_result["pooled_effect"],
            "pooled_se": meta_result["pooled_se"],
            "pooled_p": meta_result["pooled_p"],
        })
    return rows


def _build_null_summary(
    sig_id: str,
    observed_effect: float,
    null_effects: list[float],
    null_p: float,
) -> dict[str, Any]:
    """Build a null-model summary record for one signature."""
    return {
        "signature_id": sig_id,
        "observed_effect": observed_effect,
        "null_mean": float(np.mean(null_effects)) if null_effects else 0.0,
        "null_sd": float(np.std(null_effects)) if null_effects else 0.0,
        "null_draws": len(null_effects),
        "null_separation_p": null_p,
    }


def _build_cohort_overlap_summary(
    sig_manifest: pd.DataFrame,
    sig_rows: pd.DataFrame,
    cohort_data: dict[str, dict[str, Any]],
) -> pd.DataFrame:
    """For each signature x cohort, show gene overlap count and coverage."""
    rows: list[dict[str, Any]] = []
    for _, srow in sig_manifest.iterrows():
        sig_id = str(srow["signature_id"])
        genes = set(
            sig_rows.loc[sig_rows["signature_id"] == sig_id, "gene_symbol"]
            .astype(str)
            .str.upper()
            .str.strip()
        )
        for cid, cdata in cohort_data.items():
            cohort_genes = set(cdata["expr"].index.astype(str))
            overlap = genes & cohort_genes
            rows.append({
                "signature_id": sig_id,
                "cohort_id": cid,
                "signature_genes": len(genes),
                "cohort_genes": len(cohort_genes),
                "overlap_count": len(overlap),
                "coverage_fraction": len(overlap) / len(genes) if genes else 0.0,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Aggregate metrics + success rule
# ---------------------------------------------------------------------------

def _compute_auprc(scores_df: pd.DataFrame, model_name: str) -> float:
    """Compute AUPRC for a single model on the PRIMARY split.

    y_true = 1 if expected_class in DURABLE_CLASSES, else 0.
    y_score = aggregate_effect (larger = more likely durable).
    Rows with predicted_class == 'insufficient_coverage' are excluded.
    """
    subset = scores_df[
        (scores_df["model"] == model_name) & (scores_df["split"] == "primary")
    ].copy()
    subset = subset[subset["predicted_class"] != "insufficient_coverage"]
    if subset.empty:
        return 0.0
    y_true = subset["expected_class"].isin(DURABLE_CLASSES).astype(int).values
    y_score = subset["aggregate_effect"].values
    if y_true.sum() == 0 or y_true.sum() == len(y_true):
        return 0.0  # degenerate: all one class
    return float(average_precision_score(y_true, y_score))


def _compute_exact_accuracy(scores_df: pd.DataFrame, model_name: str, split: str = "primary") -> float:
    """Fraction of signatures where predicted_class == expected_class."""
    subset = scores_df[
        (scores_df["model"] == model_name) & (scores_df["split"] == split)
    ]
    if subset.empty:
        return 0.0
    correct = (subset["predicted_class"] == subset["expected_class"]).sum()
    return float(correct / len(subset))


def _compute_direction_accuracy(scores_df: pd.DataFrame, model_name: str) -> float:
    """Among expected-durable signatures, fraction predicted durable."""
    subset = scores_df[
        (scores_df["model"] == model_name)
        & (scores_df["split"] == "primary")
        & (scores_df["expected_class"].isin(DURABLE_CLASSES))
    ]
    if subset.empty:
        return 0.0
    correct = subset["predicted_class"].isin(DURABLE_CLASSES).sum()
    return float(correct / len(subset))


def _compute_confounded_rejection_accuracy(scores_df: pd.DataFrame, model_name: str) -> float:
    """Among expected-confounded, fraction predicted confounded."""
    subset = scores_df[
        (scores_df["model"] == model_name)
        & (scores_df["split"] == "primary")
        & (scores_df["expected_class"] == "confounded")
    ]
    if subset.empty:
        return 0.0
    correct = (subset["predicted_class"] == "confounded").sum()
    return float(correct / len(subset))


def _compute_blind_exact_recovery(scores_df: pd.DataFrame, model_name: str) -> float:
    """Exact accuracy on the blind split."""
    return _compute_exact_accuracy(scores_df, model_name, split="blind")


def _aggregate_metrics(
    scores_df: pd.DataFrame,
    tolerance: float,
) -> dict[str, Any]:
    """Compute all aggregate metrics and the success rule for every model."""
    per_model: dict[str, dict[str, Any]] = {}
    for model_name in MODEL_NAMES:
        auprc = _compute_auprc(scores_df, model_name)
        exact_acc = _compute_exact_accuracy(scores_df, model_name)
        dir_acc = _compute_direction_accuracy(scores_df, model_name)
        conf_rej = _compute_confounded_rejection_accuracy(scores_df, model_name)
        blind_rec = _compute_blind_exact_recovery(scores_df, model_name)
        per_model[model_name] = {
            "auprc": auprc,
            "exact_class_accuracy": exact_acc,
            "direction_accuracy": dir_acc,
            "confounded_rejection_accuracy": conf_rej,
            "blind_exact_recovery": blind_rec,
        }

    # Success rule: full_model AUPRC > overlap_only AUPRC + tolerance
    # AND secondary wins >= 2
    full = per_model["full_model"]
    overlap = per_model["overlap_only"]
    auprc_margin = full["auprc"] - overlap["auprc"]
    auprc_pass = auprc_margin > tolerance

    secondary_wins = 0
    for metric in ["exact_class_accuracy", "direction_accuracy",
                   "confounded_rejection_accuracy", "blind_exact_recovery"]:
        if full[metric] > overlap[metric]:
            secondary_wins += 1

    success = auprc_pass and secondary_wins >= 2

    return {
        "per_model": per_model,
        "success_rule": {
            "full_model_auprc": full["auprc"],
            "overlap_only_auprc": overlap["auprc"],
            "auprc_margin": auprc_margin,
            "tolerance": tolerance,
            "auprc_pass": auprc_pass,
            "secondary_wins": secondary_wins,
            "secondary_wins_required": 2,
            "success": success,
        },
    }


# ---------------------------------------------------------------------------
# Public summary
# ---------------------------------------------------------------------------

def _generate_public_summary(
    aggregate: dict[str, Any],
    scores_df: pd.DataFrame,
    start_time: float,
) -> str:
    """Generate a Markdown public summary of the benchmark run."""
    rule = aggregate["success_rule"]
    elapsed = time.time() - start_time
    lines = [
        "# Signature Durability Benchmark Results",
        "",
        f"**Run timestamp**: {now_timestamp()}",
        f"**Elapsed**: {elapsed:.1f}s",
        "",
        "## Success Rule",
        "",
        f"- full_model AUPRC: {rule['full_model_auprc']:.4f}",
        f"- overlap_only AUPRC: {rule['overlap_only_auprc']:.4f}",
        f"- AUPRC margin: {rule['auprc_margin']:.4f} (tolerance: {rule['tolerance']})",
        f"- AUPRC pass: {'YES' if rule['auprc_pass'] else 'NO'}",
        f"- Secondary wins: {rule['secondary_wins']}/{rule['secondary_wins_required']} required",
        f"- **Overall: {'PASS' if rule['success'] else 'FAIL'}**",
        "",
        "## Per-Model Metrics",
        "",
        "| Model | AUPRC | Exact Acc | Direction Acc | Conf Reject | Blind Recovery |",
        "|-------|-------|-----------|---------------|-------------|----------------|",
    ]
    for model_name in MODEL_NAMES:
        m = aggregate["per_model"][model_name]
        lines.append(
            f"| {model_name} | {m['auprc']:.4f} | {m['exact_class_accuracy']:.4f} "
            f"| {m['direction_accuracy']:.4f} | {m['confounded_rejection_accuracy']:.4f} "
            f"| {m['blind_exact_recovery']:.4f} |"
        )

    lines += [
        "",
        "## Signature Counts by Split",
        "",
    ]
    for split in ["primary", "blind"]:
        count = scores_df.loc[
            scores_df["model"] == "full_model",
            "split",
        ].eq(split).sum()
        lines.append(f"- **{split}**: {count} signatures")

    lines += [
        "",
        "## Classification Distribution (full_model, primary split)",
        "",
    ]
    primary_full = scores_df[
        (scores_df["model"] == "full_model") & (scores_df["split"] == "primary")
    ]
    for cls in sorted(ALL_CLASSES):
        n = (primary_full["predicted_class"] == cls).sum()
        lines.append(f"- {cls}: {n}")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Optional module loaders (certificates, plots)
# ---------------------------------------------------------------------------

def _try_generate_certificates(
    scores_df: pd.DataFrame,
    aggregate: dict[str, Any],
    out_path: Path,
) -> dict[str, str | None]:
    """Attempt to generate certificates if the module exists."""
    results: dict[str, str | None] = {}
    try:
        certs = importlib.import_module(".certificates", package=__package__)
    except (ImportError, ModuleNotFoundError):
        logger.info("certificates module not available; skipping certificate generation")
        return results

    cert_specs = [
        ("durability_certificate.json", "generate_durability_certificate"),
        ("platform_transfer_certificate.json", "generate_platform_transfer_certificate"),
        ("confounder_rejection_certificate.json", "generate_confounder_rejection_certificate"),
        ("coverage_certificate.json", "generate_coverage_certificate"),
    ]
    for filename, func_name in cert_specs:
        fn = getattr(certs, func_name, None)
        if fn is not None:
            try:
                cert = fn(scores_df, aggregate)
                cert_path = out_path / filename
                write_json(cert_path, cert)
                results[filename] = str(cert_path)
            except Exception as exc:
                logger.warning("Certificate %s failed: %s", filename, exc)
                results[filename] = None
        else:
            logger.info("Certificate function %s not found", func_name)
    return results


def _try_generate_plots(
    scores_df: pd.DataFrame,
    per_cohort_df: pd.DataFrame,
    null_summary_df: pd.DataFrame,
    loo_df: pd.DataFrame,
    platform_holdout_df: pd.DataFrame,
    out_path: Path,
) -> dict[str, str | None]:
    """Attempt to generate plots if the module exists."""
    results: dict[str, str | None] = {}
    try:
        plots = importlib.import_module(".plots", package=__package__)
    except (ImportError, ModuleNotFoundError):
        logger.info("plots module not available; skipping plot generation")
        return results

    plot_specs = [
        ("forest_plot.png", "forest_plot", (per_cohort_df,)),
        ("null_separation_plot.png", "null_separation_plot", (null_summary_df,)),
        ("stability_heatmap.png", "stability_heatmap", (loo_df,)),
        ("platform_transfer_panel.png", "platform_transfer_panel", (platform_holdout_df,)),
    ]
    for filename, func_name, args in plot_specs:
        fn = getattr(plots, func_name, None)
        if fn is not None:
            try:
                plot_path = out_path / filename
                fn(*args, plot_path)
                results[filename] = str(plot_path)
            except Exception as exc:
                logger.warning("Plot %s failed: %s", filename, exc)
                results[filename] = None
        else:
            logger.info("Plot function %s not found", func_name)
    return results


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_pipeline(config: SkillConfig, out_dir: str | Path) -> dict[str, Any]:
    """Execute the full benchmark pipeline.

    Parameters
    ----------
    config : SkillConfig
        Loaded configuration (via ``load_config``).
    out_dir : str | Path
        Directory where all outputs will be written.

    Returns
    -------
    dict
        The run manifest including file checksums, metrics, and success rule.
    """
    start_time = time.time()
    out_path = ensure_dir(out_dir)

    # Deterministic runtime
    seed = int(config.runtime["random_seed"])
    set_runtime_environment(seed)
    null_draws = int(config.scoring["null_draws"])
    tolerance = float(config.runtime["verification_tolerance"])

    logger.info("Loading frozen assets...")

    # Load frozen assets
    sig_manifest = read_table(config.path("signature_panel"))
    cohort_manifest = read_table(config.path("cohort_manifest"))
    sig_rows = read_table(config.path("signature_rows"))
    confounder_sets = load_confounder_sets(config.path("confounder_panel"))

    freeze_dir = config.path("freeze_dir")
    cohort_data = _load_cohort_data(cohort_manifest, freeze_dir)
    universe = _build_gene_universe(cohort_data)
    cohort_ids = list(cohort_data.keys())

    logger.info(
        "Loaded %d signatures, %d cohorts, %d genes in universe",
        len(sig_manifest), len(cohort_data), len(universe),
    )

    # Build cohort overlap summary (before scoring loop)
    cohort_overlap_df = _build_cohort_overlap_summary(sig_manifest, sig_rows, cohort_data)

    # Main scoring loop: signature x cohort
    all_classification_records: list[dict[str, Any]] = []
    all_per_cohort_records: list[dict[str, Any]] = []
    all_loo_records: list[dict[str, Any]] = []
    all_platform_holdout_records: list[dict[str, Any]] = []
    all_null_summaries: list[dict[str, Any]] = []
    normalization_audit: dict[str, dict[str, Any]] = {}

    n_sigs = len(sig_manifest)
    for sig_idx, (_, srow) in enumerate(sig_manifest.iterrows()):
        sig_id = str(srow["signature_id"])
        logger.info(
            "Scoring signature %d/%d: %s", sig_idx + 1, n_sigs, sig_id,
        )
        raw_sig = sig_rows.loc[
            sig_rows["signature_id"] == sig_id,
            ["gene_symbol", "direction", "weight"],
        ].copy()

        # Normalization audit (once per signature)
        if sig_id not in normalization_audit:
            _, audit = normalize_signature(raw_sig)
            normalization_audit[sig_id] = audit

        # Score signature in all cohorts
        scoring_result = _score_signature_across_cohorts(
            sig_id, raw_sig, cohort_data, confounder_sets,
        )
        all_per_cohort_records.extend(scoring_result["per_cohort_records"])

        # Null model
        null_effects = _compute_null_effects(
            raw_sig, universe, null_draws, seed, cohort_data,
        )

        # Build aggregate profile
        profile = _build_profile(
            scoring_result["cohort_effects"],
            scoring_result["cohort_variances"],
            scoring_result["cohort_platforms"],
            scoring_result["cohort_confounder_maxes"],
            scoring_result["coverages"],
            null_effects,
        )

        # Null summary
        all_null_summaries.append(
            _build_null_summary(
                sig_id,
                profile["aggregate_effect"],
                null_effects,
                profile["null_separation_p"],
            )
        )

        # LOO table
        all_loo_records.extend(
            _build_loo_table(
                sig_id,
                scoring_result["cohort_effects"],
                scoring_result["cohort_variances"],
                cohort_ids,
            )
        )

        # Platform holdout table
        all_platform_holdout_records.extend(
            _build_platform_holdout_table(
                sig_id,
                scoring_result["cohort_effects"],
                scoring_result["cohort_variances"],
                scoring_result["cohort_platforms"],
            )
        )

        # Classification under all models
        all_classification_records.extend(
            _classify_all_models(sig_id, srow, profile)
        )

    # -----------------------------------------------------------------------
    # Build output DataFrames
    # -----------------------------------------------------------------------
    scores_df = pd.DataFrame(all_classification_records)
    per_cohort_df = pd.DataFrame(all_per_cohort_records)
    loo_df = pd.DataFrame(all_loo_records)
    platform_holdout_df = pd.DataFrame(all_platform_holdout_records)
    null_summary_df = pd.DataFrame(all_null_summaries)

    # -----------------------------------------------------------------------
    # Aggregate metrics + success rule
    # -----------------------------------------------------------------------
    aggregate = _aggregate_metrics(scores_df, tolerance)

    logger.info(
        "Success rule: %s (AUPRC margin %.4f, secondary wins %d)",
        "PASS" if aggregate["success_rule"]["success"] else "FAIL",
        aggregate["success_rule"]["auprc_margin"],
        aggregate["success_rule"]["secondary_wins"],
    )

    # -----------------------------------------------------------------------
    # Write all outputs
    # -----------------------------------------------------------------------

    # Core data tables
    write_table(scores_df, out_path / "aggregate_durability_scores.csv")
    write_table(per_cohort_df, out_path / "per_cohort_effects.csv")
    write_table(cohort_overlap_df, out_path / "cohort_overlap_summary.csv")
    write_table(null_summary_df, out_path / "matched_null_summary.csv")
    write_table(loo_df, out_path / "leave_one_cohort_out.csv")
    write_table(platform_holdout_df, out_path / "platform_holdout_summary.csv")

    # Normalization audit
    write_json(out_path / "normalization_audit.json", normalization_audit)

    # Verification manifest
    verification = {
        "aggregate_metrics": aggregate,
        "runtime": runtime_summary(),
        "seed": seed,
        "null_draws": null_draws,
        "signature_count": len(sig_manifest),
        "cohort_count": len(cohort_data),
        "gene_universe_size": len(universe),
        "timestamp": now_timestamp(),
    }
    write_json(out_path / "verification.json", verification)

    # Benchmark protocol (copy/record what was used)
    protocol_source = config.path("benchmark_protocol")
    if protocol_source.exists():
        protocol = read_json(protocol_source)
    else:
        protocol = {
            "note": "benchmark_protocol.json not found at source; recording config",
            "scoring": config.scoring,
            "runtime": config.runtime,
            "meta_analysis": config.meta_analysis,
        }
    write_json(out_path / "benchmark_protocol.json", protocol)

    # Public summary
    summary_text = _generate_public_summary(aggregate, scores_df, start_time)
    write_text(out_path / "public_summary.md", summary_text)

    # Certificates (optional - Task 9)
    cert_results = _try_generate_certificates(scores_df, aggregate, out_path)

    # Plots (optional - Task 10)
    plot_results = _try_generate_plots(
        scores_df, per_cohort_df, null_summary_df,
        loo_df, platform_holdout_df, out_path,
    )

    # -----------------------------------------------------------------------
    # Build manifest
    # -----------------------------------------------------------------------
    elapsed = time.time() - start_time
    output_files: dict[str, dict[str, Any]] = {}
    for filename in REQUIRED_OUTPUTS:
        if filename == "manifest.json":
            # manifest.json is written as the final step; mark it as present
            output_files[filename] = {
                "path": str(out_path / filename),
                "sha256": "written_last",
                "exists": True,
            }
            continue
        filepath = out_path / filename
        if filepath.exists():
            output_files[filename] = {
                "path": str(filepath),
                "sha256": sha256_file(filepath),
                "exists": True,
            }
        else:
            output_files[filename] = {
                "path": str(filepath),
                "sha256": None,
                "exists": False,
            }

    manifest = {
        "pipeline": "signature_durability_benchmark",
        "version": "0.1.0",
        "timestamp": now_timestamp(),
        "elapsed_seconds": round(elapsed, 2),
        "config": {
            "seed": seed,
            "null_draws": null_draws,
            "tolerance": tolerance,
            "signature_count": len(sig_manifest),
            "cohort_count": len(cohort_data),
            "gene_universe_size": len(universe),
        },
        "aggregate_metrics": aggregate,
        "output_files": output_files,
        "missing_outputs": [
            f for f, info in output_files.items() if not info["exists"]
        ],
        "runtime": runtime_summary(),
    }

    write_json(out_path / "manifest.json", manifest)
    logger.info(
        "Pipeline complete in %.1fs. %d/%d outputs written.",
        elapsed,
        sum(1 for v in output_files.values() if v["exists"]),
        len(REQUIRED_OUTPUTS),
    )

    return manifest
