#!/usr/bin/env python3
"""Held-out external RNA-seq validation for the program-conditioned IFN diagnostic.

Primary panel:
  - bulk RNA-seq cohorts with explicit acute SARS-CoV-2 case/control labels
  - whole blood or upper-airway tissue
  - processed host-gene matrices available directly from GEO supplementary files

Exploratory stress test:
  - one PBMC cohort retained separately to show that cell-selected RNA-seq does not
    trivially validate every IFN-related signature.
"""

from __future__ import annotations

import csv
import gzip
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from signature_durability_benchmark.meta_analysis import guarded_hksj_random_effects_meta
from signature_durability_benchmark.scoring import score_signature_in_cohort
from signature_durability_benchmark.utils import write_json, write_text


ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "external_rnaseq"
MANIFEST = DATA_DIR / "cohort_manifest.tsv"
SIGNATURES = ROOT / "data" / "freeze" / "signatures.tsv"
GENE_INFO = DATA_DIR / "cache" / "Homo_sapiens.gene_info.gz"
OUT_DIR = ROOT / "outputs" / "canonical_v8"
OUT_JSON = OUT_DIR / "external_rnaseq_validation.json"
OUT_MD = OUT_DIR / "external_rnaseq_validation.md"

PRIMARY_SIGNATURES = [
    "hallmark_ifng_response",
    "hallmark_ifna_response",
    "schoggins_2011_irg",
    "blind_durable_ifn_composite",
    "hallmark_inflammatory_response",
    "hallmark_tnfa_nfkb",
    "hallmark_e2f_targets",
]

IFN_QUARTET = [
    "hallmark_ifng_response",
    "hallmark_ifna_response",
    "schoggins_2011_irg",
    "blind_durable_ifn_composite",
]

DISPLAY_NAMES = {
    "hallmark_ifng_response": "IFN-gamma core",
    "hallmark_ifna_response": "IFN-alpha core",
    "schoggins_2011_irg": "Schoggins 2011 IRG",
    "blind_durable_ifn_composite": "Blind IFN composite",
    "hallmark_inflammatory_response": "Inflammatory Response",
    "hallmark_tnfa_nfkb": "TNF-alpha/NF-kB",
    "hallmark_e2f_targets": "E2F Targets",
}


@dataclass(frozen=True)
class GeneMaps:
    entrez_to_symbol: dict[str, str]
    ensembl_to_symbol: dict[str, str]


def load_gene_maps(path: Path) -> GeneMaps:
    entrez_to_symbol: dict[str, str] = {}
    ensembl_to_symbol: dict[str, str] = {}
    with gzip.open(path, "rt") as handle:
        next(handle)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            entrez_id = parts[1]
            symbol = parts[2].upper()
            entrez_to_symbol[entrez_id] = symbol
            for item in parts[5].split("|"):
                if item.startswith("Ensembl:"):
                    ensembl_to_symbol[item.split(":", 1)[1].split(".", 1)[0]] = symbol
    return GeneMaps(entrez_to_symbol=entrez_to_symbol, ensembl_to_symbol=ensembl_to_symbol)


def parse_series_matrix(path: Path) -> pd.DataFrame:
    rows: list[tuple[str, list[str]]] = []
    with gzip.open(path, "rt", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            rows.append((row[0][1:], [value.strip('"') for value in row[1:]]))

    if not rows:
        raise ValueError(f"No !Sample_ rows found in {path}")

    sample_n = max(len(values) for _, values in rows)
    titles: list[str] = []
    characteristic_rows: list[list[str]] = []
    passthrough_rows: dict[str, list[str]] = {}

    for key, values in rows:
        if key == "Sample_title":
            titles = values
        elif key == "Sample_characteristics_ch1":
            characteristic_rows.append(values)
        elif key in {"Sample_source_name_ch1", "Sample_geo_accession"}:
            passthrough_rows[key[7:].lower()] = values

    records: list[dict[str, str]] = []
    for index in range(sample_n):
        record = {"title": titles[index] if index < len(titles) else f"sample_{index + 1}"}
        for key, values in passthrough_rows.items():
            record[key] = values[index] if index < len(values) else ""
        for values in characteristic_rows:
            raw = values[index] if index < len(values) else ""
            if ": " not in raw:
                continue
            key, value = raw.split(": ", 1)
            norm_key = key.strip().lower().replace(" ", "_").replace("-", "_")
            record[norm_key] = value.strip()
        records.append(record)
    return pd.DataFrame.from_records(records)


def _finalize_symbol_matrix(
    frame: pd.DataFrame,
    sample_columns: list[str],
    duplicate_agg: str,
) -> pd.DataFrame:
    frame = frame.copy()
    frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper().str.strip()
    frame = frame.loc[frame["gene_symbol"] != ""].copy()
    for column in sample_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0.0)
    if duplicate_agg == "sum":
        expr = frame.groupby("gene_symbol", as_index=True)[sample_columns].sum()
    elif duplicate_agg == "mean":
        expr = frame.groupby("gene_symbol", as_index=True)[sample_columns].mean()
    else:
        raise ValueError(f"Unsupported duplicate_agg={duplicate_agg}")
    return np.log1p(expr)


def load_wide_csv(spec: pd.Series, gene_maps: GeneMaps) -> tuple[pd.DataFrame, list[str]]:
    frame = pd.read_csv(DATA_DIR / str(spec["matrix_file"]), compression="gzip")
    gene_column = frame.columns[0]
    frame["gene_symbol"] = frame[gene_column].astype(str).map(gene_maps.entrez_to_symbol)
    sample_columns = [column for column in frame.columns if column not in {gene_column, "gene_symbol"}]
    frame = frame.dropna(subset=["gene_symbol"])
    return _finalize_symbol_matrix(frame[["gene_symbol", *sample_columns]], sample_columns, duplicate_agg="sum"), sample_columns


def load_wide_tsv(spec: pd.Series, gene_maps: GeneMaps) -> tuple[pd.DataFrame, list[str]]:
    frame = pd.read_csv(DATA_DIR / str(spec["matrix_file"]), sep="\t", compression="gzip")
    gene_column = frame.columns[0]
    if str(spec["gene_id_namespace"]) == "ensembl":
        frame["gene_symbol"] = (
            frame[gene_column].astype(str).str.split(".", n=1).str[0].map(gene_maps.ensembl_to_symbol)
        )
        frame = frame.dropna(subset=["gene_symbol"])
    else:
        frame["gene_symbol"] = frame[gene_column].astype(str)
    sample_columns = [column for column in frame.columns if column not in {gene_column, "gene_symbol"}]
    return _finalize_symbol_matrix(frame[["gene_symbol", *sample_columns]], sample_columns, duplicate_agg="sum"), sample_columns


def load_three_row_preamble_tsv(spec: pd.Series, _: GeneMaps) -> tuple[pd.DataFrame, list[str]]:
    with gzip.open(DATA_DIR / str(spec["matrix_file"]), "rt", errors="ignore", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header_row = next(reader)
        next(reader)
        next(reader)
        rows = list(reader)
    sample_columns = header_row[1:]
    frame = pd.DataFrame(rows, columns=["gene_symbol", *sample_columns])
    return _finalize_symbol_matrix(frame, sample_columns, duplicate_agg="mean"), sample_columns


def load_space_delimited_counts(spec: pd.Series, _: GeneMaps) -> tuple[pd.DataFrame, list[str]]:
    matrix_path = DATA_DIR / str(spec["matrix_file"])
    with gzip.open(matrix_path, "rt", errors="ignore") as handle:
        sample_columns = handle.readline().strip().split()
    frame = pd.read_csv(matrix_path, sep=r"\s+", compression="gzip", skiprows=1, header=None)
    frame.columns = ["gene_symbol", *sample_columns]
    return _finalize_symbol_matrix(frame, sample_columns, duplicate_agg="sum"), sample_columns


def load_long_gene_name_tsv(spec: pd.Series, _: GeneMaps) -> tuple[pd.DataFrame, list[str]]:
    frame = pd.read_csv(
        DATA_DIR / str(spec["matrix_file"]),
        sep="\t",
        compression="gzip",
        usecols=["RecordID", "gene_name", "raw_count"],
    )
    expr = frame.pivot_table(index="gene_name", columns="RecordID", values="raw_count", aggfunc="sum").fillna(0.0)
    expr.index = expr.index.astype(str).str.upper().str.strip()
    return np.log1p(expr), expr.columns.tolist()


LOADERS: dict[str, Callable[[pd.Series, GeneMaps], tuple[pd.DataFrame, list[str]]]] = {
    "wide_csv": load_wide_csv,
    "wide_tsv": load_wide_tsv,
    "wide_tsv_three_row_preamble": load_three_row_preamble_tsv,
    "space_delimited_counts": load_space_delimited_counts,
    "long_tsv_gene_name": load_long_gene_name_tsv,
}


def prepare_cohort(spec: pd.Series, gene_maps: GeneMaps) -> dict[str, object]:
    loader = LOADERS[str(spec["matrix_format"])]
    expr, sample_columns = loader(spec, gene_maps)
    meta = parse_series_matrix(DATA_DIR / str(spec["series_matrix_file"]))

    if str(spec["sample_id_source"]) == "sample_title":
        meta = meta.copy()
        meta["sample"] = meta["title"]
    elif str(spec["sample_id_source"]) == "matrix_header_row0":
        if len(meta) != len(sample_columns):
            raise ValueError(
                f"{spec['external_cohort_id']}: matrix header has {len(sample_columns)} samples but series matrix has {len(meta)} entries"
            )
        meta = meta.copy()
        meta["sample"] = sample_columns
    else:
        raise ValueError(f"Unsupported sample_id_source={spec['sample_id_source']}")

    phenotype_field = str(spec["phenotype_field"])
    case_label = str(spec["case_label"])
    control_label = str(spec["control_label"])
    if phenotype_field not in meta.columns:
        raise ValueError(f"{spec['external_cohort_id']}: phenotype field {phenotype_field} missing from series matrix metadata")

    meta = meta.loc[meta[phenotype_field].isin([case_label, control_label])].copy()
    meta["phenotype"] = meta[phenotype_field]
    missing = [sample for sample in meta["sample"].tolist() if sample not in expr.columns]
    if missing:
        raise ValueError(f"{spec['external_cohort_id']}: {len(missing)} metadata samples missing from expression matrix")
    expr = expr.loc[:, meta["sample"].tolist()]

    observed_total = int(meta.shape[0])
    observed_case = int((meta["phenotype"] == case_label).sum())
    observed_control = int((meta["phenotype"] == control_label).sum())
    for field, observed in [
        ("total_samples", observed_total),
        ("case_samples", observed_case),
        ("control_samples", observed_control),
    ]:
        expected = int(spec[field])
        if expected != observed:
            raise ValueError(
                f"{spec['external_cohort_id']}: expected {field}={expected} but observed {observed}"
            )

    return {
        "spec": spec,
        "expr": expr,
        "pheno": meta[["sample", "phenotype"]],
        "total_samples": observed_total,
        "case_samples": observed_case,
        "control_samples": observed_control,
    }


def bonferroni(p_value: float, n_tests: int) -> float:
    return min(1.0, float(p_value) * n_tests)


def main() -> None:
    manifest = pd.read_csv(MANIFEST, sep="\t")
    signatures = pd.read_csv(SIGNATURES, sep="\t")
    gene_maps = load_gene_maps(GENE_INFO)

    cohorts: dict[str, dict[str, object]] = {}
    for _, spec in manifest.iterrows():
        cohorts[str(spec["external_cohort_id"])] = prepare_cohort(spec, gene_maps)

    primary_ids = manifest.loc[manifest["analysis_group"] == "primary", "external_cohort_id"].astype(str).tolist()
    exploratory_ids = manifest.loc[manifest["analysis_group"] == "exploratory", "external_cohort_id"].astype(str).tolist()

    cohort_signature_results: dict[str, dict[str, dict[str, float | int | bool]]] = {}
    for cohort_id, cohort in cohorts.items():
        spec = cohort["spec"]
        expr = cohort["expr"]
        pheno = cohort["pheno"]
        case_label = str(spec["case_label"])
        control_label = str(spec["control_label"])
        cohort_signature_results[cohort_id] = {}
        for signature_id in PRIMARY_SIGNATURES:
            sig_df = signatures.loc[signatures["signature_id"] == signature_id, ["gene_symbol", "direction", "weight"]].copy()
            result = score_signature_in_cohort(sig_df, expr, pheno, "phenotype", case_label, control_label)
            shared_gene_count = int(sig_df["gene_symbol"].astype(str).str.upper().isin(expr.index).sum())
            cohort_signature_results[cohort_id][signature_id] = {
                "cohens_d": float(result["cohens_d"]),
                "cohens_d_var": float(result["cohens_d_var"]),
                "coverage_fraction": float(result["coverage_fraction"]),
                "case_n": int(result["case_n"]),
                "control_n": int(result["control_n"]),
                "shared_gene_count": shared_gene_count,
                "direction_consistent": bool(result["direction_consistent"]),
            }

    pooled_primary: dict[str, dict[str, float | int]] = {}
    for signature_id in PRIMARY_SIGNATURES:
        effects = [cohort_signature_results[cohort_id][signature_id]["cohens_d"] for cohort_id in primary_ids]
        variances = [cohort_signature_results[cohort_id][signature_id]["cohens_d_var"] for cohort_id in primary_ids]
        meta = guarded_hksj_random_effects_meta(effects, variances)
        sign_consistency = float(sum(effect > 0 for effect in effects) / len(effects))
        pooled_primary[signature_id] = {
            "pooled_g": float(meta["pooled_effect"]),
            "guarded_p": float(meta["pooled_p"]),
            "guarded_p_bonf_7": bonferroni(float(meta["pooled_p"]), len(PRIMARY_SIGNATURES)),
            "i_squared": float(meta["i_squared"]),
            "tau2": float(meta["tau2"]),
            "Q": float(meta["Q"]),
            "k": int(meta["k"]),
            "sign_consistency": sign_consistency,
            "mean_coverage": float(np.mean([cohort_signature_results[cohort_id][signature_id]["coverage_fraction"] for cohort_id in primary_ids])),
        }

    primary_panel = {
        "description": (
            "Primary held-out bulk RNA-seq validation panel. These cohorts were chosen to test cross-platform "
            "transfer of the IFN diagnostic on bulk whole-blood or upper-airway RNA-seq datasets with explicit "
            "acute SARS-CoV-2 case/control labels and processed matrices deposited directly in GEO."
        ),
        "selection_rules": [
            "Bulk RNA-seq supplementary matrix deposited directly in GEO",
            "Explicit two-group acute SARS-CoV-2 case/control labels in the series matrix",
            "Whole blood, leukocyte-rich whole blood, or upper-airway bulk tissue",
            "At least 25 samples after filtering to the primary two-group comparison",
        ],
        "n_cohorts": len(primary_ids),
        "n_samples": int(sum(int(cohorts[cohort_id]["total_samples"]) for cohort_id in primary_ids)),
        "cohorts": {},
        "pooled_signatures": pooled_primary,
        "ifn_focus_summary": {
            "all_four_positive_in_all_primary_cohorts": all(
                pooled_primary[signature_id]["sign_consistency"] == 1.0 for signature_id in IFN_QUARTET
            ),
            "all_four_guarded_p_bonf_7": {
                signature_id: pooled_primary[signature_id]["guarded_p_bonf_7"] for signature_id in IFN_QUARTET
            },
            "all_four_pooled_g": {signature_id: pooled_primary[signature_id]["pooled_g"] for signature_id in IFN_QUARTET},
        },
    }

    for cohort_id in primary_ids:
        spec = cohorts[cohort_id]["spec"]
        primary_panel["cohorts"][cohort_id] = {
            "geo_accession": str(spec["geo_accession"]),
            "tissue": str(spec["tissue"]),
            "design_summary": str(spec["design_summary"]),
            "total_samples": int(cohorts[cohort_id]["total_samples"]),
            "case_samples": int(cohorts[cohort_id]["case_samples"]),
            "control_samples": int(cohorts[cohort_id]["control_samples"]),
            "signatures": cohort_signature_results[cohort_id],
        }

    exploratory_panel = {
        "description": (
            "Exploratory cell-selected PBMC stress test. Retained separately from the primary bulk RNA-seq pool "
            "to show that the external extension does not trivially validate every IFN-related comparison."
        ),
        "cohorts": {},
    }
    for cohort_id in exploratory_ids:
        spec = cohorts[cohort_id]["spec"]
        exploratory_panel["cohorts"][cohort_id] = {
            "geo_accession": str(spec["geo_accession"]),
            "tissue": str(spec["tissue"]),
            "design_summary": str(spec["design_summary"]),
            "total_samples": int(cohorts[cohort_id]["total_samples"]),
            "case_samples": int(cohorts[cohort_id]["case_samples"]),
            "control_samples": int(cohorts[cohort_id]["control_samples"]),
            "signatures": cohort_signature_results[cohort_id],
            "notes": str(spec["notes"]),
        }

    payload = {
        "description": (
            "Held-out RNA-seq external validation for the program-conditioned IFN diagnostic, split into a "
            "primary bulk RNA-seq panel and an exploratory PBMC stress test."
        ),
        "primary_bulk_rnaseq_panel": primary_panel,
        "exploratory_pbmc_stress_test": exploratory_panel,
    }

    summary_lines = [
        "# External RNA-seq Validation",
        "",
        f"- Primary bulk RNA-seq panel: `{primary_panel['n_cohorts']}` cohorts, `{primary_panel['n_samples']}` samples",
        "- Primary cohorts: `GSE152641`, `GSE171110`, `GSE152075`, `GSE167000`",
        "",
        "## IFN Quartet",
        "",
    ]
    for signature_id in IFN_QUARTET:
        pooled = pooled_primary[signature_id]
        summary_lines.append(
            f"- `{DISPLAY_NAMES[signature_id]}`: pooled g=`{pooled['pooled_g']:.3f}`, "
            f"guarded p=`{pooled['guarded_p']:.4g}`, Bonferroni-7 p=`{pooled['guarded_p_bonf_7']:.4g}`, "
            f"I²=`{pooled['i_squared']:.3f}`, sign consistency=`{pooled['sign_consistency']:.3f}`"
        )

    summary_lines.extend([
        "",
        "## Comparator Signatures",
        "",
    ])
    for signature_id in ["hallmark_inflammatory_response", "hallmark_tnfa_nfkb", "hallmark_e2f_targets"]:
        pooled = pooled_primary[signature_id]
        summary_lines.append(
            f"- `{DISPLAY_NAMES[signature_id]}`: pooled g=`{pooled['pooled_g']:.3f}`, "
            f"guarded p=`{pooled['guarded_p']:.4g}`, Bonferroni-7 p=`{pooled['guarded_p_bonf_7']:.4g}`, "
            f"I²=`{pooled['i_squared']:.3f}`"
        )

    if exploratory_ids:
        summary_lines.extend([
            "",
            "## Exploratory PBMC Stress Test",
            "",
        ])
        for cohort_id in exploratory_ids:
            cohort = exploratory_panel["cohorts"][cohort_id]
            summary_lines.append(
                f"- `{cohort_id}` ({cohort['geo_accession']}): "
                f"IFN-gamma g=`{cohort['signatures']['hallmark_ifng_response']['cohens_d']:.3f}`, "
                f"IFN-alpha g=`{cohort['signatures']['hallmark_ifna_response']['cohens_d']:.3f}`, "
                f"Schoggins g=`{cohort['signatures']['schoggins_2011_irg']['cohens_d']:.3f}`, "
                f"blind IFN g=`{cohort['signatures']['blind_durable_ifn_composite']['cohens_d']:.3f}`"
            )

    write_json(OUT_JSON, payload)
    write_text(OUT_MD, "\n".join(summary_lines) + "\n")

    print("=" * 80)
    print("HELD-OUT EXTERNAL RNA-SEQ VALIDATION")
    print("=" * 80)
    for signature_id in IFN_QUARTET:
        pooled = pooled_primary[signature_id]
        print(
            f"{signature_id:28s} pooled_g={pooled['pooled_g']:.3f} "
            f"p={pooled['guarded_p']:.4g} p_bonf7={pooled['guarded_p_bonf_7']:.4g} "
            f"I2={pooled['i_squared']:.3f}"
        )
    print(f"\nSaved to {OUT_JSON}")
    print(f"Saved to {OUT_MD}")


if __name__ == "__main__":
    main()
