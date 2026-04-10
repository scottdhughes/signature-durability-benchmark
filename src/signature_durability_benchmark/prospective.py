"""Shared prospective-round declaration/evaluation helpers."""

from __future__ import annotations

import csv
import gzip
import io
import json
import re
import subprocess
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

from .meta_analysis import direction_consistency, guarded_hksj_random_effects_meta
from .scoring import score_signature_in_cohort
from .utils import ensure_dir, now_timestamp, read_json, read_table, sha256_file, write_json, write_table, write_text


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DOWNLOAD_DIR = ROOT / "data" / "prospective_holdout" / "downloads"
SIGNATURES = ROOT / "data" / "freeze" / "signatures.tsv"
GENE_INFO = ROOT / "data" / "external_rnaseq" / "cache" / "Homo_sapiens.gene_info.gz"

IFN_QUARTET = [
    "hallmark_ifng_response",
    "hallmark_ifna_response",
    "schoggins_2011_irg",
    "blind_durable_ifn_composite",
]

CONTEXT_SIGNATURES = IFN_QUARTET + [
    "hallmark_inflammatory_response",
    "hallmark_tnfa_nfkb",
    "hallmark_e2f_targets",
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

_ANNOTATION_COLUMNS = {
    "",
    "Unnamed: 0",
    "GeneID",
    "Geneid",
    "gene_id",
    "ENSEMBL_gene_ID",
    "ensembl_gene_id",
    "gene_symbol",
    "Gene.symbol",
    "gene",
    "gene_name",
    "geneid",
    "Symbol",
    "symbol",
    "description",
    "Description",
    "entrezgene",
    "uniprotswissprot",
    "refseq_mrna",
}


@dataclass(frozen=True)
class GeneMaps:
    entrez_to_symbol: dict[str, str]
    ensembl_to_symbol: dict[str, str]


def _request(url: str) -> Request:
    return Request(url, headers={"User-Agent": "signature-durability-benchmark/1.0"})


def fetch_bytes(url: str) -> bytes:
    with urlopen(_request(url), timeout=120) as response:
        return response.read()


def fetch_url(url: str, destination: Path) -> dict[str, Any]:
    ensure_dir(destination.parent)
    fetched_this_run = False
    checked_at = now_timestamp()
    if not destination.exists():
        destination.write_bytes(fetch_bytes(url))
        fetched_this_run = True
    return {
        "url": url,
        "local_path": str(destination),
        "size_bytes": destination.stat().st_size,
        "sha256": sha256_file(destination),
        "fetched_this_run": fetched_this_run,
        "checked_at_utc": checked_at,
    }


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

    passthrough_keys = {
        "Sample_source_name_ch1",
        "Sample_geo_accession",
        "Sample_description",
        "Sample_title",
    }
    for key, values in rows:
        if key == "Sample_title":
            titles = values
        elif key == "Sample_characteristics_ch1":
            characteristic_rows.append(values)
        elif key in passthrough_keys:
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


def _split_declared_values(value: str) -> list[str]:
    return [part.strip() for part in str(value).split("|") if part.strip()]


def _normalize_name(value: str) -> str:
    text = str(value).strip().lower()
    text = text.replace(".bam", "").replace(".csv", "").replace(".txt", "")
    return re.sub(r"[^a-z0-9]+", "", text)


def _resolve_metadata_samples(meta: pd.DataFrame, sample_columns: list[str] | None = None) -> pd.DataFrame:
    meta = meta.copy()
    if sample_columns is None:
        for field in ("title", "description", "geo_accession"):
            if field in meta.columns and meta[field].astype(str).str.strip().ne("").all():
                meta["sample"] = meta[field].astype(str).str.replace(r" \(RNA-seq\)$", "", regex=True)
                return meta
        raise ValueError("Unable to infer sample identifiers from series-matrix metadata")

    normalized_columns = {_normalize_name(column): column for column in sample_columns}
    for field in ("description", "title", "geo_accession"):
        if field not in meta.columns:
            continue
        mapped: list[str] = []
        success = True
        for value in meta[field].astype(str).tolist():
            normalized = _normalize_name(value)
            column = normalized_columns.get(normalized)
            if column is None:
                candidates = [real for key, real in normalized_columns.items() if normalized and (normalized in key or key in normalized)]
                candidates = sorted(set(candidates))
                if len(candidates) == 1:
                    column = candidates[0]
            if column is None:
                success = False
                break
            mapped.append(column)
        if success and len(set(mapped)) == len(mapped):
            meta["sample"] = mapped
            return meta

    if len(meta) == len(sample_columns):
        meta["sample"] = sample_columns
        return meta
    raise ValueError("Unable to align series-matrix metadata to expression-matrix sample columns")


def _build_pheno(spec: pd.Series, meta: pd.DataFrame, sample_columns: list[str] | None = None) -> pd.DataFrame:
    meta = _resolve_metadata_samples(meta, sample_columns)
    phenotype_field = str(spec["phenotype_field"])
    if phenotype_field not in meta.columns:
        raise ValueError(f"{spec['external_cohort_id']}: phenotype field {phenotype_field} missing from series-matrix metadata")

    subset_field = str(spec.get("subset_field", "")).strip()
    subset_value = str(spec.get("subset_value", "")).strip()
    if subset_field:
        if subset_field not in meta.columns:
            raise ValueError(f"{spec['external_cohort_id']}: subset field {subset_field} missing from series-matrix metadata")
        subset_values = set(_split_declared_values(subset_value)) if subset_value else set()
        if subset_values:
            meta = meta.loc[meta[subset_field].astype(str).isin(subset_values)].copy()

    case_values = set(_split_declared_values(str(spec["case_label"])))
    control_values = set(_split_declared_values(str(spec["control_label"])))
    labels = meta[phenotype_field].astype(str)
    meta = meta.loc[labels.isin(case_values | control_values)].copy()
    meta["phenotype"] = labels.map(lambda value: "case" if value in case_values else "control")

    observed_case = int((meta["phenotype"] == "case").sum())
    observed_control = int((meta["phenotype"] == "control").sum())
    if observed_case != int(spec["expected_case_n"]):
        raise ValueError(f"{spec['external_cohort_id']}: expected {spec['expected_case_n']} case samples, found {observed_case}")
    if observed_control != int(spec["expected_control_n"]):
        raise ValueError(f"{spec['external_cohort_id']}: expected {spec['expected_control_n']} control samples, found {observed_control}")
    return meta[["sample", "phenotype"]].drop_duplicates()


def _detect_gene_symbol_column(frame: pd.DataFrame, gene_maps: GeneMaps) -> tuple[str, pd.Series]:
    for column in ("gene_symbol", "Symbol", "symbol", "gene_name", "Gene.symbol"):
        if column in frame.columns:
            return column, frame[column].astype(str).str.upper().str.strip()

    for column in ("ENSEMBL_gene_ID", "ensembl_gene_id", "Geneid", "GeneID", "entrezgene", "gene_id", "Unnamed: 0", ""):
        if column not in frame.columns:
            continue
        values = frame[column].astype(str).str.strip()
        sample = values.head(250)
        if sample.str.match(r"^ENSG\d+(\.\d+)?$").mean() > 0.8:
            mapped = values.str.split(".", n=1).str[0].map(gene_maps.ensembl_to_symbol)
            return column, mapped
        if sample.str.match(r"^\d+$").mean() > 0.8:
            mapped = values.map(gene_maps.entrez_to_symbol)
            return column, mapped
        return column, values.str.upper().str.strip()

    first = frame.columns[0]
    values = frame[first].astype(str).str.strip()
    sample = values.head(250)
    if sample.str.match(r"^ENSG\d+(\.\d+)?$").mean() > 0.8:
        mapped = values.str.split(".", n=1).str[0].map(gene_maps.ensembl_to_symbol)
        return first, mapped
    return first, values.str.upper().str.strip()


def _detect_sample_columns(frame: pd.DataFrame, gene_column: str) -> list[str]:
    sample_columns: list[str] = []
    for column in frame.columns:
        if column == gene_column or column in _ANNOTATION_COLUMNS:
            continue
        series = pd.to_numeric(frame[column].head(100), errors="coerce")
        if series.notna().mean() >= 0.7:
            sample_columns.append(column)
    if not sample_columns:
        raise ValueError("Unable to detect numeric sample columns in expression matrix")
    return sample_columns


def _finalize_symbol_matrix(
    frame: pd.DataFrame,
    sample_columns: list[str],
    *,
    log_transform: bool,
    duplicate_agg: str,
) -> pd.DataFrame:
    frame = frame.copy()
    frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper().str.strip()
    frame = frame.loc[(frame["gene_symbol"] != "") & (~frame["gene_symbol"].str.startswith("__"))].copy()
    frame = frame.dropna(subset=["gene_symbol"])
    for column in sample_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0.0)
    if duplicate_agg == "sum":
        expr = frame.groupby("gene_symbol", as_index=True)[sample_columns].sum()
    elif duplicate_agg == "mean":
        expr = frame.groupby("gene_symbol", as_index=True)[sample_columns].mean()
    else:
        raise ValueError(f"Unsupported duplicate_agg={duplicate_agg}")
    return np.log1p(expr) if log_transform else expr


def _load_direct_matrix(spec: pd.Series, series_matrix_path: Path, matrix_path: Path, gene_maps: GeneMaps) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = parse_series_matrix(series_matrix_path)
    frame = pd.read_csv(matrix_path, sep=None, engine="python", compression="infer")
    gene_column, gene_symbols = _detect_gene_symbol_column(frame, gene_maps)
    sample_columns = _detect_sample_columns(frame, gene_column)
    pheno = _build_pheno(spec, meta, sample_columns)
    duplicate_agg = "sum" if str(spec["acquisition_mode"]) == "direct_combined_matrix" else "mean"
    log_transform = str(spec["acquisition_mode"]) == "direct_combined_matrix"
    expr_frame = frame[[gene_column, *sample_columns]].copy()
    expr_frame["gene_symbol"] = gene_symbols
    expr = _finalize_symbol_matrix(
        expr_frame[["gene_symbol", *sample_columns]],
        sample_columns,
        log_transform=log_transform,
        duplicate_agg=duplicate_agg,
    )
    expr = expr.loc[:, pheno["sample"].tolist()]
    return expr, pheno


def _load_raw_tar_per_sample_counts(spec: pd.Series, series_matrix_path: Path, tar_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = parse_series_matrix(series_matrix_path)
    pheno = _build_pheno(spec, meta, None)
    wanted = {_normalize_name(sample): sample for sample in pheno["sample"].tolist()}
    series_list: list[pd.Series] = []
    name_re = re.compile(r"^GSM\d+_(.+?)_raw_count")

    with tarfile.open(tar_path, "r") as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.endswith((".txt.gz", ".tsv.gz", ".csv.gz")):
                continue
            basename = Path(member.name).name
            sample = None
            match = name_re.match(basename)
            if match:
                sample = wanted.get(_normalize_name(match.group(1)))
            norm_name = _normalize_name(basename)
            if sample is None:
                for candidate_norm, original in wanted.items():
                    if candidate_norm and (candidate_norm in norm_name or norm_name in candidate_norm):
                        sample = original
                        break
            if sample is None:
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            payload = gzip.decompress(extracted.read()).decode("utf-8", "ignore")
            frame = pd.read_csv(io.StringIO(payload), sep=None, engine="python")
            if frame.shape[1] < 2:
                continue
            gene_col = frame.columns[0]
            value_col = frame.columns[1]
            series = pd.to_numeric(frame[value_col], errors="coerce").fillna(0.0)
            series.index = frame[gene_col].astype(str).str.upper().str.strip()
            series.name = sample
            series_list.append(series)

    if not series_list:
        raise ValueError(f"{spec['external_cohort_id']}: no matching sample count files extracted from {tar_path}")

    expr = pd.concat(series_list, axis=1).fillna(0.0)
    expr.index.name = "gene_symbol"
    expr = np.log1p(expr.groupby(expr.index).sum())
    expr = expr.loc[:, pheno["sample"].tolist()]
    return expr, pheno


def prepare_cohort(
    spec: pd.Series,
    gene_maps: GeneMaps,
    download_audit: list[dict[str, Any]],
    *,
    download_dir: Path = DEFAULT_DOWNLOAD_DIR,
) -> dict[str, Any]:
    cohort_dir = ensure_dir(download_dir / str(spec["geo_accession"]))
    series_matrix_path = cohort_dir / Path(str(spec["series_matrix_url"])).name
    download_audit.append(fetch_url(str(spec["series_matrix_url"]), series_matrix_path))

    acquisition_mode = str(spec["acquisition_mode"])
    matrix_path = cohort_dir / Path(str(spec["matrix_url"])).name
    download_audit.append(fetch_url(str(spec["matrix_url"]), matrix_path))

    if acquisition_mode in {"direct_combined_matrix", "direct_normalized_matrix"}:
        expr, pheno = _load_direct_matrix(spec, series_matrix_path, matrix_path, gene_maps)
    elif acquisition_mode == "raw_tar_per_sample_counts":
        expr, pheno = _load_raw_tar_per_sample_counts(spec, series_matrix_path, matrix_path)
    else:
        raise ValueError(f"Unsupported acquisition_mode={acquisition_mode}")

    return {
        "cohort_id": str(spec["external_cohort_id"]),
        "geo_accession": str(spec["geo_accession"]),
        "tissue": str(spec["tissue"]),
        "expr": expr,
        "pheno": pheno,
    }


def git_provenance(root: Path = ROOT) -> dict[str, Any]:
    commit = ""
    dirty = False
    try:
        commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
        dirty = bool(subprocess.check_output(["git", "-C", str(root), "status", "--porcelain"], text=True).strip())
    except Exception:
        pass
    return {"git_commit": commit, "git_dirty": dirty}


def verify_round_receipt(
    registry_path: Path,
    protocol_path: Path,
    receipt_path: Path,
    *,
    expected_round_id: str | None = None,
) -> dict[str, Any]:
    receipt = read_json(receipt_path)
    registry = read_table(registry_path)
    round_ids = sorted({str(value) for value in registry["round_id"].astype(str).tolist() if str(value)})
    if len(round_ids) != 1:
        raise ValueError(f"Registry must contain exactly one round_id, found {round_ids}")
    round_id = round_ids[0]

    if expected_round_id and round_id != expected_round_id:
        raise ValueError(f"Expected round_id {expected_round_id}, found {round_id}")
    if str(receipt.get("round_id", "")) != round_id:
        raise ValueError(f"Receipt round_id {receipt.get('round_id')} does not match registry round_id {round_id}")
    if str(receipt.get("verification_status", "")) not in {"OK", ""}:
        raise ValueError(f"Receipt verification_status must be OK, found {receipt.get('verification_status')}")

    registry_sha256 = sha256_file(registry_path)
    protocol_sha256 = sha256_file(protocol_path)
    if str(receipt.get("registry_sha256", "")) != registry_sha256:
        raise ValueError("Registry SHA256 does not match the declaration receipt")
    if str(receipt.get("protocol_sha256", "")) != protocol_sha256:
        raise ValueError("Protocol SHA256 does not match the declaration receipt")
    return receipt


def _guard_declared_inputs(paths: list[Path]) -> dict[Path, str]:
    return {path: sha256_file(path) for path in paths if path.is_file()}


def _check_declared_inputs_unchanged(before: dict[Path, str]) -> None:
    after = {path: sha256_file(path) for path in before}
    if after != before:
        changed = [str(path) for path, digest in before.items() if after.get(path) != digest]
        raise RuntimeError(f"Declared prospective inputs changed during evaluation: {changed}")


def build_round_summary_markdown(
    registry: pd.DataFrame,
    evaluation: dict[str, Any],
    pooled_results: dict[str, dict[str, Any]],
    *,
    title: str = "Prospective Round Evaluation",
) -> str:
    lines = [
        f"# {title}",
        "",
        f"- Round: `{evaluation['round_id']}`",
        f"- Declared at: `{evaluation['declared_at_utc']}`",
        f"- Evaluation started: `{evaluation['evaluation_started_at_utc']}`",
        f"- Evaluation finished: `{evaluation['evaluation_finished_at_utc']}`",
        f"- Registry SHA256: `{evaluation['registry_sha256']}`",
    ]
    if evaluation.get("receipt_sha256"):
        lines.append(f"- Receipt SHA256: `{evaluation['receipt_sha256']}`")
        lines.append(f"- Receipt TSA reply time: `{evaluation['receipt_tsa_reply_time']}`")
    lines.extend(
        [
            f"- Overall success: `{evaluation['prediction_summary']['overall_success']}`",
            "",
            "## Cohort-Level IFN Hits",
            "",
            "| Cohort | Tissue | IFN sign hits | Rule satisfied |",
            "|---|---|---:|---:|",
        ]
    )
    for row in evaluation["prediction_summary"]["per_cohort_hits"]:
        lines.append(
            f"| {row['external_cohort_id']} | {row['tissue']} | {row['positive_hits']}/4 | {row['all_positive']} |"
        )
    lines.extend(
        [
            "",
            "## Pooled Primary Panel",
            "",
            "| Signature | pooled g | guarded p_Bonf,4 | sign consistency | success |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for signature_id in IFN_QUARTET:
        result = pooled_results[signature_id]
        lines.append(
            f"| {DISPLAY_NAMES[signature_id]} | {result['pooled_effect']:.3f} | "
            f"{result['guarded_p_bonf_4']:.4g} | {result['sign_consistency']:.3f} | "
            f"{evaluation['prediction_summary']['per_signature_success'][signature_id]} |"
        )
    return "\n".join(lines) + "\n"


def evaluate_round(
    registry_path: Path,
    protocol_path: Path,
    *,
    out_dir: Path,
    receipt_path: Path | None = None,
    signatures_path: Path = SIGNATURES,
    gene_info_path: Path = GENE_INFO,
    download_dir: Path = DEFAULT_DOWNLOAD_DIR,
) -> dict[str, Any]:
    evaluation_started_at_utc = now_timestamp()
    registry = read_table(registry_path)
    round_ids = sorted({str(value) for value in registry["round_id"].astype(str).tolist() if str(value)})
    if len(round_ids) != 1:
        raise ValueError(f"Registry must contain exactly one round_id, found {round_ids}")
    round_id = round_ids[0]
    declared_times = sorted({str(value) for value in registry["declared_at_utc"].astype(str).tolist() if str(value)})
    if len(declared_times) != 1:
        raise ValueError(f"Registry must contain exactly one declared_at_utc value, found {declared_times}")

    declared_paths = [registry_path, protocol_path]
    receipt = None
    if receipt_path is not None:
        receipt = verify_round_receipt(registry_path, protocol_path, receipt_path, expected_round_id=round_id)
        declared_paths.extend(path for path in receipt_path.parent.rglob("*") if path.is_file())
    declared_hashes = _guard_declared_inputs(declared_paths)

    ensure_dir(out_dir)
    out_csv = out_dir / "per_cohort_effects.csv"
    out_json = out_dir / "evaluation.json"
    out_md = out_dir / "evaluation_summary.md"

    signatures = read_table(signatures_path)
    signature_panels = {
        signature_id: signatures.loc[signatures["signature_id"] == signature_id, ["gene_symbol", "direction", "weight"]].copy()
        for signature_id in CONTEXT_SIGNATURES
    }

    gene_maps = load_gene_maps(gene_info_path)
    download_audit: list[dict[str, Any]] = []
    cohort_data = [prepare_cohort(spec, gene_maps, download_audit, download_dir=download_dir) for _, spec in registry.iterrows()]

    per_cohort_records: list[dict[str, Any]] = []
    for cohort in cohort_data:
        case_label = "case"
        control_label = "control"
        for signature_id in CONTEXT_SIGNATURES:
            result = score_signature_in_cohort(
                signature_panels[signature_id],
                cohort["expr"],
                cohort["pheno"],
                phenotype_column="phenotype",
                case_label=case_label,
                control_label=control_label,
            )
            per_cohort_records.append(
                {
                    "external_cohort_id": cohort["cohort_id"],
                    "geo_accession": cohort["geo_accession"],
                    "tissue": cohort["tissue"],
                    "signature_id": signature_id,
                    "signature_label": DISPLAY_NAMES[signature_id],
                    **result,
                }
            )

    per_cohort_df = pd.DataFrame.from_records(per_cohort_records)
    write_table(per_cohort_df, out_csv)

    pooled_results: dict[str, dict[str, Any]] = {}
    for signature_id in CONTEXT_SIGNATURES:
        subset = per_cohort_df.loc[per_cohort_df["signature_id"] == signature_id].copy()
        effects = subset["cohens_d"].astype(float).tolist()
        variances = subset["cohens_d_var"].astype(float).tolist()
        meta = guarded_hksj_random_effects_meta(effects, variances)
        meta["sign_consistency"] = direction_consistency(effects)
        meta["positive_hits"] = int(sum(effect > 0 for effect in effects))
        meta["guarded_p_bonf_4"] = min(1.0, float(meta["pooled_p"]) * len(IFN_QUARTET))
        meta["guarded_p_bonf_7"] = min(1.0, float(meta["pooled_p"]) * len(CONTEXT_SIGNATURES))
        meta["cohort_effects"] = {row["external_cohort_id"]: float(row["cohens_d"]) for _, row in subset.iterrows()}
        pooled_results[signature_id] = meta

    per_cohort_hits: list[dict[str, Any]] = []
    for _, spec in registry.iterrows():
        subset = per_cohort_df.loc[
            (per_cohort_df["external_cohort_id"] == str(spec["external_cohort_id"]))
            & (per_cohort_df["signature_id"].isin(IFN_QUARTET))
        ]
        positive_hits = int((subset["cohens_d"].astype(float) > 0).sum())
        per_cohort_hits.append(
            {
                "external_cohort_id": str(spec["external_cohort_id"]),
                "geo_accession": str(spec["geo_accession"]),
                "tissue": str(spec["tissue"]),
                "positive_hits": positive_hits,
                "all_positive": positive_hits == len(IFN_QUARTET),
            }
        )

    per_signature_success = {
        signature_id: bool(
            pooled_results[signature_id]["pooled_effect"] >= 0.8
            and pooled_results[signature_id]["guarded_p_bonf_4"] < 0.05
            and pooled_results[signature_id]["sign_consistency"] == 1.0
        )
        for signature_id in IFN_QUARTET
    }
    overall_success = bool(all(row["all_positive"] for row in per_cohort_hits) and all(per_signature_success.values()))

    evaluation_finished_at_utc = now_timestamp()
    summary = {
        "round_id": round_id,
        "declared_at_utc": declared_times[0],
        "evaluation_started_at_utc": evaluation_started_at_utc,
        "evaluation_finished_at_utc": evaluation_finished_at_utc,
        "registry_sha256": sha256_file(registry_path),
        "protocol_sha256": sha256_file(protocol_path),
        "receipt_sha256": sha256_file(receipt_path) if receipt_path and receipt_path.exists() else "",
        "receipt_tsa_reply_time": str(receipt.get("tsa_reply_time", "")) if receipt else "",
        "code_provenance": git_provenance(ROOT),
        "download_audit": download_audit,
        "prediction_summary": {
            "round_id": round_id,
            "declared_at_utc": declared_times[0],
            "overall_success": overall_success,
            "per_cohort_hits": per_cohort_hits,
            "per_signature_success": per_signature_success,
            "success_rule": {
                "per_cohort": "all four IFN signatures positive (g > 0) in each primary cohort",
                "pooled": "each IFN signature pooled g >= 0.8, guarded Bonferroni-4 p < 0.05, sign consistency = 1.0",
            },
        },
        "pooled_primary_panel": pooled_results,
        "per_cohort_effects_csv": str(out_csv),
    }
    write_json(out_json, summary)
    write_text(out_md, build_round_summary_markdown(registry, summary, pooled_results))

    _check_declared_inputs_unchanged(declared_hashes)
    return summary


def legacy_v1_summary_from_evaluation(
    evaluation: dict[str, Any],
    *,
    registry_path: Path,
    protocol_path: Path,
) -> dict[str, Any]:
    return {
        "registry": {
            "path": str(registry_path),
            "sha256": evaluation["registry_sha256"],
            "protocol_path": str(protocol_path),
            "protocol_sha256": evaluation["protocol_sha256"],
            "declared_at_utc": evaluation["declared_at_utc"],
            "n_primary_cohorts": int(len(evaluation["prediction_summary"]["per_cohort_hits"])),
        },
        "download_audit": evaluation["download_audit"],
        "prediction_summary": evaluation["prediction_summary"],
        "pooled_primary_panel": evaluation["pooled_primary_panel"],
        "per_cohort_effects_csv": evaluation["per_cohort_effects_csv"],
    }


def write_legacy_v1_outputs(
    evaluation: dict[str, Any],
    *,
    registry_path: Path,
    protocol_path: Path,
    out_json: Path,
    out_md: Path,
) -> None:
    registry = read_table(registry_path)
    pooled_results = evaluation["pooled_primary_panel"]
    legacy = legacy_v1_summary_from_evaluation(evaluation, registry_path=registry_path, protocol_path=protocol_path)
    write_json(out_json, legacy)
    write_text(
        out_md,
        build_round_summary_markdown(
            registry,
            evaluation,
            pooled_results,
            title="Prospective Holdout Prediction",
        ).replace(f"- Round: `{evaluation['round_id']}`\n", "- Registry: `prediction_registry_v1.tsv`\n", 1),
    )
