#!/usr/bin/env python3
"""
Download and freeze REAL GEO expression cohorts.

Downloads series matrix files from NCBI GEO FTP, parses expression matrices
and phenotype metadata, maps probes to gene symbols, and saves frozen
gene-by-sample expression matrices and phenotype files.

Re-runnable and idempotent: skips already-downloaded/frozen cohorts.
"""

import gzip
import hashlib
import json
import logging
import re
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = PROJECT_ROOT / "data" / "raw"
FREEZE_MATRICES = PROJECT_ROOT / "data" / "freeze" / "cohort_matrices"
FREEZE_PHENO = PROJECT_ROOT / "data" / "freeze" / "cohort_phenotypes"
CONFIG_DIR = PROJECT_ROOT / "config"

for d in [DATA_RAW, FREEZE_MATRICES, FREEZE_PHENO]:
    d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Cohort registry — verified GEO datasets with correct platform IDs
# ---------------------------------------------------------------------------
COHORTS = [
    {
        "cohort_id": "sepsis_mortality_gse65682",
        "accession": "GSE65682",
        "platform_id": "GPL13667",
        "platform_type": "affymetrix",
        "tissue": "whole_blood",
        "biological_program": "inflammation",
        "pheno_key": "mortality_event_28days",
        "case_label": "died",
        "control_label": "survived",
        "notes": "Scicluna et al. ICU sepsis, 28-day mortality outcome, ~480 annotated",
    },
    {
        "cohort_id": "breast_prognosis_gse2034",
        "accession": "GSE2034",
        "platform_id": "GPL96",
        "platform_type": "affymetrix",
        "tissue": "breast_tumor",
        "biological_program": "proliferation",
        "pheno_key": "bone relapses",
        "case_label": "relapse",
        "control_label": "no_relapse",
        "notes": "Wang et al. breast cancer bone relapse, 286 samples",
    },
    {
        "cohort_id": "breast_er_gse7390",
        "accession": "GSE7390",
        "platform_id": "GPL96",
        "platform_type": "affymetrix",
        "tissue": "breast_tumor",
        "biological_program": "proliferation",
        "pheno_key": "er",
        "case_label": "er_positive",
        "control_label": "er_negative",
        "notes": "Desmedt et al. untreated breast cancer ER status, 198 samples",
    },
    {
        "cohort_id": "hypoxia_cellline_gse53012",
        "accession": "GSE53012",
        "platform_id": "GPL570",
        "platform_type": "affymetrix",
        "tissue": "cell_line",
        "biological_program": "hypoxia",
        "pheno_key": "source",
        "case_label": "hypoxia",
        "control_label": "normoxia",
        "notes": "Cancer cell lines under hypoxia vs normoxia, 27 samples",
    },
    {
        "cohort_id": "sepsis_blood_gse28750",
        "accession": "GSE28750",
        "platform_id": "GPL570",
        "platform_type": "affymetrix",
        "tissue": "whole_blood",
        "biological_program": "inflammation",
        "pheno_key": "health status",
        "case_label": "sepsis",
        "control_label": "healthy",
        "notes": "Sepsis vs healthy whole blood, 30 samples (10+20)",
    },
    {
        "cohort_id": "tb_blood_gse19491",
        "accession": "GSE19491",
        "platform_id": "GPL6947",
        "platform_type": "illumina",
        "tissue": "whole_blood",
        "biological_program": "interferon",
        "pheno_key": "illness",
        "case_label": "active_tb",
        "control_label": "healthy",
        "notes": "Berry et al. TB vs healthy blood, ~155 case+ctrl from 498 total",
    },
    {
        "cohort_id": "ipf_lung_gse47460",
        "accession": "GSE47460",
        "platform_id": "GPL14550",
        "platform_type": "agilent",
        "tissue": "lung",
        "biological_program": "emt",
        "matrix_url_override": True,
        "pheno_key": "disease state",
        "case_label": "ild",
        "control_label": "control",
        "notes": "ILD (mostly IPF) vs control lung, Agilent, 429 samples",
    },
    {
        "cohort_id": "copd_lung_gse47460",
        "accession": "GSE47460",
        "platform_id": "GPL6480",
        "platform_type": "agilent",
        "tissue": "lung",
        "biological_program": "inflammation",
        "matrix_url_override": True,
        "pheno_key": "disease state",
        "case_label": "copd",
        "control_label": "control",
        "notes": "COPD vs control lung, Agilent GPL6480, 153 samples",
    },
    {
        "cohort_id": "influenza_pbmc_gse101702",
        "accession": "GSE101702",
        "platform_id": "GPL21185",
        "platform_type": "agilent",
        "tissue": "pbmc",
        "biological_program": "interferon",
        "pheno_key": "diagnosis",
        "case_label": "influenza",
        "control_label": "healthy",
        "notes": "Influenza vs healthy PBMC, 159 samples",
    },
    {
        "cohort_id": "trauma_blood_gse36809",
        "accession": "GSE36809",
        "platform_id": "GPL570",
        "platform_type": "affymetrix",
        "tissue": "whole_blood",
        "biological_program": "inflammation",
        "pheno_key": "title",
        "case_label": "trauma",
        "control_label": "control",
        "notes": "Trauma/burns patients vs healthy controls, 857 samples",
    },
    {
        "cohort_id": "hcc_liver_gse6764",
        "accession": "GSE6764",
        "platform_id": "GPL570",
        "platform_type": "affymetrix",
        "tissue": "liver",
        "biological_program": "proliferation",
        "pheno_key": "source",
        "case_label": "hcc",
        "control_label": "normal",
        "notes": "HCC vs normal liver, 75 samples",
    },
    {
        "cohort_id": "crohn_intestine_gse112366",
        "accession": "GSE112366",
        "platform_id": "GPL13158",
        "platform_type": "affymetrix",
        "tissue": "intestine",
        "biological_program": "inflammation",
        "pheno_key": "diagnosis",
        "case_label": "crohn",
        "control_label": "control",
        "notes": "Crohn's disease vs control intestinal biopsies, 388 samples",
    },
]


# ---------------------------------------------------------------------------
# Utility: download with retries
# ---------------------------------------------------------------------------
def download_url(url: str, dest: Path, retries: int = 3, timeout: int = 300) -> bool:
    """Download a URL to a local file with retries."""
    if dest.exists() and dest.stat().st_size > 0:
        log.info(f"  Already cached: {dest.name}")
        return True
    for attempt in range(1, retries + 1):
        try:
            log.info(f"  Downloading (attempt {attempt}/{retries})")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            dest.write_bytes(data)
            log.info(f"  Downloaded {len(data):,} bytes -> {dest.name}")
            return True
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
            log.warning(f"  Attempt {attempt} failed: {e}")
            if attempt < retries:
                time.sleep(5 * attempt)
    return False


def download_text(url: str, timeout: int = 300) -> str | None:
    """Download a URL and return text."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except Exception as e:
        log.warning(f"  Text download failed: {e}")
        return None


# ---------------------------------------------------------------------------
# URL builders
# ---------------------------------------------------------------------------
def series_matrix_url(accession: str) -> str:
    prefix = accession[:-3] + "nnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}"
        f"/matrix/{accession}_series_matrix.txt.gz"
    )


def series_matrix_platform_url(accession: str, platform_id: str) -> str:
    prefix = accession[:-3] + "nnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}"
        f"/matrix/{accession}-{platform_id}_series_matrix.txt.gz"
    )


def platform_annot_url(platform_id: str) -> str:
    prefix = platform_id[:-3] + "nnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{platform_id}"
        f"/annot/{platform_id}.annot.gz"
    )


def platform_web_url(platform_id: str) -> str:
    return (
        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        f"?acc={platform_id}&targ=self&form=text&view=data"
    )


# ---------------------------------------------------------------------------
# Series matrix parser
# ---------------------------------------------------------------------------
def parse_series_matrix(gz_path: Path):
    """Parse a GEO series matrix .gz file.
    Returns (expr_df, pheno_dict, sample_ids) or (None, None, None)."""
    with gzip.open(gz_path, "rt", errors="replace") as f:
        lines = f.readlines()

    pheno: dict[str, list[list[str]]] = {}
    expr_lines: list[str] = []
    in_matrix = False
    header_line = None

    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("!Series_") or line.startswith("!series_"):
            continue
        if line.startswith("!Sample_") or line.startswith("!sample_"):
            key = line.split("\t")[0]
            vals = [v.strip('"') for v in line.split("\t")[1:]]
            pheno.setdefault(key, []).append(vals)
            continue
        if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
            in_matrix = True
            header_line = line
            continue
        if in_matrix:
            if line.startswith("!") or line == "":
                break
            expr_lines.append(line)

    if header_line is None or not expr_lines:
        return None, None, None

    sample_ids = [h.strip('"') for h in header_line.split("\t")[1:]]

    probe_ids = []
    values = []
    for el in expr_lines:
        parts = el.split("\t")
        probe_ids.append(parts[0].strip('"'))
        row = []
        for v in parts[1:]:
            v = v.strip('"').strip()
            if v in ("", "null", "nan", "NULL", "NaN"):
                row.append(np.nan)
            else:
                try:
                    row.append(float(v))
                except ValueError:
                    row.append(np.nan)
        values.append(row)

    expr_df = pd.DataFrame(values, index=probe_ids, columns=sample_ids)
    expr_df.index.name = "probe_id"
    return expr_df, pheno, sample_ids


# ---------------------------------------------------------------------------
# Platform annotation
# ---------------------------------------------------------------------------
_platform_cache: dict[str, dict[str, str]] = {}


def get_probe_to_gene(platform_id: str) -> dict[str, str]:
    """Get probe->gene mapping. Tries .annot.gz then GEO web API."""
    if platform_id in _platform_cache:
        return _platform_cache[platform_id]

    # Try .annot.gz
    annot_gz = DATA_RAW / f"{platform_id}.annot.gz"
    if download_url(platform_annot_url(platform_id), annot_gz):
        mapping = _parse_annot_gz(annot_gz)
        if mapping:
            log.info(f"  {len(mapping)} probe->gene from .annot.gz")
            _platform_cache[platform_id] = mapping
            return mapping

    # Fallback: GEO web API
    log.info(f"  Using GEO web API for {platform_id}")
    mapping = _parse_platform_web(platform_id)
    if mapping:
        log.info(f"  {len(mapping)} probe->gene from web API")
        _platform_cache[platform_id] = mapping
        return mapping

    log.warning(f"  No probe->gene mapping for {platform_id}")
    return {}


def _parse_annot_gz(gz_path: Path) -> dict[str, str]:
    """Parse .annot.gz for probe->gene."""
    probe_to_gene = {}
    with gzip.open(gz_path, "rt", errors="replace") as f:
        in_table = False
        header = None
        gene_idx = None
        probe_idx = 0
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table or not line.strip():
                continue
            parts = line.split("\t")
            if header is None:
                header = [h.strip('"').strip().upper() for h in parts]
                for i, h in enumerate(header):
                    if h in ("ID", "PROBE_ID"):
                        probe_idx = i
                # Prefer exact "GENE SYMBOL" over partial matches like "UNIGENE SYMBOL"
                for i, h in enumerate(header):
                    if h == "GENE SYMBOL" or h == "GENE_SYMBOL":
                        gene_idx = i
                        break
                if gene_idx is None:
                    # Fallback to partial match but exclude UNIGENE
                    for i, h in enumerate(header):
                        if "GENE" in h and "SYMBOL" in h and "UNIGENE" not in h:
                            gene_idx = i
                            break
                if gene_idx is None:
                    return {}
                continue
            if gene_idx is None or len(parts) <= max(probe_idx, gene_idx):
                continue
            pid = parts[probe_idx].strip('"').strip()
            raw = parts[gene_idx].strip('"').strip()
            if not raw or raw == "---":
                continue
            sym = re.split(r"\s*///?\s*", raw)[0].strip()
            if sym and sym != "---":
                probe_to_gene[pid] = sym.upper()
    return probe_to_gene


def _parse_platform_web(platform_id: str) -> dict[str, str]:
    """Download and parse platform annotation from GEO web API."""
    text = download_text(platform_web_url(platform_id))
    if not text:
        return {}

    probe_to_gene = {}
    in_table = False
    header = None
    gene_idx = None
    probe_idx = 0

    for line in text.split("\n"):
        line = line.rstrip("\r\n")
        if "!platform_table_begin" in line:
            in_table = True
            continue
        if "!platform_table_end" in line:
            break
        if not in_table or not line.strip():
            continue
        parts = line.split("\t")
        if header is None:
            header = [h.strip('"').strip().upper() for h in parts]
            for i, h in enumerate(header):
                if h in ("ID", "PROBE_ID"):
                    probe_idx = i
            # Prefer exact "GENE SYMBOL" / "GENE_SYMBOL"
            for i, h in enumerate(header):
                if h in ("GENE SYMBOL", "GENE_SYMBOL"):
                    gene_idx = i
                    break
            if gene_idx is None:
                for i, h in enumerate(header):
                    if ("GENE" in h and "SYMBOL" in h and "UNIGENE" not in h) or h == "SYMBOL":
                        gene_idx = i
                        break
            if gene_idx is None:
                log.warning(f"  No Gene Symbol column. Cols: {header[:10]}")
                return {}
            continue
        if gene_idx is None or len(parts) <= max(probe_idx, gene_idx):
            continue
        pid = parts[probe_idx].strip('"').strip()
        raw = parts[gene_idx].strip('"').strip()
        if not raw or raw in ("---", ""):
            continue
        sym = re.split(r"\s*///?\s*", raw)[0].strip()
        if " // " in sym:
            sym = sym.split(" // ")[0].strip()
        if sym and sym != "---":
            probe_to_gene[pid] = sym.upper()
    return probe_to_gene


# ---------------------------------------------------------------------------
# Phenotype extraction — per-cohort logic
# ---------------------------------------------------------------------------
def extract_phenotype(pheno: dict, cohort: dict, sample_ids: list[str]) -> pd.DataFrame:
    """Extract case/control labels from parsed GEO metadata."""
    cid = cohort["cohort_id"]
    pk = cohort["pheno_key"].lower()

    char_rows = pheno.get("!Sample_characteristics_ch1", [])
    source_rows = pheno.get("!Sample_source_name_ch1", [])
    title_rows = pheno.get("!Sample_title", [])
    n = len(sample_ids)

    labels = [None] * n

    # --- Special: title-based classification ---
    if pk == "title":
        for row_vals in title_rows:
            for i, val in enumerate(row_vals):
                if i < n:
                    labels[i] = val.strip().lower()

    # --- Special: source-based classification ---
    elif pk == "source":
        for row_vals in source_rows:
            for i, val in enumerate(row_vals):
                if i < n:
                    labels[i] = val.strip().lower()

    # --- Standard: characteristics-based ---
    else:
        for row_vals in char_rows:
            for i, val in enumerate(row_vals):
                if i >= n:
                    break
                val_lower = val.lower()
                # Match the key portion (before colon) to avoid false substring matches
                if ":" in val_lower:
                    key_part = val_lower.split(":", 1)[0].strip()
                    if pk == key_part or key_part.startswith(pk):
                        labels[i] = val_lower.split(":", 1)[1].strip()
                elif pk in val_lower:
                    labels[i] = val_lower

    # Classify
    phenotypes = []
    raw_labels_out = []
    for raw in labels:
        if raw is None:
            phenotypes.append("unknown")
            raw_labels_out.append("")
        else:
            raw_labels_out.append(raw)
            phenotypes.append(_classify(raw, cohort))

    return pd.DataFrame({
        "sample_id": sample_ids, "phenotype": phenotypes, "raw_label": raw_labels_out
    })


def _classify(raw: str, cohort: dict) -> str:
    """Map raw label to case/control/other for each cohort."""
    cid = cohort["cohort_id"]
    r = raw.strip().lower()

    # --- Sepsis mortality GSE65682 ---
    if cid == "sepsis_mortality_gse65682":
        if r == "1":
            return "died"
        if r == "0":
            return "survived"
        return "other"

    # --- Breast prognosis GSE2034 (bone relapse 0/1) ---
    if cid == "breast_prognosis_gse2034":
        if r == "1":
            return "relapse"
        if r == "0":
            return "no_relapse"
        return "other"

    # --- Breast ER GSE7390 (er: 0/1) ---
    if cid == "breast_er_gse7390":
        if r == "1":
            return "er_positive"
        if r == "0":
            return "er_negative"
        return "other"

    # --- Hypoxia cell line GSE53012 (source name) ---
    if cid == "hypoxia_cellline_gse53012":
        if "hypoxia" in r:
            return "hypoxia"
        if "control" in r:
            return "normoxia"
        return "other"

    # --- Sepsis blood GSE28750 (health status) ---
    if cid == "sepsis_blood_gse28750":
        if "sepsis" in r:
            return "sepsis"
        if "healthy" in r:
            return "healthy"
        return "other"

    # --- TB blood GSE19491 (illness: ptb, control, latent, etc.) ---
    if cid == "tb_blood_gse19491":
        if r == "ptb" or "active tb" in r or "active tuberculosis" in r:
            return "active_tb"
        if "control" in r:
            return "healthy"
        return "other"

    # --- IPF lung GSE47460 (disease state) ---
    if cid == "ipf_lung_gse47460":
        if "interstitial" in r or "ild" in r or "ipf" in r:
            return "ild"
        if "control" in r:
            return "control"
        return "other"

    # --- COPD lung GSE47460 ---
    if cid == "copd_lung_gse47460":
        if "obstructive" in r or "copd" in r:
            return "copd"
        if "control" in r:
            return "control"
        return "other"

    # --- Influenza GSE101702 ---
    if cid == "influenza_pbmc_gse101702":
        if "influenza" in r:
            return "influenza"
        if "healthy" in r or "control" in r or "no infection" in r:
            return "healthy"
        return "other"

    # --- Trauma GSE36809 (title-based) ---
    if cid == "trauma_blood_gse36809":
        if "control" in r:
            return "control"
        if "subject" in r or "blood" in r:
            return "trauma"
        return "other"

    # --- HCC liver GSE6764 (source name) ---
    if cid == "hcc_liver_gse6764":
        if "hcc" in r or "hepatocellular" in r:
            return "hcc"
        if "normal" in r:
            return "normal"
        return "other"

    # --- Crohn's GSE112366 (diagnosis) ---
    if cid == "crohn_intestine_gse112366":
        if "crohn" in r:
            return "crohn"
        if "control" in r or "normal" in r or "healthy" in r:
            return "control"
        return "other"

    return "other"


# ---------------------------------------------------------------------------
# Matrix processing
# ---------------------------------------------------------------------------
def probes_to_genes(expr_df: pd.DataFrame, probe_to_gene: dict) -> pd.DataFrame:
    mapped = [p for p in expr_df.index if p in probe_to_gene]
    if not mapped:
        return pd.DataFrame()
    df = expr_df.loc[mapped].copy()
    df["gene_symbol"] = [probe_to_gene[p] for p in df.index]
    return df.groupby("gene_symbol").mean()


def normalize_matrix(gene_df: pd.DataFrame) -> pd.DataFrame:
    if gene_df.empty:
        return gene_df
    gene_df.index = gene_df.index.str.upper()
    gene_df = gene_df.groupby(gene_df.index).mean()
    gene_df = gene_df[gene_df.isna().mean(axis=1) <= 0.5]

    # Detect if log2 transform is needed: use 75th percentile and check for
    # right-skew characteristic of raw intensities/counts
    med = gene_df.median().median()
    p75 = float(np.nanpercentile(gene_df.values, 75))
    has_negatives = float(np.nanmin(gene_df.values)) < 0

    if has_negatives:
        # Data already log-transformed (has negatives)
        log.info(f"  Median={med:.2f}, p75={p75:.1f}, has negatives -> already log-scale")
    elif p75 > 50 or med > 50:
        log.info(f"  Median={med:.1f}, p75={p75:.1f}, applying log2(x+1)")
        gene_df = np.log2(gene_df + 1)
    else:
        log.info(f"  Median={med:.2f}, p75={p75:.1f}, no log transform needed")
    return gene_df.fillna(0)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Per-cohort pipeline
# ---------------------------------------------------------------------------
def process_cohort(cohort: dict) -> dict | None:
    cid = cohort["cohort_id"]
    acc = cohort["accession"]
    pid = cohort["platform_id"]
    log.info("=" * 60)
    log.info(f"Processing {cid} ({acc}, {pid})")
    log.info("=" * 60)

    matrix_out = FREEZE_MATRICES / f"{cid}.tsv"
    pheno_out = FREEZE_PHENO / f"{cid}.tsv"

    if matrix_out.exists() and pheno_out.exists():
        log.info("  Already frozen.")
        full_expr = pd.read_csv(matrix_out, sep="\t", index_col=0)
        pheno_df = pd.read_csv(pheno_out, sep="\t")
        return {
            "cohort_id": cid, "accession": acc, "platform_id": pid,
            "platform_type": cohort["platform_type"], "tissue": cohort["tissue"],
            "biological_program": cohort["biological_program"],
            "gene_count": len(full_expr), "sample_count": len(full_expr.columns),
            "case_count": int((pheno_df["phenotype"] == cohort["case_label"]).sum()),
            "control_count": int((pheno_df["phenotype"] == cohort["control_label"]).sum()),
            "matrix_sha256": sha256_file(matrix_out),
            "phenotype_sha256": sha256_file(pheno_out),
            "download_timestamp": datetime.now(timezone.utc).isoformat(),
            "status": "frozen_previously",
        }

    # Download series matrix
    if cohort.get("matrix_url_override"):
        gz_path = DATA_RAW / f"{acc}-{pid}_series_matrix.txt.gz"
        url = series_matrix_platform_url(acc, pid)
    else:
        gz_path = DATA_RAW / f"{acc}_series_matrix.txt.gz"
        url = series_matrix_url(acc)

    if not download_url(url, gz_path):
        log.error(f"  FAILED download for {acc}")
        return None

    expr_df, pheno, sample_ids = parse_series_matrix(gz_path)
    if expr_df is None:
        log.error(f"  FAILED parse for {acc}")
        return None
    log.info(f"  Parsed: {expr_df.shape[0]} probes x {expr_df.shape[1]} samples")

    # Probe -> gene
    p2g = get_probe_to_gene(pid)
    if not p2g:
        log.error(f"  No probe mapping for {pid}")
        return None

    gene_df = probes_to_genes(expr_df, p2g)
    if gene_df.empty:
        log.error(f"  Empty gene matrix")
        return None
    log.info(f"  Mapped: {gene_df.shape[0]} genes")

    gene_df = normalize_matrix(gene_df)
    if gene_df.empty:
        log.error(f"  Empty after normalization")
        return None

    # Phenotype
    pheno_df = extract_phenotype(pheno, cohort, sample_ids)
    case_lab = cohort["case_label"]
    ctrl_lab = cohort["control_label"]
    n_case = int((pheno_df["phenotype"] == case_lab).sum())
    n_ctrl = int((pheno_df["phenotype"] == ctrl_lab).sum())
    n_other = int((~pheno_df["phenotype"].isin([case_lab, ctrl_lab])).sum())
    log.info(f"  Phenotype: {n_case} {case_lab}, {n_ctrl} {ctrl_lab}, {n_other} other")

    # Filter to case+control
    valid = pheno_df["phenotype"].isin([case_lab, ctrl_lab])
    valid_ids = pheno_df.loc[valid, "sample_id"].tolist()

    if len(valid_ids) < 10:
        log.warning(f"  <10 valid samples ({len(valid_ids)}), keeping all")
    else:
        pheno_df = pheno_df[valid].reset_index(drop=True)
        common = [s for s in valid_ids if s in gene_df.columns]
        gene_df = gene_df[common]
        n_case = int((pheno_df["phenotype"] == case_lab).sum())
        n_ctrl = int((pheno_df["phenotype"] == ctrl_lab).sum())

    log.info(f"  Final: {gene_df.shape[0]} genes x {gene_df.shape[1]} samples")

    gene_df.index.name = "gene_symbol"
    gene_df.to_csv(matrix_out, sep="\t")
    pheno_df.to_csv(pheno_out, sep="\t", index=False)
    log.info(f"  Saved: {matrix_out.name}, {pheno_out.name}")

    return {
        "cohort_id": cid, "accession": acc, "platform_id": pid,
        "platform_type": cohort["platform_type"], "tissue": cohort["tissue"],
        "biological_program": cohort["biological_program"],
        "gene_count": gene_df.shape[0], "sample_count": gene_df.shape[1],
        "case_count": n_case, "control_count": n_ctrl,
        "matrix_sha256": sha256_file(matrix_out),
        "phenotype_sha256": sha256_file(pheno_out),
        "download_timestamp": datetime.now(timezone.utc).isoformat(),
        "status": "frozen",
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("GEO Cohort Downloader — REAL DATA ONLY")
    log.info(f"Target: {len(COHORTS)} cohorts")

    audit_records = []
    successes = []
    failures = []

    for cohort in COHORTS:
        try:
            result = process_cohort(cohort)
            if result:
                audit_records.append(result)
                successes.append(cohort["cohort_id"])
            else:
                failures.append((cohort["cohort_id"], "processing failed"))
        except Exception as e:
            log.error(f"  EXCEPTION: {cohort['cohort_id']}: {e}", exc_info=True)
            failures.append((cohort["cohort_id"], str(e)))

    # Audit
    audit_path = PROJECT_ROOT / "data" / "freeze" / "freeze_audit.json"
    audit = {
        "freeze_timestamp": datetime.now(timezone.utc).isoformat(),
        "total_attempted": len(COHORTS),
        "total_success": len(successes),
        "total_failed": len(failures),
        "cohorts": audit_records,
        "failures": [{"cohort_id": c, "reason": r} for c, r in failures],
    }
    with open(audit_path, "w") as f:
        json.dump(audit, f, indent=2)
    log.info(f"\nAudit: {audit_path}")

    # Manifest
    if audit_records:
        rows = []
        for r in audit_records:
            c = next(c for c in COHORTS if c["cohort_id"] == r["cohort_id"])
            rows.append({
                "cohort_id": r["cohort_id"], "geo_accession": r["accession"],
                "platform": r["platform_type"], "platform_id": r["platform_id"],
                "tissue": r["tissue"], "biological_program": r["biological_program"],
                "phenotype_column": "phenotype",
                "case_label": c["case_label"], "control_label": c["control_label"],
                "sample_count": r["sample_count"],
            })
        pd.DataFrame(rows).to_csv(CONFIG_DIR / "cohort_manifest.tsv", sep="\t", index=False)
        log.info(f"Manifest: {CONFIG_DIR / 'cohort_manifest.tsv'}")

    # Summary
    log.info(f"\n{'='*60}\nSUMMARY\n{'='*60}")
    log.info(f"Success: {len(successes)}")
    for cid in successes:
        r = next(a for a in audit_records if a["cohort_id"] == cid)
        log.info(f"  {cid}: {r['gene_count']}g x {r['sample_count']}s "
                 f"({r['platform_type']}, {r['case_count']}c/{r['control_count']}ctrl)")
    if failures:
        log.info(f"Failed: {len(failures)}")
        for cid, reason in failures:
            log.info(f"  {cid}: {reason}")

    return 0 if len(successes) >= 8 else 1


if __name__ == "__main__":
    sys.exit(main())
