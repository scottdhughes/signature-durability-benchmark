#!/usr/bin/env python3
"""Process 10+ new GEO datasets to expand the signature-durability-benchmark.

Goal: bring every biological program to k>=6 cohorts.

Datasets:
  Inflammation (+2):
    1. GSE95233  -- Septic shock blood (GPL570, 124 samples)
    2. GSE69528  -- Sepsis melioidosis blood (GPL10558, 138 samples)
  Proliferation (+1):
    3. GSE19188  -- Lung tumor vs normal (GPL570, 156 samples)
  Interferon (+2):
    4. GSE68310  -- Influenza challenge blood (GPL10558, 880 samples)
  Hypoxia (+2):
    5. GSE36895  -- ccRCC tumor vs normal kidney (GPL570, 52 human)
    6. GSE4290   -- GBM vs non-tumor brain (GPL570, 180 samples)
  EMT (+2):
    7. GSE12548  -- EMT time series ARPE19 (GPL570, 20 samples)
    8. GSE53845  -- IPF vs control lung (GPL6480, 48 samples)
"""

import gzip
import os
import sys
import re
from pathlib import Path
import pandas as pd
import numpy as np
from collections import Counter

BASE = Path(__file__).resolve().parent.parent
RAW = str(BASE / "data" / "raw")
MATRICES_DIR = str(BASE / "data" / "freeze" / "cohort_matrices")
PHENO_DIR = str(BASE / "data" / "freeze" / "cohort_phenotypes")

# ── Marker genes for verification ────────────────────────────────────────────

MARKERS = {
    'inflammation': {
        'genes': ['IL6', 'IL1B', 'TNF', 'CXCL8', 'CCL2', 'NFKBIA', 'PTGS2', 'ICAM1', 'SOD2', 'SERPINE1'],
        'direction': 'up',  # should be UP in case vs control
    },
    'proliferation': {
        'genes': ['MCM2', 'PCNA', 'RRM2', 'TYMS', 'TK1', 'BUB1', 'CCNE1', 'CDK2', 'E2F1', 'MKI67'],
        'direction': 'up',
    },
    'interferon': {
        'genes': ['STAT1', 'IRF1', 'IFIT1', 'IFIT2', 'ISG15', 'MX1', 'OAS1', 'GBP1', 'CXCL10', 'RSAD2'],
        'direction': 'up',
    },
    'hypoxia': {
        'genes': ['VEGFA', 'SLC2A1', 'LDHA', 'PGK1', 'ENO1', 'HK2', 'BNIP3', 'CA9', 'ADM', 'NDRG1'],
        'direction': 'up',
    },
    'emt': {
        'genes': ['VIM', 'FN1', 'SNAI2', 'ZEB1', 'MMP2', 'COL1A1', 'COL1A2', 'ACTA2', 'SPARC', 'TAGLN'],
        'direction': 'up',
    },
}


# ── Helpers ──────────────────────────────────────────────────────────────────

def load_annot(gpl_file):
    """Load GPL annotation file and return probe_id -> gene_symbol mapping."""
    path = os.path.join(RAW, gpl_file)

    rows = []
    header = None
    in_table = False

    with gzip.open(path, 'rt', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.strip() == '!platform_table_begin':
                in_table = True
                continue
            if line.strip() == '!platform_table_end':
                break
            if not in_table:
                continue
            if header is None:
                header = line.split('\t')
                continue
            parts = line.split('\t')
            rows.append(parts)

    if not rows or header is None:
        print(f"  ERROR: Could not parse annotation table from {gpl_file}")
        return None

    ncols = len(header)
    rows = [r[:ncols] for r in rows]
    rows = [r + [''] * (ncols - len(r)) for r in rows]

    df = pd.DataFrame(rows, columns=header)

    sym_col = None
    for col in ['Gene Symbol', 'Gene symbol', 'GENE_SYMBOL', 'Symbol', 'gene_symbol', 'ILMN_Gene']:
        if col in df.columns:
            sym_col = col
            break

    if sym_col is None:
        print(f"  Available columns in {gpl_file}: {list(df.columns)}")
        return None

    id_col = df.columns[0]

    mapping = df[[id_col, sym_col]].copy()
    mapping = mapping[mapping[sym_col].notna()]
    mapping = mapping[mapping[sym_col] != '']
    mapping = mapping[mapping[sym_col] != '---']

    mapping[sym_col] = mapping[sym_col].str.split('///').str[0].str.strip()
    mapping[sym_col] = mapping[sym_col].str.split(' /// ').str[0].str.strip()
    mapping[sym_col] = mapping[sym_col].str.upper()

    return dict(zip(mapping[id_col], mapping[sym_col]))


def parse_series_matrix(filepath):
    """Parse a GEO series matrix file. Returns (metadata_dict, expression_df, sample_ids)."""
    metadata = {}
    data_lines = []

    opener = gzip.open if filepath.endswith('.gz') else open

    with opener(filepath, 'rt', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!'):
                key = line.split('\t')[0]
                values = line.split('\t')[1:]
                values = [v.strip('"') for v in values]
                if key in metadata:
                    if isinstance(metadata[key], list) and isinstance(metadata[key][0], list):
                        metadata[key].append(values)
                    else:
                        metadata[key] = [metadata[key], values]
                else:
                    metadata[key] = values
            elif line.startswith('"ID_REF"') or (not line.startswith('!') and '\t' in line):
                if line.strip() and not line.startswith('!'):
                    data_lines.append(line)

    if not data_lines:
        return metadata, None, None

    header = data_lines[0].split('\t')
    header = [h.strip('"') for h in header]

    rows = []
    for line in data_lines[1:]:
        if line.strip() == '' or line.startswith('!'):
            continue
        parts = line.split('\t')
        parts = [p.strip('"') for p in parts]
        rows.append(parts)

    if not rows:
        return metadata, None, None

    df = pd.DataFrame(rows, columns=header[:len(rows[0])])
    df = df.set_index(header[0])

    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    sample_ids = list(df.columns)
    return metadata, df, sample_ids


def map_probes_to_genes(expr_df, probe_to_gene):
    """Map probe IDs to gene symbols, taking max per gene."""
    expr_df = expr_df.copy()
    expr_df['gene_symbol'] = expr_df.index.map(lambda x: probe_to_gene.get(str(x), None))
    expr_df = expr_df.dropna(subset=['gene_symbol'])
    expr_df = expr_df[expr_df['gene_symbol'] != '']
    expr_df = expr_df.groupby('gene_symbol').max()
    return expr_df


def save_cohort(cohort_id, expr_df, pheno_df):
    """Save expression matrix and phenotype table in benchmark format."""
    expr_path = os.path.join(MATRICES_DIR, f"{cohort_id}.tsv")
    expr_df.index.name = 'gene_symbol'
    expr_df.to_csv(expr_path, sep='\t')

    pheno_path = os.path.join(PHENO_DIR, f"{cohort_id}.tsv")
    pheno_df.to_csv(pheno_path, sep='\t', index=False)

    n_genes = len(expr_df)
    n_samples = len(expr_df.columns)
    print(f"  Saved {cohort_id}: {n_genes} genes, {n_samples} samples")
    print(f"    Groups: {pheno_df['phenotype'].value_counts().to_dict()}")
    return n_genes, n_samples


def get_sample_titles(metadata):
    titles = metadata.get('!Sample_title', [])
    if isinstance(titles[0], list):
        titles = titles[0]
    return [t.strip('"') for t in titles]


def get_sample_geo_accessions(metadata):
    acc = metadata.get('!Sample_geo_accession', [])
    if isinstance(acc[0], list):
        acc = acc[0]
    return [a.strip('"') for a in acc]


def get_sample_source(metadata):
    src = metadata.get('!Sample_source_name_ch1', [])
    if isinstance(src[0], list):
        src = src[0]
    return [s.strip('"') for s in src]


def get_characteristics(metadata, idx=0):
    chars = metadata.get('!Sample_characteristics_ch1', [])
    if not chars:
        return []
    if isinstance(chars[0], list):
        if idx < len(chars):
            return [c.strip('"') for c in chars[idx]]
        return []
    return [c.strip('"') for c in chars]


def verify_markers(expr_df, pheno_df, program, case_label, control_label):
    """Verify marker genes show expected direction. Returns count of concordant markers."""
    marker_info = MARKERS.get(program)
    if not marker_info:
        print(f"  WARNING: No markers defined for program '{program}'")
        return 0

    case_samples = pheno_df[pheno_df['phenotype'] == case_label]['sample_id'].tolist()
    ctrl_samples = pheno_df[pheno_df['phenotype'] == control_label]['sample_id'].tolist()

    case_samples = [s for s in case_samples if s in expr_df.columns]
    ctrl_samples = [s for s in ctrl_samples if s in expr_df.columns]

    if len(case_samples) == 0 or len(ctrl_samples) == 0:
        print(f"  WARNING: No overlapping samples for marker verification")
        return 0

    concordant = 0
    checked = 0
    print(f"  Marker verification ({program}):")
    for gene in marker_info['genes']:
        if gene not in expr_df.index:
            continue
        checked += 1
        case_mean = expr_df.loc[gene, case_samples].mean()
        ctrl_mean = expr_df.loc[gene, ctrl_samples].mean()
        fc = case_mean - ctrl_mean  # log-space difference
        expected = marker_info['direction']
        is_concordant = (expected == 'up' and fc > 0) or (expected == 'down' and fc < 0)
        if is_concordant:
            concordant += 1
        status = "OK" if is_concordant else "WRONG"
        print(f"    {gene}: case={case_mean:.2f} ctrl={ctrl_mean:.2f} diff={fc:+.2f} [{status}]")

    print(f"  Concordant: {concordant}/{checked} markers")
    return concordant


# ── Dataset Processors ──────────────────────────────────────────────────────

def process_gse95233():
    """GSE95233: Septic shock blood. Inflammation program.
    Case: septic shock (Day 0/1), Control: healthy controls
    Platform: GPL570
    """
    print("\n=== Processing GSE95233 (Inflammation: Septic shock) ===")

    filepath = os.path.join(RAW, "GSE95233_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)
    titles = get_sample_titles(metadata)
    sources = get_sample_source(metadata)

    # Parse source names: "Control NNN" vs "Patient NNN, DayX"
    # Title prefixes: CS/PC = Control, SV = Survivor, NS = Non-Survivor
    pheno_rows = []
    keep = []
    for acc, src, title in zip(accessions, sources, titles):
        src_clean = src.strip().strip('"')
        title_clean = title.strip().strip('"')
        if src_clean.lower().startswith('control'):
            pheno_rows.append({'sample_id': acc, 'phenotype': 'healthy', 'raw_label': src_clean})
            keep.append(acc)
        elif src_clean.lower().startswith('patient') and 'day1' in src_clean.lower().replace(' ', ''):
            # Day 1 septic shock samples = acute inflammation
            pheno_rows.append({'sample_id': acc, 'phenotype': 'sepsis', 'raw_label': src_clean})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[c for c in keep if c in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")

    if len(pheno_df[pheno_df['phenotype'] == 'healthy']) < 10 or len(pheno_df[pheno_df['phenotype'] == 'sepsis']) < 10:
        print("  ERROR: Less than 10 samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'inflammation', 'sepsis', 'healthy')
    if n_ok < 3:
        print("  ERROR: Fewer than 3 marker genes concordant")
        return False

    save_cohort("sepsis_shock_gse95233", expr_df, pheno_df)
    return True


def process_gse69528():
    """GSE69528: Sepsis melioidosis blood. Inflammation program.
    Case: melioidosis sepsis, Control: healthy uninfected
    Platform: GPL10558 (Illumina HumanHT-12 V4.0)
    """
    print("\n=== Processing GSE69528 (Inflammation: Melioidosis sepsis) ===")

    filepath = os.path.join(RAW, "GSE69528_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL10558.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)

    # Characteristic 0: pathogens, Characteristic 1: study group
    char0 = get_characteristics(metadata, 0)
    char1 = get_characteristics(metadata, 1)

    pheno_rows = []
    keep = []
    for acc, c0, c1 in zip(accessions, char0, char1):
        c0_clean = c0.strip().strip('"').lower()
        c1_clean = c1.strip().strip('"').lower()
        if 'control' in c0_clean and 'uninfected healthy' in c1_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'healthy', 'raw_label': c1_clean})
            keep.append(acc)
        elif 'melioidosis' in c0_clean or 'burkholderia' in c0_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'sepsis', 'raw_label': c0_clean + ' | ' + c1_clean})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[c for c in keep if c in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")

    if len(pheno_df[pheno_df['phenotype'] == 'healthy']) < 10 or len(pheno_df[pheno_df['phenotype'] == 'sepsis']) < 10:
        print("  ERROR: Less than 10 samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'inflammation', 'sepsis', 'healthy')
    if n_ok < 3:
        print("  ERROR: Fewer than 3 marker genes concordant")
        return False

    save_cohort("melioidosis_blood_gse69528", expr_df, pheno_df)
    return True


def process_gse19188():
    """GSE19188: Lung cancer tumor vs normal. Proliferation program.
    Case: tumor, Control: healthy tissue
    Platform: GPL570
    """
    print("\n=== Processing GSE19188 (Proliferation: Lung tumor vs normal) ===")

    filepath = os.path.join(RAW, "GSE19188_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)

    # Char[0] = tissue type: tumor or healthy
    char0 = get_characteristics(metadata, 0)

    pheno_rows = []
    keep = []
    for acc, c0 in zip(accessions, char0):
        c0_clean = c0.strip().strip('"').lower()
        if 'tumor' in c0_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'tumor', 'raw_label': c0_clean})
            keep.append(acc)
        elif 'healthy' in c0_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'normal', 'raw_label': c0_clean})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[c for c in keep if c in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")

    if len(pheno_df[pheno_df['phenotype'] == 'tumor']) < 10 or len(pheno_df[pheno_df['phenotype'] == 'normal']) < 10:
        print("  ERROR: Less than 10 samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'proliferation', 'tumor', 'normal')
    if n_ok < 3:
        print("  ERROR: Fewer than 3 marker genes concordant")
        return False

    save_cohort("lung_tumor_gse19188", expr_df, pheno_df)
    return True


def process_gse68310():
    """GSE68310: Influenza challenge blood. Interferon program.
    Case: Day 2-4 (peak interferon response), Control: Baseline
    Platform: GPL10558 (Illumina HumanHT-12 V4.0)
    """
    print("\n=== Processing GSE68310 (Interferon: Influenza challenge) ===")

    filepath = os.path.join(RAW, "GSE68310_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL10558.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)
    titles = get_sample_titles(metadata)

    # Char[1] has timepoints, Char[3] has infection type
    char_time = get_characteristics(metadata, 1)
    char_infect = get_characteristics(metadata, 3)

    pheno_rows = []
    keep = []
    for acc, ct, ci in zip(accessions, char_time, char_infect):
        ct_clean = ct.strip().strip('"').lower()
        ci_clean = ci.strip().strip('"').lower()
        # Only include influenza A virus samples
        if 'influenza' not in ci_clean:
            continue
        if 'baseline' in ct_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'baseline', 'raw_label': ct_clean})
            keep.append(acc)
        elif any(d in ct_clean for d in ['day2', 'day 2', 'day4', 'day 4']):
            pheno_rows.append({'sample_id': acc, 'phenotype': 'infected', 'raw_label': ct_clean})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[c for c in keep if c in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")

    if len(pheno_df[pheno_df['phenotype'] == 'baseline']) < 10 or len(pheno_df[pheno_df['phenotype'] == 'infected']) < 10:
        print("  ERROR: Less than 10 samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'interferon', 'infected', 'baseline')
    if n_ok < 3:
        print("  ERROR: Fewer than 3 marker genes concordant")
        return False

    save_cohort("influenza_challenge_gse68310", expr_df, pheno_df)
    return True


def process_gse36895():
    """GSE36895: ccRCC tumor vs normal kidney. Hypoxia program.
    Case: ccRCC tumor (constitutive hypoxia due to VHL loss), Control: normal cortex
    Platform: GPL570
    """
    print("\n=== Processing GSE36895 (Hypoxia: ccRCC vs normal kidney) ===")

    filepath = os.path.join(RAW, "GSE36895_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)

    # Char[1] has tissue type
    char1 = get_characteristics(metadata, 1)

    pheno_rows = []
    keep = []
    for acc, c1 in zip(accessions, char1):
        c1_clean = c1.strip().strip('"').lower()
        if 'ccrcc' in c1_clean or 'clear-cell' in c1_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'hypoxia', 'raw_label': c1_clean})
            keep.append(acc)
        elif 'normal cortex' in c1_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'normoxia', 'raw_label': c1_clean})
            keep.append(acc)
        # Skip mouse tumorgrafts

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[c for c in keep if c in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")

    if len(pheno_df[pheno_df['phenotype'] == 'hypoxia']) < 10 or len(pheno_df[pheno_df['phenotype'] == 'normoxia']) < 10:
        print("  ERROR: Less than 10 samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'hypoxia', 'hypoxia', 'normoxia')
    if n_ok < 3:
        print("  ERROR: Fewer than 3 marker genes concordant")
        return False

    save_cohort("ccrcc_kidney_gse36895", expr_df, pheno_df)
    return True


def process_gse4290():
    """GSE4290: Glioblastoma vs non-tumor brain. Hypoxia program.
    Case: GBM (glioblastoma grade 4, highly hypoxic), Control: non-tumor (epilepsy)
    Platform: GPL570
    """
    print("\n=== Processing GSE4290 (Hypoxia: GBM vs non-tumor brain) ===")

    filepath = os.path.join(RAW, "GSE4290_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)

    # Char[0] has histopathological diagnosis
    char0 = get_characteristics(metadata, 0)

    pheno_rows = []
    keep = []
    for acc, c0 in zip(accessions, char0):
        c0_clean = c0.strip().strip('"').lower()
        if 'glioblastoma' in c0_clean and 'grade 4' in c0_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'hypoxia', 'raw_label': c0_clean})
            keep.append(acc)
        elif 'non-tumor' in c0_clean or 'epilepsy' in c0_clean or c0_clean == 'non tumor':
            pheno_rows.append({'sample_id': acc, 'phenotype': 'normoxia', 'raw_label': c0_clean})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[c for c in keep if c in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")

    if len(pheno_df[pheno_df['phenotype'] == 'hypoxia']) < 10 or len(pheno_df[pheno_df['phenotype'] == 'normoxia']) < 10:
        print("  ERROR: Less than 10 samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'hypoxia', 'hypoxia', 'normoxia')
    if n_ok < 3:
        print("  ERROR: Fewer than 3 marker genes concordant")
        return False

    save_cohort("gbm_brain_gse4290", expr_df, pheno_df)
    return True


def process_gse12548():
    """GSE12548: EMT time series in ARPE19 cells. EMT program.
    Case: late TGF-beta + TNF treatment (16h, 24h, 42h, 60h)
    Control: untreated (0h)
    Platform: GPL570
    """
    print("\n=== Processing GSE12548 (EMT: ARPE19 TGF-beta time series) ===")

    filepath = os.path.join(RAW, "GSE12548_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)
    titles = get_sample_titles(metadata)

    print(f"  Total samples: {len(accessions)}")
    print(f"  Titles: {titles}")

    # Parse titles: "ARPE-19 0h-0", "ARPE-19 42h-2"
    pheno_rows = []
    keep = []
    for acc, title in zip(accessions, titles):
        m = re.search(r'(\d+)h', title)
        if m:
            hour = int(m.group(1))
            if hour == 0:
                pheno_rows.append({'sample_id': acc, 'phenotype': 'epithelial', 'raw_label': title})
                keep.append(acc)
            elif hour >= 16:  # 16h, 24h, 42h, 60h = late EMT
                pheno_rows.append({'sample_id': acc, 'phenotype': 'mesenchymal', 'raw_label': title})
                keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[c for c in keep if c in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")

    # Check minimum per group (may be small, but we need at least 3)
    epi_count = len(pheno_df[pheno_df['phenotype'] == 'epithelial'])
    mes_count = len(pheno_df[pheno_df['phenotype'] == 'mesenchymal'])
    if epi_count < 3 or mes_count < 3:
        print(f"  ERROR: Too few samples per group (epithelial={epi_count}, mesenchymal={mes_count})")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'emt', 'mesenchymal', 'epithelial')
    if n_ok < 3:
        print("  ERROR: Fewer than 3 marker genes concordant")
        return False

    save_cohort("emt_arpe19_gse12548", expr_df, pheno_df)
    return True


def process_gse53845():
    """GSE53845: IPF lung vs control. EMT program.
    Case: IPF (idiopathic pulmonary fibrosis - strong EMT phenotype)
    Control: healthy lung
    Platform: GPL6480 (Agilent)
    """
    print("\n=== Processing GSE53845 (EMT: IPF vs control lung) ===")

    filepath = os.path.join(RAW, "GSE53845_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL6480.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)

    # Char[2] has diagnosis: IPF or Control
    char2 = get_characteristics(metadata, 2)

    pheno_rows = []
    keep = []
    for acc, c2 in zip(accessions, char2):
        c2_clean = c2.strip().strip('"').lower()
        if 'ipf' in c2_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'mesenchymal', 'raw_label': c2_clean})
            keep.append(acc)
        elif 'control' in c2_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'epithelial', 'raw_label': c2_clean})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[c for c in keep if c in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")

    if len(pheno_df[pheno_df['phenotype'] == 'mesenchymal']) < 10 or len(pheno_df[pheno_df['phenotype'] == 'epithelial']) < 3:
        print("  ERROR: Too few samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'emt', 'mesenchymal', 'epithelial')
    if n_ok < 3:
        print("  ERROR: Fewer than 3 marker genes concordant")
        return False

    save_cohort("ipf_lung_gse53845", expr_df, pheno_df)
    return True


def process_gse111368():
    """GSE111368: Severe influenza blood. Interferon program.
    Case: H1N1 influenza at T1 (early infection, no bacterial co-infection)
    Control: Healthy controls (HC)
    Platform: GPL10558 (Illumina HumanHT-12 V4.0)
    """
    print("\n=== Processing GSE111368 (Interferon: Severe influenza) ===")

    filepath = os.path.join(RAW, "GSE111368_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL10558.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)

    # Char[1] = flu_type, Char[2] = timepoint, Char[8] = bacterial_status
    char_flu = get_characteristics(metadata, 1)
    char_time = get_characteristics(metadata, 2)
    char_bact = get_characteristics(metadata, 8)

    pheno_rows = []
    keep = []
    for acc, cf, ct, cb in zip(accessions, char_flu, char_time, char_bact):
        cf_clean = cf.strip().strip('"').lower()
        ct_clean = ct.strip().strip('"').lower()
        cb_clean = cb.strip().strip('"').lower()

        if 'hc' in ct_clean:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'healthy', 'raw_label': 'HC'})
            keep.append(acc)
        elif 'h1n1' in cf_clean and 't1' in ct_clean and 'yes' not in cb_clean:
            # H1N1, first timepoint, no bacterial co-infection
            pheno_rows.append({'sample_id': acc, 'phenotype': 'infected', 'raw_label': f'{cf_clean} {ct_clean}'})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[c for c in keep if c in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")

    if len(pheno_df[pheno_df['phenotype'] == 'healthy']) < 10 or len(pheno_df[pheno_df['phenotype'] == 'infected']) < 10:
        print("  ERROR: Less than 10 samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'interferon', 'infected', 'healthy')
    if n_ok < 3:
        print("  ERROR: Fewer than 3 marker genes concordant")
        return False

    save_cohort("influenza_severe_gse111368", expr_df, pheno_df)
    return True


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    results = {}

    # Inflammation
    results['GSE95233 (inflammation)'] = process_gse95233()
    results['GSE69528 (inflammation)'] = process_gse69528()

    # Proliferation
    results['GSE19188 (proliferation)'] = process_gse19188()

    # Interferon
    results['GSE68310 (interferon)'] = process_gse68310()
    results['GSE111368 (interferon)'] = process_gse111368()

    # Hypoxia
    results['GSE36895 (hypoxia)'] = process_gse36895()
    results['GSE4290 (hypoxia)'] = process_gse4290()

    # EMT
    results['GSE12548 (emt)'] = process_gse12548()
    results['GSE53845 (emt)'] = process_gse53845()

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, success in results.items():
        status = "OK" if success else "FAILED"
        print(f"  {name}: {status}")

    succeeded = sum(1 for v in results.values() if v)
    failed = sum(1 for v in results.values() if not v)
    print(f"\n  Total: {succeeded} succeeded, {failed} failed")


if __name__ == '__main__':
    main()
