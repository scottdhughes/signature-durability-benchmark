#!/usr/bin/env python3
"""Process 8 new GEO datasets for the signature-durability-benchmark.

Datasets:
  Interferon:
    1. GSE34205 — RSV bronchiolitis blood (GPL570)
    2. GSE73072 — Viral respiratory challenge (GPL14604, entrez IDs)
  Proliferation:
    3. GSE3494  — Breast cancer Uppsala (GPL96)
    4. GSE1456  — Breast cancer Stockholm (GPL96)
  Hypoxia:
    5. GSE4086  — Hypoxia in P493-6 lymphoma cells (GPL570)
    6. GSE47533 — Hypoxia time course MCF7 (GPL6884)
  EMT:
    7. GSE17708 — TGF-β EMT A549 time course (GPL570)
    8. GSE24202 — EMT in HMLE cells (GPL3921)
"""

import gzip
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(__file__).resolve().parent.parent
RAW = str(BASE / "data" / "raw")
MATRICES_DIR = str(BASE / "data" / "freeze" / "cohort_matrices")
PHENO_DIR = str(BASE / "data" / "freeze" / "cohort_phenotypes")

# ── Helpers ──────────────────────────────────────────────────────────────────

def load_annot(gpl_file):
    """Load GPL annotation file and return probe_id -> gene_symbol mapping."""
    path = os.path.join(RAW, gpl_file)

    # Read lines, skipping header metadata until we find the table
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

    # Ensure consistent column count
    ncols = len(header)
    rows = [r[:ncols] for r in rows]  # truncate extra cols
    rows = [r + [''] * (ncols - len(r)) for r in rows]  # pad short rows

    df = pd.DataFrame(rows, columns=header)

    # Find gene symbol column
    sym_col = None
    for col in ['Gene Symbol', 'Gene symbol', 'GENE_SYMBOL', 'Symbol', 'gene_symbol', 'ILMN_Gene']:
        if col in df.columns:
            sym_col = col
            break

    if sym_col is None:
        print(f"  Available columns in {gpl_file}: {list(df.columns)}")
        return None

    # Find probe ID column
    id_col = df.columns[0]  # Usually 'ID'

    mapping = df[[id_col, sym_col]].copy()
    mapping = mapping[mapping[sym_col].notna()]
    mapping = mapping[mapping[sym_col] != '']
    mapping = mapping[mapping[sym_col] != '---']

    # Handle multi-gene probes: take first symbol
    mapping[sym_col] = mapping[sym_col].str.split('///').str[0].str.strip()
    mapping[sym_col] = mapping[sym_col].str.split(' /// ').str[0].str.strip()
    mapping[sym_col] = mapping[sym_col].str.upper()

    return dict(zip(mapping[id_col], mapping[sym_col]))


def load_entrez_mapping():
    """Load NCBI gene_info to get entrez_id -> gene_symbol mapping."""
    path = "/tmp/Homo_sapiens.gene_info.gz"
    df = pd.read_csv(path, sep='\t', usecols=['GeneID', 'Symbol'], dtype={'GeneID': str})
    df['Symbol'] = df['Symbol'].str.upper()
    return dict(zip(df['GeneID'], df['Symbol']))


def parse_series_matrix(filepath):
    """Parse a GEO series matrix file. Returns (metadata_dict, expression_df, sample_ids)."""
    metadata = {}
    data_lines = []
    in_data = False

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
            elif line.startswith('"ID_REF"') or (not line.startswith('!') and '\t' in line and not line.startswith('!')):
                if line.strip() and not line.startswith('!'):
                    data_lines.append(line)

    if not data_lines:
        return metadata, None, None

    # Parse header
    header = data_lines[0].split('\t')
    header = [h.strip('"') for h in header]

    # Parse data
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

    # Convert to float
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    sample_ids = list(df.columns)
    return metadata, df, sample_ids


def map_probes_to_genes(expr_df, probe_to_gene):
    """Map probe IDs to gene symbols, taking max per gene."""
    # Add gene symbol column
    expr_df = expr_df.copy()
    expr_df['gene_symbol'] = expr_df.index.map(lambda x: probe_to_gene.get(str(x), None))

    # Drop probes without gene mapping
    expr_df = expr_df.dropna(subset=['gene_symbol'])
    expr_df = expr_df[expr_df['gene_symbol'] != '']

    # Group by gene symbol and take max
    expr_df = expr_df.groupby('gene_symbol').max()

    return expr_df


def save_cohort(cohort_id, expr_df, pheno_df):
    """Save expression matrix and phenotype table in benchmark format."""
    # Expression matrix
    expr_path = os.path.join(MATRICES_DIR, f"{cohort_id}.tsv")
    expr_df.index.name = 'gene_symbol'
    expr_df.to_csv(expr_path, sep='\t')

    # Phenotype table
    pheno_path = os.path.join(PHENO_DIR, f"{cohort_id}.tsv")
    pheno_df.to_csv(pheno_path, sep='\t', index=False)

    n_genes = len(expr_df)
    n_samples = len(expr_df.columns)
    case_count = (pheno_df['phenotype'] == pheno_df['phenotype'].unique()[0]).sum()
    ctrl_count = n_samples - case_count

    print(f"  Saved {cohort_id}: {n_genes} genes, {n_samples} samples")
    print(f"    Case: {pheno_df['phenotype'].value_counts().to_dict()}")
    return n_genes, n_samples


def get_sample_titles(metadata):
    """Extract sample titles from metadata."""
    titles = metadata.get('!Sample_title', [])
    if isinstance(titles[0], list):
        titles = titles[0]
    return [t.strip('"') for t in titles]


def get_sample_geo_accessions(metadata):
    """Extract sample GEO accessions."""
    acc = metadata.get('!Sample_geo_accession', [])
    if isinstance(acc[0], list):
        acc = acc[0]
    return [a.strip('"') for a in acc]


def get_sample_characteristics(metadata, idx=0):
    """Extract sample characteristics at a given index."""
    chars = metadata.get('!Sample_characteristics_ch1', [])
    if not chars:
        return []
    if isinstance(chars[0], list):
        if idx < len(chars):
            return [c.strip('"') for c in chars[idx]]
        return []
    return [c.strip('"') for c in chars]


def get_sample_source(metadata):
    """Extract sample source names."""
    src = metadata.get('!Sample_source_name_ch1', [])
    if isinstance(src[0], list):
        src = src[0]
    return [s.strip('"') for s in src]


# ── Dataset Processors ──────────────────────────────────────────────────────

def process_gse34205():
    """GSE34205: RSV bronchiolitis blood. Interferon program.
    Case: RSV infection, Control: healthy
    Platform: GPL570 (Affymetrix HG-U133 Plus 2.0)
    """
    print("\n=== Processing GSE34205 (Interferon: RSV bronchiolitis) ===")

    filepath = os.path.join(RAW, "GSE34205_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    # Get sample titles for classification
    titles = get_sample_titles(metadata)
    accessions = get_sample_geo_accessions(metadata)

    # Build phenotype: RSV = case, healthy = control
    pheno_rows = []
    keep_samples = []
    for acc, title in zip(accessions, titles):
        title_lower = title.lower()
        if 'rsv' in title_lower and 'healthy' not in title_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'rsv_infection', 'raw_label': title})
            keep_samples.append(acc)
        elif 'healthy' in title_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'healthy', 'raw_label': title})
            keep_samples.append(acc)
        # Skip influenza samples for a cleaner RSV vs healthy comparison

    pheno_df = pd.DataFrame(pheno_rows)

    # Filter expression to keep only selected samples
    expr_df = expr_df[[c for c in keep_samples if c in expr_df.columns]]

    # Map probes to genes
    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        print("  ERROR: Could not load GPL570 annotation")
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    # Align
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    save_cohort("rsv_blood_gse34205", expr_df, pheno_df)
    return True


def process_gse73072():
    """GSE73072: Viral respiratory challenge. Interferon program.
    Case: peak infection timepoint, Control: baseline (pre-inoculation)
    Platform: GPL14604 (custom Affymetrix with Entrez Gene IDs)
    """
    print("\n=== Processing GSE73072 (Interferon: Viral challenge) ===")

    filepath = os.path.join(RAW, "GSE73072_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    # Get sample titles
    titles = get_sample_titles(metadata)
    accessions = get_sample_geo_accessions(metadata)

    # Map Entrez Gene IDs to symbols
    entrez_map = load_entrez_mapping()

    # The probe IDs are like "10000_at" - extract entrez ID
    def probe_to_symbol(probe_id):
        eid = str(probe_id).replace('_at', '')
        return entrez_map.get(eid, None)

    expr_df_mapped = expr_df.copy()
    expr_df_mapped['gene_symbol'] = expr_df_mapped.index.map(probe_to_symbol)
    expr_df_mapped = expr_df_mapped.dropna(subset=['gene_symbol'])
    expr_df_mapped = expr_df_mapped[expr_df_mapped['gene_symbol'] != '']
    expr_df_mapped['gene_symbol'] = expr_df_mapped['gene_symbol'].str.upper()
    expr_df_mapped = expr_df_mapped.groupby('gene_symbol').max()
    expr_df = expr_df_mapped

    print(f"  Total samples: {len(accessions)}")
    print(f"  Sample titles (first 5): {titles[:5]}")

    # Parse hour from title: "H1N1 DEE3 Subject 1, Hour -21"
    # Use Hour <= 0 as baseline, Hour 45-96 as peak infection
    import re
    pheno_rows = []
    keep_samples = []

    for acc, title in zip(accessions, titles):
        m = re.search(r'Hour\s+(-?\d+)', title)
        if m:
            hour = int(m.group(1))
            if hour <= 0:
                pheno_rows.append({'sample_id': acc, 'phenotype': 'baseline', 'raw_label': title})
                keep_samples.append(acc)
            elif 45 <= hour <= 96:
                pheno_rows.append({'sample_id': acc, 'phenotype': 'infected', 'raw_label': title})
                keep_samples.append(acc)
            # Skip intermediate timepoints for cleaner comparison

    if len(pheno_rows) < 10:
        print(f"  ERROR: Only {len(pheno_rows)} samples classified")
        return False

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[c for c in keep_samples if c in expr_df.columns]]

    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    save_cohort("viral_challenge_gse73072", expr_df, pheno_df)
    return True


def process_gse3494():
    """GSE3494: Breast cancer Uppsala. Proliferation program.
    Since grade is not in metadata, use p53 status or ER status as proliferation proxy.
    Fall back to using all samples split by median of a proliferation marker.
    Platform: GPL96 (Affymetrix HG-U133A)
    """
    print("\n=== Processing GSE3494 (Proliferation: breast Uppsala) ===")

    filepath = os.path.join(RAW, "GSE3494-GPL96_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    # Load GPL96 annotation
    probe_map = load_annot("GPL96.annot.gz")
    if probe_map is None:
        print("  ERROR: Could not load GPL96 annotation")
        return False

    # Map probes to genes
    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)
    titles = get_sample_titles(metadata)

    # GSE3494 is the Miller et al. breast cancer dataset
    # Grade is not in GEO metadata. Use MKI67 expression as proliferation proxy
    # MKI67 is the gold-standard proliferation marker
    if 'MKI67' in expr_df.index:
        mki67 = expr_df.loc['MKI67']
        median_val = mki67.median()

        pheno_rows = []
        for acc in expr_df.columns:
            val = mki67[acc]
            if val >= median_val:
                pheno_rows.append({'sample_id': acc, 'phenotype': 'high_proliferation', 'raw_label': f'MKI67={val:.2f}'})
            else:
                pheno_rows.append({'sample_id': acc, 'phenotype': 'low_proliferation', 'raw_label': f'MKI67={val:.2f}'})

        pheno_df = pd.DataFrame(pheno_rows)
        print("  NOTE: GSE3494 is a legacy exploratory proliferation cohort and is not part of the active 35-cohort freeze; not saving into data/freeze/")
        return True
    else:
        print("  WARNING: MKI67 not found, trying TOP2A...")
        if 'TOP2A' in expr_df.index:
            marker = expr_df.loc['TOP2A']
            median_val = marker.median()
            pheno_rows = []
            for acc in expr_df.columns:
                val = marker[acc]
                if val >= median_val:
                    pheno_rows.append({'sample_id': acc, 'phenotype': 'high_proliferation', 'raw_label': f'TOP2A={val:.2f}'})
                else:
                    pheno_rows.append({'sample_id': acc, 'phenotype': 'low_proliferation', 'raw_label': f'TOP2A={val:.2f}'})
            pheno_df = pd.DataFrame(pheno_rows)
            print("  NOTE: GSE3494 is a legacy exploratory proliferation cohort and is not part of the active 35-cohort freeze; not saving into data/freeze/")
            return True
        print("  ERROR: No proliferation marker found")
        return False


def process_gse1456():
    """GSE1456: Breast cancer Stockholm. Proliferation program.
    Has RELAPSE in characteristics. Use relapse=case, no_relapse=control.
    Platform: GPL96 (Affymetrix HG-U133A)
    """
    print("\n=== Processing GSE1456 (Proliferation: breast Stockholm) ===")

    filepath = os.path.join(RAW, "GSE1456-GPL96_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    # Load GPL96 annotation
    probe_map = load_annot("GPL96.annot.gz")
    if probe_map is None:
        print("  ERROR: Could not load GPL96 annotation")
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)

    # Get RELAPSE characteristic
    # Find which characteristic line has RELAPSE
    all_chars = metadata.get('!Sample_characteristics_ch1', [])
    relapse_chars = None

    if isinstance(all_chars[0], list):
        for char_row in all_chars:
            if any('RELAPSE:' in str(c) for c in char_row):
                relapse_chars = [c.strip('"') for c in char_row]
                break
    else:
        if any('RELAPSE:' in str(c) for c in all_chars):
            relapse_chars = [c.strip('"') for c in all_chars]

    if relapse_chars is None:
        print("  ERROR: Could not find RELAPSE characteristic")
        return False

    pheno_rows = []
    for acc, char in zip(accessions, relapse_chars):
        if 'RELAPSE: 1' in char:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'relapse', 'raw_label': '1'})
        elif 'RELAPSE: 0' in char:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'no_relapse', 'raw_label': '0'})

    pheno_df = pd.DataFrame(pheno_rows)

    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    save_cohort("breast_relapse_gse1456", expr_df, pheno_df)
    return True


def process_gse4086():
    """GSE4086: Hypoxia in P493-6 lymphoma cells. Hypoxia program.
    Case: hypoxia, Control: non-hypoxia
    Platform: GPL570
    Note: Only 4 samples (2 hypoxia, 2 normoxia) — too small.
    """
    print("\n=== Processing GSE4086 (Hypoxia: P493-6 cells) ===")

    filepath = os.path.join(RAW, "GSE4086_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    titles = get_sample_titles(metadata)
    accessions = get_sample_geo_accessions(metadata)

    print(f"  Total samples: {len(accessions)}")
    print(f"  Titles: {titles}")

    # Only 4 samples - too small for benchmark (need >=10 per group)
    if len(accessions) < 10:
        print(f"  WARNING: Only {len(accessions)} samples — below minimum threshold of 10 per group")
        print("  Skipping GSE4086 — will use GSE47533 instead for hypoxia")
        return False

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    pheno_rows = []
    for acc, title in zip(accessions, titles):
        if 'hypoxia' in title.lower() and 'non' not in title.lower():
            pheno_rows.append({'sample_id': acc, 'phenotype': 'hypoxia', 'raw_label': title})
        else:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'normoxia', 'raw_label': title})

    pheno_df = pd.DataFrame(pheno_rows)

    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    save_cohort("hypoxia_p493_gse4086", expr_df, pheno_df)
    return True


def process_gse3188():
    """GSE3188: HIF knockdown in MCF7 cells under hypoxia. Hypoxia program.
    Case: oligofectamine control (intact hypoxia response)
    Control: HIF1+HIF2 siRNA (ablated hypoxia response)
    Platform: GPL570
    """
    print("\n=== Processing GSE3188 (Hypoxia: MCF7 HIF knockdown) ===")

    filepath = os.path.join(RAW, "GSE3188-GPL570_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    titles = get_sample_titles(metadata)
    accessions = get_sample_geo_accessions(metadata)

    print(f"  Total samples: {len(accessions)}")
    print(f"  Titles: {titles}")

    if len(accessions) < 6:
        print(f"  WARNING: Only {len(accessions)} samples — very small dataset")

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    # OF = oligofectamine (control reagent, hypoxia response intact) = hypoxia case
    # HIF12 = HIF1+HIF2 knockdown = normoxia-like (ablated hypoxia response) = control
    pheno_rows = []
    for acc, title in zip(accessions, titles):
        if 'OF' in title and 'HIF' not in title:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'hypoxia', 'raw_label': title})
        elif 'HIF12' in title or 'HIF1_' in title or 'HIF2_' in title:
            # Include all HIF knockdowns as "normoxia-like"
            pheno_rows.append({'sample_id': acc, 'phenotype': 'normoxia', 'raw_label': title})

    pheno_df = pd.DataFrame(pheno_rows)

    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    save_cohort("hypoxia_mcf7_gse3188", expr_df, pheno_df)
    return True


def process_gse47533():
    """GSE47533: Hypoxia time course in MCF7 cells. Hypoxia program.
    Case: hypoxia (16h, 32h, 48h), Control: normoxia
    Platform: GPL6884 (Illumina HumanWG-6 v3.0)
    """
    print("\n=== Processing GSE47533 (Hypoxia: MCF7 time course) ===")

    filepath = os.path.join(RAW, "GSE47533_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    titles = get_sample_titles(metadata)
    accessions = get_sample_geo_accessions(metadata)

    print(f"  Total samples: {len(accessions)}")
    print(f"  Titles: {titles}")

    # Filter to mRNA samples only (not miRNA)
    # Load annotation
    probe_map = load_annot("GPL6884.annot.gz")
    if probe_map is None:
        print("  ERROR: Could not load GPL6884 annotation")
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    pheno_rows = []
    for acc, title in zip(accessions, titles):
        title_lower = title.lower()
        if '[mrna]' not in title_lower:
            continue  # Skip non-mRNA samples
        if 'normoxia' in title_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'normoxia', 'raw_label': title})
        elif 'hypoxia' in title_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'hypoxia', 'raw_label': title})

    pheno_df = pd.DataFrame(pheno_rows)

    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    save_cohort("hypoxia_timecourse_gse47533", expr_df, pheno_df)
    return True


def process_gse18494():
    """GSE18494: Hypoxia in HepG2, U87, MDA-MB231 cells. Hypoxia program.
    Case: hypoxia (4h, 8h, 12h), Control: normoxia (0h)
    Platform: GPL9419 (Custom Affymetrix with RefSeq probe IDs)
    36 samples total (3 cell lines x 4 timepoints x 3 replicates)
    """
    print("\n=== Processing GSE18494 (Hypoxia: multi-cell type) ===")

    filepath = os.path.join(RAW, "GSE18494_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    titles = get_sample_titles(metadata)
    accessions = get_sample_geo_accessions(metadata)

    print(f"  Total samples: {len(accessions)}")

    # Build RefSeq -> Gene Symbol mapping from HGNC
    import csv
    refseq_map = {}
    with open('/tmp/hgnc_complete_set.txt', 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            symbol = str(row['symbol']).upper()
            refseq = row.get('refseq_accession', '')
            if refseq and refseq != '':
                for acc in str(refseq).split('|'):
                    acc = acc.strip()
                    if acc.startswith('NM_'):
                        refseq_map[acc] = symbol

    # Map probe IDs (NM_000014_at) -> RefSeq (NM_000014) -> Gene Symbol
    def probe_to_symbol(probe_id):
        refseq = str(probe_id).replace('_at', '')
        return refseq_map.get(refseq, None)

    expr_df_mapped = expr_df.copy()
    expr_df_mapped['gene_symbol'] = expr_df_mapped.index.map(probe_to_symbol)
    expr_df_mapped = expr_df_mapped.dropna(subset=['gene_symbol'])
    expr_df_mapped = expr_df_mapped[expr_df_mapped['gene_symbol'] != '']
    expr_df_mapped['gene_symbol'] = expr_df_mapped['gene_symbol'].str.upper()
    expr_df_mapped = expr_df_mapped.groupby('gene_symbol').max()
    expr_df = expr_df_mapped

    # Parse titles: "HepG2 hypoxia (0.5%O2) treatment 0h, biological rep1"
    # 0h = normoxia, 4h/8h/12h = hypoxia
    pheno_rows = []
    for acc, title in zip(accessions, titles):
        if 'treatment 0h' in title:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'normoxia', 'raw_label': title})
        elif any(f'treatment {t}h' in title for t in ['4', '8', '12']):
            pheno_rows.append({'sample_id': acc, 'phenotype': 'hypoxia', 'raw_label': title})

    pheno_df = pd.DataFrame(pheno_rows)

    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    save_cohort("hypoxia_multicell_gse18494", expr_df, pheno_df)
    return True


def process_gse17708():
    """GSE17708: TGF-β EMT in A549 cells. EMT program.
    Case: TGF-β treated (24h and 72h = late timepoints), Control: untreated (0h)
    Platform: GPL570
    """
    print("\n=== Processing GSE17708 (EMT: TGF-β A549) ===")

    filepath = os.path.join(RAW, "GSE17708_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    titles = get_sample_titles(metadata)
    accessions = get_sample_geo_accessions(metadata)

    print(f"  Total samples: {len(accessions)}")

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    # Case: late TGF-β treatment (8h, 16h, 24h, 72h), Control: untreated
    # Using 8h+ because EMT transcriptional program is well-engaged by 8h
    pheno_rows = []
    for acc, title in zip(accessions, titles):
        title_lower = title.lower()
        if 'untreated' in title_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'epithelial', 'raw_label': title})
        elif any(x in title_lower for x in ['72 h', '24 h', '16 h', '8 h']):
            pheno_rows.append({'sample_id': acc, 'phenotype': 'mesenchymal', 'raw_label': title})
        # Skip very early timepoints (0.5h, 1h, 2h, 4h)

    pheno_df = pd.DataFrame(pheno_rows)

    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    save_cohort("emt_tgfb_gse17708", expr_df, pheno_df)
    return True


def process_gse24202():
    """GSE24202: EMT in HMLE cells. EMT program.
    Case: EMT-inducing (TGFb, Twist, Gsc, Snail, shEcad)
    Control: control vectors (pWZL, shGFP)
    Platform: GPL3921 (Affymetrix HT_HG-U133A)
    """
    print("\n=== Processing GSE24202 (EMT: HMLE cells) ===")

    filepath = os.path.join(RAW, "GSE24202_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    titles = get_sample_titles(metadata)
    accessions = get_sample_geo_accessions(metadata)

    print(f"  Total samples: {len(accessions)}")
    print(f"  Titles: {titles}")

    # Load GPL3921 annotation
    probe_map = load_annot("GPL3921.annot.gz")
    if probe_map is None:
        print("  ERROR: Could not load GPL3921 annotation")
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    # EMT-inducing: TGFb, Twist, Gsc, Snail, shEcad
    # Control: pWZL, shGFP
    pheno_rows = []
    for acc, title in zip(accessions, titles):
        title_lower = title.lower()
        if any(x in title_lower for x in ['pwzl', 'shgfp']):
            pheno_rows.append({'sample_id': acc, 'phenotype': 'epithelial', 'raw_label': title})
        elif any(x in title_lower for x in ['tgfb', 'twist', 'gsc', 'snail', 'shecad']):
            pheno_rows.append({'sample_id': acc, 'phenotype': 'mesenchymal', 'raw_label': title})

    pheno_df = pd.DataFrame(pheno_rows)

    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    save_cohort("emt_hmle_gse24202", expr_df, pheno_df)
    return True


def process_gse43495():
    """GSE43495: EMT in HMLE cells with Twist/Snail/Slug. EMT program.
    Case: EMT transcription factors (Twist, Snail, Slug)
    Control: untreated + empty vector
    Platform: GPL6883 (Illumina HumanRef-8 v3.0)
    """
    print("\n=== Processing GSE43495 (EMT: HMLE Twist/Snail/Slug) ===")

    filepath = os.path.join(RAW, "GSE43495_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)

    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    titles = get_sample_titles(metadata)
    accessions = get_sample_geo_accessions(metadata)

    print(f"  Total samples: {len(accessions)}")
    print(f"  Titles: {titles}")

    # Load GPL6883 annotation
    probe_map = load_annot("GPL6883.annot.gz")
    if probe_map is None:
        print("  ERROR: Could not load GPL6883 annotation")
        return False

    expr_df = map_probes_to_genes(expr_df, probe_map)

    pheno_rows = []
    for acc, title in zip(accessions, titles):
        title_lower = title.lower()
        if 'untreated' in title_lower or 'empty vector' in title_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'epithelial', 'raw_label': title})
        elif any(x in title_lower for x in ['twist', 'snail', 'slug']):
            pheno_rows.append({'sample_id': acc, 'phenotype': 'mesenchymal', 'raw_label': title})

    pheno_df = pd.DataFrame(pheno_rows)

    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    save_cohort("emt_mammary_gse43495", expr_df, pheno_df)
    return True


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    results = {}

    # Interferon
    results['GSE34205'] = process_gse34205()
    results['GSE73072'] = process_gse73072()

    # Proliferation
    results['GSE3494'] = process_gse3494()
    results['GSE1456'] = process_gse1456()

    # Hypoxia
    results['GSE4086'] = process_gse4086()
    results['GSE3188'] = process_gse3188()
    results['GSE47533'] = process_gse47533()
    results['GSE18494'] = process_gse18494()

    # EMT
    results['GSE17708'] = process_gse17708()
    results['GSE24202'] = process_gse24202()
    results['GSE43495'] = process_gse43495()

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for gse, success in results.items():
        status = "OK" if success else "FAILED"
        print(f"  {gse}: {status}")

if __name__ == '__main__':
    main()
