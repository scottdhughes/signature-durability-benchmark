#!/usr/bin/env python3
"""Process 5 new IFN-program GEO cohorts to expand the interferon panel
from k=6 to k=11, covering tissue-resident IFN (psoriasis skin, SLE PBMC/blood)
and additional viral infection (dengue).

Expansion targets:
  Tissue IFN (+2 SLE, +2 psoriasis):
    1. GSE50772  -- SLE PBMCs vs controls (GPL570, 81 samples)
    2. GSE49454  -- SLE whole blood vs healthy (GPL10558, 177 samples)
    3. GSE13355  -- Psoriasis lesional vs normal skin (GPL570, 180 samples)
    4. GSE14905  -- Psoriasis lesional vs normal skin (GPL570, 82 samples)
  Viral IFN (+1):
    5. GSE51808  -- Dengue acute vs healthy blood (GPL13158, 56 samples)

All use the "interferon" biological_program assignment. Validation: marker gene
verification (STAT1, IRF1, IFIT1, IFIT2, ISG15, MX1, OAS1, GBP1, CXCL10, RSAD2).
"""

import gzip
import os
import sys
import pandas as pd
import numpy as np

# Import shared helpers from the existing expansion script
sys.path.insert(0, os.path.dirname(__file__))
from process_expansion_cohorts import (
    BASE, RAW, MATRICES_DIR, PHENO_DIR,
    MARKERS,
    load_annot,
    parse_series_matrix,
    map_probes_to_genes,
    save_cohort,
    get_sample_titles,
    get_sample_geo_accessions,
    get_sample_source,
    get_characteristics,
    verify_markers,
)


# ── Dataset Processors ──────────────────────────────────────────────────────

def process_gse50772():
    """GSE50772: SLE PBMCs vs controls. Interferon program.
    Platform: GPL570 (Affymetrix HG-U133 Plus 2.0)
    Characteristics: 'disease status: SLE' vs 'disease status: Control'
    """
    print("\n=== Processing GSE50772 (Interferon: SLE PBMC) ===")
    filepath = os.path.join(RAW, "GSE50772_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)
    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False
    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)
    chars = get_characteristics(metadata, idx=0)  # disease status

    pheno_rows = []
    keep = []
    for acc, c in zip(accessions, chars):
        c_lower = c.lower()
        if "disease status: sle" in c_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'sle', 'raw_label': c})
            keep.append(acc)
        elif "disease status: control" in c_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'healthy', 'raw_label': c})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[col for col in keep if col in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")
    if (pheno_df['phenotype'] == 'sle').sum() < 10 or (pheno_df['phenotype'] == 'healthy').sum() < 5:
        print("  ERROR: Insufficient samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'interferon', 'sle', 'healthy')
    if n_ok < 3:
        print("  WARNING: Fewer than 3 marker genes concordant, but saving anyway for sensitivity test")
    save_cohort("sle_pbmc_gse50772", expr_df, pheno_df)
    return True


def process_gse49454():
    """GSE49454: SLE whole blood vs healthy. Interferon program.
    Platform: GPL10558 (Illumina HumanHT-12 V4.0)
    Characteristics: 'group: SLE' vs 'group: Healthy control of SLE'
    """
    print("\n=== Processing GSE49454 (Interferon: SLE whole blood) ===")
    filepath = os.path.join(RAW, "GSE49454_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)
    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL10558.annot.gz")
    if probe_map is None:
        return False
    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)
    chars = get_characteristics(metadata, idx=0)  # group

    pheno_rows = []
    keep = []
    for acc, c in zip(accessions, chars):
        c_lower = c.lower()
        if "group: sle" in c_lower and "healthy" not in c_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'sle', 'raw_label': c})
            keep.append(acc)
        elif "group: healthy" in c_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'healthy', 'raw_label': c})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[col for col in keep if col in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")
    if (pheno_df['phenotype'] == 'sle').sum() < 10 or (pheno_df['phenotype'] == 'healthy').sum() < 5:
        print("  ERROR: Insufficient samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'interferon', 'sle', 'healthy')
    if n_ok < 3:
        print("  WARNING: Fewer than 3 marker genes concordant")
    save_cohort("sle_blood_gse49454", expr_df, pheno_df)
    return True


def process_gse13355():
    """GSE13355: Psoriasis lesional vs normal skin. Interferon program.
    Platform: GPL570
    Characteristics: 'normal skin from controls' / 'involved skin from cases' / 'uninvolved skin from cases'
    We use NN (normal controls) vs PP (involved / lesional) only.
    """
    print("\n=== Processing GSE13355 (Interferon: Psoriasis lesional skin) ===")
    filepath = os.path.join(RAW, "GSE13355_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)
    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False
    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)
    chars = get_characteristics(metadata, idx=0)

    pheno_rows = []
    keep = []
    for acc, c in zip(accessions, chars):
        c_lower = c.lower()
        if "normal skin from controls" in c_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'healthy', 'raw_label': c})
            keep.append(acc)
        elif "involved skin from cases" in c_lower and "uninvolved" not in c_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'psoriasis', 'raw_label': c})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[col for col in keep if col in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")
    if (pheno_df['phenotype'] == 'psoriasis').sum() < 10 or (pheno_df['phenotype'] == 'healthy').sum() < 10:
        print("  ERROR: Insufficient samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'interferon', 'psoriasis', 'healthy')
    if n_ok < 3:
        print("  WARNING: Fewer than 3 marker genes concordant")
    save_cohort("psoriasis_skin_gse13355", expr_df, pheno_df)
    return True


def process_gse14905():
    """GSE14905: Psoriasis lesional vs normal skin. Interferon program.
    Platform: GPL570
    Characteristics: 'skin: normal' vs 'skin: lesional' (keep only these two classes,
    exclude non-lesional from cases to keep the contrast clean).
    """
    print("\n=== Processing GSE14905 (Interferon: Psoriasis skin) ===")
    filepath = os.path.join(RAW, "GSE14905_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)
    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL570.annot.gz")
    if probe_map is None:
        return False
    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)
    chars0 = get_characteristics(metadata, idx=0)  # "skin: normal" or "disease: psoriasis"
    chars1 = get_characteristics(metadata, idx=1)  # "" or "skin: lesional"/"skin: non-lesional"

    pheno_rows = []
    keep = []
    for acc, c0, c1 in zip(accessions, chars0, chars1 if chars1 else [''] * len(accessions)):
        c0_lower = c0.lower() if c0 else ""
        c1_lower = c1.lower() if c1 else ""
        if "skin: normal" in c0_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'healthy', 'raw_label': c0})
            keep.append(acc)
        elif "disease: psoriasis" in c0_lower and "skin: lesional" in c1_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'psoriasis', 'raw_label': f"{c0} | {c1}"})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[col for col in keep if col in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")
    if (pheno_df['phenotype'] == 'psoriasis').sum() < 10 or (pheno_df['phenotype'] == 'healthy').sum() < 10:
        print("  ERROR: Insufficient samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'interferon', 'psoriasis', 'healthy')
    if n_ok < 3:
        print("  WARNING: Fewer than 3 marker genes concordant")
    save_cohort("psoriasis_skin_gse14905", expr_df, pheno_df)
    return True


def process_gse51808():
    """GSE51808: Dengue acute infection vs healthy. Interferon program.
    Platform: GPL13158 (Affymetrix HT-HG-U133A)
    Characteristics (idx 2): 'infection: DENV' vs 'infection: control'
    Characteristics (idx 3): 'status: DF'/'DHF' (acute) vs 'status: control'/'convalescent'
    We use acute (DF + DHF) vs healthy control only (exclude convalescent).
    """
    print("\n=== Processing GSE51808 (Interferon: Dengue acute infection) ===")
    filepath = os.path.join(RAW, "GSE51808_series_matrix.txt.gz")
    metadata, expr_df, sample_ids = parse_series_matrix(filepath)
    if expr_df is None:
        print("  ERROR: Could not parse expression data")
        return False

    probe_map = load_annot("GPL13158.annot.gz")
    if probe_map is None:
        return False
    expr_df = map_probes_to_genes(expr_df, probe_map)

    accessions = get_sample_geo_accessions(metadata)
    chars_inf = get_characteristics(metadata, idx=1)   # 'infection: DENV' / 'infection: control'
    chars_status = get_characteristics(metadata, idx=2)  # 'status: DF'/'DHF'/'convalescent'/'control'

    pheno_rows = []
    keep = []
    for acc, inf, status in zip(accessions, chars_inf, chars_status):
        inf_lower = inf.lower() if inf else ""
        st_lower = status.lower() if status else ""
        if "infection: control" in inf_lower and "status: control" in st_lower:
            pheno_rows.append({'sample_id': acc, 'phenotype': 'healthy', 'raw_label': f"{inf}|{status}"})
            keep.append(acc)
        elif "infection: denv" in inf_lower and ("status: df" in st_lower or "status: dhf" in st_lower):
            pheno_rows.append({'sample_id': acc, 'phenotype': 'dengue', 'raw_label': f"{inf}|{status}"})
            keep.append(acc)

    pheno_df = pd.DataFrame(pheno_rows)
    expr_df = expr_df[[col for col in keep if col in expr_df.columns]]
    common = [s for s in pheno_df['sample_id'] if s in expr_df.columns]
    expr_df = expr_df[common]
    pheno_df = pheno_df[pheno_df['sample_id'].isin(common)]

    print(f"  Samples: {pheno_df['phenotype'].value_counts().to_dict()}")
    if (pheno_df['phenotype'] == 'dengue').sum() < 10 or (pheno_df['phenotype'] == 'healthy').sum() < 5:
        print("  ERROR: Insufficient samples per group")
        return False

    n_ok = verify_markers(expr_df, pheno_df, 'interferon', 'dengue', 'healthy')
    if n_ok < 3:
        print("  WARNING: Fewer than 3 marker genes concordant")
    save_cohort("dengue_blood_gse51808", expr_df, pheno_df)
    return True


if __name__ == "__main__":
    os.makedirs(MATRICES_DIR, exist_ok=True)
    os.makedirs(PHENO_DIR, exist_ok=True)

    processors = [
        process_gse50772,
        process_gse49454,
        process_gse13355,
        process_gse14905,
        process_gse51808,
    ]

    results = {}
    for fn in processors:
        try:
            ok = fn()
            results[fn.__name__] = ok
        except Exception as e:
            print(f"  EXCEPTION in {fn.__name__}: {e}")
            import traceback
            traceback.print_exc()
            results[fn.__name__] = False

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, ok in results.items():
        print(f"  {name}: {'OK' if ok else 'FAILED'}")
    total_ok = sum(1 for v in results.values() if v)
    print(f"\nProcessed: {total_ok}/{len(results)} cohorts")
