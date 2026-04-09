---
name: signature-durability-benchmark
description: Score human gene signatures against 35 frozen GEO cohorts with program-conditioned meta-analysis and orthogonal interferon panel validation.
allowed-tools: Bash(uv *, python *, python3 *, ls *, test *, shasum *, tectonic *)
requires_python: "3.12.x"
package_manager: uv
repo_root: .
canonical_output_dir: outputs/canonical_v8
---

# Signature Durability Benchmark (Expanded Panel, v8)

This skill scores published gene signatures against 35 frozen real GEO expression cohorts (5,922 samples, 14 microarray platforms) covering five biological programs: interferon (k=11), inflammation (k=7), hypoxia (k=6), EMT (k=6), and proliferation (k=5). The expanded interferon arm (11 cohorts) spans viral infection (k=7), autoimmunity (SLE, psoriasis; k=4), and three tissues (whole blood, PBMC, skin).

The benchmark answers: is a transcriptomic signature reproducible **within** its native biological program, and does it correctly fail **outside** it? We use Cochran's Q decomposition (Q_within / Q_between), DerSimonian-Laird random-effects meta-analysis with Hartung-Knapp-Sidik-Jonkman guarded inference, 10,000-iteration permutation testing, and leave-one-program-out Q_B leverage diagnostics.

The interferon panel includes an **orthogonal Schoggins 2011 IRG panel** (76 genes from viral overexpression screens, 25% Hallmark overlap) and a **41-gene strictly-unique** Schoggins subset with zero overlap against the 9 marker genes used for expansion cohort admission — a direct test for curation-level circularity.

## Runtime Expectations

- Platform: CPU-only
- Python: 3.12.x
- Package manager: uv
- Offline after initial clone (all GEO data pre-frozen in `data/freeze/cohort_matrices/` and `data/freeze/cohort_phenotypes/`)
- Typical wall-clock for full rerun: ~5 minutes on a 2023 laptop

## Step 1: Install the Locked Environment

```bash
uv sync --frozen
```

## Step 2: Compute Per-Cohort Effects (Expanded Panel)

Scores all 35 cohorts × 30 signatures (1,050 pairs) using rank-based z-score normalization and Hedges' g with the small-sample correction. Writes `outputs/canonical_v8/per_cohort_effects.csv`.

```bash
uv run python scripts/compute_expanded_effects.py
```

Success condition: `outputs/canonical_v8/per_cohort_effects.csv` contains 1,050 rows (35 cohorts × 30 signatures).

## Step 3: Run All Meta-Analyses on the Expanded Panel

Runs I² decomposition, within-program DL+HKSJ meta-analysis, 10,000-iteration permutation test, and leave-one-program-out Q_B stability for 9 target signatures (Hallmark IFN-γ / IFN-α, Hallmark inflammatory / TNF-α-NFκB, Hallmark hypoxia, Hallmark E2F targets, Hallmark EMT, Schoggins 2011 IRG, blind IFN composite).

```bash
uv run python scripts/rerun_all_expanded.py
```

Success condition: all four of these JSONs exist in `outputs/canonical_v8/`:
- `i2_decomposition_expanded.json`
- `hartung_knapp_expanded.json`
- `permutation_validation_expanded.json`
- `lopo_cross_validation_expanded.json`

## Step 4: Orthogonal Schoggins Anti-Circularity Test

Scores the 41-gene strictly-unique Schoggins subset (zero overlap with Hallmark IFN-γ, Hallmark IFN-α, and the 9 marker genes used for cohort admission: STAT1, IRF1, IFIT1, IFIT2, ISG15, MX1, OAS1, GBP1, CXCL10, RSAD2). If this panel still produces a large within-IFN effect, curation-level circularity is falsified.

```bash
uv run python scripts/strictly_unique_schoggins.py
```

Success condition: `outputs/canonical_v8/strictly_unique_schoggins.json` shows:
- `expansion_marker_overlap`: 0
- `within_ifn_meta.pooled_g`: ≥ 1.2
- `within_ifn_meta.p_bonf`: < 0.001

## Step 5: Within-IFN Heterogeneity Meta-Regression

Decomposes residual within-IFN I² into contributions from tissue (blood/PBMC vs skin), etiology (viral vs autoimmune), and platform (Affymetrix vs Illumina). The combined tissue × etiology × platform model should explain ≥ 80% of within-IFN Q_total.

```bash
uv run python scripts/within_ifn_metaregression.py
```

Success condition: `outputs/canonical_v8/within_ifn_metaregression.json` shows `combined_r2.R2 > 0.8`.

## Step 6: Confirm Required Artifacts

Required files in `outputs/canonical_v8/`:
- `per_cohort_effects.csv` (1,050 rows)
- `i2_decomposition_expanded.json`
- `hartung_knapp_expanded.json`
- `permutation_validation_expanded.json`
- `lopo_cross_validation_expanded.json`
- `strictly_unique_schoggins.json`
- `within_ifn_metaregression.json`

## Expected Headline Results

| Signature | Within-IFN g | HKSJ-guarded p_Bonf (9 tests) | Permutation p | Outside-IFN g (NS) |
|-----------|-------------:|-------------------------------:|--------------:|-------------------:|
| Hallmark IFN-γ | +1.383 | 0.0049 | 0.0008 | +0.138 |
| Hallmark IFN-α | +1.458 | 0.0003 | 0.0004 | +0.219 |
| Schoggins 2011 IRG | +1.393 | 0.0006 | 0.0008 | +0.241 |
| Strictly-unique Schoggins (41g) | +1.247 | 0.00088 | — | +0.241 |
| Blind IFN composite | +1.402 | 0.0010 | 0.0007 | +0.196 |

## Scope Rules

- Human bulk transcriptomic signatures only
- No live data fetching in scored path
- Frozen GEO cohorts from real public data
- Blind panel is held out from all threshold-tuning decisions
- Source leakage between signature sources and cohort sources is forbidden
- LOPO Q_B is a leverage/sensitivity diagnostic, not an anti-circularity test (the tautology: hiding the driving program mechanically reduces Q_B). The permutation test and strictly-unique Schoggins are the primary anti-circularity evidence.
