---
name: signature-durability-benchmark
description: Score and triage human gene signatures against 35 frozen GEO cohorts with program-conditioned meta-analysis, held-out validation, external RNA-seq transfer tests, metadata-first prospective holdout auditing, an externally timestamped pending second prospective round, and orthogonal interferon panel checks.
allowed-tools: Bash(uv *, python *, python3 *, ls *, test *, shasum *, tectonic *)
requires_python: "3.12.x"
package_manager: uv
repo_root: .
canonical_output_dir: outputs/canonical_v8
---

# Signature Durability Benchmark (Expanded Panel, v8)

This skill scores and triages gene signatures against 35 frozen real GEO expression cohorts (5,922 samples, 14 microarray platforms) covering five biological programs: interferon (k=11), inflammation (k=7), hypoxia (k=6), EMT (k=6), and proliferation (k=5). The expanded interferon arm (11 cohorts) spans viral infection (k=7), autoimmunity (SLE, psoriasis; k=4), and three tissues (whole blood, PBMC, skin).

The benchmark answers: is a transcriptomic signature reproducible **within** its native biological program, and does it correctly fail **outside** it? We use Cochran's Q decomposition (Q_within / Q_between), DerSimonian-Laird random-effects meta-analysis with Hartung-Knapp-Sidik-Jonkman guarded inference, 10,000-iteration permutation testing, leave-one-program-out Q_B leverage diagnostics, held-out validation, a held-out external RNA-seq transfer layer, a deterministic second non-IFN breadth case, a metadata-first prospective holdout registry, an externally timestamped pending second prospective round, a generic locked-round evaluator for future readouts, a reusable `triage` interface for new signatures, and a provenance audit that verifies the active scored panel uses real GEO cohorts while quarantining the broader benchmark's explicit synthetic control signatures from the paper-facing target set.

The interferon panel includes an **orthogonal Schoggins 2011 IRG panel** (76 genes from viral overexpression screens, 25% overlap with the frozen Hallmark IFN-γ core) and a **41-gene strictly-unique** Schoggins subset with zero overlap against the 10 marker genes used for expansion cohort admission — a direct test against exact admission-marker reuse.

## Runtime Expectations

- Platform: CPU-only
- Python: 3.12.x
- Package manager: uv
- Offline after initial clone for the frozen benchmark; the v1 prospective holdout step fetches predeclared external GEO files into `data/prospective_holdout/downloads/`, and verifying the shipped v2 timestamp receipt is offline once those receipt files are present
- Typical wall-clock for full rerun: ~5 minutes on a 2023 laptop

## Step 1: Install the Locked Environment

```bash
uv sync --frozen
```

## Step 2: Compute Per-Cohort Effects (Expanded Panel)

Scores all 35 cohorts × 30 signatures (1,050 pairs) using per-gene z-score normalization across samples within cohort and Hedges' g with the small-sample correction. Writes `outputs/canonical_v8/per_cohort_effects.csv`.

```bash
uv run python scripts/compute_expanded_effects.py
```

Success condition: `outputs/canonical_v8/per_cohort_effects.csv` contains 1,050 rows (35 cohorts × 30 signatures).

## Step 3: Run All Meta-Analyses on the Expanded Panel

Runs I² decomposition, within-program DL+HKSJ meta-analysis, 10,000-iteration permutation test, and leave-one-program-out Q_B stability for the 9 paper target signatures enumerated in `config/paper_target_signatures.tsv` (Hallmark IFN-γ / IFN-α cores, Hallmark inflammatory / TNF-α-NFκB, Hallmark hypoxia, Hallmark E2F targets, Hallmark EMT, Schoggins 2011 IRG, blind IFN composite).

```bash
uv run python scripts/rerun_all_expanded.py
```

Success condition: all four of these JSONs exist in `outputs/canonical_v8/`:
- `i2_decomposition_expanded.json`
- `hartung_knapp_expanded.json`
- `permutation_validation_expanded.json`
- `lopo_cross_validation_expanded.json`

## Step 4: Orthogonal Schoggins Anti-Circularity Test

Scores the 41-gene strictly-unique Schoggins subset (zero overlap with the frozen Hallmark IFN-γ and IFN-α cores, and the 10 marker genes used for cohort admission: STAT1, IRF1, IFIT1, IFIT2, ISG15, MX1, OAS1, GBP1, CXCL10, RSAD2). If this panel still produces a large within-IFN effect, exact admission-marker reuse is unlikely to be the driver.

```bash
uv run python scripts/strictly_unique_schoggins.py
```

Success condition: `outputs/canonical_v8/strictly_unique_schoggins.json` shows:
- `expansion_marker_overlap`: 0
- `within_ifn_meta.p_bonf`: < 0.0056

## Step 5: Within-IFN Heterogeneity Decomposition

Descriptively decomposes residual within-IFN I² for the IFN-γ core into contributions from tissue (blood/PBMC vs skin), etiology (viral vs autoimmune), and platform (Affymetrix vs Illumina). The combined tissue × etiology × platform partition should explain ≥ 80% of within-IFN Q_total.

```bash
uv run python scripts/within_ifn_metaregression.py
```

Success condition: `outputs/canonical_v8/within_ifn_metaregression.json` shows `combined_r2.R2 > 0.8`.

## Step 6: Held-Out Validation

Evaluates leave-one-cohort-out prediction and split-half replication within each signature's home program.

```bash
uv run python scripts/external_validation.py
```

Success condition: `outputs/canonical_v8/external_validation.json` exists and reports:
- IFN quartet leave-one-out accuracy: 1.0
- IFN quartet split-half sign agreement: 1.0

## Step 7: Failure-Mode Analysis

Summarizes ambiguous comparator signatures where the diagnostic is better interpreted as exposing coarse program labels than a broken signature.

```bash
uv run python scripts/failure_mode_analysis.py
```

Success condition: `outputs/canonical_v8/failure_mode_analysis.json` exists and shows:
- inflammatory largest LOPO drop when hiding `proliferation`
- TNF-α/NFκB largest LOPO drop when hiding `proliferation`
- E2F `within_proliferation_effect_span.span > 3.0`

## Step 8: External RNA-seq Transfer Validation

Runs the held-out cross-platform extension on four primary bulk RNA-seq cohorts (GSE152641, GSE171110, GSE152075, GSE167000; 719 samples total) and reports one additional PBMC cohort (GSE152418) as an exploratory stress test.

```bash
uv run python scripts/external_rnaseq_validation.py
```

Success condition: `outputs/canonical_v8/external_rnaseq_validation.json` exists and shows:
- `primary_bulk_rnaseq_panel.n_cohorts`: 4
- `primary_bulk_rnaseq_panel.ifn_focus_summary.all_four_positive_in_all_primary_cohorts`: `true`
- `primary_bulk_rnaseq_panel.pooled_signatures.hallmark_ifng_response.guarded_p_bonf_7`: < 0.05
- `primary_bulk_rnaseq_panel.pooled_signatures.hallmark_ifna_response.guarded_p_bonf_7`: < 0.05

## Step 9: Prospective Holdout Audit

Runs the metadata-first held-out prediction round defined in `data/prospective_holdout/prediction_registry_v1.tsv`. This step audits a real predeclared registry rather than requiring the forecast itself to succeed.

```bash
uv run python scripts/prospective_holdout_prediction.py
```

Success condition: `outputs/canonical_v8/prospective_holdout_validation.json` exists and shows:
- `registry.sha256` matches the current SHA256 of `data/prospective_holdout/prediction_registry_v1.tsv`
- `download_audit` contains the newly fetched source files
- `per_cohort_effects_csv` points to `outputs/canonical_v8/prospective_holdout_per_cohort_effects.csv`

## Optional Step 10: Verify the Externally Timestamped v2 Declaration

Verifies or, if absent, creates the RFC3161 timestamp receipt for the fresh v2 prospective registry. This is intentionally separate from the scored v1 round.

```bash
uv run python -m signature_durability_benchmark.cli declare-prospective-round \
  --registry data/prospective_holdout/prediction_registry_v2.tsv \
  --protocol data/prospective_holdout/PREDICTION_PROTOCOL_v2.md
```

Success condition: `data/prospective_holdout/external_timestamps/prospective_holdout_v2/declaration_receipt.json` exists and shows:
- `round_id`: `prospective_holdout_v2`
- `verification_status`: `OK`
- `status`: either `created_new_receipt` or `verified_existing_receipt`

## Optional Step 11: Build the External Hypoxia Breadth Layer

Scores the exact-perturbation external hypoxia cohort selected under `data/external_hypoxia/SEARCH_PROTOCOL.md`.

```bash
uv run python scripts/external_hypoxia_validation.py
```

Success condition: `outputs/canonical_v8/external_hypoxia_validation.json` exists and shows:
- `chosen_signature_support.signature_id`: `hallmark_hypoxia`
- `chosen_signature_support.rank_among_scored_signatures`: `1`
- `chosen_signature_support.positive_direction`: `true`

## Optional Step 12: Build the Deterministic Second-Case Study

Ranks Hallmark Hypoxia versus Hallmark EMT from frozen artifacts only, reruns `triage` for the selected case, and records the supportive external layer.

```bash
uv run python scripts/generalization_case_study.py
```

Success condition: `outputs/canonical_v8/generalization_case_study.json` exists and shows:
- `selected_case.signature_id`: `hallmark_hypoxia`
- `selection_rule_satisfied`: `true`
- `selected_case.triage_best_program_matches_home`: `true`

## Optional Step 13: Generate the Workflow Figure

Builds the high-level method figure for the paper and repo.

```bash
uv run python scripts/generate_diagnostic_workflow_figure.py
```

Success condition:
- `paper/figure_workflow.pdf` exists
- `paper/figure_workflow.png` exists

## Optional Step 14: Score a Locked Future Prospective Round

This command is the future scoring path for v2 or later rounds. It must only be run against an already declared registry/protocol/receipt bundle.

```bash
uv run python -m signature_durability_benchmark.cli prospective-round-evaluate \
  --registry data/prospective_holdout/prediction_registry_v2.tsv \
  --protocol data/prospective_holdout/PREDICTION_PROTOCOL_v2.md \
  --receipt data/prospective_holdout/external_timestamps/prospective_holdout_v2/declaration_receipt.json \
  --out outputs/canonical_v8/prospective_rounds/prospective_holdout_v2
```

Success condition: if executed, the output directory contains:
- `evaluation.json`
- `evaluation_summary.md`
- `per_cohort_effects.csv`

## Optional Step 15: Build the Release-Ready Declaration Archive Bundle

Prepares the local declaration bundle for later GitHub Release / Zenodo publication.

```bash
uv run python scripts/build_archive_release_bundle.py
```

Success condition: `submission/archive_bundles/prospective_holdout_v2_declaration/` exists and contains:
- `prediction_registry_v2.tsv`
- `PREDICTION_PROTOCOL_v2.md`
- `prospective_holdout_v2/declaration_receipt.json`
- `CHECKSUMS.sha256`

## Optional Step 16: Build the Rescued Signature Portability Case

Runs `triage` on the fixed 6-gene Ayers IFN-gamma-related profile to show how the diagnostic can rescue a published external signature that would otherwise look heterogeneous in aggregate.

```bash
uv run python scripts/rescued_signature_case_study.py
```

Success condition: `outputs/canonical_v8/rescued_signature_case_study.json` exists and shows:
- `signature.name`: `Ayers IFN-gamma-related 6-gene profile`
- `triage.inferred_program`: `interferon`
- `triage.full_model_class`: `mixed`
- `triage.within_program_class`: `durable`

## Optional Step 17: Run the Provenance Audit

Verifies that the active scored cohort panel is the real 35-cohort GEO freeze, that the paper-facing target signatures are all non-synthetic, that no unexpected `synthetic` / `stub` / `mock` placeholders leak into the runtime or paper surface, and that the broader synthetic controls remain documented benchmark controls rather than headline evidence.

```bash
uv run python scripts/provenance_audit.py
```

Success condition: `outputs/canonical_v8/provenance_audit.json` exists and shows:
- `active_cohort_panel.active_cohorts`: `35`
- `active_cohort_panel.active_sample_sum`: `5922`
- `signature_panel.paper_target_all_non_synthetic`: `true`
- `keyword_audit.unexpected_hits`: `[]`

## Step 18: Confirm Required Artifacts

Required files in `outputs/canonical_v8/`:
- `per_cohort_effects.csv` (1,050 rows)
- `i2_decomposition_expanded.json`
- `hartung_knapp_expanded.json`
- `permutation_validation_expanded.json`
- `lopo_cross_validation_expanded.json`
- `strictly_unique_schoggins.json`
- `within_ifn_metaregression.json`
- `external_validation.json`
- `external_rnaseq_validation.json`
- `external_hypoxia_validation.json`
- `failure_mode_analysis.json`
- `generalization_case_study.json`
- `prospective_holdout_validation.json`
- `prospective_holdout_per_cohort_effects.csv`
- `rescued_signature_case_study.json`
- `provenance_audit.json`

## Expected Headline Results

| Signature | Within-IFN g | HKSJ-guarded p_Bonf (9 tests) | Permutation p | Outside-IFN g (NS) |
|-----------|-------------:|-------------------------------:|--------------:|-------------------:|
| Hallmark IFN-γ core | +1.383 | 0.0049 | 0.0008 | +0.138 |
| Hallmark IFN-α core | +1.458 | 0.0003 | 0.0004 | +0.219 |
| Schoggins 2011 IRG | +1.393 | 0.0006 | 0.0008 | +0.241 |
| Strictly-unique Schoggins (41g) | +1.247 | 0.00088 | — | +0.241 |
| Blind IFN composite | +1.402 | 0.0010 | 0.0007 | +0.196 |

Expected external RNA-seq pooled results in `external_rnaseq_validation.json` (4 primary bulk RNA-seq cohorts, 7 tested signatures):

| Signature | External RNA-seq pooled g | Guarded p_Bonf,7 | I² |
|-----------|--------------------------:|-----------------:|---:|
| Hallmark IFN-γ core | +0.922 | 0.022 | 0.000 |
| Hallmark IFN-α core | +1.193 | 0.011 | 0.000 |
| Schoggins 2011 IRG | +0.996 | 0.018 | 0.000 |
| Blind IFN composite | +1.177 | 0.020 | 0.244 |

Expected prospective audit structure in `prospective_holdout_validation.json`:

- `prediction_summary.round_id`: `prospective_holdout_v1`
- `prediction_summary.registry_sha256`: exact SHA256 of the frozen registry file
- `prediction_summary.per_cohort_hits`: 3 rows, one per predeclared cohort
- `pooled_primary_panel`: pooled IFN quartet plus comparator summaries for the predeclared panel

Expected externally timestamped v2 declaration receipt in `data/prospective_holdout/external_timestamps/prospective_holdout_v2/declaration_receipt.json`:

- `round_id`: `prospective_holdout_v2`
- `tsa_url`: `https://freetsa.org/tsr`
- `verification_status`: `OK`
- `cohorts`: 3 fresh predeclared rows (yellow fever PBMC, influenza blood, RSV nasal challenge)

## Scope Rules

- Human bulk transcriptomic signatures only
- No live data fetching in scored path
- Frozen GEO cohorts from real public data
- Blind panel is held out from all threshold-tuning decisions
- Source leakage between signature sources and cohort sources is forbidden
- LOPO Q_B is a leverage/sensitivity diagnostic, not an anti-circularity test (the tautology: hiding the driving program mechanically reduces Q_B). The permutation test and strictly-unique Schoggins are the primary anti-circularity evidence.

## Optional: Triage a New Signature

To diagnose an arbitrary input signature against the frozen panel, prepare a TSV/CSV with `gene_symbol` and optional `direction` / `weight` columns, then run:

```bash
uv run python -m signature_durability_benchmark.cli triage \
  --config config/benchmark_config.yaml \
  --input my_signature.tsv \
  --out outputs/my_signature
```

Success condition: `outputs/my_signature/diagnostic_summary.md` names a best-supported program and reports within-program versus outside-program pooled effects plus the permutation p-value for program structure.
