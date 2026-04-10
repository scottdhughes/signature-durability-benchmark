# Signature Durability Benchmark

Executable benchmark and **program-conditioned diagnostic** for transcriptomic signatures across 35 frozen GEO cohorts (5,922 samples, 14 microarray platforms). The current validation paper uses interferon as the strongest reference case, but the package is meant to answer the broader question: when a signature looks heterogeneous across cohorts, is the signature broken, or is the cohort context mismatched?

The diagnostic does three things:

- scores a signature across the frozen panel
- compares **within-program durability** against **outside-program mismatch**
- flags **failure modes** where the label system itself appears biologically coarse

## Canonical Reproduction

```bash
uv sync --frozen
uv run python scripts/compute_expanded_effects.py
uv run python scripts/rerun_all_expanded.py
uv run python scripts/strictly_unique_schoggins.py
uv run python scripts/within_ifn_metaregression.py
uv run python scripts/external_validation.py
uv run python scripts/external_rnaseq_validation.py
uv run python scripts/external_hypoxia_validation.py
uv run python scripts/failure_mode_analysis.py
uv run python scripts/generalization_case_study.py
uv run python scripts/rescued_signature_case_study.py
uv run python scripts/prospective_holdout_prediction.py
uv run python scripts/generate_diagnostic_workflow_figure.py
uv run python scripts/build_archive_release_bundle.py
```

Canonical outputs land in `outputs/canonical_v8/`.

## Prospective Rounds

- `v1`: scored, mixed, and preserved as the honest evaluated result in `outputs/canonical_v8/prospective_holdout_validation.json`
- `v2`: fresh, externally timestamped, and intentionally pending evaluation in `data/prospective_holdout/prediction_registry_v2.tsv`

Verify the RFC3161 declaration receipt for v2 with:

```bash
uv run python -m signature_durability_benchmark.cli declare-prospective-round \
  --registry data/prospective_holdout/prediction_registry_v2.tsv \
  --protocol data/prospective_holdout/PREDICTION_PROTOCOL_v2.md
```

The declaration bundle lives in `data/prospective_holdout/external_timestamps/prospective_holdout_v2/` and binds the registry plus protocol hashes to a third-party timestamp.

When the paper is ready for the locked v2 readout, score it with:

```bash
uv run python -m signature_durability_benchmark.cli prospective-round-evaluate \
  --registry data/prospective_holdout/prediction_registry_v2.tsv \
  --protocol data/prospective_holdout/PREDICTION_PROTOCOL_v2.md \
  --receipt data/prospective_holdout/external_timestamps/prospective_holdout_v2/declaration_receipt.json \
  --out outputs/canonical_v8/prospective_rounds/prospective_holdout_v2
```

The release-ready declaration bundle for GitHub Releases / Zenodo is built locally at `submission/archive_bundles/prospective_holdout_v2_declaration/`.

## Triage a New Signature

Provide a TSV or CSV containing `gene_symbol` and optional `direction` / `weight` columns:

```bash
uv run python -m signature_durability_benchmark.cli triage \
  --config config/benchmark_config.yaml \
  --input my_signature.tsv \
  --out outputs/my_signature
```

The triage run writes:

- `diagnostic.json`
- `diagnostic_summary.md`
- `per_cohort_effects.csv`

The summary reports the best-supported program, within-program and outside-program pooled effects, permutation support for program structure, and full-model versus within-program class labels.

## Breadth Beyond IFN

The repo now carries one deterministic second breadth case in addition to the primary IFN validation:

- Hallmark Hypoxia is selected over Hallmark EMT by a frozen ranking rule using current benchmark artifacts only
- `outputs/canonical_v8/generalization_case_study.json` records the full ranking, fresh triage report, home-program cohorts, and largest foreign-program effects
- `outputs/canonical_v8/external_hypoxia_validation.json` adds a supportive exact-perturbation external cohort (GSE179885), where Hallmark Hypoxia ranks first with Hedges' g = +5.11

This breadth case is intentionally not as clean as IFN. The method recovers hypoxia as the best-supported program while still showing that foreign-program bleed-through keeps the full-model result mixed.

## Rescued Signature Case

The repo also includes a portability demo on a published external signature rather than a benchmark-native panel:

- `data/case_studies/ayers_ifng6_signature.tsv` encodes the 6-gene Ayers IFN-gamma-related profile (IDO1, CXCL10, CXCL9, HLA-DRA, STAT1, IFNG)
- `scripts/rescued_signature_case_study.py` reruns `triage` on that fixed signature against the frozen 35-cohort panel
- `outputs/canonical_v8/rescued_signature_case_study.json` and `.md` record the result: aggregate `mixed`, best-supported program `interferon`, within-program class `durable`, within-IFN positive, and outside-IFN near null

This is intentionally framed as a portability demonstration, not a reanalysis of the original oncology response cohorts. The point is that a published IFN-gamma response signature that would look heterogeneous in the aggregate is recoverable as context-matched once triaged against the frozen panel.

## Workflow Figure

`paper/figure_workflow.pdf` and `paper/figure_workflow.png` now provide the high-level packaged view of the method:

- input signature
- frozen-panel scoring
- within/outside separation
- permutation check
- durable-versus-failure-mode branch
- held-out external transfer
- prospective v1 scored challenge
- prospective v2 declared future challenge

## Current Scope

- 35 frozen GEO cohorts across interferon, inflammation, proliferation, hypoxia, and EMT
- 30 frozen benchmark signatures scored against all 35 cohorts (1,050 cohort-signature pairs)
- Primary IFN validation on 11 IFN-engaged cohorts spanning viral infection, SLE, and psoriasis
- Orthogonal validation with the 76-gene Schoggins 2011 IRG panel
- Strictly-disjoint 41-gene Schoggins holdout with zero overlap against the 10 IFN cohort-admission markers
- Held-out validation via leave-one-cohort-out and split-half replication
- A rescued-signature portability case study using the published Ayers IFN-gamma-related 6-gene profile
- Held-out external bulk RNA-seq validation across 4 GEO cohorts (719 samples) plus an exploratory PBMC stress test
- Supportive external exact-perturbation hypoxia validation on GSE179885
- Metadata-first prospective v1 holdout registry plus auditable evaluation outputs
- A second fresh prospective v2 registry that is externally timestamped and kept separate from v1 pending evaluation
- Generic locked-round evaluator for scoring future prospective registries without changing declared inputs
- Failure-mode analysis for ambiguous comparator programs
- Release metadata for public archival publication (`LICENSE`, `CITATION.cff`, `.zenodo.json`)
- CPU-only, deterministic, and offline after clone
