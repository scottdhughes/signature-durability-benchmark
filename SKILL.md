---
name: signature-durability-benchmark
description: Score human gene signatures against frozen real GEO cohorts to determine cross-cohort transcriptomic durability with self-verification and confounder rejection.
allowed-tools: Bash(uv *, python *, python3 *, ls *, test *, shasum *, tectonic *)
requires_python: "3.12.x"
package_manager: uv
repo_root: .
canonical_output_dir: outputs/canonical
---

# Signature Durability Benchmark

This skill scores published gene signatures against 22 frozen real GEO expression cohorts (4,730 samples, 3 microarray platforms) to determine whether each signature is durable, brittle, mixed, confounded, or insufficiently covered across independent cohorts. The full model is compared against 4 baselines with a success rule.

## Runtime Expectations

- Platform: CPU-only
- Python: 3.12.x
- Package manager: uv
- Offline after initial clone (all GEO data pre-frozen)

## Step 1: Install the Locked Environment

```bash
uv sync --frozen
```

## Step 2: Build Freeze (Validate Frozen Assets)

```bash
uv run --frozen --no-sync signature-durability-benchmark build-freeze --config config/benchmark_config.yaml --out data/freeze
```

Success condition: freeze_audit.json shows valid=true

## Step 3: Run the Canonical Benchmark

```bash
uv run --frozen --no-sync signature-durability-benchmark run --config config/benchmark_config.yaml --out outputs/canonical
```

Success condition: outputs/canonical/manifest.json exists

## Step 4: Verify the Run

```bash
uv run --frozen --no-sync signature-durability-benchmark verify --config config/benchmark_config.yaml --run-dir outputs/canonical
```

Success condition: verification status is passed

## Step 5: Confirm Required Artifacts

Required files in outputs/canonical/:
- manifest.json
- normalization_audit.json
- cohort_overlap_summary.csv
- per_cohort_effects.csv
- aggregate_durability_scores.csv
- matched_null_summary.csv
- leave_one_cohort_out.csv
- platform_holdout_summary.csv
- durability_certificate.json
- platform_transfer_certificate.json
- confounder_rejection_certificate.json
- coverage_certificate.json
- benchmark_protocol.json
- verification.json
- public_summary.md
- within_program_durability.csv
- forest_plot.png
- null_separation_plot.png
- stability_heatmap.png
- platform_transfer_panel.png

## Scope Rules

- Human bulk transcriptomic signatures only
- No live data fetching in scored path
- Frozen GEO cohorts from real public data
- Blind panel never influences thresholds
- Source leakage between signature sources and cohort sources is forbidden
