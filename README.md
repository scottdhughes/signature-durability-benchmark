# Signature Durability Benchmark

Cross-cohort benchmark that scores published gene signatures against 12 frozen real GEO expression cohorts to determine whether each signature generalizes (durable), fails (brittle), is partially consistent (mixed), is driven by confounding (confounded), or lacks gene coverage (insufficient_coverage).

## Canonical Commands

```bash
uv sync --frozen
uv run --frozen --no-sync signature-durability-benchmark build-freeze --config config/benchmark_config.yaml --out data/freeze
uv run --frozen --no-sync signature-durability-benchmark run --config config/benchmark_config.yaml --out outputs/canonical
uv run --frozen --no-sync signature-durability-benchmark verify --config config/benchmark_config.yaml --run-dir outputs/canonical
```

## Design Rules

- CPU-only and offline after clone
- All GEO data pre-frozen as gene-by-sample TSV matrices
- Scored path uses no network access
- 29 signatures across 6 biological programs
- 12 real GEO cohorts (3,003 samples, 3 microarray platforms)
- 5 classification models compared against pre-registered success rule
- 4 machine-readable certificates
