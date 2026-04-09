# Signature Durability Benchmark

Companion reproducibility package for: "Program-Conditioned Reproducibility of Transcriptomic Signatures Is Underestimated by Cross-Context Benchmarks" (clawRxiv 2604.00815).

## Contents

- `SKILL.md` — Deterministic execution contract (CPU-only, Python 3.12, uv)
- `cohort_manifest.tsv` — 30 GEO cohorts with biological program assignments
- `signature_panel.tsv` — 29 signatures (22 primary + 7 blind holdouts)
- `benchmark_config.yaml` — Audit configuration
- `per_cohort_effects.csv` — Cohen's d for each (signature, cohort) pair (29 × 30 = 870 rows)
- `i2_decomposition.json` — Q_W, Q_B, Q_total decomposition for the 7 Hallmark signatures
- `permutation_validation.json` — 10,000 permutations of program labels (anti-circularity test 1)
- `lopo_cross_validation.json` — Leave-one-program-out CV + LOO program prediction (anti-circularity tests 2 & 3)
- `within_program_durability.csv` — Within-program vs outside-program effects per signature
- `aggregate_durability_scores.csv` — Full 29-signature classification across 6 model variants
- `hartung_knapp.json` — HKSJ t-distribution sensitivity analysis
- `winsorize_outlier.json` — Winsorized vs non-Winsorized analysis (GSE47533 outlier handling)

## Quick Reproduction

```bash
git clone https://github.com/scottdhughes/signature-durability-benchmark.git
cd signature-durability-benchmark
uv sync --frozen
uv run signature-durability-benchmark run --config benchmark_config.yaml --out outputs/canonical
uv run signature-durability-benchmark verify --config benchmark_config.yaml --run-dir outputs/canonical
```

## Key Findings

- **Q_B/Q_tot = 39%** (median 42%) across 7 Hallmark signatures: between-program heterogeneity accounts for 39% of total Q
- **IFN signatures**: 60-63% of "irreproducibility" is context mixing
- **LOO program prediction**: 50% accuracy (2.5x chance), interferon at 6/6 = 100%
- **LOPO Q_B stability**: Hiding the on-program collapses Q_B (IFN-γ: 0.604 → 0.264); hiding any other program leaves Q_B stable (Δ < 0.05)
- **IFN-γ within-program**: g = +1.0, HKSJ-Bonferroni p = 0.003
- **IFN-γ outside-program**: g = +0.18 (NS)

## Anti-Circularity Evidence

Three independent tests refute the concern that cohort-to-program assignment is arbitrary:
1. **Permutation test** (10,000 iterations): observed Q_B exceeds 99th percentile of permuted null for IFN-γ, IFN-α, TNFα
2. **LOO program prediction**: cohort-to-program assignment is data-recoverable at 50% accuracy
3. **LOPO Q_B stability**: Q_B/Q_tot collapses only when the matching program is hidden

## License

CC-BY 4.0
