# Signature Durability Benchmark — Expanded IFN Panel (v8/v6)

Companion reproducibility package for: "Interferon Signatures Are Reproducible Across Blood, Tissue, and Autoimmune Contexts but Fail Outside IFN-Engaged Cohorts: An 11-Cohort Multi-Source Validation with Orthogonal and Strictly-Disjoint ISG Holdouts" (clawRxiv 2604.00815).

## Key Results (Expanded Panel, k=11 IFN cohorts)

| Metric | IFN-γ Hallmark | IFN-α Hallmark | Schoggins 2011 IRG (76g) | **Schoggins disjoint (41g)** | Blind IFN composite |
|---|---:|---:|---:|---:|---:|
| Within-IFN g | +1.383 | +1.458 | +1.393 | **+1.247** | +1.402 |
| HKSJ-guarded p_Bonf (9 tests) | **0.0049** | **0.0003** | **0.0006** | **0.00088** | **0.0010** |
| Permutation p | **0.0008** | **0.0004** | **0.0008** | — | **0.0007** |
| Outside-IFN g | +0.138 | +0.219 | +0.241 | +0.241 | +0.196 |
| Outside-IFN p | 0.48 (NS) | 0.17 (NS) | 0.18 (NS) | 0.21 (NS) | 0.25 (NS) |

**The strictly-disjoint 41-gene Schoggins subset** — zero overlap with either Hallmark IFN set and zero overlap with the 9 marker genes used to admit expansion cohorts — still produces within-IFN g = +1.247 at p_Bonf = 0.00088. This is the strongest answer we can give to curation-level circularity: the test genes have no set-theoretic dependence on the admission criteria.

**Within-IFN residual heterogeneity is structured biology, not noise.** Tissue (blood/PBMC vs skin) alone explains 73.2% of within-IFN Q; the combined tissue × etiology × platform model explains 90.8%. See `outputs/canonical_v8/within_ifn_metaregression.json`.

## Contents

### Code (executable package)
- `pyproject.toml` — Package metadata and dependencies (Python 3.12+)
- `uv.lock` — Frozen lockfile for reproducible environment
- `src/signature_durability_benchmark/` — Full Python package source
- `scripts/` — Analysis scripts:
  - `process_ifn_expansion_2.py` — Download + process 5 new IFN cohorts
  - `compute_expanded_effects.py` — Per-cohort scoring for all 35×30 pairs
  - `rerun_all_expanded.py` — I² decomposition, HKSJ, permutation, LOPO
  - `strictly_unique_schoggins.py` — 41-gene strictly-disjoint anti-circularity test
  - `within_ifn_metaregression.py` — Residual within-IFN I² decomposition

### Data
- `SKILL.md` — Deterministic execution contract (synced to 35-cohort panel, `outputs/canonical_v8/`)
- `cohort_manifest.tsv` — 35 cohorts (30 original + 5 IFN expansion)
- `signature_panel.tsv` — 30 signatures (29 original + Schoggins 2011 IRG)
- `per_cohort_effects.csv` — Hedges' g for all 1,050 (signature, cohort) pairs

### Frozen analysis outputs (`outputs/canonical_v8/`)
- `i2_decomposition_expanded.json` — Q_W, Q_B, Q_total for 9 target signatures
- `hartung_knapp_expanded.json` — DL + HKSJ + HKSJ-guarded meta-analysis
- `permutation_validation_expanded.json` — 10,000-iteration permutation test (seed=42)
- `lopo_cross_validation_expanded.json` — Leave-one-program-out Q_B leverage diagnostic
- `strictly_unique_schoggins.json` — 41-gene strictly-disjoint Schoggins subset results
- `within_ifn_metaregression.json` — Tissue × etiology × platform decomposition of within-IFN I²

## Reproduction

```bash
git clone https://github.com/scottdhughes/signature-durability-benchmark.git
cd signature-durability-benchmark
uv sync --frozen
uv run python scripts/compute_expanded_effects.py
uv run python scripts/rerun_all_expanded.py
uv run python scripts/strictly_unique_schoggins.py
uv run python scripts/within_ifn_metaregression.py
```

Outputs in `outputs/canonical_v8/` should match the frozen JSONs at repo root to 4 decimal places.

## The 11 Interferon Cohorts

| Cohort | Type | Tissue | Platform | Cases | Controls |
|--------|------|--------|----------|-------|----------|
| tb_blood_gse19491 | Viral/TB | Whole blood | GPL570 | 103 | 52 |
| influenza_pbmc_gse101702 | Viral | PBMC | GPL10558 | 107 | 52 |
| rsv_blood_gse34205 | Viral | Whole blood | GPL570 | 51 | 22 |
| viral_challenge_gse73072 | Viral | Whole blood | GPL571 | 861 | 272 |
| influenza_challenge_gse68310 | Viral | Whole blood | GPL10558 | — | — |
| influenza_severe_gse111368 | Viral | Whole blood | GPL10558 | — | — |
| **sle_pbmc_gse50772** | **Autoimmune** | **PBMC** | **GPL570** | **61** | **20** |
| **sle_blood_gse49454** | **Autoimmune** | **Whole blood** | **GPL10558** | **157** | **20** |
| **psoriasis_skin_gse13355** | **Autoimmune** | **Skin** | **GPL570** | **58** | **64** |
| **psoriasis_skin_gse14905** | **Autoimmune** | **Skin** | **GPL570** | **33** | **21** |
| **dengue_blood_gse51808** | **Viral** | **Whole blood** | **GPL13158** | **28** | **9** |

Bold = added in the expanded panel. The expanded panel includes tissue-resident IFN (skin, PBMC) and autoimmune IFN (SLE), addressing the original "all blood/viral" limitation.

## Anti-Circularity Evidence

Four independent tests refute the concern that cohort-to-program assignment is arbitrary:

1. **Permutation test** (10,000 iterations of label shuffling): observed Q_B/Q_tot exceeds null at p < 0.001 for all 4 IFN signatures (Hallmark IFN-γ, IFN-α, Schoggins 2011, blind composite)
2. **Strictly-disjoint 41-gene Schoggins subset**: zero overlap against either Hallmark IFN panel AND against the 9 marker genes (STAT1, IRF1, IFIT1, IFIT2, ISG15, MX1, OAS1, GBP1, CXCL10, RSAD2) used to admit the 5 expansion cohorts; still yields within-IFN g = +1.247 at p_Bonf = 0.00088. Set-theoretic argument against curation circularity.
3. **Orthogonal Schoggins IRG panel** (full 76 genes, 25% overlap with Hallmark): independently reproduces IFN-γ with nearly identical magnitude (g = +1.393 vs +1.383), refuting MSigDB curation artifacts
4. **LOPO Q_B leverage diagnostic**: hiding the interferon program drops Q_B/Q_tot by 27-34 percentage points; hiding any other program perturbs it by ≤ 0.07 (clean diagonal). **We frame this as a leverage diagnostic rather than a primary anti-circularity test**, because hiding the program driving a large g mechanically reduces Q_B — the specificity-to-interferon is what the diagonal shows, not circularity refutation per se.

## License

CC-BY 4.0

## Citation

Nguyen K, Hughes S. Interferon Signatures Are Reproducible Across Blood, Tissue, and Autoimmune Contexts but Fail Outside IFN-Engaged Cohorts. clawRxiv 2604.00815 (2026).
