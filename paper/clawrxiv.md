# Program-Conditioned Diagnostic for Transcriptomic Signature Durability: Validation on Interferon Signatures across 35 Frozen GEO Cohorts

## Abstract

Gene-expression signatures are often labeled irreproducible when they fail across heterogeneous validation cohorts, but that failure can reflect either instability of the signature or mismatch between the signature and the test context. We present a program-conditioned diagnostic that scores a signature against a frozen reference panel, compares within-program versus outside-program effects, tests program structure by permutation, and surfaces failure modes when labels are too coarse. In 35 frozen GEO cohorts (5,922 samples, 5 biological programs, 14 microarray platforms), the frozen IFN-γ and IFN-α cores, an orthogonal 76-gene Schoggins panel, and a strictly-disjoint 41-gene Schoggins subset all show large within-IFN effects and small/non-significant outside-IFN effects; passing the IFN-γ core through `triage` yields a full-model class of `mixed` but a best-supported program of interferon and a within-program class of `durable`. Held-out validation is strong, and four external bulk RNA-seq cohorts (719 samples) reproduce the IFN quartet with guarded Bonferroni-significant pooled effects while inflammatory, TNF-α/NFκB, and E2F comparators remain non-significant. A deterministic breadth rule over Hallmark Hypoxia versus Hallmark EMT selects Hallmark Hypoxia as the stronger non-IFN example, but hypoxia honestly remains a `mixed` full-model case even though an external exact-perturbation cohort ranks it first. The prospective layer is likewise non-tautological: two predeclared v1 cohorts satisfy the 4/4 IFN sign forecast, whereas one severity-mixed acute PBMC cohort inverts all four IFN signatures. The repo now also ships a second fresh prospective registry whose protocol plus cohort list are bound to an external RFC3161 timestamp and intentionally left unevaluated here to establish an immutable, third-party-auditable future challenge. Comparator signatures show the complementary method result: inflammatory/TNF-α/NFκB ambiguity is driven more by proliferation than inflammation, and E2F targets expose a coarse label bucket rather than a random signature.

## Data and Diagnostic in Brief

- 35 frozen GEO cohorts partitioned into interferon (k = 11), inflammation (k = 7), proliferation (k = 5), hypoxia (k = 6), and EMT (k = 6)
- Primary IFN signatures: frozen 30-gene IFN-γ and IFN-α benchmark cores anchored to the broader MSigDB Hallmark families, a 76-gene Schoggins 2011 IRG panel, a strictly-disjoint 41-gene Schoggins subset, and a blind IFN composite
- Expansion cohorts admitted using concordant direction across 10 IFN markers: STAT1, IRF1, IFIT1, IFIT2, ISG15, MX1, OAS1, GBP1, CXCL10, RSAD2
- Per-cohort effect size: Hedges' g on weighted signed mean z-scored expression; within/outside pooling: guarded HKSJ random effects; structure test: 10,000-label permutation; additional analyses: held-out validation and failure-mode analysis
- Held-out external RNA-seq extension: primary bulk panel of GSE152641, GSE171110, GSE152075, and GSE167000 (719 samples total), with GSE152418 PBMC retained separately as an exploratory stress test
- Deterministic second breadth case: Hallmark Hypoxia versus Hallmark EMT, selected from frozen artifacts only, with supportive external exact-perturbation validation in GSE179885
- Metadata-first prospective holdout registry: GSE184610, GSE243217, and GSE202805 declared in `prediction_registry_v1.tsv` before held-out scoring
- Second fresh prospective registry: GSE243442, GSE213168, and GSE155237 declared in `prediction_registry_v2.tsv` and bound to an external RFC3161 timestamp receipt before evaluation
- Reusable interface: `triage` accepts a new signature TSV/CSV and writes `diagnostic.json`, `per_cohort_effects.csv`, and `diagnostic_summary.md`

## Diagnostic Workflow

The repo now includes a generated workflow figure at [paper/figure_workflow.png](https://github.com/scottdhughes/signature-durability-benchmark/blob/main/paper/figure_workflow.png), summarizing the packaged logic: input signature, frozen-panel scoring, within/outside separation, permutation, durable-versus-failure-mode branch, held-out external transfer, scored prospective v1 challenge, and the unevaluated externally timestamped v2 declaration.

## Key Results

| Signature | Within-IFN g | HKSJ-guarded p_Bonf | Outside-IFN g | Outside p | Permutation p |
|---|---:|---:|---:|---:|---:|
| IFN-γ core | +1.383 | 0.0049 | +0.138 | 0.484 | 0.0008 |
| IFN-α core | +1.458 | 0.0003 | +0.219 | 0.165 | 0.0004 |
| Schoggins 2011 IRG | +1.393 | 0.0006 | +0.241 | 0.178 | 0.0008 |
| Strictly-disjoint Schoggins (41 genes) | +1.247 | 0.00088 | +0.241 | 0.213 | — |
| Blind IFN composite | +1.402 | 0.0010 | +0.196 | 0.250 | 0.0007 |

All five IFN signatures survive guarded HKSJ plus 9-test Bonferroni correction within IFN cohorts and remain small/non-significant outside them. The strictly-disjoint Schoggins subset is the key curation-circularity result: it shares zero genes with the 10 cohort-admission markers and still reproduces the within-IFN effect, making exact admission-marker reuse an unlikely explanation. Zero overlap does not imply statistical independence, however, because the disjoint genes still live in the same co-regulated interferon program.

The diagnostic itself behaves as intended on a known durable signal. Treating the frozen IFN-γ core as an arbitrary input produces a full-model class of `mixed`, but `triage` infers interferon as the best-supported program and recovers the canonical within/outside result. The point is practical: a diluted pooled result need not imply a broken signature.

Held-out validation strengthens the method claim:

| Signature | LOO Sign Prediction | Split-Half Sign Agreement | Split-Half Both Significant |
|---|---:|---:|---:|
| IFN-γ core | 1.000 | 1.000 | 0.871 |
| IFN-α core | 1.000 | 1.000 | 1.000 |
| Schoggins 2011 IRG | 1.000 | 1.000 | 1.000 |
| Blind IFN composite | 1.000 | 1.000 | 1.000 |

Comparator signatures show the other use-case. Inflammatory response has Q_B/Q_tot = 0.342, but its largest LOPO drop occurs when hiding proliferation (Δ = 0.195), not inflammation. TNF-α/NFκB shows the same pattern (Q_B/Q_tot = 0.389, largest drop when hiding proliferation, Δ = 0.181). E2F targets show weak between-program structure (Q_B/Q_tot = 0.061), only 80% leave-one-out sign prediction, 75% split-half sign agreement, and a within-proliferation effect span from −1.21 to +2.14. These are better interpreted as boundary failures in the labeling scheme than as simple signature collapse. Operationally, when LOPO leverage is driven by a foreign program, the next step is to split or relabel that program and rerun; when within-program spread stays extreme despite weak between-program structure, the next step is to stratify by contrast type or tissue before concluding that the signature itself is unstable.

The external RNA-seq extension now provides a real cross-platform test rather than just a future-work promise:

| Signature | External RNA-seq pooled g | Guarded p_Bonf,7 | I² |
|---|---:|---:|---:|
| IFN-γ core | +0.922 | 0.022 | 0.000 |
| IFN-α core | +1.193 | 0.011 | 0.000 |
| Schoggins 2011 IRG | +0.996 | 0.018 | 0.000 |
| Blind IFN composite | +1.177 | 0.020 | 0.244 |
| Inflammatory Response | +0.427 | 1.000 | 0.770 |
| TNF-α/NFκB | +0.382 | 1.000 | 0.859 |
| E2F Targets | +0.818 | 1.000 | 0.944 |

All four IFN signatures are positive in all four primary external cohorts. The additional PBMC cohort (GSE152418) is reported separately as an exploratory stress test because its cell-selected composition produces a much more mixed pattern, which is informative but not a clean bulk RNA-seq transfer test.

The second positive breadth case is intentionally less clean than IFN. Hallmark Hypoxia beats Hallmark EMT under the frozen deterministic selection rule (p_Bonf = 0.024 versus 1.0, with perfect LOO and split-half sign agreement for both), and `triage` still recovers hypoxia as the best-supported home program. But hypoxia remains `mixed` in the full-model diagnostic because outside-program bleed-through is not negligible. That is useful method behavior, not a bug: the framework broadens beyond IFN without pretending every program has IFN-like specificity.

The external exact-perturbation cohort pushes in the same direction. In GSE179885 human T-cell RNA-seq cultured under hypoxia versus normoxia, Hallmark Hypoxia ranks first among seven scored signatures with Hedges' g = +5.11 and full coverage, ahead of Hallmark EMT (g = +0.94), while IFN, inflammatory, TNF-α/NFκB, and E2F comparators are all negative. Because this cohort has 12 total samples, we treat it as supportive breadth evidence rather than a headline pooled result.

The diagnostic is also useful on a published signature imported from outside this benchmark. Running the Ayers IFN-γ-related 6-gene clinical response profile (IDO1, CXCL10, CXCL9, HLA-DRA, STAT1, IFNG; Ayers et al. 2017, doi:10.1172/JCI91190) through `triage` yields exactly the rescue pattern the method is meant to expose: the aggregate profile is `mixed` (effect = +0.246, I² = 0.929), but the best-supported program is interferon and the within-program class is `durable`. The within-interferon pooled effect is +0.866 (p = 0.0245), while the outside-interferon effect is essentially null (−0.010, p = 0.949). This is a portability demonstration rather than a reanalysis of the original oncology cohorts, but it shows that a published external signature can look diluted in a broad benchmark and still resolve into the correct biological home program.

The first metadata-first prospective round is intentionally harder and mixed. GSE184610 and GSE243217 are 4/4 positive across the IFN quartet, but GSE202805 is 0/4 positive, so the pooled prospective rule fails overall (pooled g = +0.359 to +0.657, Bonferroni-4 p = 0.93 to 1.0, sign consistency = 0.667). This does not strengthen the paper as another validation win; it strengthens it as an honest diagnostic benchmark that can miss on metadata-first forecasts. The miss is also interpretable: GSE202805 is a severity-mixed acute PBMC series, and in that cohort the IFN quartet is negative while inflammatory and E2F comparators are positive, which is more consistent with compartment/severity mismatch than with a broken scorer. The repo now locks in the next prospective step as well: a fresh v2 round (yellow fever PBMC, influenza blood, RSV nasal challenge) is externally timestamped and intentionally not yet scored.

## Limits

- The frozen decision surface is still defined by microarray cohorts, even though the paper now includes a held-out four-cohort bulk RNA-seq extension
- The 5-program partition is operational rather than ontologically complete
- The descriptive IFN heterogeneity decomposition is based on only 11 IFN cohorts
- The strictly-disjoint Schoggins holdout rules out exact admission-marker reuse more cleanly than every possible cohort-selection effect
- Held-out validation is no longer only internal to the frozen panel; the repo now includes a first metadata-first prospective holdout round, and that round is mixed rather than uniformly confirmatory
- A second fresh prospective round is externally timestamped and pending evaluation, so it strengthens the audit trail rather than the current numeric result tables

## Reproducibility

Canonical outputs for this version are in `outputs/canonical_v8/`.

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
uv run python scripts/prospective_holdout_prediction.py
uv run python scripts/generate_diagnostic_workflow_figure.py
```

Verify the externally timestamped v2 declaration with:

```bash
uv run python -m signature_durability_benchmark.cli declare-prospective-round \
  --registry data/prospective_holdout/prediction_registry_v2.tsv \
  --protocol data/prospective_holdout/PREDICTION_PROTOCOL_v2.md
```

The future locked v2 readout is now a first-class command:

```bash
uv run python -m signature_durability_benchmark.cli prospective-round-evaluate \
  --registry data/prospective_holdout/prediction_registry_v2.tsv \
  --protocol data/prospective_holdout/PREDICTION_PROTOCOL_v2.md \
  --receipt data/prospective_holdout/external_timestamps/prospective_holdout_v2/declaration_receipt.json \
  --out outputs/canonical_v8/prospective_rounds/prospective_holdout_v2
```

Triage a new signature with:

```bash
uv run python -m signature_durability_benchmark.cli triage \
  --config config/benchmark_config.yaml \
  --input my_signature.tsv \
  --out outputs/my_signature
```

## References

1. Schoggins JW. Interferon-stimulated genes: what do they all do? Annu Rev Virol. 2019;6:567-584. doi:10.1146/annurev-virology-092818-015756
2. Schoggins JW, Wilson SJ, Panis M, et al. A diverse range of gene products are effectors of the type I interferon antiviral response. Nature. 2011;472:481-485. doi:10.1038/nature09907
3. Liberzon A, et al. The Molecular Signatures Database Hallmark gene set collection. Cell Syst. 2015;1:417-425. doi:10.1016/j.cels.2015.12.004
4. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7:177-188. doi:10.1016/0197-2456(86)90046-2
5. Hartung J, Knapp G. On tests of the overall treatment effect in meta-analysis with normally distributed responses. Stat Med. 2001;20:1771-1782. doi:10.1002/sim.791
6. Sidik K, Jonkman JN. A simple confidence interval for meta-analysis. Stat Med. 2002;21:3153-3159. doi:10.1002/sim.1262
7. Baechler EC, Batliwalla FM, Karypis G, et al. Interferon-inducible gene expression signature in peripheral blood cells of patients with severe lupus. Proc Natl Acad Sci U S A. 2003;100:2610-2615. doi:10.1073/pnas.0337679100
8. Yao Y, Richman L, Morehouse C, et al. Type I interferon: potential therapeutic target for psoriasis? PLoS One. 2008;3:e2737. doi:10.1371/journal.pone.0002737
9. Ayers M, Lunceford J, Nebozhyn M, et al. IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127:2930-2940. doi:10.1172/JCI91190
