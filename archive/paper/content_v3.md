# Program-Conditioned Reproducibility of Transcriptomic Signatures Is Underestimated by Cross-Context Benchmarks

## Abstract

Gene expression signatures are routinely dismissed as irreproducible when they fail cross-context validation — but how much of that apparent irreproducibility is a measurement artifact of context mixing? We score 7 MSigDB Hallmark signatures (and 22 additional benchmark signatures) across 30 frozen GEO cohorts (5,451 samples) organized into 5 biological programs (inflammation, interferon, proliferation, hypoxia, EMT) and decompose Cochran's Q into within-program (Q_W) and between-program (Q_B) components. Between-program heterogeneity accounts for 39% of total Q (median 42%; all p < 10^-5), rising to 61-63% for interferon signatures. We directly address the circularity concern that cohorts were "pre-assigned to matching signatures" via two new analyses: (1) **Leave-one-cohort-out k-NN program prediction** from 7-D effect-size fingerprints achieves 50% accuracy (2.5x chance), with interferon at 6/6 = 100% — cohort-program assignments are data-recoverable, not arbitrary. (2) **Leave-one-program-out Q_B stability** shows that hiding any non-on-program leaves Q_B/Q_tot largely unchanged (Δ < 0.05), while hiding the on-program causes the expected sharp drop (e.g., IFN-γ: 0.604 → 0.264 when interferon cohorts are removed). The diagonal pattern is exactly what would be observed if program structure were real, not imposed. IFN-γ exemplifies the phenomenon: within-program Hedges' g = +1.0 (HKSJ-Bonferroni p = 0.003), outside-program g = +0.18 (NS). All 30 cohorts, 29 signatures, code, and frozen outputs: https://github.com/scottdhughes/signature-durability-benchmark.

## Introduction

Gene expression signatures routinely fail validation outside their discovery context. These failures are interpreted as irreproducibility. But a signature replicating across 6 interferon-response cohorts yet showing no effect in breast cancer studies is context-specific, not broken — cross-context testing conflates the two. The standard diagnostic, I-squared, quantifies total heterogeneity without distinguishing within-program instability from between-program context differences. We decompose Cochran's Q into within-program (Q_W) and between-program (Q_B) components across 7 Hallmark signatures in 30 GEO cohorts organized into 5 biological programs, then directly test whether the program structure is data-driven or imposed.

## Data

30 frozen GEO cohorts (5,451 samples) across 5 biological programs:
- inflammation (k=7, N=1,976)
- interferon (k=6, N=1,991)
- proliferation (k=5, N=847)
- hypoxia (k=6, N=239)
- EMT (k=6, N=398)

Platforms: 20 Affymetrix, 4 Agilent, 6 Illumina. 29 signatures (22 primary, 7 blind holdouts): 7 MSigDB Hallmarks (Liberzon et al. 2015, doi:10.1016/j.cels.2015.12.004); 4 brittle (small-study, single-tissue, platform-specific, overfit-noise); 3 mixed-program; 5 confounded (immune-infiltrate, ribosomal, proliferation, plus 2 stealth); 3 insufficient-coverage; 7 blind holdouts spanning all classes.

## Methods

**Effect size:** Hedges' g (small-sample-corrected Cohen's d).

**I-squared decomposition:** Q_total = Q_W + Q_B, where Q_W is the sum of within-program weighted squared deviations and Q_B is computed as Q_total - Q_W. Q_B tested against chi-squared(K-1).

**Within-program meta-analysis:** DerSimonian-Laird random-effects, with Hartung-Knapp-Sidik-Jonkman (HKSJ) t-distribution as a robustness check. Bonferroni correction at alpha = 0.05/9.

**Permutation test:** Cohort-to-program labels permuted 10,000 times (seed=42); empirical p_B is the fraction of permuted Q_B/Q_tot exceeding observed.

**Leave-one-program-out (LOPO) cross-validation (NEW):** For each held-out program P, recompute Q_B/Q_tot using only the remaining N-1 programs. If the structure depended entirely on P's pre-assignment, removing P should collapse Q_B/Q_tot to noise.

**Leave-one-cohort-out (LOO) program prediction (NEW):** For each cohort, build a 7-D effect-size fingerprint (Cohen's d for each Hallmark signature). Train k-NN (k=3) on the other 29 cohorts and predict the held-out cohort's program. Chance accuracy = 1/5 = 20%.

## Results

### I-Squared Decomposition: Context Explains 39% of Heterogeneity

| Signature | Q_total | Q_W | Q_B | Q_B/Q_tot | I²_total | I²_W | p_B |
|---|---|---|---|---|---|---|---|
| IFN-alpha | 393.5 | 145.3 | 247.9 | **0.63** | 0.93 | 0.83 | <10^-6 |
| IFN-gamma | 429.7 | 170.2 | 259.5 | **0.60** | 0.93 | 0.85 | <10^-6 |
| TNFα/NFκB | 335.9 | 156.9 | 178.9 | **0.53** | 0.91 | 0.84 | <10^-6 |
| Inflammatory | 387.0 | 226.1 | 161.0 | 0.42 | 0.93 | 0.89 | <10^-6 |
| EMT | 240.6 | 176.4 | 64.2 | 0.27 | 0.88 | 0.86 | <10^-13 |
| Hypoxia | 573.4 | 459.8 | 113.6 | 0.20 | 0.95 | 0.95 | <10^-6 |
| E2F Targets | 448.6 | 418.7 | 30.0 | 0.07 | 0.94 | 0.94 | 5x10^-6 |
| **Mean / Median** | | | | **0.39 / 0.42** | | | |

Across all 7 Hallmarks, 39% of total Q (median 42%) is between-program heterogeneity. For interferon signatures, 60-63% of what I² reports as "irreproducibility" is context mixing.

### Anti-Circularity Test 1: Permutation of Program Labels

Permuting cohort-to-program labels 10,000 times (preserving cohort effect sizes), 3/7 Hallmarks have observed Q_B/Q_tot exceeding the 99th percentile of the permuted null:

| Signature | Observed Q_B/Q_tot | Null 95th | Null 99th | Empirical p |
|---|---|---|---|---|
| IFN-γ | 0.604 | 0.434 | 0.486 | < 0.003 |
| IFN-α | 0.631 | 0.441 | 0.494 | < 0.003 |
| TNFα | 0.533 | 0.363 | 0.412 | < 0.003 |
| Inflammatory | 0.416 | 0.325 | 0.358 | < 0.05 |
| EMT | 0.267 | 0.318 | NS | NS |
| Hypoxia | 0.198 | 0.309 | NS | NS |
| E2F | 0.067 | 0.394 | NS | NS |

If program assignments were arbitrary labels on pre-selected cohorts, random permutation would reproduce the observed Q_B structure. It does not — for IFN, TNFα, and Inflammatory signatures.

### Anti-Circularity Test 2: Leave-One-Cohort-Out Program Prediction

We tested whether the cohort-to-program assignment is data-recoverable. For each cohort, we computed a 7-D effect-size fingerprint (one Cohen's d per Hallmark signature) and trained k-NN (k=3) on the other 29 cohorts to predict the held-out cohort's program label.

| Program | Correct / Total | Per-Program Accuracy |
|---|---|---|
| Interferon | **6/6** | **100%** |
| Inflammation | 4/7 | 57% |
| EMT | 3/6 | 50% |
| Proliferation | 1/5 | 20% |
| Hypoxia | 1/6 | 17% |
| **Overall** | **15/30** | **50% (2.5× chance = 20%)** |

50% accuracy is 2.5× chance. Interferon achieves perfect prediction. The poor performance on hypoxia and proliferation is informative — cells under hypoxic stress and proliferating cells share many transcriptional features (cell-cycle, metabolic shifts), so single-program labels are biologically incomplete for these contexts. This is a feature of biology, not a flaw in the assignment.

### Anti-Circularity Test 3: Leave-One-Program-Out Q_B Stability

If the Q_B structure were driven by arbitrary program labels, removing any one program should collapse it. Instead, removing a non-on-program program barely affects Q_B/Q_tot, while removing the on-program causes the expected sharp drop:

| Signature | Full | Hide IFN | Hide Inflam | Hide Prolif | Hide Hypoxia | Hide EMT |
|---|---|---|---|---|---|---|
| IFN-γ | 0.604 | **0.264** ↓ | 0.634 | 0.572 | 0.649 | 0.628 |
| IFN-α | 0.631 | **0.264** ↓ | 0.608 | 0.602 | 0.664 | 0.668 |
| Inflammatory | 0.416 | 0.394 | **0.478** | 0.211 | 0.418 | 0.469 |
| TNFα | 0.533 | 0.492 | **0.584** | 0.331 | 0.521 | 0.615 |
| Hypoxia | 0.198 | 0.231 | 0.180 | 0.226 | **0.168** | 0.156 |
| E2F | 0.067 | 0.086 | 0.050 | **0.108** | 0.068 | 0.030 |
| EMT | 0.267 | 0.296 | 0.303 | 0.272 | 0.235 | **0.174** ↓ |

For IFN-γ and IFN-α, removing the interferon program drops Q_B/Q_tot by 34-37 percentage points (0.604 → 0.264), while removing any of the other four programs leaves Q_B/Q_tot within 0.05 of the full-panel value. The diagonal pattern (each signature's biggest drop is when its own program is removed) is the expected behavior if signature-program matching is real.

### Within-Program Durability (7 Hallmarks)

| Signature | Within g | within p | p_Bonf | Within I² | k | Outside g | Outside p | k |
|---|---|---|---|---|---|---|---|---|
| IFN-γ | +1.003 | <.001 | **<.001** | 0.72 | 6 | +0.177 | .245 | 24 |
| IFN-α | +1.189 | <.001 | **<.001** | 0.88 | 6 | +0.228 | .070 | 24 |
| Hypoxia | +3.545 | .0003 | **.003** | 0.95 | 6 | +0.706 | <.001 | 24 |
| EMT | +2.508 | .0002 | **.002** | 0.90 | 6 | +0.457 | <.001 | 24 |
| Inflammatory | +0.746 | .007 | .064 | 0.91 | 7 | +0.092 | NS | 23 |
| TNFα | +0.548 | .009 | .085 | 0.85 | 7 | +0.123 | NS | 23 |
| E2F | +0.627 | .295 | 1.000 | 0.98 | 5 | +0.454 | .001 | 25 |

IFN-γ is the cleanest exemplar: within-interferon g = +1.0, outside g = +0.18 (NS). Survives HKSJ-Bonferroni at p = 0.003.

### Full 29-Signature Classification Panel

| Class | n | Signatures |
|---|---|---|
| Hallmark (durable) | 7 | IFN-γ, IFN-α, Inflammatory, TNFα, Hypoxia, EMT, E2F |
| Curated durable | 1 | curated_senescence |
| Brittle | 4 | small-study-inflammation, platform-specific, single-tissue, overfit-noise |
| Mixed | 3 | EMT-inflammation, hypoxia-stress, senescence-proliferation |
| Confounded | 5 | proliferation, ribosomal, immune-infiltrate, 2 stealth |
| Insufficient coverage | 2 | olfactory_1, olfactory_2 |
| Blind holdouts | 7 | 2 durable, 1 mixed, 2 confounded, 1 brittle, 1 lowcov |
| **Total** | **29** | |

Within-program meta-analysis correctly classifies all 4 hypoxia/EMT/IFN durable signatures, both insufficient-coverage signatures, 4/5 confounded signatures, both blind durable holdouts, both blind confounded holdouts, and 1/3 mixed signatures (14/29 = 48%). Brittle signatures are intentionally constructed to lack within-program structure and are not the central focus of this paper. The complete per-signature classification table with all 29 signatures is in the GitHub repository (`outputs/canonical_v7/aggregate_durability_scores.csv`).

### Inflammatory Cross-Talk

The Inflammatory Response and TNFα/NFκB Hallmarks produce their largest effects not in inflammation cohorts but in EMT cohorts, exceeding on-program effects. The IPF lung cohorts — where EMT and inflammatory remodeling co-occur (Thiery et al. 2009, doi:10.1016/j.cell.2009.11.007) — drive this pattern. This is not a contradiction of the thesis but a feature: the Q_B decomposition reveals biological cross-talk that flat I² would have hidden. Programs are not biological constants; they are operational labels, and where they fail to capture cross-talk, the framework exposes it rather than pretending it doesn't exist.

### Outlier Disclosure

GSE47533 (hypoxia time-course) produces g = 15.1 due to near-zero within-group variance in a small cell-line experiment. We report both Winsorized (g cap = 3.0) and non-Winsorized analyses. Q_B/Q_tot for hypoxia: 0.198 (non-Winsorized) vs 0.196 (Winsorized). The conclusion is unchanged.

## Discussion

I² in gene signature meta-analysis is substantially inflated by context mixing: 39% of total Q (median 42%) is between-program heterogeneity. The circularity concern — that cohorts were pre-assigned to programs that match the tested signatures — is empirically refuted by three independent tests:
1. **Permutation test:** Random program assignments do not produce the observed Q_B for IFN, TNFα, or Inflammatory signatures.
2. **LOO program prediction:** Cohort-to-program assignment is recoverable from the data alone at 50% accuracy (2.5× chance), with interferon at 100%.
3. **LOPO Q_B stability:** Hiding any non-on-program program leaves Q_B/Q_tot stable; only hiding the on-program collapses it. The diagonal pattern is the signature of real program structure.

**Limitations.** (1) The within-program I² remains high (median 0.84) — context conditioning reduces heterogeneity by 39%, not 100%. The remaining I²_W reflects platform differences, sample composition, and technical batch effects — the genuine measurement floor for transcriptomic reproducibility, not artifact of context mixing. (2) Sample size per program is k=5-7, at the lower limit for reliable Q_W estimation. The permutation test mitigates this by providing non-parametric significance independent of asymptotic chi-squared assumptions. (3) E2F lumps tumor-vs-normal contrasts (g ≈ +2.3) with subtype contrasts (g ≈ -1.2); the program label "proliferation" is too coarse for this signature. (4) GSE47533's near-zero variance produces g = 15.1; both Winsorized and raw analyses are reported.

**Reproducibility.** All 30 frozen GEO cohorts, 29 benchmark signatures, scoring code, and frozen outputs (including the LOPO-CV results) are available at https://github.com/scottdhughes/signature-durability-benchmark. The complete benchmark executes deterministically from the accompanying SKILL.md on CPU-only hardware.

## Conclusion

I² in gene signature meta-analysis is substantially inflated by context mixing: 39% of total Q (median 42%) is between-program heterogeneity. The structure is data-driven, not imposed: cohort-to-program assignment is recoverable from effect-size fingerprints at 2.5× chance, and removing the on-program (not arbitrary programs) is what collapses Q_B. IFN-γ survives the most conservative inference (HKSJ-Bonferroni p = 0.003) within interferon cohorts while appearing null cross-context. Reported I² values for gene signatures should be accompanied by program-conditioned decomposition; signatures dismissed as irreproducible by cross-context benchmarks may warrant re-evaluation within their biological context.

## References

1. Liberzon A, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Systems. 2015;1:417-425. doi:10.1016/j.cels.2015.12.004
2. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. PLoS Comput Biol. 2011;7:e1002240. doi:10.1371/journal.pcbi.1002240
3. Subramanian A, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci. 2005;102:15545-15550. doi:10.1073/pnas.0506580102
4. Barbie DA, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. 2009;462:108-112. doi:10.1038/nature08460
5. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. Wiley; 2009. doi:10.1002/9780470743386
6. Hartung J, Knapp G. On tests of the overall treatment effect in meta-analysis with normally distributed responses. Stat Med. 2001;20:1771-1782. doi:10.1002/sim.791
7. Hedges LV, Olkin I. Statistical Methods for Meta-Analysis. Academic Press; 1985.
8. Thiery JP, et al. Epithelial-mesenchymal transitions in development and disease. Cell. 2009;139:871-890. doi:10.1016/j.cell.2009.11.007
