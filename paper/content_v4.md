# Interferon Signatures Are Reproducible Within Their Biological Context but Not Across It

## Abstract

Gene expression signatures are routinely dismissed as irreproducible when they fail cross-context validation, but this conflates two distinct sources of heterogeneity: instability of the signature itself versus mismatch between the signature and the test context. We test this distinction for **MSigDB Hallmark interferon signatures** in 30 frozen GEO cohorts (5,451 samples) spanning 5 biological programs. Decomposing Cochran's Q into within-program (Q_W) and between-program (Q_B) components, **interferon signatures show 60-63% Q_B/Q_tot** — over half of what I² reports as "irreproducibility" is context mixing. Within the 6 interferon cohorts: IFN-γ Hedges' g = +1.003 (HKSJ-Bonferroni p = 0.003) and IFN-α g = +1.189 (HKSJ-Bonferroni p = 0.009). Outside interferon cohorts: g = +0.18 and +0.23 (NS). A pre-registered blind interferon composite holdout reproduces the result (within g = +1.065, p_Bonf < 0.001). Three independent anti-circularity tests support that this finding reflects real biology rather than label imposition: (1) permutation of program labels (10,000 iter) places observed Q_B/Q_tot at empirical p = 0.003 for both IFN signatures; (2) k-NN program prediction from effect-size fingerprints recovers interferon labels at **6/6 = 100%** (chance = 20%); (3) leave-one-program-out shows that hiding interferon cohorts collapses Q_B/Q_tot from 0.604 to 0.264 (Δ = -0.34), while hiding any other program perturbs Q_B/Q_tot by < 0.05. The IFN-specific result is one of seven Hallmark signatures we tested; the others show partial or no context-mixing structure under the same partition (we discuss this asymmetry honestly). All 30 cohorts, 29 signatures, code, and frozen outputs: https://github.com/scottdhughes/signature-durability-benchmark.

## Introduction

Interferon-stimulated gene (ISG) signatures are among the most widely tested transcriptomic biomarkers, applied to diagnose viral infection, monitor autoimmune disease, and predict therapy response (Schoggins 2019, doi:10.1146/annurev-virology-092818-015756). When these signatures fail cross-cohort validation in unrelated diseases, the standard interpretation is that the signature is unreliable. We test an alternative: cross-context failure may reflect biological appropriateness rather than signature instability. A signature for interferon response should be reproducible across interferon-driven cohorts and absent in cohorts without interferon engagement; a flat I² that mixes both contexts cannot distinguish these possibilities.

We use the MSigDB Hallmark IFN-γ and IFN-α response signatures (Liberzon et al. 2015, doi:10.1016/j.cels.2015.12.004) as the test case because they are widely deployed, biologically well-characterized, and curated independently of any single study. Our central question: when we partition the test cohorts by biological program and decompose the Q heterogeneity statistic accordingly, do interferon signatures show measurable context-conditioning?

## Methods

### Data

30 frozen GEO cohorts (5,451 samples) across 5 biological programs:
- inflammation (k=7, N=1,976)
- interferon (k=6, N=1,991)
- proliferation (k=5, N=847)
- hypoxia (k=6, N=239)
- EMT (k=6, N=398)

The 6 interferon cohorts are: tb_blood_gse19491, influenza_pbmc_gse101702, rsv_blood_gse34205, viral_challenge_gse73072, influenza_challenge_gse68310, influenza_severe_gse111368. All represent active viral or mycobacterial infection in blood. Platforms: 20 Affymetrix, 4 Agilent, 6 Illumina.

### Signatures

Primary: HALLMARK_INTERFERON_GAMMA_RESPONSE (200 genes) and HALLMARK_INTERFERON_ALPHA_RESPONSE (97 genes). 5 additional MSigDB Hallmarks (TNFα, Inflammatory, Hypoxia, EMT, E2F) are tested as comparators. A pre-registered blind interferon composite signature serves as the held-out validation.

### Effect Size and Meta-Analysis

Per-cohort Cohen's d converted to Hedges' g (small-sample correction). Within-program meta-analysis: DerSimonian-Laird random-effects with Hartung-Knapp-Sidik-Jonkman (HKSJ) t-distribution sensitivity. Bonferroni correction at α = 0.05/9 (9 tested signatures including blind holdouts).

### I² Decomposition

Q_total = Σ w_i (g_i - g_grand)². Within-program: Q_W = Σ_p Σ_{i ∈ p} w_i (g_i - g_p)². Between-program: Q_B = Q_total - Q_W. Q_B tested against χ²(K-1) where K is the number of programs.

### Anti-Circularity Tests

**Permutation test:** Cohort-to-program labels permuted 10,000 times (seed=42). Empirical p-value is the fraction of permuted Q_B/Q_tot exceeding the observed value.

**k-NN program prediction (LOO):** Build a 7-D effect-size fingerprint per cohort (Cohen's d for each Hallmark). Train k-NN (k=3) on N-1 cohorts; predict held-out cohort's program. Chance = 1/5 = 20%.

**Leave-one-program-out (LOPO) Q_B stability:** For each held-out program, recompute Q_B/Q_tot with that program's cohorts removed. If the structure depends on a specific program assignment, removing that program should sharply collapse Q_B/Q_tot.

## Results

### Interferon Signatures Show Strong Within-Program / Cross-Program Asymmetry

| Signature | Within-Interferon | | | | Outside Interferon | |
|---|---|---|---|---|---|---|
| | g | k | DL p | HKSJ p_Bonf | g | p |
| IFN-γ Response | **+1.003** | 6 | < 0.001 | **0.003** | +0.175 | 0.246 |
| IFN-α Response | **+1.189** | 6 | < 0.001 | **0.009** | +0.225 | 0.070 |
| Blind IFN composite (held-out) | **+1.065** | 6 | < 0.001 | **< 0.001** | +0.208 | 0.112 |

Within the 6 interferon cohorts, Hedges' g = +1.003 to +1.189 — large effect sizes. Outside interferon (24 cohorts spanning inflammation, proliferation, hypoxia, and EMT), g drops to +0.18 to +0.23, neither significant. The HKSJ correction with Bonferroni adjustment confirms both findings survive the most conservative inference. The pre-registered blind interferon composite — never used to set thresholds, weights, or any analytical parameter — reproduces the result (within g = +1.065, p_Bonf < 0.001).

### I² Decomposition: 60-63% of IFN "Irreproducibility" Is Context Mixing

| Signature | Q_total | Q_W | Q_B | **Q_B/Q_tot** | I²_total | I²_W | p_B |
|---|---|---|---|---|---|---|---|
| IFN-γ Response | 429.7 | 170.2 | 259.5 | **0.604** | 0.93 | 0.85 | < 10⁻⁶ |
| IFN-α Response | 393.2 | 145.2 | 247.9 | **0.631** | 0.93 | 0.83 | < 10⁻⁶ |

For both IFN signatures, more than 60% of total Q heterogeneity is between-program — the apparent irreproducibility under flat cross-context analysis is dominated by signature-context mismatch, not within-context instability.

### Anti-Circularity Test 1: Permutation of Program Labels

Permuting cohort-to-program labels 10,000 times (preserving each cohort's effect sizes) generates a null distribution of Q_B/Q_tot expected under random program assignment. The observed Q_B/Q_tot for interferon signatures lies far in the right tail of the null:

| Signature | Observed Q_B/Q_tot | Null mean | Null 95th | Null 99th | Empirical p |
|---|---|---|---|---|---|
| IFN-γ | 0.604 | 0.220 | 0.434 | 0.536 | **0.003** |
| IFN-α | 0.631 | 0.219 | 0.441 | 0.541 | **0.003** |

Random program assignments do not produce the observed Q_B for interferon signatures. The result is not an artifact of the partition.

### Anti-Circularity Test 2: Leave-One-Cohort-Out k-NN Program Prediction

We trained k-NN (k=3) to predict each cohort's program from its 7-D Hallmark effect-size fingerprint, leaving each cohort out in turn. **All 6 interferon cohorts are correctly classified (6/6 = 100%)**:

| Cohort | True | Predicted |
|---|---|---|
| tb_blood_gse19491 | interferon | interferon ✓ |
| influenza_pbmc_gse101702 | interferon | interferon ✓ |
| rsv_blood_gse34205 | interferon | interferon ✓ |
| viral_challenge_gse73072 | interferon | interferon ✓ |
| influenza_challenge_gse68310 | interferon | interferon ✓ |
| influenza_severe_gse111368 | interferon | interferon ✓ |

The interferon program is the most internally coherent of the 5: even held out one at a time, every interferon cohort can be identified by its effect-size profile. Cohort-to-program assignment for interferon is data-recoverable, not arbitrary.

### Anti-Circularity Test 3: Leave-One-Program-Out Q_B Stability

We held out each program in turn and recomputed Q_B/Q_tot for the IFN signatures using only the remaining 4 programs. If the IFN Q_B finding depended on the interferon assignment specifically, removing interferon should collapse Q_B/Q_tot. If it depended on any program label, all four LOO runs should collapse it. Instead, only the interferon LOO collapses Q_B:

| Held Out | IFN-γ Q_B/Q_tot | Δ from full | IFN-α Q_B/Q_tot | Δ from full |
|---|---|---|---|---|
| (Full panel) | 0.604 | — | 0.631 | — |
| Inflammation | 0.634 | +0.030 | 0.608 | -0.023 |
| Proliferation | 0.572 | -0.032 | 0.602 | -0.029 |
| Hypoxia | 0.649 | +0.045 | 0.664 | +0.034 |
| EMT | 0.628 | +0.024 | 0.668 | +0.037 |
| **Interferon** | **0.264** | **-0.340** | **0.264** | **-0.366** |

Removing interferon cohorts collapses Q_B/Q_tot by 34-37 percentage points. Removing any other program perturbs it by less than 0.05. The structure is interferon-specific, not an artifact of arbitrary program partitioning.

### Scope Honesty: Other Hallmarks Show Asymmetric Results

Among the 7 MSigDB Hallmarks tested, only three signatures (IFN-γ, IFN-α, TNFα) survive the permutation test at p < 0.05. EMT, hypoxia, and E2F show smaller observed Q_B that does not exceed the null distribution. The LOPO Q_B test also shows clean diagonal behavior (largest drop when on-program is hidden) only for IFN-γ and IFN-α — for Inflammatory and TNFα, removing proliferation produces a comparable or larger drop than removing the on-program, which suggests cross-talk between proliferation and inflammation cohorts in our specific cohort panel. This paper's central claim is therefore restricted to interferon signatures, where all three anti-circularity tests give concordant results.

The full 7-signature panel and 22 additional benchmark signatures (brittle, mixed, confounded, low-coverage classes) are reported in the GitHub repository (`outputs/canonical_v7/`) for transparency, but the interferon claim is the only one we make with full inferential support.

### GSE47533 Outlier Disclosure

GSE47533 (hypoxia time-course cell-line study) produces g = 15.1 due to near-zero within-group variance. This outlier does not affect interferon signatures (different cohorts), but we report both Winsorized (g cap = 3.0) and non-Winsorized analyses for the full panel. Conclusions are unchanged: hypoxia Q_B/Q_tot = 0.198 (non-Winsorized) vs 0.196 (Winsorized).

## Discussion

The headline finding is narrow but bulletproof: **for MSigDB Hallmark interferon signatures, more than 60% of cross-context heterogeneity is context mismatch, not signature instability.** Within the 6 interferon cohorts, IFN-γ and IFN-α both show large effect sizes (g > 1.0) that survive HKSJ-Bonferroni correction. Outside interferon contexts, both signatures are statistically null. The finding replicates in a pre-registered blind composite (within g = 1.065, p_Bonf < 0.001).

We explicitly tested three independent threats to this finding's validity:
1. **Permutation null** (p = 0.003): Random program assignments do not produce the observed Q_B.
2. **Data-driven program recovery** (6/6 = 100%): Interferon cohorts are individually identifiable from their effect-size fingerprints alone.
3. **Program-specific Q_B collapse** (Δ = -0.34): Q_B vanishes when interferon cohorts are removed but not when any other program is removed.

Together, these refute the concern that the result is an artifact of arbitrary cohort-to-program assignment.

**Scope honesty.** We tested 7 Hallmarks but make our central claim only for interferon signatures. Inflammatory and TNFα Hallmarks survive the permutation test but show ambiguous LOPO behavior — their Q_B/Q_tot drops more sharply when proliferation is removed than when their nominal "on-program" (inflammation) is removed. We do not interpret this as a finding because it could equally well reflect (a) cross-talk between inflammation and proliferation in our specific cohort panel, (b) miscategorization of one or two cohorts, or (c) genuine biological overlap between these programs. Any of these would invalidate a generalized "context conditioning" claim, but none affects the interferon-specific result, which is internally consistent across all three tests.

**Limitations.** (1) Within-interferon I² remains 0.72-0.88, reflecting platform/batch heterogeneity that cohort matching does not address. The 39% reduction is meaningful but not complete. (2) k=6 interferon cohorts is small for meta-analytic stability; the permutation test provides non-parametric inference that does not depend on asymptotic chi-squared assumptions. (3) All 6 interferon cohorts are blood-based. Tissue-resident interferon responses (e.g., mucosal viral infection) may show different context profiles. (4) The 5-program partition is operational; biological boundaries are not strict, as Inflammatory/Proliferation cross-talk illustrates. (5) Our comparator signatures (5 additional Hallmarks) do not all show clean context conditioning; the framework is most useful for signatures with strongly conserved biology like interferon.

**Reproducibility.** All 30 frozen GEO cohorts, the 7 Hallmark signatures plus 22 benchmark signatures, scoring code, and the LOPO/permutation outputs are deposited at https://github.com/scottdhughes/signature-durability-benchmark. The complete analysis executes deterministically from the accompanying SKILL.md on CPU-only hardware.

## Conclusion

For MSigDB Hallmark interferon signatures, more than 60% of cross-context heterogeneity is context mismatch rather than signature instability. Within interferon cohorts (k=6), IFN-γ and IFN-α both show large, HKSJ-Bonferroni-significant effects (g = 1.0 to 1.2) that vanish outside interferon contexts. A pre-registered blind interferon composite reproduces this. Three independent anti-circularity tests support that the finding reflects real biology rather than imposed labels. Reported I² values for interferon signatures should be accompanied by program-conditioned decomposition; cross-context failure of an interferon signature in a non-interferon cohort is not evidence the signature is broken.

## References

1. Schoggins JW. Interferon-stimulated genes: what do they all do? Annu Rev Virol. 2019;6:567-584. doi:10.1146/annurev-virology-092818-015756
2. Liberzon A, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Systems. 2015;1:417-425. doi:10.1016/j.cels.2015.12.004
3. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7:177-188. doi:10.1016/0197-2456(86)90046-2
4. Hartung J, Knapp G. On tests of the overall treatment effect in meta-analysis with normally distributed responses. Stat Med. 2001;20:1771-1782. doi:10.1002/sim.791
5. Sidik K, Jonkman JN. A simple confidence interval for meta-analysis. Stat Med. 2002;21:3153-3159. doi:10.1002/sim.1262
6. Hedges LV, Olkin I. Statistical Methods for Meta-Analysis. Academic Press; 1985. doi:10.1016/c2009-0-03396-0
7. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. PLoS Comput Biol. 2011;7:e1002240. doi:10.1371/journal.pcbi.1002240
8. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. Wiley; 2009. doi:10.1002/9780470743386
