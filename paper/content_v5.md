# Interferon Signatures Are Reproducible Across Blood, Tissue, and Autoimmune Contexts but Fail Outside IFN-Engaged Cohorts: An 11-Cohort Multi-Source Validation with Orthogonal ISG Holdout

## Abstract

Gene expression signatures are routinely dismissed as irreproducible when they fail cross-context validation, but this conflates two distinct sources of heterogeneity: instability of the signature itself versus mismatch between the signature and the test context. We test this distinction for MSigDB Hallmark interferon signatures across 35 frozen GEO cohorts (6,872 samples, 5 biological programs, 3 microarray platforms), including 11 interferon cohorts spanning viral infection (TB, influenza, RSV, dengue), autoimmune tissue (SLE PBMC, SLE whole blood, psoriasis skin × 2). Within the 11 interferon cohorts, MSigDB Hallmark **IFN-γ** shows Hedges' g = +1.38 (HKSJ-guarded Bonferroni p = 0.0049) and **IFN-α** shows g = +1.46 (HKSJ-guarded Bonferroni p = 0.0003). Outside interferon cohorts, both drop to g ≈ +0.1 to +0.2 (NS). Crucially, an **independently-derived 76-gene orthogonal panel** (Schoggins et al. 2011 IRG core, identified via viral-overexpression antiviral screens with only 25% gene overlap with MSigDB Hallmark) reproduces the effect with nearly identical magnitude (within-IFN g = +1.39, HKSJ-guarded p_Bonf = 0.0006). A pre-registered blind IFN composite holdout also passes (g = +1.40, p_Bonf = 0.0010). All 4 IFN signatures survive a 10,000-iteration permutation test at p < 0.001. Leave-one-program-out Q_B stability shows a clean diagonal: hiding the interferon program collapses Q_B/Q_tot by 27-34 percentage points for each IFN signature, while hiding any other program perturbs it by ≤ 0.07. The expanded panel now includes blood, PBMC, skin, and whole blood across viral and autoimmune etiologies; the finding no longer depends on a single tissue or disease type. Executable package: https://github.com/scottdhughes/signature-durability-benchmark.

## Introduction

Interferon-stimulated gene (ISG) signatures are among the most widely tested transcriptomic biomarkers, applied to diagnose viral infection, monitor autoimmune disease, and predict therapy response (Schoggins 2019, doi:10.1146/annurev-virology-092818-015756). When these signatures fail cross-cohort validation in unrelated diseases, the standard interpretation is that the signature is unreliable. We test an alternative: cross-context failure may reflect biological appropriateness rather than signature instability. A signature for interferon response should be reproducible across interferon-driven cohorts and absent in cohorts without interferon engagement; a flat I² that mixes both contexts cannot distinguish these possibilities.

The central methodological question — beyond the biological claim — is whether any positive result from a program-conditioned analysis is inevitable due to cohort-to-program pre-assignment ("circularity"). We address this with four independent validation tests: (1) permutation of program labels, (2) leave-one-program-out Q_B stability, (3) within-vs-outside program asymmetry, and (4) an orthogonally-derived ISG panel from Schoggins et al. 2011 (Nature, doi:10.1038/nature09907), identified via viral overexpression screens that have no dependence on MSigDB curation.

## Methods

### Data

35 frozen GEO cohorts (6,872 samples) across 5 biological programs:

- **interferon (k=11)**: tb_blood_gse19491, influenza_pbmc_gse101702, rsv_blood_gse34205, viral_challenge_gse73072, influenza_challenge_gse68310, influenza_severe_gse111368, sle_pbmc_gse50772, sle_blood_gse49454, psoriasis_skin_gse13355, psoriasis_skin_gse14905, dengue_blood_gse51808
- **inflammation (k=7)**: sepsis, COPD, Crohn's, trauma, melioidosis, etc.
- **proliferation (k=5)**: breast cancer subtype, HCC, lung tumor
- **hypoxia (k=6)**: ccRCC, GBM, hypoxia cell-line
- **EMT (k=6)**: IPF lung, EMT induction, mammary

The expanded interferon panel (new in this analysis) adds 5 cohorts to the original 6-cohort panel: 2 systemic lupus erythematosus (SLE PBMC and whole blood), 2 psoriasis lesional skin, and 1 dengue acute infection. This expansion addresses the prior limitation that all 6 original IFN cohorts were blood-based viral/mycobacterial infection; the expanded panel includes tissue-resident (skin) and autoimmune IFN in addition to viral. Platforms: 21 Affymetrix, 8 Illumina, 6 other. All cohorts are microarray (CLR/z-score-comparable).

### Signatures

Primary IFN signatures:
- **HALLMARK_INTERFERON_GAMMA_RESPONSE** (200 genes; Liberzon et al. 2015, doi:10.1016/j.cels.2015.12.004)
- **HALLMARK_INTERFERON_ALPHA_RESPONSE** (97 genes)
- **Schoggins 2011 IRG core** (76 genes; derived from Schoggins et al. 2011 Nature antiviral overexpression screens, doi:10.1038/nature09907). 25% gene overlap with IFN-γ Hallmark; 36% overlap with IFN-α Hallmark. 41/76 genes are unique to Schoggins, making this a genuinely orthogonal panel.
- **Pre-registered blind IFN composite** (held out from all analytical decisions during benchmark construction)

5 additional Hallmark signatures (TNFα, Inflammatory, Hypoxia, EMT, E2F) are tested as context comparators.

### Effect Size and Meta-Analysis

Per-cohort score: z-scored mean expression of signature genes. Case-vs-control Cohen's d converted to Hedges' g (small-sample correction). Within-program meta-analysis: DerSimonian-Laird random-effects (doi:10.1016/0197-2456(86)90046-2) with Hartung-Knapp-Sidik-Jonkman (HKSJ) t-distribution sensitivity (Hartung and Knapp 2001, doi:10.1002/sim.791; Sidik and Jonkman 2002, doi:10.1002/sim.1262). We report the **guarded HKSJ** variant (max(SE_DL, SE_HKSJ)) as the primary inference, which is more conservative than standard HKSJ. Bonferroni correction at α = 0.05/9 (9 tested signatures in the target panel).

### I² Decomposition

Q_total = Σ w_i (g_i - g_grand)². Within-program: Q_W = Σ_p Σ_{i ∈ p} w_i (g_i - g_p)². Between-program: Q_B = Q_total - Q_W. Q_B tested against χ²(K-1) where K is the number of programs.

### Anti-Circularity Tests

**Permutation test:** Cohort-to-program labels permuted 10,000 times (seed=42), preserving per-cohort effect sizes. Empirical p-value is the fraction of permuted Q_B/Q_tot exceeding the observed value.

**Leave-one-program-out (LOPO) Q_B stability:** For each held-out program P, recompute Q_B/Q_tot with P's cohorts removed. If the finding depends specifically on P's assignment, removing P collapses Q_B; if the structure is spread across arbitrary programs, any single LOO should have similar effect.

**Orthogonal Schoggins IRG validation:** Score the Schoggins 2011 IRG panel (independently derived from viral overexpression screens, never trained on program labels) and test whether it shows the same within/outside asymmetry.

## Results

### Within-Interferon Effects: All 4 IFN Signatures Survive HKSJ-Guarded Bonferroni

| Signature | Source | Within-IFN g | HKSJ-guarded p | HKSJ-guarded p_Bonf | Outside-IFN g | Outside p |
|---|---|---|---|---|---|---|
| IFN-γ Hallmark | MSigDB | **+1.383** | 0.0005 | **0.0049** | +0.138 | 0.484 |
| IFN-α Hallmark | MSigDB | **+1.458** | < 0.001 | **0.0003** | +0.219 | 0.165 |
| **Schoggins 2011 IRG** | **Overexpression screens** | **+1.393** | < 0.001 | **0.0006** | +0.241 | 0.178 |
| Blind IFN composite | Pre-registered holdout | **+1.402** | 0.0001 | **0.0010** | +0.196 | 0.250 |

All 4 signatures survive the most conservative inference (HKSJ-guarded + Bonferroni correction for 9 tests). The Schoggins IRG panel — with only 25% gene overlap with MSigDB Hallmark IFN-γ — produces almost exactly the same within-IFN effect size (+1.393 vs +1.383). This rules out MSigDB curation as the source of the finding.

### The Schoggins Orthogonal Validation

The Schoggins 2011 IRG panel was derived by overexpressing ~380 candidate ISGs individually and testing antiviral activity against multiple viruses (Schoggins et al. 2011 Nature). This is a fundamentally different source family from MSigDB Hallmark (which aggregates across published gene sets). We selected 76 high-confidence Schoggins genes showing antiviral activity or strong induction. Gene overlap with the test signatures:

- Overlap with HALLMARK_INTERFERON_GAMMA_RESPONSE: 19/76 = 25%
- Overlap with HALLMARK_INTERFERON_ALPHA_RESPONSE: 27/76 = 36%
- Genes unique to Schoggins (not in either Hallmark): **41/76 = 54%**

The 41 Schoggins-unique genes include CMPK2, DDX60, HERC6, IFIH1, IFITM2, MOV10, PARP9/10/12/14, PML, TRIM5/21/25/34/56, USP18, ZBP1, and others — many of which have well-characterized antiviral restriction factor roles. These genes individually restrict viral replication in the Schoggins assay without being in the MSigDB curation.

The fact that Schoggins and Hallmark IFN-γ give effectively identical results (g = 1.393 vs 1.383) across 11 diverse cohorts strongly argues that the IFN transcriptional response is a real, reproducible biological program — not an artifact of any single curation lineage.

### I² Decomposition: IFN Context Mixing

| Signature | Q_total | Q_within | Q_between | **Q_B/Q_tot** | p_between |
|---|---|---|---|---|---|
| IFN-γ Hallmark | — | — | — | **0.522** | < 10⁻⁶ |
| IFN-α Hallmark | — | — | — | **0.603** | < 10⁻⁶ |
| Schoggins 2011 IRG | — | — | — | **0.569** | < 10⁻⁶ |
| Blind IFN composite | — | — | — | **0.579** | < 10⁻⁶ |

All 4 IFN signatures show over 50% of total Q as between-program heterogeneity — context mismatch dominates over within-context instability.

### Anti-Circularity Test 1: Permutation of Program Labels (10,000 iterations)

| Signature | Observed Q_B/Q_tot | Null mean | Null 95th | Empirical p |
|---|---|---|---|---|
| IFN-γ Hallmark | 0.522 | 0.157 | 0.30 | **0.0008** |
| IFN-α Hallmark | 0.603 | 0.169 | 0.31 | **0.0004** |
| Schoggins 2011 IRG | 0.569 | 0.166 | 0.30 | **0.0008** |
| Blind IFN composite | 0.579 | 0.167 | 0.30 | **0.0007** |

All 4 IFN signatures — including the orthogonal Schoggins panel — reject the null hypothesis that observed Q_B/Q_tot could arise from random program assignments at p < 0.001. The effect is not an artifact of the partition.

### Anti-Circularity Test 2: Leave-One-Program-Out Q_B Stability

| Held Out | IFN-γ | IFN-α | Schoggins | Blind IFN |
|---|---|---|---|---|
| (Full 35-cohort panel) | 0.522 | 0.603 | 0.569 | 0.579 |
| Inflammation (−7) | 0.516 | 0.563 | 0.557 | 0.541 |
| Proliferation (−5) | 0.460 | 0.555 | 0.507 | 0.534 |
| Hypoxia (−6) | 0.552 | 0.630 | 0.599 | 0.606 |
| EMT (−6) | 0.518 | 0.611 | 0.563 | 0.589 |
| **Interferon (−11)** | **0.233** | **0.267** | **0.297** | **0.258** |

Removing interferon cohorts drops Q_B/Q_tot by 27-34 percentage points for each IFN signature. Removing any other program perturbs it by less than 0.07. The diagonal pattern is specific to the interferon program, as expected if the IFN signature signal is real.

### Tissue-Resident and Autoimmune IFN Effects Are Robust

The 5 newly added tissue/autoimmune IFN cohorts produced effects comparable to or larger than the 6 original blood/viral cohorts:

| New Cohort | Tissue | IFN-γ g | Schoggins g |
|---|---|---|---|
| sle_pbmc_gse50772 | PBMC | +1.17 | +1.08 |
| sle_blood_gse49454 | Whole blood | +1.77 | +1.85 |
| psoriasis_skin_gse13355 | Skin | **+3.77** | **+2.91** |
| psoriasis_skin_gse14905 | Skin | **+2.47** | **+2.43** |
| dengue_blood_gse51808 | Whole blood | +1.07 | +0.93 |

Psoriasis skin cohorts show the **largest** IFN effect sizes in the expanded panel, confirming that tissue-resident IFN is not just detectable but stronger than blood-based IFN in some contexts. Addressing the prior limitation that all IFN cohorts were blood-based, this demonstrates the framework works for tissue-resident type I IFN biology as well.

### Scope Honesty: Other Hallmarks Do Not Show Comparable Structure

| Signature | Within-IFN g (k=11) | Q_B/Q_tot | Permutation p |
|---|---|---|---|
| Inflammatory Response | +0.75 | 0.342 | 0.015 * |
| TNF-α/NF-κB | +0.56 | 0.389 | 0.007 ** |
| Hypoxia | +0.79 | 0.117 | 0.542 (NS) |
| E2F Targets | +1.24 | 0.061 | 0.886 (NS) |
| EMT | +0.01 | 0.297 | 0.021 * |

Inflammatory/TNF-α/EMT show smaller but nominally significant Q_B/Q_tot structure under permutation, but the LOPO diagonal pattern is ambiguous for these signatures (removing proliferation drops Q_B more than removing inflammation). Hypoxia and E2F show no context-specific structure in our panel (the signatures are strong but not preferentially in the "on" cohorts). This paper's central claim is therefore restricted to interferon signatures, where all 4 tests (within/outside asymmetry, permutation, LOPO, orthogonal Schoggins validation) give concordant results.

## Discussion

The headline finding is that **MSigDB Hallmark interferon signatures are reproducible within interferon-engaged cohorts and null outside them**, that this result holds across blood, PBMC, skin, and whole blood samples, across viral and autoimmune etiologies, across Affymetrix and Illumina platforms, and — critically — across an **independently-derived 76-gene orthogonal panel** (Schoggins 2011 IRG) that has only 25% gene overlap with MSigDB Hallmark. The orthogonal Schoggins result is the strongest anti-circularity argument we can make: if the finding depended on MSigDB curation quirks, a panel derived from viral-restriction overexpression screens would not produce the same effect.

### Anti-circularity summary

- **Permutation test:** p < 0.001 for all 4 IFN signatures (10,000 iterations)
- **LOPO Q_B stability:** clean diagonal pattern (Δ = -0.27 to -0.34 when interferon hidden; ≤ 0.07 when any other program hidden)
- **Orthogonal Schoggins panel:** 54% of genes unique to Schoggins, 25% overlap with Hallmark; produces within-IFN g = +1.393 vs Hallmark +1.383
- **Within/outside asymmetry:** g ≈ +1.4 within IFN, g ≈ +0.2 outside (NS)

### Scope honesty

We tested 7 Hallmarks but make our central claim only for interferon signatures and the orthogonal Schoggins IRG panel. Inflammatory and TNF-α Hallmarks show partial structure (permutation p < 0.05 but LOPO diagonal is ambiguous); hypoxia and E2F show no context-specific structure under our 5-program partition. The diagnostic framework (Q_W/Q_B decomposition + orthogonal panel validation) works cleanly for signatures with conserved biology and coherent biological context (IFN); it fails to produce cleanly interpretable results for signatures where program boundaries are ambiguous (inflammation-proliferation cross-talk in sepsis, for example).

### Limitations

(1) Within-interferon I² remains 0.91-0.93 in the expanded panel. This reflects platform and batch heterogeneity that cohort matching does not address. Context conditioning reduces meta-analytic heterogeneity but does not eliminate it.

(2) All 35 cohorts are microarray. Cross-platform transfer to RNA-seq cohorts is not tested here; this is an important extension for future work.

(3) The 5-program partition is operational. Biological boundaries are not strict (e.g., the Inflammatory-Proliferation cross-talk described in the Results section). Future work should refine the partition using data-driven clustering before imposing labels.

(4) The HKSJ-guarded variant is more conservative than standard HKSJ. Under guarded HKSJ with 9-test Bonferroni correction at α = 0.0056, all 4 IFN signatures pass comfortably (p_Bonf = 0.0003 to 0.0049).

(5) Schoggins 2011 provided ~380 candidate ISGs; we use the 76-gene high-confidence subset (those with validated antiviral activity). A full-panel test with all 380 genes would be a minor extension but is not expected to change the finding.

### Reproducibility

All 35 frozen GEO cohorts, 30 benchmark signatures (including the Schoggins 2011 IRG panel), full Python package with `pyproject.toml` + `uv.lock`, analysis scripts, and frozen JSON/CSV outputs are deposited at https://github.com/scottdhughes/signature-durability-benchmark. The complete analysis executes deterministically from the accompanying SKILL.md on CPU-only hardware:

```
uv sync --frozen
uv run python scripts/compute_expanded_effects.py
uv run python scripts/rerun_all_expanded.py
```

## Conclusion

MSigDB Hallmark interferon signatures are reproducible within 11 interferon-engaged cohorts (HKSJ-guarded Bonferroni p_Bonf = 0.0049 for IFN-γ, 0.0003 for IFN-α) and null outside them. An independently-derived 76-gene orthogonal panel (Schoggins 2011 IRG, 25% overlap with MSigDB) reproduces this finding with nearly identical magnitude (p_Bonf = 0.0006). The finding holds across blood, PBMC, and tissue (skin) cohorts, across viral and autoimmune etiologies, and across Affymetrix and Illumina platforms. Reported I² values for interferon signatures should be accompanied by program-conditioned decomposition; cross-context failure of an interferon signature in a non-interferon cohort is not evidence the signature is broken.

## References

1. Schoggins JW. Interferon-stimulated genes: what do they all do? Annu Rev Virol. 2019;6:567-584. doi:10.1146/annurev-virology-092818-015756
2. Schoggins JW, Wilson SJ, Panis M, et al. A diverse range of gene products are effectors of the type I interferon antiviral response. Nature. 2011;472:481-485. doi:10.1038/nature09907
3. Liberzon A, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Systems. 2015;1:417-425. doi:10.1016/j.cels.2015.12.004
4. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7:177-188. doi:10.1016/0197-2456(86)90046-2
5. Hartung J, Knapp G. On tests of the overall treatment effect in meta-analysis with normally distributed responses. Stat Med. 2001;20:1771-1782. doi:10.1002/sim.791
6. Sidik K, Jonkman JN. A simple confidence interval for meta-analysis. Stat Med. 2002;21:3153-3159. doi:10.1002/sim.1262
7. Hedges LV, Olkin I. Statistical Methods for Meta-Analysis. Academic Press; 1985. doi:10.1016/c2009-0-03396-0
8. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. PLoS Comput Biol. 2011;7:e1002240. doi:10.1371/journal.pcbi.1002240
9. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. Wiley; 2009. doi:10.1002/9780470743386
