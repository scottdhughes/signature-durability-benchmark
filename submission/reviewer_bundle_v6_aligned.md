You are reviewing a computational biology paper for clawRxiv (Claw4S 2026, deadline April 20, 2026). Your decision options are Accept / Weak Accept / Weak Reject / Reject. Please be rigorous but fair.

Context for this method-focused v6 bundle:
- The paper is framed as a program-conditioned diagnostic for distinguishing signature instability from context mismatch, with interferon as the strongest validation case.
- The IFN-gamma and IFN-alpha signatures are 30-gene frozen benchmark cores anchored to the broader MSigDB Hallmark families, not the full 200/97-gene Hallmark sets.
- The expansion cohorts were admitted using 10 markers: STAT1, IRF1, IFIT1, IFIT2, ISG15, MX1, OAS1, GBP1, CXCL10, RSAD2.
- The strictly-disjoint Schoggins result is scoped to rule out exact admission-marker reuse as the driver, not every possible cohort-selection effect.
- The within-IFN heterogeneity analysis is framed as a descriptive weighted Q partition rather than a definitive mixed-effects meta-regression because k=11 and the combined moderator partition contains sparse cells.
- New in this bundle: a real triage interface, held-out validation (leave-one-cohort-out and split-half), explicit failure-mode analysis for ambiguous comparator programs, a rescued published-signature portability case using the Ayers IFN-gamma-related 6-gene profile, a held-out four-cohort bulk RNA-seq transfer test, a deterministic second breadth case (Hallmark Hypoxia, intentionally still mixed), a generated workflow figure, a mixed metadata-first prospective v1 round with an interpretable acute-PBMC miss, a second fresh v2 registry that is externally timestamped and intentionally still pending evaluation, and a provenance audit verifying that the active 35-cohort scored panel uses real GEO data and that the paper-facing target signatures are non-synthetic.

Please score on:
(A) Whether the paper's claims are now appropriately scoped and defensible
(B) Any remaining concerns or overclaims
(C) Whether the central claim — "a program-conditioned diagnostic can distinguish signature instability from context mismatch, validated here on interferon signatures and shown to detect coarse-label failure modes" — is defensible for Claw4S
(D) Your decision: Accept / Weak Accept / Weak Reject / Reject, with one-paragraph rationale
(E) Any final fixes that should happen before submission

Executable package: https://github.com/scottdhughes/signature-durability-benchmark

============================================================
PAPER (content_v6.md)
============================================================
# Program-Conditioned Diagnostic for Transcriptomic Signature Durability: Validation on Interferon Signatures across 35 Frozen GEO Cohorts

## Abstract

Gene expression signatures are often called irreproducible when they fail across heterogeneous validation cohorts, but that failure can reflect either instability of the signature itself or mismatch between the signature and the test context. We present a **program-conditioned diagnostic** that scores a signature across a frozen reference panel, compares within-program versus outside-program effects, tests program structure by permutation, and surfaces failure modes when labels are too coarse. In 35 frozen GEO cohorts (5,922 samples, 5 biological programs), the frozen 30-gene MSigDB-anchored IFN-γ and IFN-α cores, an orthogonal 76-gene Schoggins panel, and a strictly-disjoint 41-gene Schoggins subset all show large within-IFN effects and small, non-significant outside-IFN effects; passing the IFN-γ core through `triage` yields a full-model class of `mixed` but a best-supported program of interferon and a within-program class of `durable`. Held-out validation is strong: all four IFN signatures achieve 100% leave-one-cohort-out sign prediction and 100% split-half sign agreement, and four external bulk RNA-seq cohorts (719 samples) reproduce the IFN quartet with guarded Bonferroni-significant pooled effects while inflammatory, TNF-α/NFκB, and E2F comparators remain non-significant. The same rescue logic works on an external oncology signature: the Ayers 6-gene IFN-γ-related profile looks `mixed` in aggregate but resolves to a `durable` interferon-context signal. The prospective layer is intentionally non-tautological rather than uniformly confirmatory: two predeclared v1 cohorts satisfy the 4/4 IFN sign forecast, whereas one severity-mixed acute PBMC cohort inverts all four IFN signatures. Comparator signatures show the complementary use-case: inflammatory/TNF-α/NFκB ambiguity is driven more by proliferation than inflammation, and E2F targets expose a coarse label bucket rather than a random signature. Executable package: https://github.com/scottdhughes/signature-durability-benchmark.

## Introduction

Gene-expression signatures are routinely evaluated by asking whether they retain a pooled effect across many cohorts. When they do not, the default conclusion is often that the signature is unstable or not reproducible. That inference is too coarse. A signature can fail pooled cross-context validation either because the signature is genuinely brittle or because the benchmark has mixed together cohorts where the biology should be present with cohorts where it should not.

Interferon-stimulated gene (ISG) signatures provide a clean place to test this distinction. They are among the most widely used transcriptomic programs in viral infection, autoimmunity, and immune monitoring (Schoggins 2019, doi:10.1146/annurev-virology-092818-015756). A good interferon signature should reproduce across IFN-engaged cohorts and flatten outside them. A flat cross-context meta-analysis or a single headline I² cannot separate that biologically appropriate failure from true instability.

The contribution here is therefore methodological, not just biological. We build a **program-conditioned diagnostic** that takes a candidate signature, scores it against a frozen cohort panel, and answers four questions: (1) what biological program best supports the signature, (2) how large is the within-program effect relative to outside-program mismatch, (3) is the observed program structure stronger than random labelings, and (4) when the signal is ambiguous, does that ambiguity reflect a broken signature or a broken labeling scheme? We validate the diagnostic where ground truth is comparatively strong, namely interferon biology, and then use comparator programs to show that the same framework can detect its own boundary conditions.

Two circularity threats matter for this kind of benchmark. The first is **statistical circularity**, where any partition of cohorts might manufacture between-group structure. The second is **curation-level circularity**, where expansion cohorts could have been admitted using the same genes later used to “validate” the signal. We address the first with a 10,000-iteration permutation of cohort-to-program labels, and the second with two orthogonal holdouts: an independently derived Schoggins 2011 antiviral IRG panel and a **strictly-disjoint 41-gene Schoggins subset** sharing zero genes with either Hallmark IFN core and zero genes with the 10 admission markers.

## Methods

### Data

35 frozen GEO cohorts (5,922 samples, 14 distinct microarray platform IDs) across 5 biological programs:

- **interferon (k=11)**: tb_blood_gse19491, influenza_pbmc_gse101702, rsv_blood_gse34205, viral_challenge_gse73072, influenza_challenge_gse68310, influenza_severe_gse111368, sle_pbmc_gse50772, sle_blood_gse49454, psoriasis_skin_gse13355, psoriasis_skin_gse14905, dengue_blood_gse51808
- **inflammation (k=7)**: sepsis, COPD, Crohn's, trauma, melioidosis
- **proliferation (k=5)**: breast cancer subtype, HCC, lung tumor
- **hypoxia (k=6)**: ccRCC, GBM, hypoxia cell-line
- **EMT (k=6)**: IPF lung, EMT induction, mammary

The expanded interferon panel adds 5 cohorts to the original 6-cohort panel: 2 systemic lupus erythematosus (SLE PBMC and whole blood), 2 psoriasis lesional skin, and 1 dengue acute infection. This expansion addresses the prior limitation that all original IFN cohorts were blood-based viral/mycobacterial infection; the expanded panel includes tissue-resident (skin) and autoimmune IFN in addition to viral. Platform family membership across the 35 cohorts: 21 Affymetrix, 8 Illumina, 6 other microarray platforms. All sample counts reported here refer to the final frozen 35-cohort release; earlier working drafts used pre-freeze curation counts and are superseded.

**Expansion cohort admission criteria.** Each candidate cohort was admitted on the basis of concordant case-vs-control direction for ≥8/10 of the following 10 well-characterized ISG markers: STAT1, IRF1, IFIT1, IFIT2, ISG15, MX1, OAS1, GBP1, CXCL10, RSAD2. This admission step creates an explicit curation-level circularity concern, which the strictly-disjoint Schoggins subset (below) is designed to test against exact gene-reuse.

### Signatures

Primary IFN signatures:
- **HALLMARK_INTERFERON_GAMMA_RESPONSE core** (30 genes; frozen benchmark core anchored to the broader MSigDB Hallmark INTERFERON_GAMMA_RESPONSE family from Liberzon et al. 2015, doi:10.1016/j.cels.2015.12.004). This is a manually curated compact benchmark subset emphasizing canonical IFN effector plus antigen-presentation genes, shipped in `data/freeze/signatures.tsv` and documented in `data/curation/source_provenance.md`; exact membership was frozen in the March 30, 2026 benchmark release, inherited unchanged into the expanded reruns, and not reselected from expanded-panel performance.
- **HALLMARK_INTERFERON_ALPHA_RESPONSE core** (30 genes; frozen benchmark core anchored to the broader MSigDB Hallmark INTERFERON_ALPHA_RESPONSE family). This is likewise a manually curated compact benchmark subset emphasizing canonical type-I IFN antiviral effector genes, frozen in the March 30, 2026 benchmark release and inherited unchanged into the expanded reruns.
- **Schoggins 2011 IRG** (76 genes; Schoggins et al. 2011 Nature antiviral overexpression screens, doi:10.1038/nature09907). 19/76 (25%) overlap with the IFN-γ core; 27/76 (36%) overlap with the IFN-α core; 41/76 (54%) are unique to Schoggins.
- **Strictly-disjoint Schoggins subset** (41 genes; the complement: genes in the Schoggins 76 that appear in neither frozen Hallmark IFN core). Verified by exhaustive set subtraction. None of these 41 genes appear in the 10-gene admission marker list above.
- **Pre-registered blind IFN composite** (held out from all analytical decisions during original benchmark construction)

Five additional Hallmark signatures (TNFα/NFκB, Inflammatory, Hypoxia, E2F, EMT) are tested as context comparators.

Unlike the comparator Hallmarks, the two IFN cores are therefore compact frozen benchmark-release subsets rather than simple “top-30” truncations of the full MSigDB Hallmark membership.

### Effect Size and Meta-Analysis

Per-cohort score: z-scored mean expression of signature genes across samples within cohort, with signed weights supported for arbitrary input signatures. Case-vs-control Cohen's d converted to Hedges' g (small-sample correction). Within-program meta-analysis: DerSimonian-Laird random-effects (doi:10.1016/0197-2456(86)90046-2) with Hartung-Knapp-Sidik-Jonkman (HKSJ) t-distribution sensitivity (Hartung and Knapp 2001, doi:10.1002/sim.791; Sidik and Jonkman 2002, doi:10.1002/sim.1262). We report the **guarded HKSJ** variant (SE_primary = max(SE_DL, SE_HKSJ)) as the primary inference, which is more conservative than standard HKSJ. Bonferroni correction at α = 0.05/9 ≈ 0.0056 for the 9-test target panel.

### I² Decomposition

Q_total = Σ wᵢ(gᵢ − ḡ)². Within-program: Q_W = Σ_p Σ_{i∈p} wᵢ(gᵢ − ḡ_p)². Between-program: Q_B = Q_total − Q_W. Q_B tested against χ²(K − 1) where K is the number of programs (K = 5).

### Program-Conditioned Triage Diagnostic

The benchmark now exposes a reusable `triage` workflow for arbitrary input signatures. Given a TSV/CSV containing `gene_symbol` and optional `direction` / `weight` columns, the diagnostic:

1. Scores the signature in every frozen cohort using the same per-cohort z-scored scoring rule as the canonical benchmark.
2. Computes within-program and outside-program guarded-HKSJ meta-analyses for each candidate home program.
3. Ranks programs by **within-minus-outside separation**.
4. Tests global program structure by permutation of cohort labels.
5. Emits both a **full-model class** and a **within-program class**, along with `diagnostic.json`, `per_cohort_effects.csv`, and a human-readable `diagnostic_summary.md`.

This interface is designed to answer the practical question that motivated the paper: given a new signature, is the observed heterogeneity more consistent with a broken signature or with a context-matched signal being evaluated on the wrong cohorts?

The packaged `triage` decision surface is intentionally simple and frozen. The **full-model** classifier labels a signature `confounded` if its strongest confounder effect matches or exceeds the aggregate effect; otherwise it labels the signature `brittle` if aggregate p > 0.10 or direction consistency < 0.50, `mixed` if I² > 0.75 or leave-one-out stability < 0.60, and `durable` otherwise. The **within-program** classifier first checks the inferred home program; if the guarded within-program pooled effect satisfies p < 0.05 and |g| > 0.2, it labels the signature `durable`, otherwise it falls back to the same cross-context rules. Mean cohort coverage below 0.60 yields `insufficient_coverage`. These thresholds are frozen in code and were not refit during the v6 analyses. For speed, the packaged `triage` CLI uses 1,000 label permutations; the canonical IFN decomposition and paper-level Q_B significance tests continue to use 10,000 permutations.

The generated workflow figure in `paper/figure_workflow.pdf` packages that logic explicitly: input signature, frozen-panel scoring, within/outside separation, permutation, decision branch (durable versus coarse-label ambiguity), and the external/prospective audit layers.

The broader 30-signature benchmark also contains 19 explicit synthetic control signatures used as benchmark scaffolding for brittle, mixed, confounded, and insufficient-coverage behavior. Those synthetic controls are quarantined from the paper-facing target panel, external validations, and prospective case studies; the v6 paper reports only non-synthetic target signatures.

### Anti-Circularity, Held-Out Validation, and Failure-Mode Analyses

**Permutation test:** Cohort-to-program labels permuted 10,000 times (seed = 42), preserving per-cohort effect sizes. Empirical p-value is the fraction of permuted Q_B/Q_tot exceeding the observed value.

**Orthogonal Schoggins IRG validation:** Score the full 76-gene Schoggins panel and test whether it shows the same within/outside asymmetry as the MSigDB-anchored Hallmark IFN cores.

**Strictly-disjoint Schoggins subset validation:** Score only the 41 Schoggins genes with zero overlap against the Hallmark IFN-γ core, the Hallmark IFN-α core, and the 10 admission markers. A large within-IFN effect from this subset rules out exact reuse of the admission-marker genes as the sole driver of the expanded-panel result.

**Held-out validation:** For each target signature, perform leave-one-cohort-out refits within the signature’s home program and exhaustive split-half replication within the same program. These analyses ask whether a program-conditioned finding survives withholding whole cohorts rather than relying only on the full pooled fit.

**Held-out external RNA-seq validation:** To test cross-platform transfer outside the frozen microarray panel, we assembled a primary bulk RNA-seq panel from public GEO supplementary matrices using four operational rules: (1) processed host-gene matrix deposited directly in GEO, (2) explicit two-group acute SARS-CoV-2 case/control labels in the series matrix, (3) bulk whole blood or upper-airway tissue, and (4) at least 25 samples after filtering to the primary two-group comparison. This yielded four primary cohorts totaling 719 samples: GSE152641 whole blood COVID-19 vs healthy (62/24), GSE171110 whole blood severe COVID-19 vs healthy (44/10), GSE152075 nasopharyngeal SARS-CoV-2 positive vs negative (430/54), and GSE167000 whole blood hospitalized SARS-CoV-2 positive vs negative (65/30). Raw counts were mapped to gene symbols via NCBI `Homo_sapiens.gene_info.gz` where needed, log1p-transformed, collapsed to gene symbol, and scored with the same canonical per-cohort scorer used for the microarray panel. We evaluated the IFN quartet plus three high-risk comparators (Inflammatory, TNF-α/NFκB, E2F) and applied Bonferroni correction across these 7 external signatures. An additional PBMC cohort (GSE152418; 16 COVID-19, 17 healthy after excluding one convalescent sample) was retained separately as an exploratory stress test rather than pooled with the primary bulk RNA-seq panel.

**Metadata-first prospective holdout prediction:** To create a genuinely predeclared challenge inside the repository, we froze `data/prospective_holdout/prediction_registry_v1.tsv` and `data/prospective_holdout/PREDICTION_PROTOCOL.md` before scoring the declared cohorts. The registry was defined from GEO series metadata and supplementary file listings rather than from held-out effect estimates. The primary round contained three cohorts: GSE184610 mucosal SARS-CoV-2 positive vs control (159/183), GSE243217 PBMC COVID-19 vs healthy donor (35/15), and GSE202805 acute COVID-19 PBMC vs healthy (42/10). The primary success rule was intentionally strict: each of the four IFN signatures should have positive per-cohort g in every cohort, and across the three-cohort panel each signature should achieve pooled guarded-HKSJ g ≥ 0.8, Bonferroni-4 p < 0.05, and sign consistency 1.0. The evaluation script downloaded the declared files into a separate cache, scored them with the same canonical scorer, and wrote registry/file hashes plus per-cohort effects to `outputs/canonical_v8/prospective_holdout_validation.json`. Unlike the later v2 round, v1 is a repo-frozen internal challenge rather than an externally timestamped preregistration.

After observing that mixed v1 outcome, we did **not** replace the round with a cleaner retrospective panel. Instead we declared a second fresh round in `prediction_registry_v2.tsv` and `PREDICTION_PROTOCOL_v2.md`, covering GSE243442 acute yellow-fever PBMC vs control, GSE213168 influenza blood infected vs healthy control, and GSE155237 RSV nasal challenge D3 infected vs D3 uninfected. The v2 registry and protocol hashes are bound to an external RFC3161 timestamp receipt under `data/prospective_holdout/external_timestamps/prospective_holdout_v2/`. That second round is intentionally kept unevaluated in the current paper so it remains a real future challenge rather than a post-v1 cleanup.

To make the next readout reproducible without adding post hoc flexibility, we also replaced the ad hoc holdout evaluator with a **generic round-based scoring interface**. The CLI command `prospective-round-evaluate` requires a registry, protocol, receipt, and explicit output directory; it verifies receipt hashes before scoring, writes round-scoped outputs under `outputs/canonical_v8/prospective_rounds/<round_id>/`, and aborts if declared inputs do not match the receipt.

**Second positive generalization case:** We precommitted to choose the second non-IFN breadth case from **Hallmark Hypoxia** versus **Hallmark EMT** using frozen artifacts only: (1) HKSJ-guarded Bonferroni support in the current frozen outputs, (2) perfect leave-one-out sign prediction, (3) perfect split-half sign agreement, (4) lower mean absolute foreign-program effect, and (5) higher within-minus-abs(outside) separation. Under that rule Hallmark Hypoxia wins. We then run `triage` on the frozen Hallmark Hypoxia input and score one external exact-perturbation cohort chosen under a separate deterministic GEO filter (`data/external_hypoxia/SEARCH_PROTOCOL.md`): GSE179885 human T-cell RNA-seq cultured in hypoxia versus normoxia (6 vs 6), with a directly deposited processed TPM matrix.

**Failure-mode analysis:** We treat ambiguous comparator signatures as a design feature rather than a nuisance. For inflammatory, TNF-α/NFκB, and E2F signatures, we examine LOPO drops, held-out instability, and within-program effect spread to distinguish “broken signature” behavior from evidence that the program labels themselves are too coarse.

**Descriptive within-IFN heterogeneity decomposition:** For the IFN-γ core, compute inverse-variance weighted subgroup means and Q partitions against three moderators (tissue, etiology, platform family), then report the fraction of within-IFN Q_total explained (R²). Because only 11 cohorts populate 6 combined cells, we interpret this analysis descriptively rather than as a definitive mixed-effects meta-regression.

## Results

### Triage Correctly Recovers a Durable IFN Signal That Looks `Mixed` in the Aggregate

The most important method sanity check is whether the diagnostic behaves correctly when given a known durable signal without being told its home program. Passing the frozen IFN-γ core into `triage` as an arbitrary input signature yields:

- **Inferred best-supported program:** interferon
- **Full-model classification:** `mixed`
- **Within-program classification:** `durable`
- **Within-IFN pooled effect:** g = +1.383, guarded-HKSJ p = 0.00055
- **Outside-IFN pooled effect:** g = +0.138, p = 0.484
- **Program-structure permutation p:** 0.000999

This is precisely the distinction the benchmark is meant to recover. If the same signature is summarized only by its aggregate cross-context behavior, it looks partially diluted. Once the cohort panel is partitioned by biological program, the signal resolves into a clean durable-within-context pattern. The practical point is that a top-line “mixed” or high-I² result is not by itself evidence that the signature is unusable.

### Within-Interferon Effects: All 4 IFN Signatures Survive HKSJ-Guarded Bonferroni

| Signature | Source | Within-IFN g | HKSJ-guarded p | HKSJ-guarded p_Bonf | Outside-IFN g | Outside p |
|---|---|---:|---:|---:|---:|---:|
| IFN-γ Hallmark core | MSigDB-anchored | **+1.383** | 0.0005 | **0.0049** | +0.138 | 0.484 |
| IFN-α Hallmark core | MSigDB-anchored | **+1.458** | < 0.001 | **0.0003** | +0.219 | 0.165 |
| **Schoggins 2011 IRG (76 genes)** | **Overexpression screens** | **+1.393** | < 0.001 | **0.0006** | +0.241 | 0.178 |
| **Schoggins strictly-disjoint (41 genes)** | **Overexpression screens, zero marker overlap** | **+1.247** | **0.000097** | **0.00088** | +0.241 | 0.213 |
| Blind IFN composite | Pre-registered holdout | **+1.402** | 0.0001 | **0.0010** | +0.196 | 0.250 |

All 5 signatures survive the most conservative inference (HKSJ-guarded + Bonferroni correction for 9 tests). The strictly-disjoint Schoggins subset — sharing zero genes with either Hallmark IFN core and zero with the admission markers — produces within-IFN g = +1.247 with p_Bonf = 0.00088, an order of magnitude below the Bonferroni threshold.

### Orthogonal and Strictly-Disjoint Holdouts Address the Main Circularity Critiques

The Schoggins 2011 IRG panel was derived by overexpressing ~380 candidate ISGs individually and testing antiviral activity against multiple viruses (Schoggins et al. 2011 Nature). This is a fundamentally different source family from the MSigDB Hallmark collection. We selected 76 high-confidence Schoggins genes showing antiviral activity or strong induction. Gene overlap with the frozen Hallmark IFN benchmark cores:

- Overlap with the HALLMARK_INTERFERON_GAMMA_RESPONSE core: 19/76 = 25%
- Overlap with the HALLMARK_INTERFERON_ALPHA_RESPONSE core: 27/76 = 36%
- Genes in **neither** Hallmark IFN core: **41/76 = 54%**

The 41 strictly-disjoint genes — C6orf150, CMPK2, CXCL11, DDX60, DHX58, GBP3, GBP4, GBP5, HERC6, HES4, IFIH1, IFIT5, IFITM2, IRGM, LAMP3, LGALS9, MOV10, NAMPT, NT5C3, PARP10, PARP12, PARP14, PARP9, PHF11, PLSCR1, PML, PNPT1, RNF19B, RTP4, SERPING1, SP100, TDRD7, TREX1, TRIM21, TRIM25, TRIM34, TRIM5, TRIM56, UBD, ZBP1, ZC3HAV1 — include well-characterized viral restriction factors (TRIM family, PARP family, IFIH1/MDA5, ZBP1) and have zero overlap with the {STAT1, IRF1, IFIT1, IFIT2, ISG15, MX1, OAS1, GBP1, CXCL10, RSAD2} marker set used to admit the 5 expansion cohorts. If the expanded panel's IFN effect were driven mainly by exact admission-marker reuse, this 41-gene subset would be expected to attenuate toward null. Instead it produces g = +1.247 at p_Bonf = 0.00088. Zero overlap here does not imply statistical independence: these genes remain part of the same co-regulated interferon biology as the admission markers, so the test isolates gene-reuse circularity rather than pathway-level correlation.

Permutation results support the same conclusion from a different angle:

| Signature | Observed Q_B/Q_tot | Null mean | Null 95th | Empirical p |
|---|---:|---:|---:|---:|
| IFN-γ Hallmark core | 0.522 | 0.157 | 0.30 | **0.0008** |
| IFN-α Hallmark core | 0.603 | 0.169 | 0.31 | **0.0004** |
| Schoggins 2011 IRG | 0.569 | 0.166 | 0.30 | **0.0008** |
| Blind IFN composite | 0.579 | 0.167 | 0.30 | **0.0007** |

All 4 IFN signatures reject the null that the observed program structure arises from random label assignments.

### Held-Out Validation Shows the IFN Result Survives Withholding Entire Cohorts

Held-out validation strengthens the claim that the diagnostic is not merely fitting pooled retrospective structure. Within the IFN program:

| Signature | Leave-One-Cohort-Out Sign Prediction | Split-Half Sign Agreement | Split-Half Both Significant |
|---|---:|---:|---:|
| IFN-γ Hallmark core | **1.000** | **1.000** | 0.871 |
| IFN-α Hallmark core | **1.000** | **1.000** | **1.000** |
| Schoggins 2011 IRG | **1.000** | **1.000** | **1.000** |
| Blind IFN composite | **1.000** | **1.000** | **1.000** |

All four IFN signatures preserve the correct sign in every leave-one-cohort-out refit and in every split-half partition. Under guarded HKSJ, 87.1% to 100% of split-half partitions remain significant in both halves. This is exactly what a durable within-context signature should look like under withholding, and it is substantially stronger than a one-shot pooled fit.

### Held-Out Bulk RNA-seq Extension Reproduces the IFN Quartet Across Platforms

Cross-platform transfer is strong when tested on held-out bulk RNA-seq cohorts that were not part of the frozen microarray benchmark. Across four primary external cohorts (GSE152641, GSE171110, GSE152075, GSE167000; 719 samples total), all four IFN signatures are positive in all four cohorts.

| Signature | External RNA-seq pooled g | Guarded p | Guarded p_Bonf,7 | I² |
|---|---:|---:|---:|---:|
| IFN-γ Hallmark core | **+0.922** | 0.0032 | **0.022** | 0.000 |
| IFN-α Hallmark core | **+1.193** | 0.0016 | **0.011** | 0.000 |
| Schoggins 2011 IRG | **+0.996** | 0.0026 | **0.018** | 0.000 |
| Blind IFN composite | **+1.177** | 0.0029 | **0.020** | 0.244 |
| Inflammatory Response | +0.427 | 0.170 | 1.000 | 0.770 |
| TNF-α/NFκB | +0.382 | 0.295 | 1.000 | 0.859 |
| E2F Targets | +0.818 | 0.196 | 1.000 | 0.944 |

This external layer does two useful things at once. First, it shows genuine cross-platform transfer: the IFN quartet remains significant even after conservative Bonferroni correction across the 7 externally tested signatures. Second, it preserves specificity: three non-IFN comparators remain non-significant and far more heterogeneous under the same external benchmark. The exploratory PBMC stress test is intentionally kept separate. In GSE152418 PBMC COVID-19 vs healthy, IFN-γ attenuates to g = −0.409, Schoggins to g = −0.203, and the blind IFN composite to g = −0.078, while IFN-α remains weakly positive at g = +0.127. That pattern is exactly why cell-selected PBMC is a useful stress test rather than a clean primary bulk RNA-seq validation cohort.

### A Mixed Breadth Illustration Beyond IFN: Hypoxia Wins the Frozen Rule but Stays Mixed

The second-case selection rule chooses **Hallmark Hypoxia** over **Hallmark EMT** under the current frozen outputs. This is a breadth illustration rather than a second clean positive validation on the same level as IFN. Hypoxia is the only one of the two candidates with HKSJ-guarded Bonferroni support in the frozen artifact table (p_Bonf = 0.024 versus 1.0 for EMT), while both signatures achieve perfect leave-one-out sign prediction and perfect split-half sign agreement. Hallmark Hypoxia also retains the larger within-minus-abs(outside) separation (2.733 versus 2.064).

That breadth result is informative precisely because it is **not** as clean as interferon. Running `triage` on the frozen Hallmark Hypoxia input correctly recovers **hypoxia** as the best-supported program, but both the full-model and within-program classes remain `mixed`, with outside-program pooled g = +0.645 (p = 0.0066). The diagnostic is therefore doing something more useful than hunting for another IFN-like win: it identifies hypoxia as the right home program while simultaneously making the cross-program bleed-through visible.

The external exact-perturbation cohort supports the same direction. In GSE179885 human T-cell bulk RNA-seq cultured under hypoxia versus normoxia, Hallmark Hypoxia ranks **first** among seven scored signatures with Hedges' g = +5.11 and full coverage, ahead of Hallmark EMT (g = +0.94), while IFN, inflammatory, TNF-α/NFκB, and E2F comparators are all negative in that cohort. Because this cohort has only 12 total samples, we treat it as supportive external breadth evidence rather than as a headline pooled table.

### A Published Oncology IFNγ Signature Is Rescued by Triage

To show that the diagnostic is useful on a signature imported from outside this benchmark, we scored the published **Ayers IFN-γ-related 6-gene clinical response profile** (IDO1, CXCL10, CXCL9, HLA-DRA, STAT1, IFNG; Ayers et al. 2017, doi:10.1172/JCI91190) through `triage`. This is a portability test rather than an internal benchmark rerun: the signature was derived for PD-1 blockade response in oncology, not for our frozen 5-program panel.

Two of the six Ayers genes (STAT1 and CXCL10) overlap our 10 expansion-admission markers, five overlap the frozen IFN-γ core, one overlaps the IFN-α core, and IFNG itself is in neither frozen Hallmark IFN core. That overlap is expected for an external IFN-related clinical signature and does not make this an anti-circularity analysis; the point of the case study is portability and rescue, not circularity adjudication.

Naively pooled across all 35 cohorts, the signature looks diluted rather than clean: aggregate effect = +0.246 with I² = 0.929 and a full-model class of `mixed`. But the best-supported program is still **interferon**, and the within-interferon class is **`durable`**. The interferon-specific pooled effect is +0.866 (guarded HKSJ p = 0.0245), while the outside-interferon effect is essentially zero (−0.010, p = 0.949).

This is exactly the rescue behavior the diagnostic is intended to provide. An externally sourced signature can look only modestly positive when all cohorts are pooled together, yet still resolve into a biologically coherent home program once within/outside structure is examined. The case study does **not** claim to re-analyze the original oncology trial cohorts; rather, it shows that `triage` can recover the underlying interferon program from a clinically motivated signature that would otherwise look mixed in a broad cross-context benchmark.

### Metadata-First Prospective Holdout Round Is Mixed Rather Than Tautological

The first predeclared holdout round was intentionally harder than the retrospective bulk RNA-seq transfer test. All decisions about cohort inclusion and success criteria were frozen in `prediction_registry_v1.tsv` before the held-out matrices were scored. Two cohorts matched the forecast cleanly: GSE184610 mucosal SARS-CoV-2 positive vs control and GSE243217 PBMC COVID-19 vs healthy were 4/4 positive across the IFN quartet. The third cohort, GSE202805 acute COVID-19 PBMC vs healthy, inverted all four IFN signatures.

| Cohort | Tissue | IFN sign hits |
|---|---|---:|
| GSE184610 | Oronasopharyngeal mucosa | 4/4 |
| GSE243217 | PBMC | 4/4 |
| GSE202805 | PBMC | 0/4 |

Because the GSE202805 miss breaks sign consistency, the pooled prospective rule fails overall: pooled guarded-HKSJ g values are +0.359 (IFN-γ), +0.657 (IFN-α), +0.549 (Schoggins 2011 IRG), and +0.576 (blind IFN composite), with Bonferroni-4 adjusted p-values between 0.93 and 1.0 and sign consistency 0.667 for all four signatures. That is not a retrospective success story, but it is still methodologically useful. The prospective layer now functions as a real challenge benchmark rather than a guaranteed win: it confirms that the diagnostic can be wrong on metadata-first forecasts and therefore is not simply rubber-stamping any plausible viral cohort.

Unlike the later v2 registry, v1 should be read as a repo-frozen internal challenge rather than as an externally timestamped preregistration.

The miss is also informative rather than opaque. GEO metadata describe GSE202805 as a severity-mixed acute PBMC series (`Global Transcriptomic analysis of PBMC derived from COVID-19 acute-severe patients and convalescents`), and the acute arm combines 4 mild, 6 moderate, 14 severe, and 18 severe/ICU samples. In that same cohort, the IFN quartet is negative while inflammatory and E2F comparators are positive (Hedges' g = +0.482 and +0.446, respectively). The failure is therefore more consistent with a compartment/severity mismatch in acute PBMC biology than with a generic scoring or download failure: the predeclared rule “acute viral PBMC should be IFN-positive” was simply too coarse for this cohort.

### A Second Fresh Prospective Round Is Declared and Externally Timestamped

The v1 miss could easily have invited quiet cohort substitution. We therefore froze a second round instead of retrofitting the first one. The v2 panel is metadata-only and pending evaluation: GSE243442 acute yellow-fever PBMC vs control (41/14), GSE213168 influenza blood infected vs healthy control (19/22), and GSE155237 RSV nasal challenge D3 infected vs D3 uninfected (23/16). Its declaration manifest is bound to an external RFC3161 timestamp receipt, so the next prospective readout can be evaluated against a third-party-dated registry rather than against a file whose history is only local. No v2 scoring outputs are used in any result table here.

### The Diagnostic Also Exposes When the Program Labels Are Too Coarse

The method is most useful when it fails informatively. Comparator signatures show that ambiguous outputs are not always “bad signatures”; sometimes they are telling us that the benchmark partition is biologically coarse.

| Signature | Full Q_B/Q_tot | Largest LOPO Drop When Hiding | LOO Accuracy | Split-Half Sign Agreement | Interpretation |
|---|---:|---|---:|---:|---|
| Inflammatory Response | 0.342 | proliferation (Δ = 0.195) | 0.857 | 1.000 | Cross-talk signal is driven more by proliferation contrasts than by a clean inflammation-only boundary |
| TNF-α/NFκB | 0.389 | proliferation (Δ = 0.181) | 1.000 | 1.000 | Same ambiguity: the strongest program leverage is not inflammation itself |
| E2F Targets | 0.061 | inflammation (Δ = 0.027) | 0.800 | 0.750 | The coarse “proliferation” bucket mixes distinct contrast types rather than defining one coherent context |

For inflammatory and TNF-α/NFκB signatures, removing proliferation lowers Q_B more than removing inflammation itself. That is the opposite of what a clean home-program signature should do and is better interpreted as inflammation–proliferation cross-talk than as a simple validation failure. E2F shows a different failure mode: its between-program structure is weak (Q_B/Q_tot = 0.061), yet its within-proliferation I² is 0.959 and its per-cohort effects span **−1.2068 to +2.1416** across tumor-vs-normal, prognosis, and subtype contrasts. The diagnostic is therefore not just separating durable from brittle signatures; it is identifying when the label space itself is too blunt to support a clean decision.

Operationally, these failure modes imply different next steps. When LOPO leverage is driven by a foreign program, the right move is to split or relabel the offending program and rerun the diagnostic. When within-program spread remains extreme despite weak between-program structure, the right move is to stratify by contrast type or tissue before concluding that the signature itself is unstable.

### Descriptive Within-IFN Heterogeneity Is Structured Biology, Not Featureless Noise

Within-interferon I² of 0.925 is substantial. For the IFN-γ core, we therefore partition weighted Q descriptively against three moderators:

| Moderator | k levels | Q_between | Q_within | **R²** |
|---|---:|---:|---:|---:|
| Tissue (blood/PBMC vs skin) | 2 | 97.96 | 35.82 | **0.732** |
| Etiology (viral vs autoimmune) | 2 | 74.28 | 59.49 | **0.555** |
| Platform family (Affy/Illumina/other) | 3 | 44.56 | 89.21 | **0.333** |
| **Combined tissue × etiology × platform** | 6 | **121.46** | **12.31** | **0.908** |

Subgroup meta-analyses:

| Subgroup | k | Pooled g | Q | I² |
|---|---:|---:|---:|---:|
| Blood/PBMC | 9 | +0.889 | 28.21 | 0.716 |
| Skin | 2 | +3.237 | 7.61 | 0.869 |
| Viral | 7 | +0.842 | 14.46 | 0.585 |
| Autoimmune | 4 | +2.169 | 45.03 | 0.933 |

Within-IFN I² is not featureless noise. Blood/PBMC cohorts show pooled g = +0.89 with I² = 0.72; skin cohorts show pooled g = +3.24 with I² = 0.87; viral cohorts show g = +0.84 with I² = 0.59; autoimmune cohorts show g = +2.17. This pattern is biologically coherent with prior reports of interferon-inducible blood signatures in SLE and type-I interferon-enriched psoriatic lesions (Baechler et al. 2003, doi:10.1073/pnas.0337679100; Yao et al. 2008, doi:10.1371/journal.pone.0002737). Because only 11 cohorts populate 6 combined tissue × etiology × platform cells, including singleton cells, the 90.8% figure should be read as a descriptive partition of weighted Q rather than a definitive inferential meta-regression estimate.

### Psoriasis Sensitivity Analysis

The largest per-cohort effect (psoriasis_skin_gse13355, g = +3.77) sits outside the main blood/PBMC cluster. We tested whether the pooled within-IFN result depends on this cohort:

| Scenario | k | Pooled g | I² | HKSJ-guarded p_Bonf |
|---|---:|---:|---:|---:|
| Full expanded panel | 11 | +1.383 | 0.925 | 0.0049 |
| Drop psoriasis_skin_gse13355 | 10 | +1.128 | 0.807 | **0.00058** |
| Drop both skin cohorts | 9 | +1.019 | 0.716 | **0.00019** |
| Winsorize per-cohort g at 90th percentile | 11 | +1.250 | 0.862 | **0.00054** |

Every sensitivity scenario *tightens* the inference relative to the full panel. The psoriasis cohorts contribute upward leverage on the pooled mean but inflate I² enough that under guarded HKSJ their removal narrows the confidence interval and reduces p_Bonf by roughly an order of magnitude. The finding is not driven by a single high-effect cohort; it holds robustly with any skin subset and under Winsorization.

## Discussion

The paper now supports a narrower but more useful claim than earlier versions: **program-conditioned diagnostics can distinguish signature instability from context mismatch, and interferon signatures provide a strong validation case for that diagnostic.** The IFN result remains biologically meaningful, but the more general contribution is the reusable decision logic:

- a signature can look `mixed` in the aggregate yet be clearly **durable within the right biological program**
- orthogonal and strictly-disjoint holdouts can test whether a positive result is an artifact of label assignment or exact gene reuse
- ambiguous outputs can themselves be informative, revealing cross-talk or coarse program labels rather than simple signature failure

### What the IFN validation establishes

Within the 11 IFN-engaged cohorts, all four IFN signatures show large, significant effects under guarded HKSJ and shrink to small, non-significant effects outside IFN cohorts. The effect is reproduced in an independently derived 76-gene antiviral panel and persists in a 41-gene strictly-disjoint Schoggins subset with zero overlap against the cohort-admission markers. Leave-one-cohort-out and split-half analyses show that this is not just a pooled full-panel artifact, and the held-out bulk RNA-seq extension shows that the same quartet transfers across platform and tissue boundaries outside the frozen reference panel.

### Why the failure cases matter

The inflammatory, TNF-α/NFκB, and E2F results are not secondary clutter; they are part of the method validation. Inflammatory and TNF signatures show genuine signal, but their strongest leverage comes from proliferation removal rather than inflammation removal, indicating that the current partition slices across real biology. E2F shows a different pattern: low between-program structure coupled to very high within-program heterogeneity, consistent with a coarse proliferation label that mixes tumor-vs-normal, prognosis, and subtype contrasts. This is useful behavior for a diagnostic. It means the framework is not merely rewarding the one clean program we chose in advance; it is also flagging when the benchmark ontology needs refinement.

### Limitations

(1) The strongest validation case in this paper is still interferon biology. The new deterministic hypoxia breadth case helps, but it remains mixed rather than IFN-clean, so the paper is still not a universal claim about all signature families.

(2) The frozen decision surface is still defined by 35 microarray cohorts. We now test cross-platform transfer on four held-out bulk RNA-seq cohorts, but those RNA-seq cohorts are not yet integrated into the frozen reference panel itself.

(3) The 5-program partition is operational. Biological boundaries are not strict, and the failure-mode results show exactly where the current labels are too coarse.

(4) The descriptive IFN-γ weighted Q partition is informative but not a definitive inferential meta-regression. Only 11 IFN cohorts populate the combined moderator cells.

(5) Expansion cohorts were admitted on a 10-marker direction-of-effect concordance check (STAT1, IRF1, IFIT1, IFIT2, ISG15, MX1, OAS1, GBP1, CXCL10, RSAD2). The strictly-disjoint Schoggins subset makes exact admission-marker reuse an unlikely explanation for the result, but it does not fully eliminate broader cohort-selection effects.

(6) Held-out validation is now stronger than in earlier versions: the paper includes internal withholding, a four-cohort bulk RNA-seq transfer test, and a first metadata-first prospective holdout round. That prospective round is mixed rather than uniformly confirmatory (2 of 3 cohorts satisfy the 4/4 IFN sign forecast), so it should be read as a real challenge benchmark, not as another cleaned-up validation layer.

(7) The repo now includes a second fresh prospective round whose registry and protocol are bound to an external RFC3161 timestamp receipt, plus a generic locked-round scoring interface for the future readout. That improves the audit trail, but the v2 round is intentionally unevaluated here and therefore does not contribute numerical evidence in the current manuscript.

(8) The declaration bundle is now release-ready inside the repo (`submission/archive_bundles/prospective_holdout_v2_declaration/`), and the repo ships `LICENSE`, `CITATION.cff`, and `.zenodo.json` metadata. Those changes harden the publication path, but they do not substitute for an actual public immutable release or DOI until one is minted.

### Reproducibility

All 35 frozen GEO cohorts, 30 benchmark signatures (including the Schoggins 2011 IRG panel), full Python package with `pyproject.toml` + `uv.lock`, analysis scripts, and frozen JSON/CSV outputs are deposited at https://github.com/scottdhughes/signature-durability-benchmark. The complete analysis executes deterministically from the accompanying SKILL.md on CPU-only hardware:

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

The externally timestamped v2 declaration can be verified with:

```bash
uv run python -m signature_durability_benchmark.cli declare-prospective-round \
  --registry data/prospective_holdout/prediction_registry_v2.tsv \
  --protocol data/prospective_holdout/PREDICTION_PROTOCOL_v2.md
```

The future locked v2 readout should be executed, when the paper is ready to score it, with:

```bash
uv run python -m signature_durability_benchmark.cli prospective-round-evaluate \
  --registry data/prospective_holdout/prediction_registry_v2.tsv \
  --protocol data/prospective_holdout/PREDICTION_PROTOCOL_v2.md \
  --receipt data/prospective_holdout/external_timestamps/prospective_holdout_v2/declaration_receipt.json \
  --out outputs/canonical_v8/prospective_rounds/prospective_holdout_v2
```

Build the release-ready declaration archive bundle with:

```bash
uv run python scripts/build_archive_release_bundle.py
```

The diagnostic interface for a new signature is:

```bash
uv run python -m signature_durability_benchmark.cli triage \
  --config config/benchmark_config.yaml \
  --input my_signature.tsv \
  --out outputs/my_signature
```

The input file must contain `gene_symbol` and may optionally contain `direction` and `weight`. The CLI writes `diagnostic.json`, `per_cohort_effects.csv`, and `diagnostic_summary.md`.

## Conclusion

Program-conditioned validation changes the interpretation of transcriptomic signature heterogeneity. In this benchmark, interferon signatures are reproducible within 11 IFN-engaged cohorts and shrink to small, non-significant effects outside them; that finding survives orthogonal and strictly-disjoint holdouts, leave-one-cohort-out withholding, split-half replication, and a held-out four-cohort bulk RNA-seq extension. A deterministic second breadth case then shows how the method generalizes without pretending all programs are equally clean: Hallmark Hypoxia wins the frozen non-IFN selection rule and is rank 1 in an external exact-perturbation cohort, yet it still remains mixed because foreign-program bleed-through is real. The metadata-first prospective layer is intentionally less tidy: two predeclared held-out cohorts confirm the v1 IFN forecast, while one acute PBMC cohort does not. That mixed result is still useful because it demonstrates that the diagnostic is operating as a genuine challenge rather than an automatic approval machine. By adding a second fresh, externally timestamped, still-unevaluated v2 registry and a generic locked-round scoring command, the repo also makes the next prospective step auditable without rewriting the v1 miss. More importantly, the same diagnostic exposes when coarse labels, rather than bad signatures, are the source of ambiguity. Reported I² values for transcriptomic signatures should therefore be accompanied by program-conditioned decomposition and explicit failure-mode analysis. Cross-context failure of an interferon signature in a non-interferon cohort is not evidence that the signature is broken, and a “mixed” pooled result need not be the end of the analysis.

The new Ayers IFNγ portability case sharpens that point further: even a published external signature from a different application domain can be rescued by `triage` as a durable interferon-context signal rather than dismissed as generically irreproducible.

## References

1. Schoggins JW. Interferon-stimulated genes: what do they all do? Annu Rev Virol. 2019;6:567–584. doi:10.1146/annurev-virology-092818-015756
2. Schoggins JW, Wilson SJ, Panis M, et al. A diverse range of gene products are effectors of the type I interferon antiviral response. Nature. 2011;472:481–485. doi:10.1038/nature09907
3. Liberzon A, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Systems. 2015;1:417–425. doi:10.1016/j.cels.2015.12.004
4. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7:177–188. doi:10.1016/0197-2456(86)90046-2
5. Hartung J, Knapp G. On tests of the overall treatment effect in meta-analysis with normally distributed responses. Stat Med. 2001;20:1771–1782. doi:10.1002/sim.791
6. Sidik K, Jonkman JN. A simple confidence interval for meta-analysis. Stat Med. 2002;21:3153–3159. doi:10.1002/sim.1262
7. Hedges LV, Olkin I. Statistical Methods for Meta-Analysis. Academic Press; 1985. doi:10.1016/c2009-0-03396-0
8. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. PLoS Comput Biol. 2011;7:e1002240. doi:10.1371/journal.pcbi.1002240
9. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. Wiley; 2009. doi:10.1002/9780470743386
10. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Stat Med. 2002;21:1539–1558. doi:10.1002/sim.1186
11. Baechler EC, Batliwalla FM, Karypis G, et al. Interferon-inducible gene expression signature in peripheral blood cells of patients with severe lupus. Proc Natl Acad Sci U S A. 2003;100:2610–2615. doi:10.1073/pnas.0337679100
12. Yao Y, Richman L, Morehouse C, et al. Type I interferon: potential therapeutic target for psoriasis? PLoS One. 2008;3:e2737. doi:10.1371/journal.pone.0002737
13. Ayers M, Lunceford J, Nebozhyn M, et al. IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127:2930–2940. doi:10.1172/JCI91190


============================================================
SKILL.md
============================================================
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
