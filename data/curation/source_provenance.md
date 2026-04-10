# Source Provenance

This benchmark uses a curated panel of gene signatures anchored to primary literature and public gene-set databases. Each entry below records the anchor source, extraction rule, and freeze rationale. Signature source families are deliberately non-overlapping with GEO cohort source families to prevent source leakage.

The conceptual backbone is the MSigDB Hallmark gene-set collection (Liberzon et al., "The Molecular Signatures Database Hallmark Gene Set Collection," Cell Systems 2015), supplemented by curated senescence literature and synthetic negative-control signatures designed to test specific failure modes.

## Durable Signatures

### hallmark_ifng_response
Anchor family: MSigDB Hallmark INTERFERON_GAMMA_RESPONSE. Resource location: MSigDB v2024.2.Hs, gene set M5913. Extraction rule: frozen 30-gene benchmark release definition anchored to the broader Hallmark family and inherited unchanged from the original benchmark release, all direction "up" with uniform weight 1.0 (unsigned enrichment set). Exact membership is enumerated in `data/freeze/signatures.tsv`; this release definition was frozen before the expanded 35-cohort reruns and was not reselected from expanded-panel performance. Key genes include STAT1, IRF1, GBP1/2, CXCL9/10, IDO1, OAS1/2/3, ISG15, and MHC-II components. Frozen because IFNg response is among the most robustly replicated immune signatures across cancer, infection, and autoimmune cohorts.

### hallmark_ifna_response
Anchor family: MSigDB Hallmark INTERFERON_ALPHA_RESPONSE. Resource location: MSigDB v2024.2.Hs, gene set M5911. Extraction rule: frozen 30-gene benchmark release definition anchored to the broader Hallmark family and inherited unchanged from the original benchmark release, direction "up" with uniform weight 1.0. Exact membership is enumerated in `data/freeze/signatures.tsv`; this release definition was frozen before the expanded reruns and was not reselected from expanded-panel performance. Key genes include IFIT1/2/3, ISG15, MX1/2, OAS1/2/3, RSAD2, IFI44/44L, STAT1/2, IRF7/9. Frozen because type-I interferon response is a canonical antiviral program with consistent cross-platform reproducibility.

### hallmark_inflammatory_response
Anchor family: MSigDB Hallmark INFLAMMATORY_RESPONSE. Resource location: MSigDB v2024.2.Hs, gene set M5932. Extraction rule: top 30 genes from the ~200-gene Hallmark set, direction "up" with uniform weight 1.0. Key genes include IL6, IL1B, TNF, CXCL8, CCL2/3/4, PTGS2, ICAM1, VCAM1, NF-kB pathway members. Frozen because acute inflammatory response is well-validated across sepsis, autoimmune, and cancer inflammation cohorts.

### hallmark_hypoxia
Anchor family: MSigDB Hallmark HYPOXIA. Resource location: MSigDB v2024.2.Hs, gene set M5891. Extraction rule: top 30 genes from the ~200-gene Hallmark set, direction "up" with uniform weight 1.0. Key genes include VEGFA, SLC2A1, LDHA, glycolytic enzymes (PGK1, ENO1, ALDOA, PKM, HK2), HIF1A, CA9, LOX. Frozen because HIF1A-regulated metabolic adaptation is consistently activated across solid tumors and ischemic tissues.

### hallmark_e2f_targets
Anchor family: MSigDB Hallmark E2F_TARGETS. Resource location: MSigDB v2024.2.Hs, gene set M5925. Extraction rule: top 30 genes from the ~200-gene Hallmark set, direction "up" with uniform weight 1.0. Key genes include MCM2-7, PCNA, RRM1/2, TYMS, TK1, replication fork components, and mitotic regulators. Frozen because E2F/cell-cycle genes are the most consistently reproducible proliferation signature across cancer types.

### hallmark_emt
Anchor family: MSigDB Hallmark EPITHELIAL_MESENCHYMAL_TRANSITION. Resource location: MSigDB v2024.2.Hs, gene set M5930. Extraction rule: top 30 genes from the ~200-gene Hallmark set, direction "up" with uniform weight 1.0. Key genes include VIM, CDH2, FN1, SNAI1/2, ZEB1/2, TWIST1, collagens, MMPs, ACTA2, SPARC, POSTN. Frozen because EMT/mesenchymal markers are robustly reproducible across fibrosis, metastasis, and wound-healing contexts.

### hallmark_tnfa_nfkb
Anchor family: MSigDB Hallmark TNFA_SIGNALING_VIA_NFKB. Resource location: MSigDB v2024.2.Hs, gene set M5890. Extraction rule: top 30 genes from the ~200-gene Hallmark set, direction "up" with uniform weight 1.0. Key genes include NFKBIA, TNFAIP3, BIRC3, TRAF1, RELB, immediate-early transcription factors (JUNB, ATF3, EGR1, KLF6), and NF-kB target chemokines. Frozen because TNFa/NF-kB immediate-early response is a canonical and highly reproducible inflammatory program.

### curated_senescence
Anchor family: cellular senescence consensus literature. Resource location: Coppe et al. "Senescence-Associated Secretory Phenotypes Reveal Cell-Nonautonomous Functions of Oncogenic RAS and the p53 Tumor Suppressor" (PLoS Biology 2008); Hernandez-Segura et al. "Hallmarks of Cellular Senescence" (Trends in Cell Biology 2018); Casella et al. "Transcriptome signature of cellular senescence" (NAR 2019). Extraction rule: 20 consensus senescence genes including CDK inhibitors (CDKN2A, CDKN1A, CDKN2B), tumor suppressors (TP53, RB1), SASP cytokines (IL6, CXCL8, CCL2), SASP proteases (MMP3, SERPINE1), senescence markers (GLB1, LMNB1-down, H2AFX, HMGA1/2), and apoptosis regulators (BCL2L1, BAX). LMNB1 is direction "down" per Freund et al. "Lamin B1 loss is a senescence-associated biomarker" (MBoC 2012). Frozen because this curated set captures the consensus senescence phenotype validated across fibroblast, epithelial, and in vivo aging models.

## Brittle Signatures

### brittle_small_study_inflammation
Anchor family: synthetic brittle control mimicking underpowered single-study extraction. Extraction rule: random subset of 15 inflammatory-response genes (IL1B, CXCL2, CCL4, PTGS2, ICAM1, SOD2, PLAUR, CD44, CSF2, IL1RN, CXCL3, SELE, TIMP1, BCL3, TRAF1) plus 5 housekeeping/non-specific genes (ACTB, B2M, HPRT1, YWHAZ, SDHA) that add noise without biological coherence. Frozen because this pattern (partial pathway + housekeeping noise) represents a common failure mode in small-N microarray studies.

### brittle_platform_specific
Anchor family: synthetic brittle control for platform-dependent expression artifacts. Extraction rule: 10 ribosomal protein genes (RPL13A, RPL7, RPL3, RPL11, RPL5, RPS4X, RPS27A, RPS3, RPS8, RPS18) plus 8 mitochondrial-encoded genes (MT-ND1, MT-ND4, MT-CO1, MT-CO2, MT-ATP6, MT-CYB, MT-RNR1, MT-RNR2) plus 2 long non-coding RNAs (MALAT1, NEAT1). Frozen because ribosomal/mitochondrial gene ratios vary dramatically between array and RNA-seq platforms and across library preparation methods, producing strong batch effects rather than biology.

### brittle_single_tissue
Anchor family: synthetic brittle control for tissue-specificity failure. Extraction rule: 13 liver-specific genes (ALB, APOA1, APOB, APOC3, CYP3A4, CYP2E1, CYP1A2, HP, TF, FGA, FGB, FGG, SERPINA1) plus 2 acute-phase proteins (SAA1, CRP) plus 5 generic immune markers (CD3D, CD8A, IL6, TNF, PTPRC). Frozen because liver-expressed genes are absent or negligible in non-hepatic tissues; mixing them with ubiquitous immune markers creates an incoherent signature that only works in liver cohorts.

### brittle_overfit_noise
Anchor family: synthetic brittle control for overfit/noise selection. Extraction rule: 15 genes selected for large genomic size or high somatic mutation frequency (ANKRD36, ZNF518A, KIAA1549, OBSCN, SYNE1, TTN, DNAH5, MUC16, CSMD3, RYR2, LRP1B, USH2A, PKHD1, ZFHX4, HMCN1) with random direction assignments. Frozen because these genes share no transcriptional program or regulatory coherence; any apparent cross-cohort enrichment would be a false positive.

## Mixed Signatures

### mixed_emt_inflammation
Anchor family: synthetic mixed signature combining EMT and inflammatory biology. Extraction rule: 10 EMT markers (VIM, FN1, SNAI2, ZEB1, COL1A1, COL3A1, ACTA2, SPARC, FAP, LOXL2) plus 10 inflammatory cytokines/effectors (IL6, IL1B, TNF, CXCL8, CCL2, PTGS2, MMP9, MMP3, SERPINE1, TIMP1). Frozen because this combination replicates in fibrosis and wound healing (where EMT and inflammation co-occur) but fails in pure sepsis or immune-activation contexts where EMT is absent.

### mixed_hypoxia_stress
Anchor family: synthetic mixed signature combining hypoxia and stress/UPR response. Extraction rule: 10 hypoxia genes (VEGFA, SLC2A1, LDHA, CA9, BNIP3, ADM, NDRG1, LOX, HIF1A, PDK1) plus 10 stress/UPR genes (HSPA1A, HSPA1B, HSP90AA1, DDIT3, ATF4, XBP1, GADD45A, CDKN1A, FOS, JUN). Frozen because hypoxia + stress co-occur in solid tumors but dissociate in hematologic malignancies and non-hypoxic inflammatory settings, producing inconsistent cross-cohort scores.

### mixed_senescence_proliferation
Anchor family: synthetic mixed signature combining senescence arrest and residual proliferation markers. Extraction rule: 10 senescence markers (CDKN2A, CDKN1A, TP53, GLB1, SERPINE1, IL6, CXCL8, IGFBP3, RB1, GADD45A) plus 10 proliferation genes (MKI67, TOP2A, PCNA, CCNB1, CDK1, CCNA2, BUB1, AURKA, PLK1, CDC20). Frozen because senescence and proliferation are biologically contradictory (CDKN2A/1A arrest vs. CDK1/CCNB1 cycling); co-elevation occurs only in specific contexts (oncogene-induced senescence bypass, senescence-adjacent proliferating cells in tumors) and produces unreliable cross-cohort behavior.

## Confounded Signatures

### confounded_proliferation
Anchor family: synthetic negative control — pure proliferation/mitotic biology. Extraction rule: 20 core cell-cycle/mitotic genes (MKI67, TOP2A, PCNA, CCNB1, CDK1, CCNA2, BUB1, AURKA, PLK1, CDC20, CCNB2, CDK2, MCM2, BIRC5, CENPE, KIF11, TTK, NEK2, MELK, FOXM1), all direction "up" with uniform weight. Frozen because pure proliferation signatures dominate any growth-vs-quiescence contrast but reflect nuisance biology (sample proliferative fraction) rather than specific pathway activation. Confounder panel overlap with proliferation_cell_cycle set should trigger rejection.

### confounded_ribosomal
Anchor family: synthetic negative control — ribosomal protein housekeeping genes. Extraction rule: 20 ribosomal protein genes (10 RPL + 10 RPS family members), all direction "up" with uniform weight. Frozen because ribosomal protein expression tracks cell size, translation rate, and sample quality rather than specific biological programs. Confounder panel overlap with batch_ribosomal set should trigger rejection.

### confounded_immune_infiltrate
Anchor family: synthetic negative control — immune cell-type deconvolution markers. Extraction rule: 20 immune lineage markers spanning T cells (CD3D, CD3E, CD8A, CD8B, CD4, FOXP3), B cells (MS4A1, CD19), monocytes/macrophages (CD14, CD68, CSF1R, ITGAM, ITGAX), NK cells (NKG7, GZMB, GZMA, PRF1, KLRD1), and pan-leukocyte (PTPRC, FCGR3A). Frozen because these markers track immune cell fraction in bulk tissue rather than pathway activation within cells. Confounder panel overlap with immune_infiltration set should trigger rejection.

### stealth_confounded_hypoxia_prolif
Anchor family: synthetic stealth confounded — real hypoxia signal contaminated by proliferation confounder. Extraction rule: 15 core hypoxia genes (VEGFA, SLC2A1, LDHA, CA9, BNIP3, ADM, NDRG1, LOX, HIF1A, PDK1, ANGPTL4, P4HA1, EGLN3, DDIT4, SLC2A3) at weight 1.0 plus 5 proliferation markers (MKI67, TOP2A, CDK1, CCNB1, AURKA) at weight 0.8. Frozen because simple enrichment models see a strong hypoxia signal and score this as durable, but the confounder overlap check detects that 25% of the signature (by gene count) overlaps with the proliferation_cell_cycle confounder set. This is the primary mechanism by which the full model demonstrates value over enrichment-only baselines.

### stealth_confounded_inflam_immune
Anchor family: synthetic stealth confounded — real inflammatory signal contaminated by immune infiltration confounder. Extraction rule: 15 core inflammatory genes (IL6, IL1B, TNF, CXCL8, CCL2, CCL3, PTGS2, ICAM1, NFKBIA, CXCL1, VCAM1, MMP9, SERPINE1, SOD2, TNFAIP3) at weight 1.0 plus 5 immune infiltration markers (CD3D, CD8A, PTPRC, NKG7, GZMB) at weight 0.8. Frozen because the inflammatory core is genuinely durable, but the immune-fraction markers cause the signature to track immune cell proportion in bulk tissues rather than intracellular inflammatory activation. Confounder check catches the immune_infiltration overlap.

## Insufficient Coverage Sentinels

### lowcov_olfactory_1
Anchor family: out-of-universe sentinel — olfactory receptor genes. Extraction rule: 15 olfactory receptor genes (OR1A1, OR2W1, OR7D4, OR10A1, OR51E1, OR6A2, OR2T8, OR13C9, OR5K1, OR2B6, OR4D1, OR8D1, OR11H1, OR14A16, OR52N1). Frozen because olfactory receptors are expressed exclusively in olfactory epithelium and are absent from standard expression arrays and bulk tissue RNA-seq. Coverage threshold should reject these outright.

### lowcov_olfactory_2
Anchor family: out-of-universe sentinel — taste receptor and olfactory receptor genes. Extraction rule: 5 taste receptor genes (TAS2R38, TAS2R16, TAS2R46, TAS1R1, TAS1R2) plus 10 olfactory receptor genes (OR3A1, OR1D2, OR5A1, OR2J3, OR10G4, OR12D3, OR56A4, OR52I2, OR4E2, OR8B8). Frozen because taste/olfactory receptors have negligible expression in standard cohort tissues, providing an independent coverage rejection sentinel.

## Blind Panel

### blind_durable_ifn_composite
Anchor family: composite interferon signature drawn from the intersection of IFNg and IFNa hallmark cores. Extraction rule: 25 genes combining the most-replicated interferon-stimulated genes from both type-I and type-II response programs (STAT1, IRF1, IFIT1/2/3, ISG15, MX1/2, OAS1/2/3, RSAD2, GBP1/2, CXCL9/10, IDO1, IFI44/44L, IRF7, IFITM1, BST2, DDX58, HERC5, USP18). Held out from threshold tuning. Frozen because it is expected to replicate robustly.

### blind_durable_hypoxia_core
Anchor family: core HIF1A-target hypoxia signature. Extraction rule: 20 genes restricted to the most consistently HIF1A-regulated targets across solid tumor and ischemia cohorts (VEGFA, SLC2A1, LDHA, PGK1, ENO1, ALDOA, PKM, HK2, CA9, ADM, NDRG1, BNIP3, LOX, HIF1A, EGLN1, PDK1, PFKFB3, ANGPTL4, P4HA1, DDIT4). Held out from threshold tuning. Frozen because it is expected to replicate robustly.

### blind_brittle_random
Anchor family: synthetic blind brittle control — neuronal/receptor genes with no coherent expression program in standard bulk-tissue cohorts. Extraction rule: 15 genes (ABCA7, GPR37, SLC17A7, CHRNA7, GABRA1, HTR2A, GRIN2B, DRD2, SCN1A, KCNJ11, TRPM8, PCDH15, CDH23, TENM4, SLITRK6) with random direction assignments. Held out from threshold tuning. Frozen because no biological coherence is expected.

### blind_mixed_emt_stress
Anchor family: synthetic blind mixed signature combining EMT markers with stress/UPR genes. Extraction rule: 10 EMT genes (VIM, FN1, SNAI1, ZEB2, COL1A2, COL5A1, ACTA2, TAGLN, POSTN, TGFBI) plus 10 stress genes (HSPA1A, HSPA1B, HSP90AA1, DDIT3, ATF4, XBP1, FOS, JUN, GADD45A, CDKN1A). Held out from threshold tuning. Frozen because EMT + stress co-occur in wound healing and fibrosis but dissociate in other contexts.

### blind_confounded_prolif_inflam
Anchor family: synthetic blind confounded control — proliferation plus inflammation. Extraction rule: 10 proliferation genes (MKI67, TOP2A, PCNA, CCNB1, CDK1, CCNA2, BUB1, AURKA, PLK1, CDC20) plus 10 inflammatory genes (IL6, IL1B, TNF, CXCL8, CCL2, PTGS2, ICAM1, NFKBIA, SOD2, TNFAIP3). Held out from threshold tuning. Frozen because confounder overlap with both proliferation_cell_cycle and generic_inflammation sets should trigger rejection.

### blind_confounded_stealth_emt_immune
Anchor family: synthetic blind stealth confounded — EMT signal plus immune infiltration markers. Extraction rule: 15 EMT genes (VIM, FN1, SNAI2, ZEB1, COL1A1, COL3A1, ACTA2, SPARC, FAP, LOXL2, POSTN, TGFBI, CTGF, SERPINE1, THY1) at weight 1.0 plus 5 immune markers (CD3D, CD8A, PTPRC, NKG7, GZMB) at weight 0.8. Held out from threshold tuning. Frozen because the EMT core is genuine but the immune-fraction markers contaminate the signature in bulk tissue.

### blind_lowcov_keratin
Anchor family: out-of-universe sentinel — hair keratin genes. Extraction rule: 15 hair-specific keratin genes (KRT31-38, KRT81-86) absent from standard expression platforms and bulk non-skin tissues. Held out from threshold tuning. Frozen because hair keratins provide an independent out-of-universe coverage rejection test.

## Confounder Panel Sources

### proliferation_cell_cycle
Source: consensus cell-cycle markers from Whitfield et al. "Identification of genes periodically expressed in the human cell cycle" (MBoC 2002) and MSigDB Hallmark E2F/G2M sets. Genes: MKI67, TOP2A, PCNA, CCNB1, CDK1, CCNA2, BUB1, AURKA, PLK1, CDC20.

### generic_inflammation
Source: MSigDB Hallmark INFLAMMATORY_RESPONSE and TNFA_SIGNALING_VIA_NFKB core overlap. Genes: IL6, IL1B, TNF, CXCL8, CCL2, NFKB1, PTGS2, ICAM1, VCAM1, SELE.

### stress_response
Source: consensus heat-shock and UPR markers from Vihervaara et al. "HSF1 at a glance" (J Cell Sci 2014) and ER stress literature. Genes: HSPA1A, HSPA1B, HSP90AA1, DDIT3, ATF4, XBP1, GADD45A, CDKN1A, FOS, JUN.

### batch_ribosomal
Source: ribosomal protein genes commonly flagged as batch-effect drivers in microarray QC literature (McCall et al., "Frozen robust multiarray analysis," Biostatistics 2010). Genes: RPS4X, RPL13A, RPS27A, RPL7, RPS3, RPL3, RPS8, RPL11, RPS18, RPL5.

### immune_infiltration
Source: CIBERSORTx immune cell-type markers (Newman et al., "Determining cell type abundance and expression from bulk tissues with digital cytometry," Nature Biotechnology 2019). Genes: PTPRC, CD3D, CD8A, CD4, MS4A1, CD14, ITGAM, FCGR3A, NKG7, GZMB.
