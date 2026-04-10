# Second Positive Generalization Case

- Selected signature: `hallmark_hypoxia`
- Home program: `hypoxia`
- Ranking expectation: `hallmark_hypoxia`
- Selection rule satisfied: `True`

## Candidate Ranking

| Rank | Signature | p_Bonf<0.05 | LOO=1.0 | Split=1.0 | mean abs foreign d | within-minus-abs(outside) gap |
|---:|---|---:|---:|---:|---:|---:|
| 1 | hallmark_hypoxia | True | True | True | 1.102 | 2.733 |
| 2 | hallmark_emt | False | True | True | 0.501 | 2.064 |

## Selected Case Checks

- Inferred program matches home program: `True`
- Within-program class: `mixed`
- Full-model class: `mixed`
- Positive statistically supported pooled effect (frozen artifact): `True`
- Leave-one-out sign prediction: `1.000`
- Split-half sign agreement: `1.000`

## Home-Program Cohorts

| Cohort | Hedges' g | Contrast |
|---|---:|---|
| hypoxia_timecourse_gse47533 | +10.095 | hypoxia vs normoxia |
| ccrcc_kidney_gse36895 | +6.659 | hypoxia vs normoxia |
| hypoxia_cellline_gse53012 | +3.521 | hypoxia vs normoxia |
| hypoxia_multicell_gse18494 | +1.539 | hypoxia vs normoxia |
| hypoxia_mcf7_gse3188 | +1.257 | hypoxia vs normoxia |
| gbm_brain_gse4290 | +0.438 | hypoxia vs normoxia |

## Largest Foreign-Program Effects

| Cohort | Foreign program | Hedges' g | Contrast |
|---|---|---:|---|
| emt_mammary_gse43495 | emt | -4.933 | mesenchymal vs epithelial |
| melioidosis_blood_gse69528 | inflammation | +2.839 | sepsis vs healthy |
| sepsis_shock_gse95233 | inflammation | +2.684 | sepsis vs healthy |
| sepsis_blood_gse28750 | inflammation | +2.407 | sepsis vs healthy |
| psoriasis_skin_gse13355 | interferon | +1.973 | psoriasis vs healthy |
| lung_tumor_gse19188 | proliferation | +1.766 | tumor vs normal |

## External Support

- Cohort: `gse179885_tcells_hypoxia_vs_normoxia` (GSE179885)
- Chosen-signature external effect: `+5.113`
- Rank among externally scored signatures: `1`
- Positive direction: `True`

## Triage Interpretation

# Diagnostic Report: hallmark_hypoxia

- Input: `hallmark_hypoxia.tsv`
- Inferred best-supported program: `hypoxia`
- Full-model classification: `mixed`
- Within-program classification: `mixed`
- Mean coverage: `0.960`
- Aggregate effect: `0.492` (p=`1.048e-48`)
- Program-structure permutation p: `0.5175`

hallmark_hypoxia shows real signal but does not reduce cleanly to one program. The aggregate effect is 0.492 with I²=0.941, suggesting either cross-context activation or overly broad program boundaries.

## Program Ranking

- `hypoxia`: separation=2.733, within=3.378 (p=0.05725), outside=0.645 (p=0.006584)
- `inflammation`: separation=0.707, within=1.448 (p=0.01303), outside=0.741 (p=0.02944)
- `interferon`: separation=-0.178, within=0.791 (p=0.002678), outside=0.968 (p=0.03133)
- `proliferation`: separation=-0.650, within=0.344 (p=0.506), outside=0.994 (p=0.002776)
- `emt`: separation=-2.888, within=-0.738 (p=0.3641), outside=1.150 (p=0.0002243)
