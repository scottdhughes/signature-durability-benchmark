# Prospective Holdout Prediction

- Registry: `prediction_registry_v1.tsv`
- Declared at: `2026-04-10T05:09:01Z`
- Evaluation started: `2026-04-10T13:09:50.239288+00:00`
- Evaluation finished: `2026-04-10T13:09:55.071706+00:00`
- Registry SHA256: `0b758bbf83fb1784276679878d6d564db7dbe2c24e5cf622981510b38af848bf`
- Overall success: `False`

## Cohort-Level IFN Hits

| Cohort | Tissue | IFN sign hits | Rule satisfied |
|---|---|---:|---:|
| gse184610_mucosa_sarscov2_positive_vs_control | oronasopharyngeal_mucosa | 4/4 | True |
| gse243217_pbmc_covid19_vs_healthy | pbmc | 4/4 | True |
| gse202805_pbmc_acute_covid19_vs_healthy | pbmc | 0/4 | False |

## Pooled Primary Panel

| Signature | pooled g | guarded p_Bonf,4 | sign consistency | success |
|---|---:|---:|---:|---:|
| IFN-gamma core | 0.359 | 1 | 0.667 | False |
| IFN-alpha core | 0.657 | 0.942 | 0.667 | False |
| Schoggins 2011 IRG | 0.549 | 1 | 0.667 | False |
| Blind IFN composite | 0.576 | 0.9307 | 0.667 | False |
