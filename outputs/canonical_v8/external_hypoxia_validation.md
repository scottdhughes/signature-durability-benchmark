# External Hypoxia Validation

- Cohort: `gse179885_tcells_hypoxia_vs_normoxia` (GSE179885)
- Tissue: `t_cells`
- Contrast: Human T-cell bulk RNA-seq cultured in hypoxia versus normoxia
- Search protocol: `data/external_hypoxia/SEARCH_PROTOCOL.md`

## Ranked Signature Effects

| Rank | Signature | Hedges' g | Coverage |
|---:|---|---:|---:|
| 1 | hallmark_hypoxia | +5.113 | 1.000 |
| 2 | hallmark_emt | +0.940 | 0.867 |
| 3 | IFN-alpha core | -1.236 | 0.967 |
| 4 | E2F Targets | -1.249 | 0.967 |
| 5 | TNF-alpha/NF-kB | -1.365 | 0.900 |
| 6 | IFN-gamma core | -1.635 | 0.967 |
| 7 | Inflammatory Response | -2.290 | 0.933 |

This is an exact-perturbation external layer rather than a main-table pooled result because the cohort has 12 total samples, below the preferred >=20-sample rule.
