# External RNA-seq Validation

- Primary bulk RNA-seq panel: `4` cohorts, `719` samples
- Primary cohorts: `GSE152641`, `GSE171110`, `GSE152075`, `GSE167000`

## IFN Quartet

- `IFN-gamma core`: pooled g=`0.922`, guarded p=`0.003196`, Bonferroni-7 p=`0.02237`, I²=`0.000`, sign consistency=`1.000`
- `IFN-alpha core`: pooled g=`1.193`, guarded p=`0.001582`, Bonferroni-7 p=`0.01107`, I²=`0.000`, sign consistency=`1.000`
- `Schoggins 2011 IRG`: pooled g=`0.996`, guarded p=`0.002605`, Bonferroni-7 p=`0.01824`, I²=`0.000`, sign consistency=`1.000`
- `Blind IFN composite`: pooled g=`1.177`, guarded p=`0.002912`, Bonferroni-7 p=`0.02038`, I²=`0.244`, sign consistency=`1.000`

## Comparator Signatures

- `Inflammatory Response`: pooled g=`0.427`, guarded p=`0.17`, Bonferroni-7 p=`1`, I²=`0.770`
- `TNF-alpha/NF-kB`: pooled g=`0.382`, guarded p=`0.2949`, Bonferroni-7 p=`1`, I²=`0.859`
- `E2F Targets`: pooled g=`0.818`, guarded p=`0.1964`, Bonferroni-7 p=`1`, I²=`0.944`

## Exploratory PBMC Stress Test

- `gse152418_pbmc_covid_vs_healthy` (GSE152418): IFN-gamma g=`-0.409`, IFN-alpha g=`0.127`, Schoggins g=`-0.203`, blind IFN g=`-0.078`
