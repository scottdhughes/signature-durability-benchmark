# Provenance Audit

- Status: `pass_with_notes`
- Active scored cohort panel: `35` cohorts, `5922` samples
- Frozen signatures: `30` IDs
- Synthetic control signatures in broader benchmark: `19`

## Findings

- Config and frozen cohort manifests are synchronized at 35 cohorts / 5,922 samples.
- All paper-facing target signatures are biologically grounded, non-synthetic panels.
- The broader 30-signature benchmark still contains 19 explicit synthetic control signatures; they are benchmark controls, not headline evidence.
- No unexpected synthetic/stub/mock keywords were found in the runtime code, paper, or skill surface.

## Real Cohort Checks

- Config vs frozen cohort manifest synchronized: `True`
- Phenotype row counts match manifest sample counts: `True`
- Matrix column counts match manifest sample counts: `True`
- Missing matrix assets: `[]`
- Missing phenotype assets: `[]`
- Unused frozen assets: `[]`

## Signature Provenance Checks

- Config and frozen signature manifests synchronized: `True`
- Paper target signatures all manifested: `True`
- Paper target signatures all non-synthetic: `True`
- Missing from config signature manifest: `[]`
- Missing from frozen signature manifest: `[]`

## Output Checks

- `per_cohort_effects.csv` rows: `1050` / expected `1050`
- Unique cohorts in scored output: `35`
- Unique signatures in scored output: `30`
- Per-cohort rows match expected grid: `True`

## External Layer Checks

- External RNA-seq manifest rows: `5`
- External hypoxia manifest rows: `1`
- Prospective v1 rows: `3`
- Prospective v2 rows: `3`
- Prospective v2 receipt verification: `OK`

## Keyword Audit

- Unexpected synthetic/stub/mock hits outside allowed provenance files: `0`
