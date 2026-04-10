# Prospective Holdout Prediction Protocol (v2)

Freeze time (UTC): `2026-04-10T05:32:41Z`

This protocol defines a second metadata-first held-out prediction round that is
kept explicitly separate from the mixed evaluated v1 result. The purpose of v2
is to prevent quiet cohort swapping after observing the v1 miss: the fresh
cohorts, success rule, and declaration manifest are frozen before any v2
signature scoring.

This round is externally timestamped through an RFC3161 time-stamp authority.
The declaration receipt is stored under
`data/prospective_holdout/external_timestamps/prospective_holdout_v2/`.

## What Was Allowed Before Declaration

- GEO series titles and summaries
- GEO series sample annotations from the series-matrix metadata
- GEO supplementary file listings, file names, and archive structure
- Cohort-level facts needed to define the contrast:
  - tissue or cell source
  - case/control labels
  - timepoint or subset restriction if present
  - expected sample counts
  - file acquisition mode

## What Was Not Part of the Declaration

- signature scoring on the v2 held-out cohorts
- pooled or per-cohort Hedges' g values
- permutation or meta-analysis results on the v2 held-out cohorts
- any claim that v2 improves or rescues the mixed v1 result

## Primary Prediction

Primary panel: the three cohorts listed in `prediction_registry_v2.tsv`.

Primary IFN quartet:

- `hallmark_ifng_response`
- `hallmark_ifna_response`
- `schoggins_2011_irg`
- `blind_durable_ifn_composite`

Primary round design:

- GSE243442 PBMC acute yellow fever versus controls
- GSE213168 whole-blood influenza infection versus healthy controls
- GSE155237 nasal RSV challenge, D3 infected versus D3 uninfected

Success criteria:

1. Each of the four IFN signatures should have a positive per-cohort effect
   (`g > 0`) in each primary held-out cohort.
2. Across the three-cohort primary panel, each IFN signature should have:
   - pooled guarded-HKSJ `g >= 0.8`
   - Bonferroni-4 adjusted `p < 0.05`
   - sign consistency `= 1.0`

## Status Rules

- v2 is **declared and pending evaluation**.
- v1 remains the only scored prospective round in this submission.
- No v2 result files should exist in `outputs/canonical_v8/` at declaration time.

## Notes

- This round is not presented as a result in the paper; it is an auditable
  future challenge.
- The external timestamp binds the declaration manifest, which in turn binds the
  registry and protocol hashes.
