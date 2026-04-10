# Prospective Holdout Prediction Protocol (v1)

Freeze time (UTC): `2026-04-10T05:09:01Z`

This directory defines a metadata-first held-out prediction round for the IFN
diagnostic. The goal is to separate prediction from evaluation as cleanly as
possible inside the repository itself.

## What Was Allowed Before Declaration

- GEO series titles and summaries
- GEO series sample annotations from the series-matrix metadata
- GEO supplementary file listings, file names, and archive structure
- Cohort-level facts needed to define the contrast:
  - tissue
  - case/control labels
  - expected sample counts
  - file acquisition mode

## What Was Not Part of the Declaration

- signature scoring on the held-out cohorts
- pooled or per-cohort Hedges' g values
- permutation or meta-analysis results on the held-out cohorts
- any claim that comparator signatures must replicate

## Primary Prediction

Primary panel: the three cohorts listed in
`prediction_registry_v1.tsv`.

Primary IFN quartet:

- `hallmark_ifng_response`
- `hallmark_ifna_response`
- `schoggins_2011_irg`
- `blind_durable_ifn_composite`

Success criteria:

1. Each of the four IFN signatures should have a positive per-cohort effect
   (`g > 0`) in each primary held-out cohort.
2. Across the three-cohort primary panel, each IFN signature should have:
   - pooled guarded-HKSJ `g >= 0.8`
   - Bonferroni-4 adjusted `p < 0.05`
   - sign consistency `= 1.0`

## Secondary Context

- Comparator signatures may be scored for context, but they are not part of the
  primary prospective claim.
- The evaluation script must report the SHA256 of the frozen registry and the
  download audit for all newly acquired files.

## Notes

- This protocol is repository-local rather than externally preregistered.
- The intended claim is still narrower than a clinical forecast:
  predictions are declared from GEO metadata before held-out cohort scoring.
