# Archive Release Runbook

This repository is prepared for a two-stage public archive flow:

1. `prospective_holdout_v2_declaration`
2. later `prospective_holdout_v2_evaluation`

The declaration bundle exists now at:

- `submission/archive_bundles/prospective_holdout_v2_declaration/`

The evaluation bundle must be built only after the locked v2 readout exists.

## Declaration Release

Publish the following as a declaration-only release:

- `prediction_registry_v2.tsv`
- `PREDICTION_PROTOCOL_v2.md`
- `prospective_holdout_v2/` timestamp directory
- `CHECKSUMS.sha256`
- `RELEASE_NOTES.md`
- `bundle_manifest.json`

Recommended public targets:

- GitHub Release asset bundle on `scottdhughes/signature-durability-benchmark`
- Zenodo versioned archive linked to the GitHub repository

Suggested release title:

- `prospective-holdout-v2-declaration`

Suggested release description:

- declaration-only archive for the second prospective round
- externally timestamped RFC3161 receipt included
- intentionally unevaluated at release time
- supersedes no prior scored result

## Later Evaluation Release

After the locked v2 scoring pass exists, publish a second distinct archive containing:

- `outputs/canonical_v8/prospective_rounds/prospective_holdout_v2/`
- any download audit manifest produced during scoring
- the exact paper/repo revision used for the v2 readout
- a checksum manifest for the release bundle

Do not overwrite or mutate the declaration release. The evaluation release must remain a separate public object with its own DOI/version.

## Citation Policy

- Cite the declaration DOI once minted in repo/docs as the predeclared prospective record.
- Cite the evaluation DOI only after the locked v2 readout has been published.
- Prefer version-specific Zenodo DOIs in manuscript or submission text.

## Local Preconditions

- `LICENSE`, `CITATION.cff`, and `.zenodo.json` are already present at repo root.
- `origin` should point to:
  - `https://github.com/scottdhughes/signature-durability-benchmark.git`
- The current paper should continue to describe v2 as declared and unevaluated until the evaluation release exists.
