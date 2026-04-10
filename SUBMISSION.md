# Submission Notes

Paper title: "Program-Conditioned Diagnostic for Transcriptomic Signature Durability: Validation on Interferon Signatures across 35 Frozen GEO Cohorts"

Authors: Karen Nguyen, Scott Hughes, Claw (corresponding co-author)
Agent: Longevist (@longevist)

## Core submission state

- The current paper already includes:
  - the program-conditioned diagnostic framing
  - IFN validation across the frozen 35-cohort panel
  - a rescued published-signature portability case using the Ayers IFN-gamma-related 6-gene profile
  - held-out bulk RNA-seq IFN transfer
  - a deterministic second breadth case (`hallmark_hypoxia`)
  - a generated workflow figure (`paper/figure_workflow.pdf`)
  - a scored prospective v1 challenge round
  - an externally timestamped but intentionally unscored prospective v2 declaration
- Do **not** score v2 before a separate locked readout pass. The current submission should keep the paper wording as "declared and unevaluated".

## Submission prep

1. Verify canonical run outputs in `outputs/canonical_v8/`
2. Verify rescued-signature portability case in `outputs/canonical_v8/rescued_signature_case_study.json`
3. Verify held-out external RNA-seq validation in `outputs/canonical_v8/external_rnaseq_validation.json`
4. Verify deterministic breadth case outputs:
   - `outputs/canonical_v8/generalization_case_study.json`
   - `outputs/canonical_v8/external_hypoxia_validation.json`
5. Rebuild workflow figure:
   - `paper/figure_workflow.pdf`
   - `paper/figure_workflow.png`
6. Build the declaration-only archive bundle for v2:
   - `submission/archive_bundles/prospective_holdout_v2_declaration/`
7. Build clawRxiv payload from the current paper and SKILL
8. Rebuild conference PDF from `paper/main.tex`
9. Submit via API

## Public archive hardening

- Canonical repo target:
  - `https://github.com/scottdhughes/signature-durability-benchmark`
- Metadata shipped locally:
  - `LICENSE`
  - `CITATION.cff`
  - `.zenodo.json`
- Declaration-only public archive to publish now:
  - `prediction_registry_v2.tsv`
  - `PREDICTION_PROTOCOL_v2.md`
  - `data/prospective_holdout/external_timestamps/prospective_holdout_v2/`
  - `submission/archive_bundles/prospective_holdout_v2_declaration/CHECKSUMS.sha256`
- Publish targets:
  - GitHub Release asset bundle
  - Zenodo versioned archive
- Local runbook:
  - `submission/ARCHIVE_RELEASES.md`
- Cite the declaration DOI once minted.
- The later v2 evaluation bundle must be published as a separate release after the locked scoring pass exists.
