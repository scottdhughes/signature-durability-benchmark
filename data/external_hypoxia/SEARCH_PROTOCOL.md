# External Hypoxia Validation Selection Protocol

This directory defines the deterministic external-validation layer for the
second non-IFN positive generalization case.

## Goal

Add one human bulk RNA-seq hypoxia-vs-normoxia cohort that is:

- not already in the frozen benchmark
- directly downloadable from the GEO series record
- scored from a processed matrix rather than reconstructed from raw FASTQ
- tied to a binary hypoxia contrast that matches the benchmark semantics

## Selection Rule

Candidate cohorts must satisfy all of the following:

1. `Homo sapiens`
2. GEO series record with `Expression profiling by high throughput sequencing`
3. explicit hypoxia versus normoxia contrast in the series metadata
4. processed matrix deposited directly on the GEO series record
5. not already present in the frozen microarray benchmark

Preferred candidates also have at least 20 total samples. If no such cohort is
available with a direct processed matrix and a clean binary contrast, the
fallback is a smaller exact perturbation cohort, which is then treated as a
supportive external layer rather than a headline pooled result.

## Selected Cohort

- `GSE179885`
- human T-cell bulk RNA-seq
- 6 hypoxia samples vs 6 normoxia samples
- processed TPM matrix available directly from the GEO series record

This cohort satisfies the exact-perturbation requirement and provides a clean
external transfer test for the hypoxia signature, while remaining visibly
separate from the main in-panel frozen benchmark.
