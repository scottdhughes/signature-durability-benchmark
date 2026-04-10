#!/usr/bin/env python3
"""External hypoxia transfer layer for the second positive generalization case."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from signature_durability_benchmark.prospective import DISPLAY_NAMES, GENE_INFO, prepare_cohort, load_gene_maps
from signature_durability_benchmark.scoring import score_signature_in_cohort
from signature_durability_benchmark.utils import read_table, repo_relpath, write_json, write_text


ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "data" / "external_hypoxia" / "cohort_manifest.tsv"
SEARCH_PROTOCOL = ROOT / "data" / "external_hypoxia" / "SEARCH_PROTOCOL.md"
SIGNATURES = ROOT / "data" / "freeze" / "signatures.tsv"
DOWNLOAD_DIR = ROOT / "data" / "external_hypoxia" / "downloads"
OUT_DIR = ROOT / "outputs" / "canonical_v8"
OUT_JSON = OUT_DIR / "external_hypoxia_validation.json"
OUT_MD = OUT_DIR / "external_hypoxia_validation.md"

SIGNATURE_IDS = [
    "hallmark_hypoxia",
    "hallmark_emt",
    "hallmark_inflammatory_response",
    "hallmark_tnfa_nfkb",
    "hallmark_e2f_targets",
    "hallmark_ifng_response",
    "hallmark_ifna_response",
]


def main() -> None:
    manifest = read_table(MANIFEST)
    signatures = read_table(SIGNATURES)
    gene_maps = load_gene_maps(GENE_INFO)
    download_audit: list[dict[str, object]] = []

    if manifest.shape[0] != 1:
        raise ValueError("external_hypoxia cohort_manifest.tsv must contain exactly one cohort row")
    spec = manifest.iloc[0]
    cohort = prepare_cohort(spec, gene_maps, download_audit, download_dir=DOWNLOAD_DIR)

    results = []
    for signature_id in SIGNATURE_IDS:
        signature = signatures.loc[
            signatures["signature_id"] == signature_id,
            ["gene_symbol", "direction", "weight"],
        ].copy()
        scored = score_signature_in_cohort(
            signature,
            cohort["expr"],
            cohort["pheno"],
            phenotype_column="phenotype",
            case_label="case",
            control_label="control",
        )
        results.append(
            {
                "signature_id": signature_id,
                "signature_label": DISPLAY_NAMES.get(signature_id, signature_id),
                **scored,
            }
        )

    results = sorted(results, key=lambda row: (-float(row["cohens_d"]), row["signature_id"]))
    chosen = next(row for row in results if row["signature_id"] == "hallmark_hypoxia")
    summary = {
        "description": (
            "External hypoxia exact-perturbation layer selected under the deterministic "
            "external_hypoxia SEARCH_PROTOCOL. This cohort is supportive rather than "
            "headline because it contains 12 total samples."
        ),
        "search_protocol_path": repo_relpath(SEARCH_PROTOCOL),
        "download_audit": download_audit,
        "cohort": {
            "external_cohort_id": str(spec["external_cohort_id"]),
            "geo_accession": str(spec["geo_accession"]),
            "contrast_description": str(spec["contrast_description"]),
            "tissue": str(spec["tissue"]),
            "expected_case_n": int(spec["expected_case_n"]),
            "expected_control_n": int(spec["expected_control_n"]),
        },
        "signature_results": results,
        "chosen_signature_support": {
            "signature_id": "hallmark_hypoxia",
            "effect_size_g": float(chosen["cohens_d"]),
            "coverage_fraction": float(chosen["coverage_fraction"]),
            "positive_direction": bool(float(chosen["cohens_d"]) > 0),
            "rank_among_scored_signatures": next(i for i, row in enumerate(results, start=1) if row["signature_id"] == "hallmark_hypoxia"),
        },
    }
    write_json(OUT_JSON, summary)

    lines = [
        "# External Hypoxia Validation",
        "",
        f"- Cohort: `{spec['external_cohort_id']}` ({spec['geo_accession']})",
        f"- Tissue: `{spec['tissue']}`",
        f"- Contrast: {spec['contrast_description']}",
        f"- Search protocol: `{SEARCH_PROTOCOL.relative_to(ROOT)}`",
        "",
        "## Ranked Signature Effects",
        "",
        "| Rank | Signature | Hedges' g | Coverage |",
        "|---:|---|---:|---:|",
    ]
    for rank, row in enumerate(results, start=1):
        lines.append(
            f"| {rank} | {row['signature_label']} | {float(row['cohens_d']):+.3f} | {float(row['coverage_fraction']):.3f} |"
        )
    lines.extend(
        [
            "",
            "This is an exact-perturbation external layer rather than a main-table pooled result "
            "because the cohort has 12 total samples, below the preferred >=20-sample rule.",
            "",
        ]
    )
    write_text(OUT_MD, "\n".join(lines))
    print(json.dumps(summary["chosen_signature_support"], indent=2))


if __name__ == "__main__":
    main()
