#!/usr/bin/env python3
"""Run a literature-grounded rescue case study through triage."""

from __future__ import annotations

import json
from pathlib import Path

import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from signature_durability_benchmark.config import load_config
from signature_durability_benchmark.diagnostic import run_triage
from signature_durability_benchmark.utils import ensure_dir, read_json, write_json, write_text


CONFIG = ROOT / "config" / "benchmark_config.yaml"
SIGNATURE = ROOT / "data" / "case_studies" / "ayers_ifng6_signature.tsv"
OUT_DIR = ROOT / "outputs" / "canonical_v8"
CASE_DIR = OUT_DIR / "rescued_signature_case_study"
TRIAGE_DIR = CASE_DIR / "triage_ayers_ifng6"
OUT_JSON = OUT_DIR / "rescued_signature_case_study.json"
OUT_MD = OUT_DIR / "rescued_signature_case_study.md"


def render_markdown(payload: dict[str, object]) -> str:
    headline = payload["headline_metrics"]
    triage = payload["triage_summary"].strip()
    source = payload["source_signature"]
    return "\n".join(
        [
            "# Rescued Signature Case Study",
            "",
            "## Source Signature",
            "",
            f"- Signature: `{source['name']}`",
            f"- Source paper: {source['citation']}",
            f"- DOI: `{source['doi']}`",
            f"- Genes: `{', '.join(source['genes'])}`",
            "",
            "## Why This Case",
            "",
            payload["motivation"],
            "",
            "## Headline Metrics",
            "",
            f"- Inferred best-supported program: `{headline['inferred_program']}`",
            f"- Full-model classification: `{headline['full_model_class']}`",
            f"- Within-program classification: `{headline['within_program_class']}`",
            f"- Aggregate effect: `{headline['aggregate_effect']:+.3f}`",
            f"- Aggregate I²: `{headline['aggregate_i_squared']:.3f}`",
            f"- Interferon within effect: `{headline['within_interferon_effect']:+.3f}` (p=`{headline['within_interferon_p']:.4g}`)",
            f"- Outside-interferon effect: `{headline['outside_interferon_effect']:+.3f}` (p=`{headline['outside_interferon_p']:.4g}`)",
            f"- Program-structure permutation p: `{headline['permutation_p']:.4g}`",
            "",
            "## Interpretation",
            "",
            payload["interpretation"],
            "",
            "## Triage Output",
            "",
            triage,
            "",
        ]
    )


def main() -> None:
    ensure_dir(CASE_DIR)
    config = load_config(CONFIG)
    run_triage(config, SIGNATURE, TRIAGE_DIR)
    diagnostic = read_json(TRIAGE_DIR / "diagnostic.json")
    triage_summary = (TRIAGE_DIR / "diagnostic_summary.md").read_text(encoding="utf-8")

    payload = {
        "case_id": "ayers_ifng6_rescue",
        "source_signature": {
            "name": "Ayers IFN-gamma-related 6-gene profile",
            "citation": (
                "Ayers M, Lunceford J, Nebozhyn M, et al. "
                "IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade. "
                "J Clin Invest. 2017;127:2930-2940."
            ),
            "doi": "10.1172/JCI91190",
            "genes": ["IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG"],
        },
        "motivation": (
            "This published oncology response signature comes from outside the frozen 5-program panel. "
            "Taken naively across all 35 cohorts it looks mixed, but triage can still recover the underlying "
            "interferon home program and the within-context durable signal."
        ),
        "headline_metrics": {
            "inferred_program": diagnostic["inferred_program"]["program"],
            "full_model_class": diagnostic["classification"]["full_model"],
            "within_program_class": diagnostic["classification"]["within_program"],
            "aggregate_effect": diagnostic["aggregate_profile"]["aggregate_effect"],
            "aggregate_i_squared": diagnostic["aggregate_profile"]["i_squared"],
            "within_interferon_effect": diagnostic["inferred_program"]["within_effect"],
            "within_interferon_p": diagnostic["inferred_program"]["within_p"],
            "outside_interferon_effect": diagnostic["inferred_program"]["outside_effect"],
            "outside_interferon_p": diagnostic["inferred_program"]["outside_p"],
            "permutation_p": diagnostic["permutation_program_structure"]["p_value"],
        },
        "interpretation": (
            "This is the rescue pattern the method is meant to expose: a signature imported from a different "
            "application domain is not cleanly one-program in the aggregate (`mixed`, aggregate effect "
            f"{diagnostic['aggregate_profile']['aggregate_effect']:+.3f}, I²={diagnostic['aggregate_profile']['i_squared']:.3f}), "
            "but the best-supported program is interferon, the within-interferon effect is positive and "
            f"guarded-significant ({diagnostic['inferred_program']['within_effect']:+.3f}, p={diagnostic['inferred_program']['within_p']:.4g}), "
            "and the outside-interferon effect is essentially null. In other words, the cross-context dilution "
            "does not imply that the signature is broken; it indicates context mismatch."
        ),
        "triage_output_dir": str(TRIAGE_DIR),
        "triage_summary": triage_summary,
    }
    write_json(OUT_JSON, payload)
    write_text(OUT_MD, render_markdown(payload))
    print(json.dumps({"case_id": payload["case_id"], "inferred_program": payload["headline_metrics"]["inferred_program"]}, indent=2))


if __name__ == "__main__":
    main()
