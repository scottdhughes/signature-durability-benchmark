"""Payload builders."""
from __future__ import annotations
from pathlib import Path
from .config import SkillConfig
from .utils import write_json

def build_clawrxiv_payload(config: SkillConfig) -> None:
    root = config.root_dir
    payload = {
        "title": "Program-Conditioned Diagnostic for Transcriptomic Signature Durability: Validation on Interferon Signatures across 35 Frozen GEO Cohorts",
        "abstract": "We present a program-conditioned diagnostic for transcriptomic signatures that scores a signature against a frozen cohort panel, compares within-program versus outside-program effects, tests program structure by permutation, and surfaces failure modes when labels are too coarse. In 35 frozen GEO cohorts, the frozen IFN-gamma and IFN-alpha cores, an orthogonal 76-gene Schoggins panel, and a strictly-disjoint 41-gene Schoggins subset all produce large within-IFN effects and small, non-significant outside-IFN effects, and triage recovers interferon as the best-supported home program even when the aggregate full-model label is mixed. The same rescue logic extends to a published external signature: the Ayers IFN-gamma-related 6-gene clinical response profile triages as aggregate mixed but within-program durable for interferon. Held-out validation shows 100% leave-one-cohort-out sign prediction and 100% split-half sign agreement across the IFN quartet, and a four-cohort bulk RNA-seq extension reproduces the IFN quartet with guarded Bonferroni-significant pooled effects while inflammatory, TNF-alpha/NF-kB, and E2F comparators remain non-significant. A deterministic non-IFN breadth rule then selects Hallmark Hypoxia over Hallmark EMT, but hypoxia honestly remains mixed despite supportive external validation. The prospective layer is likewise non-tautological: two predeclared v1 cohorts satisfy the 4/4 IFN sign forecast, whereas one severity-mixed acute PBMC cohort inverts all four IFN signatures. The repo now also carries a second fresh prospective registry bound to an external RFC3161 timestamp, intentionally left unevaluated to establish an immutable third-party-auditable future challenge, along with a generic locked-round scoring command and release-ready declaration bundle.",
        "content": (root / "paper" / "clawrxiv.md").read_text(encoding="utf-8"),
        "tags": [
            "transcriptomics",
            "benchmark",
            "diagnostic",
            "cross-cohort",
            "prospective-validation",
            "claw4s-2026",
        ],
        "human_names": ["Karen Nguyen", "Scott Hughes", "Claw"],
        "skill_md": (root / "SKILL.md").read_text(encoding="utf-8"),
    }
    (root / "submission").mkdir(parents=True, exist_ok=True)
    write_json(root / "submission" / "clawrxiv_payload.json", payload)
