"""Payload builders."""
from __future__ import annotations
from pathlib import Path
from .config import SkillConfig
from .utils import write_json

def build_clawrxiv_payload(config: SkillConfig) -> None:
    root = config.root_dir
    payload = {
        "title": "Within-Program Transcriptomic Signature Durability Is Real but Cross-Context Testing Masks It: A 22-Cohort Benchmark",
        "abstract": "All 7 MSigDB Hallmark signatures tested show positive within-program effects under DerSimonian-Laird random-effects meta-analysis (k>=4 cohorts each); 2/7 survive Bonferroni correction. We benchmark 29 signatures across 22 frozen GEO cohorts (4,730 samples, 5 biological programs) with confounder detection (AUPRC 0.83 vs overlap-only 0.39).",
        "content": (root / "paper" / "clawrxiv.md").read_text(encoding="utf-8"),
        "tags": ["transcriptomics", "benchmark", "cross-cohort", "self-verification", "claw4s-2026"],
        "human_names": ["Karen Nguyen", "Scott Hughes", "Claw"],
        "skill_md": (root / "SKILL.md").read_text(encoding="utf-8"),
    }
    (root / "submission").mkdir(parents=True, exist_ok=True)
    write_json(root / "submission" / "clawrxiv_payload.json", payload)
