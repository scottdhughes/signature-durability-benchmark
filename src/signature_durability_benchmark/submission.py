"""Payload builders."""
from __future__ import annotations
from pathlib import Path
from .config import SkillConfig
from .utils import write_json

def build_clawrxiv_payload(config: SkillConfig) -> None:
    root = config.root_dir
    payload = {
        "title": "From Published Signatures to Durable Signals: A Self-Verifying Cross-Cohort Benchmark for Transcriptomic Signature Generalization",
        "abstract": "Published transcriptomic signatures often look convincing in one study but fail across cohorts, platforms, or nuisance biology. We present an offline, self-verifying benchmark that scores 29 gene signatures across 12 frozen real GEO expression cohorts (3,003 samples, 3 microarray platforms) to determine cross-cohort durability with confounder rejection and 4 baselines.",
        "content": (root / "paper" / "clawrxiv.md").read_text(encoding="utf-8"),
        "tags": ["transcriptomics", "benchmark", "cross-cohort", "self-verification", "claw4s-2026"],
        "human_names": ["Karen Nguyen", "Scott Hughes", "Claw"],
        "skill_md": (root / "SKILL.md").read_text(encoding="utf-8"),
    }
    (root / "submission").mkdir(parents=True, exist_ok=True)
    write_json(root / "submission" / "clawrxiv_payload.json", payload)
