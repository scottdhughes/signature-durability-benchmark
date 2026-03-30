"""Paper generation from benchmark results."""
from __future__ import annotations
from pathlib import Path
from typing import Any
from .config import SkillConfig
from .utils import read_json, ensure_dir, write_text

def build_paper(config: SkillConfig, run_dir: str | Path, out_dir: str | Path) -> None:
    run_path = Path(run_dir)
    out_path = ensure_dir(out_dir)
    manifest = read_json(run_path / "manifest.json")
    rule = manifest.get("aggregate_metrics", {}).get("success_rule", {})

    results_tex = f"""\\section{{Frozen Results}}
The final frozen run compared the full model against the overlap-only baseline.
The full model achieved AUPRC {rule.get('full_model_auprc', 0):.4f} versus
overlap-only {rule.get('overlap_only_auprc', 0):.4f} (margin {rule.get('auprc_margin', 0):.4f}).
Secondary wins: {rule.get('secondary_wins', 0)}.
Success rule: {'passed' if rule.get('success') else 'failed'}.
"""
    write_text(out_path / "generated_results.tex", results_tex)
    print(f"Paper artifacts written to {out_path}")
