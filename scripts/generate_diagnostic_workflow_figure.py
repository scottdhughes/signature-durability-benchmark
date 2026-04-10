#!/usr/bin/env python3
"""Generate the high-level diagnostic workflow figure."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from signature_durability_benchmark.workflow_figure import build_workflow_figure_data, render_workflow_figure


PDF = ROOT / "paper" / "figure_workflow.pdf"
PNG = ROOT / "paper" / "figure_workflow.png"


def main() -> None:
    data = build_workflow_figure_data(ROOT)
    render_workflow_figure(data, PDF, PNG)
    print(json.dumps({"pdf": str(PDF), "png": str(PNG)}, indent=2))


if __name__ == "__main__":
    main()
