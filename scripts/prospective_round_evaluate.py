#!/usr/bin/env python3
"""Evaluate a locked prospective round against an existing declaration receipt."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from signature_durability_benchmark.prospective import evaluate_round


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    summary = evaluate_round(
        Path(args.registry).resolve(),
        Path(args.protocol).resolve(),
        receipt_path=Path(args.receipt).resolve(),
        out_dir=Path(args.out).resolve(),
    )
    print(json.dumps({"round_id": summary["round_id"], "overall_success": summary["prediction_summary"]["overall_success"], "out": args.out}, indent=2))


if __name__ == "__main__":
    main()
