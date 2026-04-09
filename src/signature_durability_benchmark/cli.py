"""Command-line interface."""
from __future__ import annotations
import argparse
import sys

def main() -> None:
    parser = argparse.ArgumentParser(prog="signature-durability-benchmark")
    sub = parser.add_subparsers(dest="command")

    # build-freeze
    p = sub.add_parser("build-freeze", help="Validate and freeze benchmark assets")
    p.add_argument("--config", required=True)
    p.add_argument("--out", required=True)

    # run
    p = sub.add_parser("run", help="Run the canonical benchmark")
    p.add_argument("--config", required=True)
    p.add_argument("--out", required=True)

    # verify
    p = sub.add_parser("verify", help="Verify a benchmark run")
    p.add_argument("--config", required=True)
    p.add_argument("--run-dir", required=True)

    # triage
    p = sub.add_parser("triage", help="Score an arbitrary new signature")
    p.add_argument("--config", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)

    # build-paper
    p = sub.add_parser("build-paper", help="Generate paper from benchmark results")
    p.add_argument("--config", required=True)
    p.add_argument("--run-dir", required=True)
    p.add_argument("--out", required=True)

    # build-clawrxiv-payload
    p = sub.add_parser("build-clawrxiv-payload", help="Build clawRxiv submission payload")
    p.add_argument("--config", required=True)

    # build-freeze-data
    p = sub.add_parser("build-freeze-data", help="Download and freeze GEO cohort data")
    p.add_argument("--config", required=True)
    p.add_argument("--out", required=True)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(1)

    from .config import load_config
    config = load_config(args.config)

    if args.command == "build-freeze":
        from .freeze import build_freeze
        build_freeze(config, args.out)
    elif args.command == "run":
        from .benchmark import run_pipeline
        run_pipeline(config, args.out)
    elif args.command == "verify":
        from .verify import run_verification
        run_verification(config, args.run_dir)
    elif args.command == "triage":
        from .benchmark import run_pipeline
        # triage is a simplified version - to be implemented
        print("Triage not yet implemented")
    elif args.command == "build-paper":
        from .paper_build import build_paper
        build_paper(config, args.run_dir, args.out)
    elif args.command == "build-clawrxiv-payload":
        from .submission import build_clawrxiv_payload
        build_clawrxiv_payload(config)
    elif args.command == "build-freeze-data":
        print("Run: uv run python scripts/download_geo_cohorts.py")

    print(f"{args.command} complete.")

if __name__ == "__main__":
    main()
