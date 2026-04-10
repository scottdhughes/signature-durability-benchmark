"""Command-line interface."""
from __future__ import annotations
import argparse
import runpy
import sys
from pathlib import Path

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
    p.add_argument("--program", required=False)

    # build-paper
    p = sub.add_parser("build-paper", help="Generate paper from benchmark results")
    p.add_argument("--config", required=True)
    p.add_argument("--run-dir", required=True)
    p.add_argument("--out", required=True)

    # build-clawrxiv-payload
    p = sub.add_parser("build-clawrxiv-payload", help="Build clawRxiv submission payload")
    p.add_argument("--config", required=True)

    # prospective-holdout
    p = sub.add_parser("prospective-holdout", help="Run the metadata-first prospective holdout evaluation")
    p.add_argument("--config", required=True)

    # declare-prospective-round
    p = sub.add_parser(
        "declare-prospective-round",
        help="Create or verify an externally timestamped prospective-round declaration",
    )
    p.add_argument("--registry", required=True)
    p.add_argument("--protocol", required=True)
    p.add_argument("--round-id", required=False)
    p.add_argument("--force", action="store_true")

    # prospective-round-evaluate
    p = sub.add_parser(
        "prospective-round-evaluate",
        help="Evaluate a locked prospective round against an existing declaration receipt",
    )
    p.add_argument("--registry", required=True)
    p.add_argument("--protocol", required=True)
    p.add_argument("--receipt", required=True)
    p.add_argument("--out", required=True)

    # build-freeze-data
    p = sub.add_parser("build-freeze-data", help="Download and freeze GEO cohort data")
    p.add_argument("--config", required=True)
    p.add_argument("--out", required=True)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(1)

    config = None
    if hasattr(args, "config"):
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
        from .diagnostic import run_triage
        run_triage(config, args.input, args.out, declared_program=args.program)
    elif args.command == "build-paper":
        from .paper_build import build_paper
        build_paper(config, args.run_dir, args.out)
    elif args.command == "build-clawrxiv-payload":
        from .submission import build_clawrxiv_payload
        build_clawrxiv_payload(config)
    elif args.command == "prospective-holdout":
        script_path = Path(__file__).resolve().parents[2] / "scripts" / "prospective_holdout_prediction.py"
        runpy.run_path(str(script_path), run_name="__main__")
    elif args.command == "declare-prospective-round":
        script_path = Path(__file__).resolve().parents[2] / "scripts" / "declare_prospective_round.py"
        saved_argv = sys.argv[:]
        try:
            sys.argv = [
                str(script_path),
                "--registry",
                args.registry,
                "--protocol",
                args.protocol,
            ]
            if args.round_id:
                sys.argv.extend(["--round-id", args.round_id])
            if args.force:
                sys.argv.append("--force")
            runpy.run_path(str(script_path), run_name="__main__")
        finally:
            sys.argv = saved_argv
    elif args.command == "prospective-round-evaluate":
        script_path = Path(__file__).resolve().parents[2] / "scripts" / "prospective_round_evaluate.py"
        saved_argv = sys.argv[:]
        try:
            sys.argv = [
                str(script_path),
                "--registry",
                args.registry,
                "--protocol",
                args.protocol,
                "--receipt",
                args.receipt,
                "--out",
                args.out,
            ]
            runpy.run_path(str(script_path), run_name="__main__")
        finally:
            sys.argv = saved_argv
    elif args.command == "build-freeze-data":
        print("Run: uv run python scripts/download_geo_cohorts.py")

    print(f"{args.command} complete.")

if __name__ == "__main__":
    main()
