#!/usr/bin/env python3
"""Create or verify an externally timestamped prospective-round declaration."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any
from urllib.request import Request, urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from signature_durability_benchmark.utils import ensure_dir, read_table, sha256_file, write_json, write_text


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TSA_URL = "https://freetsa.org/tsr"
DEFAULT_TSA_CA_URL = "https://freetsa.org/files/cacert.pem"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def fetch_bytes(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": "signature-durability-benchmark/1.0"})
    with urlopen(request, timeout=120) as response:
        return response.read()


def run_command(args: list[str]) -> str:
    completed = subprocess.run(args, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        message = (completed.stdout + completed.stderr).strip()
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(args)}\n{message}")
    return completed.stdout + completed.stderr


def build_manifest(registry_path: Path, protocol_path: Path, registry: Any) -> dict[str, Any]:
    round_id = str(registry["round_id"].iloc[0])
    declared_at_utc = str(registry["declared_at_utc"].iloc[0])
    preview = []
    for _, row in registry.iterrows():
        preview.append(
            {
                "external_cohort_id": str(row["external_cohort_id"]),
                "geo_accession": str(row["geo_accession"]),
                "tissue": str(row["tissue"]),
                "contrast_description": str(row["contrast_description"]),
                "expected_case_n": int(row["expected_case_n"]),
                "expected_control_n": int(row["expected_control_n"]),
                "acquisition_mode": str(row["acquisition_mode"]),
                "predicted_ifn_direction": str(row["predicted_ifn_direction"]),
            }
        )

    return {
        "round_id": round_id,
        "declared_at_utc": declared_at_utc,
        "status": "declared_pending_evaluation",
        "declaration_scope": (
            "Fresh cohorts are declared from metadata and supplementary-file listings only. "
            "No held-out effect-size scoring, meta-analysis, or signature-level evaluation is included "
            "in this declaration bundle."
        ),
        "registry": {
            "path": str(registry_path.relative_to(ROOT)),
            "sha256": sha256_file(registry_path),
            "n_rows": int(registry.shape[0]),
            "n_primary_rows": int((registry["prediction_group"].astype(str) == "primary").sum()),
        },
        "protocol": {
            "path": str(protocol_path.relative_to(ROOT)),
            "sha256": sha256_file(protocol_path),
        },
        "cohorts": preview,
    }


def parse_timestamp(reply_text: str) -> str:
    match = re.search(r"^Time stamp:\s+(.+)$", reply_text, flags=re.MULTILINE)
    if not match:
        return ""
    return match.group(1).strip()


def render_summary(receipt: dict[str, Any]) -> str:
    lines = [
        f"# External Timestamp Receipt: {receipt['round_id']}",
        "",
        f"- Status: `{receipt['status']}`",
        f"- Declared at (registry): `{receipt['declared_at_utc']}`",
        f"- RFC3161 TSA: [{receipt['tsa_url']}]({receipt['tsa_url']})",
        f"- TSA reply time: `{receipt['tsa_reply_time']}`",
        f"- Manifest SHA256: `{receipt['manifest_sha256']}`",
        f"- Registry SHA256: `{receipt['registry_sha256']}`",
        f"- Protocol SHA256: `{receipt['protocol_sha256']}`",
        f"- Verification: `{receipt['verification_status']}`",
        "",
        "## Declared Cohorts",
        "",
        "| Cohort | GEO | Tissue | Case n | Control n |",
        "|---|---|---|---:|---:|",
    ]
    for cohort in receipt["cohorts"]:
        lines.append(
            f"| {cohort['external_cohort_id']} | {cohort['geo_accession']} | {cohort['tissue']} | "
            f"{cohort['expected_case_n']} | {cohort['expected_control_n']} |"
        )
    lines.extend(
        [
            "",
            "This bundle is intentionally separate from the mixed evaluated v1 round. "
            "No v2 scoring outputs were generated when this receipt was minted.",
            "",
        ]
    )
    return "\n".join(lines)


def declare_round(
    registry_path: Path,
    protocol_path: Path,
    expected_round_id: str | None,
    force: bool,
    tsa_url: str,
    tsa_ca_url: str,
) -> dict[str, Any]:
    registry = read_table(registry_path)
    if "round_id" not in registry.columns or "declared_at_utc" not in registry.columns:
        raise ValueError("Registry must contain 'round_id' and 'declared_at_utc' columns")

    round_ids = sorted({str(value) for value in registry["round_id"].astype(str).tolist() if str(value)})
    if len(round_ids) != 1:
        raise ValueError(f"Registry must contain exactly one round_id, found {round_ids}")
    round_id = round_ids[0]
    if expected_round_id and expected_round_id != round_id:
        raise ValueError(f"Expected round_id {expected_round_id}, found {round_id}")

    declared_times = sorted({str(value) for value in registry["declared_at_utc"].astype(str).tolist() if str(value)})
    if len(declared_times) != 1:
        raise ValueError(f"Registry must contain exactly one declared_at_utc value, found {declared_times}")

    artifact_dir = ensure_dir(ROOT / "data" / "prospective_holdout" / "external_timestamps" / round_id)
    manifest_path = artifact_dir / "declaration_manifest.json"
    manifest_hash_path = artifact_dir / "declaration_manifest.sha256"
    query_path = artifact_dir / "declaration_request.tsq"
    response_path = artifact_dir / "declaration_response.tsr"
    ca_path = artifact_dir / "freetsa_cacert.pem"
    verify_path = artifact_dir / "verification.txt"
    reply_text_path = artifact_dir / "tsa_reply.txt"
    summary_path = artifact_dir / "README.md"
    receipt_path = artifact_dir / "declaration_receipt.json"

    manifest_payload = build_manifest(registry_path, protocol_path, registry)
    manifest_text = json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n"
    manifest_changed = manifest_path.exists() and manifest_path.read_text(encoding="utf-8") != manifest_text
    if manifest_changed and not force:
        raise ValueError(
            f"Existing manifest for {round_id} does not match the current registry/protocol. "
            "Re-run with --force to regenerate the declaration bundle."
        )

    write_text(manifest_path, manifest_text)
    manifest_sha256 = sha256_bytes(manifest_text.encode("utf-8"))
    write_text(manifest_hash_path, f"{manifest_sha256}  {manifest_path.name}\n")

    if force or not ca_path.exists():
        ca_path.write_bytes(fetch_bytes(tsa_ca_url))

    run_command(
        [
            "openssl",
            "ts",
            "-query",
            "-data",
            str(manifest_path),
            "-sha256",
            "-cert",
            "-no_nonce",
            "-out",
            str(query_path),
        ]
    )

    status = "verified_existing_receipt"
    if force or not response_path.exists() or manifest_changed:
        request = Request(
            tsa_url,
            data=query_path.read_bytes(),
            headers={
                "Content-Type": "application/timestamp-query",
                "User-Agent": "signature-durability-benchmark/1.0",
            },
        )
        with urlopen(request, timeout=120) as response:
            response_path.write_bytes(response.read())
        status = "created_new_receipt"

    verification_text = run_command(
        [
            "openssl",
            "ts",
            "-verify",
            "-in",
            str(response_path),
            "-queryfile",
            str(query_path),
            "-CAfile",
            str(ca_path),
        ]
    )
    write_text(verify_path, verification_text)

    reply_text = run_command(["openssl", "ts", "-reply", "-in", str(response_path), "-text"])
    write_text(reply_text_path, reply_text)

    receipt = {
        "round_id": round_id,
        "status": status,
        "declared_at_utc": declared_times[0],
        "registry_path": str(registry_path.relative_to(ROOT)),
        "protocol_path": str(protocol_path.relative_to(ROOT)),
        "registry_sha256": manifest_payload["registry"]["sha256"],
        "protocol_sha256": manifest_payload["protocol"]["sha256"],
        "manifest_path": str(manifest_path.relative_to(ROOT)),
        "manifest_sha256": manifest_sha256,
        "query_path": str(query_path.relative_to(ROOT)),
        "query_sha256": sha256_file(query_path),
        "response_path": str(response_path.relative_to(ROOT)),
        "response_sha256": sha256_file(response_path),
        "tsa_url": tsa_url,
        "tsa_ca_url": tsa_ca_url,
        "tsa_reply_time": parse_timestamp(reply_text),
        "verification_status": "OK" if "Verification: OK" in verification_text else "UNKNOWN",
        "cohorts": manifest_payload["cohorts"],
    }
    write_json(receipt_path, receipt)
    write_text(summary_path, render_summary(receipt))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--round-id")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--tsa-url", default=DEFAULT_TSA_URL)
    parser.add_argument("--tsa-ca-url", default=DEFAULT_TSA_CA_URL)
    args = parser.parse_args()

    receipt = declare_round(
        registry_path=Path(args.registry).resolve(),
        protocol_path=Path(args.protocol).resolve(),
        expected_round_id=args.round_id,
        force=args.force,
        tsa_url=args.tsa_url,
        tsa_ca_url=args.tsa_ca_url,
    )
    print(
        json.dumps(
            {
                "round_id": receipt["round_id"],
                "status": receipt["status"],
                "receipt": str(
                    Path("data")
                    / "prospective_holdout"
                    / "external_timestamps"
                    / receipt["round_id"]
                    / "declaration_receipt.json"
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
