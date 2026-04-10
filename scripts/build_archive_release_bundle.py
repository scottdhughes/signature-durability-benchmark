#!/usr/bin/env python3
"""Build a release-ready declaration bundle for the v2 prospective archive."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
BUNDLE_ROOT = ROOT / "submission" / "archive_bundles"
DECLARATION_DIR = BUNDLE_ROOT / "prospective_holdout_v2_declaration"
REGISTRY = ROOT / "data" / "prospective_holdout" / "prediction_registry_v2.tsv"
PROTOCOL = ROOT / "data" / "prospective_holdout" / "PREDICTION_PROTOCOL_v2.md"
TIMESTAMP_DIR = ROOT / "data" / "prospective_holdout" / "external_timestamps" / "prospective_holdout_v2"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def copy_tree(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def main() -> None:
    if DECLARATION_DIR.exists():
        shutil.rmtree(DECLARATION_DIR)
    DECLARATION_DIR.mkdir(parents=True, exist_ok=True)

    shutil.copy2(REGISTRY, DECLARATION_DIR / REGISTRY.name)
    shutil.copy2(PROTOCOL, DECLARATION_DIR / PROTOCOL.name)
    copy_tree(TIMESTAMP_DIR, DECLARATION_DIR / TIMESTAMP_DIR.name)

    release_note = "\n".join(
        [
            "# Release Notes: prospective_holdout_v2 declaration",
            "",
            "This bundle is the declaration-only archive for the second prospective round.",
            "It contains the frozen registry, frozen protocol, and the externally timestamped RFC3161 receipt bundle.",
            "",
            "Status: declared and unevaluated.",
            "",
            "Intended publication targets:",
            "- GitHub Release asset bundle",
            "- Zenodo versioned archive",
            "",
            "The scored evaluation bundle must be published separately after the locked v2 readout exists.",
            "",
        ]
    )
    (DECLARATION_DIR / "RELEASE_NOTES.md").write_text(release_note, encoding="utf-8")

    checksum_lines = []
    for path in sorted(DECLARATION_DIR.rglob("*")):
        if not path.is_file() or path.name == "CHECKSUMS.sha256":
            continue
        rel = path.relative_to(DECLARATION_DIR)
        checksum_lines.append(f"{sha256_file(path)}  {rel.as_posix()}")
    (DECLARATION_DIR / "CHECKSUMS.sha256").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")

    manifest = {
        "bundle": "prospective_holdout_v2_declaration",
        "status": "declared_pending_evaluation",
        "files": sorted(
            str(path.relative_to(DECLARATION_DIR))
            for path in DECLARATION_DIR.rglob("*")
            if path.is_file() and path.name != "bundle_manifest.json"
        ),
    }
    (DECLARATION_DIR / "bundle_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"bundle_dir": str(DECLARATION_DIR)}, indent=2))


if __name__ == "__main__":
    main()
