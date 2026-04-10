"""General helpers."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def ensure_dir(path: str | Path) -> Path:
    target = Path(path)
    target.mkdir(parents=True, exist_ok=True)
    return target


def read_json(path: str | Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path: str | Path, payload: Any) -> None:
    Path(path).write_text(json.dumps(json_ready(payload), indent=2, sort_keys=True), encoding="utf-8")


def write_text(path: str | Path, text: str) -> None:
    Path(path).write_text(text, encoding="utf-8")


def read_table(path: str | Path) -> pd.DataFrame:
    file_path = Path(path)
    separator = "\t" if file_path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(file_path, sep=separator).fillna("")


def write_table(frame: pd.DataFrame, path: str | Path) -> None:
    file_path = Path(path)
    separator = "\t" if file_path.suffix.lower() in {".tsv", ".txt"} else ","
    frame.to_csv(file_path, sep=separator, index=False, quoting=csv.QUOTE_MINIMAL)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repo_relpath(path: str | Path) -> str:
    candidate = Path(path)
    if not candidate.is_absolute():
        return str(candidate)
    try:
        return str(candidate.relative_to(ROOT))
    except ValueError:
        return str(candidate)


def set_runtime_environment(seed: int) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)


def runtime_summary() -> dict[str, str]:
    return {
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "pythonhashseed": os.environ.get("PYTHONHASHSEED", ""),
    }


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_ready(item) for item in value]
    if isinstance(value, tuple):
        return [json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if math.isnan(float(value)) else float(value)
    if isinstance(value, float):
        return None if math.isnan(value) else value
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    return value


def now_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def stable_name_seed(name: str) -> int:
    return sum((index + 1) * ord(char) for index, char in enumerate(name))


def ensure_columns(frame: pd.DataFrame, required: list[str] | tuple[str, ...], label: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} missing required columns: {missing}")


def markdown_has_heading(text: str, heading: str) -> bool:
    return f"## {heading}" in text or f"# {heading}" in text


def parse_float(value: Any, default: float = 0.0) -> float:
    if value in ("", None):
        return default
    return float(value)


def format_rate(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "0/0"
    return f"{numerator}/{denominator} ({numerator / denominator:.3f})"
