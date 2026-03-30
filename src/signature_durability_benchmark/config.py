"""Configuration loading."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml

@dataclass(frozen=True)
class SkillConfig:
    root_dir: Path
    raw: dict[str, Any]

    @property
    def title(self) -> str:
        return str(self.raw["title"])

    @property
    def abstract(self) -> str:
        return str(self.raw["abstract"])

    @property
    def contract_version(self) -> str:
        return str(self.raw["contract_version"])

    @property
    def runtime(self) -> dict[str, Any]:
        return dict(self.raw["runtime"])

    @property
    def paths(self) -> dict[str, str]:
        return dict(self.raw["paths"])

    @property
    def outputs(self) -> dict[str, str]:
        return dict(self.raw["outputs"])

    @property
    def sources(self) -> dict[str, dict[str, Any]]:
        return dict(self.raw["sources"])

    @property
    def freeze(self) -> dict[str, Any]:
        return dict(self.raw["freeze"])

    @property
    def scoring(self) -> dict[str, Any]:
        return dict(self.raw["scoring"])

    @property
    def meta_analysis(self) -> dict[str, Any]:
        return dict(self.raw["meta_analysis"])

    def path(self, key: str) -> Path:
        return self.root_dir / self.paths[key]

def load_config(path: str | Path) -> SkillConfig:
    config_path = Path(path).resolve()
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    return SkillConfig(root_dir=config_path.parent.parent, raw=raw)
