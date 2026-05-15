"""
GeneScreen 1.0 - 配置管理

用于存储用户可配置路径（如数据库目录）
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Any


CONFIG_DIR = Path.home() / ".genescreen"
CONFIG_PATH = CONFIG_DIR / "config.json"
DEFAULT_DB_DIR = CONFIG_DIR / "genomes"
DEFAULT_OUTPUT_DIR_NAME = "GeneScreenOutput"


def _load_config() -> Dict[str, Any]:
    if not CONFIG_PATH.exists():
        return {}
    try:
        with CONFIG_PATH.open("r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}


def _save_config(data: Dict[str, Any]) -> None:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    with CONFIG_PATH.open("w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=True, indent=2)


def get_db_dir() -> Path:
    data = _load_config()
    db_dir = data.get("db_dir")
    if db_dir:
        return Path(db_dir)
    return DEFAULT_DB_DIR


def set_db_dir(path: str | Path) -> None:
    data = _load_config()
    data["db_dir"] = str(Path(path))
    _save_config(data)


def get_db_path() -> Path:
    return get_db_dir() / "genescreen.db"


def get_output_dir() -> Path:
    data = _load_config()
    output_dir = data.get("output_dir")
    if output_dir:
        return Path(output_dir)
    db_dir = get_db_dir()
    return db_dir.parent / DEFAULT_OUTPUT_DIR_NAME


def set_output_dir(path: str | Path) -> None:
    data = _load_config()
    data["output_dir"] = str(Path(path))
    _save_config(data)
