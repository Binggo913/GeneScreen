"""
GeneScreen 1.0 - Gene ID utilities

Extract gene IDs from annotation files or read prebuilt ID lists.
"""
from __future__ import annotations

import gzip
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, List

from utils import run_subprocess


# 文件后缀常量
JSON_LIST_SUFFIX = ".gene_ids.json"   # 旧 JSON 格式（兼容）
TXT_LIST_SUFFIX = ".gene_ids.txt"     # 新 TXT 格式（当前使用）
# 兼容旧代码
GENE_ID_LIST_SUFFIX = JSON_LIST_SUFFIX
LEGACY_LIST_SUFFIX = TXT_LIST_SUFFIX


def read_gene_id_list(path: str | Path) -> List[str]:
    path = str(path)
    lower = path.lower()
    if lower.endswith(".json"):
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, list):
            return [str(item) for item in data if str(item).strip()]
        if isinstance(data, dict) and isinstance(data.get("ids"), list):
            return [str(item) for item in data["ids"] if str(item).strip()]
        return []
    ids = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            value = line.strip()
            if value and not value.startswith("#"):
                ids.append(value)
    return ids


def write_gene_id_list(ids: Iterable[str], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".json":
        with path.open("w", encoding="utf-8") as f:
            json.dump(list(ids), f, ensure_ascii=True, indent=2)
        return
    with path.open("w", encoding="utf-8") as f:
        for gid in ids:
            f.write(f"{gid}\n")


def count_gene_id_list(path: str | Path) -> int:
    path = str(path)
    if path.lower().endswith(".json"):
        return len(read_gene_id_list(path))
    count = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                count += 1
    return count


def load_gene_ids(path: str | Path) -> List[str]:
    path = str(path)
    lower = path.lower()
    if lower.endswith(".json") or lower.endswith(GENE_ID_LIST_SUFFIX) or lower.endswith(LEGACY_LIST_SUFFIX) or lower.endswith(".ids") or lower.endswith(".list"):
        return read_gene_id_list(path)
    return extract_gene_ids_from_annotation(path)


def extract_gene_ids_from_annotation(path: str | Path) -> List[str]:
    gene_set = set()
    fallback_set = set()
    path = str(path)
    if not path.lower().endswith(".gz"):
        gene_lines = _rg_gene_lines(path)
        if gene_lines:
            for line in gene_lines:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                attrs = parts[8]
                gid = _extract_attr(attrs, ("gene_id", "ID", "Name", "locus_tag"))
                if gid:
                    gene_set.add(gid)
            if gene_set:
                return sorted(gene_set)
    open_func = gzip.open if path.lower().endswith(".gz") else open
    with open_func(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2].lower()
            attrs = parts[8]
            if feature in {"gene", "pseudogene"}:
                gid = _extract_attr(attrs, ("gene_id", "ID", "Name", "locus_tag"))
                if gid:
                    gene_set.add(gid)
            else:
                gid = _extract_attr(attrs, ("gene_id", "ID", "Name", "locus_tag"))
                if gid:
                    fallback_set.add(gid)
    ids = gene_set or fallback_set
    return sorted(ids)


def _rg_gene_lines(path: str) -> List[str]:
    if shutil.which("rg") is None:
        return []
    try:
        result = run_subprocess(
            [
                "rg",
                "--no-filename",
                "--no-line-number",
                "-e",
                r"\tgene\t",
                "-e",
                r"\tpseudogene\t",
                path,
            ],
            timeout=120
        )
        if result.returncode in (0, 1) and result.stdout:
            return result.stdout.splitlines()
    except Exception:
        return []
    return []

def _extract_attr(attrs: str, keys: tuple) -> str:
    for key in keys:
        match = re.search(rf"(?:^|;)\s*{re.escape(key)}=([^;]+)", attrs)
        if match:
            return match.group(1).strip().strip('"')
    for key in keys:
        match = re.search(rf"\b{re.escape(key)}\s+\"([^\"]+)\"", attrs)
        if match:
            return match.group(1).strip()
    for key in keys:
        match = re.search(rf"\b{re.escape(key)}\s+([^;]+)", attrs)
        if match:
            return match.group(1).strip().strip('"')
    return ""
