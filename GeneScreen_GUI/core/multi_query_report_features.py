#!/usr/bin/env python3
"""Feature metadata helpers for multi-query report assets."""

import os
from typing import Any, Dict, List, Optional


def parse_gff_attributes(attributes: str) -> Dict[str, str]:
    parsed: Dict[str, str] = {}
    field = ""
    quoted = False
    for char in f"{attributes};":
        if char == '"':
            quoted = not quoted
            field += char
            continue
        if char == ";" and not quoted:
            if field:
                if "=" in field:
                    key, value = field.split("=", 1)
                elif " " in field:
                    key, value = field.split(" ", 1)
                else:
                    key, value = field, ""
                parsed[key.strip()] = value.strip().strip('"')
            field = ""
        else:
            field += char
    return parsed


def _clean_gene_label(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    for prefix in ("gene:", "transcript:", "mrna:", "mRNA:"):
        if text.startswith(prefix):
            return text[len(prefix):]
    return text


def _gene_label(attributes: Dict[str, str], fallback: str) -> str:
    for key in ("gene_id", "Name", "name", "locus_tag", "ID"):
        label = _clean_gene_label(attributes.get(key))
        if label:
            return label
    return _clean_gene_label(fallback) or "gene"


def read_gene_features(
    source_gff: Optional[str],
    source_chr: Optional[Any] = None,
    region_start: Optional[int] = None,
    region_end: Optional[int] = None,
    relative: bool = False,
) -> List[Dict[str, Any]]:
    """Read gene features overlapping an optional source region.

    Returned coordinates use the same coordinate system as the report track:
    absolute coordinates for query genome tracks, and relative coordinates when
    the source GFF is converted into a ref-derived sequence.
    """
    if not source_gff or not os.path.exists(source_gff):
        return []

    def chrom_key(value: Any) -> str:
        text = str(value or "").split()[0]
        lowered = text.lower()
        if lowered.startswith("chr"):
            lowered = lowered[3:]
        return lowered.lstrip("0") or lowered

    def chrom_matches(left: Any, right: Any) -> bool:
        return str(left or "") == str(right or "") or chrom_key(left) == chrom_key(right)

    features: List[Dict[str, Any]] = []
    with open(source_gff, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2].lower() != "gene":
                continue
            if source_chr is not None and not chrom_matches(parts[0], source_chr):
                continue
            try:
                feat_start = int(parts[3])
                feat_end = int(parts[4])
            except ValueError:
                continue
            if feat_end < feat_start:
                feat_start, feat_end = feat_end, feat_start
            if region_start is not None and region_end is not None:
                left, right = sorted((int(region_start), int(region_end)))
                if feat_end < left or feat_start > right:
                    continue
                if relative:
                    display_start = max(1, feat_start - left + 1)
                    display_end = max(1, feat_end - left + 1)
                else:
                    display_start, display_end = feat_start, feat_end
            else:
                display_start, display_end = feat_start, feat_end
            attributes = parse_gff_attributes(parts[8])
            features.append(
                {
                    "id": _clean_gene_label(attributes.get("ID")) or _gene_label(attributes, parts[0]),
                    "label": _gene_label(attributes, parts[0]),
                    "seqid": parts[0],
                    "start": display_start,
                    "end": display_end,
                    "source_start": feat_start,
                    "source_end": feat_end,
                    "strand": parts[6],
                }
            )
    features.sort(key=lambda item: (int(item.get("start") or 0), int(item.get("end") or 0), str(item.get("label") or "")))
    return features
