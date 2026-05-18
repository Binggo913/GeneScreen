#!/usr/bin/env python3
"""
Shared report data helpers for the 1 ref + N query result model.

Phase 1 uses this model for the current two-genome workflow (N=1) while
keeping the legacy HTML report available.
"""

import json
import os
import re
import shutil
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple


def sanitize_path_segment(value: Any, fallback: str = "item") -> str:
    text = str(value or "").strip()
    text = re.sub(r'[<>:"/\\|?*\s]+', "_", text)
    text = text.strip("._")
    return text or fallback


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def copy_if_exists(src: Optional[str], dst: str) -> Optional[str]:
    if not src or not os.path.exists(src):
        return None
    ensure_dir(os.path.dirname(dst))
    if os.path.abspath(src) != os.path.abspath(dst):
        shutil.copy2(src, dst)
    return dst


def relpath(path: Optional[str], start: str) -> Optional[str]:
    if not path:
        return None
    try:
        return os.path.relpath(path, start)
    except ValueError:
        return path


def read_fasta_length(fasta_file: Optional[str]) -> int:
    if not fasta_file or not os.path.exists(fasta_file):
        return 0
    length = 0
    with open(fasta_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(">"):
                length += len(line.strip())
    return length


def parse_coords_candidates(coords_file: Optional[str], ref_length: int = 0) -> List[Dict[str, Any]]:
    """Build query-genome candidate summaries from the BLAST-derived coords file."""
    if not coords_file or not os.path.exists(coords_file):
        return []

    grouped: Dict[str, Dict[str, Any]] = {}
    with open(coords_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("[") or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            try:
                s1, e1 = int(parts[0]), int(parts[1])
                s2, e2 = int(parts[2]), int(parts[3])
                len1 = int(float(parts[4]))
                len2 = int(float(parts[5]))
                identity = float(parts[6])
            except ValueError:
                continue

            query_chr = parts[7].split()[0]
            ref_seq = parts[8].split()[0]
            key = query_chr
            item = grouped.setdefault(
                key,
                {
                    "query_chr": query_chr,
                    "ref_seq": ref_seq,
                    "blocks": [],
                    "total_aln_len": 0,
                    "identity_weighted_sum": 0.0,
                    "ref_aligned_len": 0,
                    "query_start": min(s1, e1),
                    "query_end": max(s1, e1),
                    "strand_votes": {"+": 0, "-": 0},
                },
            )
            block_len = max(len1, len2)
            strand = "+" if e1 >= s1 else "-"
            item["blocks"].append(
                {
                    "query_start": s1,
                    "query_end": e1,
                    "ref_start": s2,
                    "ref_end": e2,
                    "query_aln_len": len1,
                    "ref_aln_len": len2,
                    "identity": identity,
                    "query_chr": query_chr,
                    "ref_seq": ref_seq,
                    "strand": strand,
                }
            )
            item["total_aln_len"] += block_len
            item["identity_weighted_sum"] += identity * block_len
            item["ref_aligned_len"] += len2
            item["query_start"] = min(item["query_start"], s1, e1)
            item["query_end"] = max(item["query_end"], s1, e1)
            item["strand_votes"][strand] += block_len

    candidates: List[Dict[str, Any]] = []
    for item in grouped.values():
        total = item["total_aln_len"]
        avg_identity = item["identity_weighted_sum"] / total if total else 0.0
        coverage = item["ref_aligned_len"] / ref_length if ref_length else 0.0
        strand = "+" if item["strand_votes"]["+"] >= item["strand_votes"]["-"] else "-"
        candidates.append(
            {
                "query_chr": item["query_chr"],
                "ref_seq": item["ref_seq"],
                "query_start": item["query_start"],
                "query_end": item["query_end"],
                "strand": strand,
                "total_aln_len": total,
                "identity": round(avg_identity, 4),
                "coverage": round(coverage, 6),
                "block_count": len(item["blocks"]),
                "blocks": item["blocks"],
            }
        )

    candidates.sort(key=lambda c: (-c["total_aln_len"], -c["identity"], -c["coverage"], c["query_chr"]))
    for index, candidate in enumerate(candidates, start=1):
        candidate["candidate_id"] = f"q1_{index}"
        candidate["rank"] = index
        candidate["is_best"] = index == 1
    return candidates


def parse_variant_stats(snps_file: Optional[str]) -> Dict[str, int]:
    stats = {"SNP": 0, "INS": 0, "DEL": 0, "INDEL": 0, "total": 0}
    if not snps_file or not os.path.exists(snps_file):
        return stats
    with open(snps_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("[") or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            var_type = parts[4].upper()
            if var_type in stats:
                stats[var_type] += 1
            if var_type in ("INS", "DEL", "INDEL"):
                stats["INDEL"] += 1
            stats["total"] += 1
    return stats


def write_pair_hl(snps_file: Optional[str], output_path: str) -> Optional[str]:
    """Write compatibility .hl markers using the current ref/query semantics."""
    ensure_dir(os.path.dirname(output_path))
    with open(output_path, "w", encoding="utf-8") as handle:
        if not snps_file or not os.path.exists(snps_file):
            return output_path
        with open(snps_file, "r", encoding="utf-8", errors="ignore") as snps:
            for line in snps:
                line = line.strip()
                if not line or line.startswith("[") or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                try:
                    query_genome_pos = int(parts[0])
                    ref_source_pos = int(parts[3])
                except ValueError:
                    continue
                var_type = parts[4].upper()
                query_chr = parts[5].split()[0] if len(parts) > 5 else "query"
                ref_seq = parts[6].split()[0] if len(parts) > 6 else "ref"
                color = "orange" if var_type == "SNP" else "blue"
                if var_type == "SNP":
                    handle.write(f"{ref_seq}\t{ref_source_pos - 1}\t{ref_source_pos}\t{color}\n")
                    handle.write(f"{query_chr}\t{query_genome_pos - 1}\t{query_genome_pos}\t{color}\n")
                elif var_type == "INS":
                    handle.write(f"{query_chr}\t{query_genome_pos - 1}\t{query_genome_pos}\t{color}\n")
                elif var_type == "DEL":
                    handle.write(f"{ref_seq}\t{ref_source_pos - 1}\t{ref_source_pos}\t{color}\n")
    return output_path


def _mode_ref_query_names(mode: str, result: Dict[str, Any], ref_name: str, qry_name: Optional[str]) -> Tuple[str, str]:
    result_id = result.get("id", "input")
    if mode == "sequence":
        return result_id, ref_name or "query_genome"
    return ref_name or "ref", qry_name or "query"


def build_single_query_payload(
    result: Dict[str, Any],
    output_dir: str,
    mode: str,
    ref_name: str = "",
    qry_name: Optional[str] = "",
    identity: float = 90,
    min_aln_len: int = 100,
    genome_files: Optional[Dict[str, Any]] = None,
    legacy_report: Optional[str] = None,
    extra_files: Optional[List[Tuple[str, Optional[str], str]]] = None,
) -> Dict[str, Any]:
    genome_files = genome_files or {}
    extra_files = extra_files or []
    result_id = result.get("id", "input")
    safe_id = sanitize_path_segment(result_id, "input")
    ref_entry_name, query_entry_name = _mode_ref_query_names(mode, result, ref_name, qry_name)
    query_safe = sanitize_path_segment(query_entry_name, "query")
    pair_id = f"ref__{query_safe}"

    ref_dir = ensure_dir(os.path.join(output_dir, "ref"))
    query_dir = ensure_dir(os.path.join(output_dir, "queries", query_safe))
    report_dir = ensure_dir(os.path.join(output_dir, "report"))

    ref_fasta_copy = copy_if_exists(result.get("fasta"), os.path.join(ref_dir, f"{safe_id}.fasta"))
    ref_gff_copy = copy_if_exists(result.get("gff"), os.path.join(ref_dir, f"{safe_id}.gff3"))

    pair_prefix = os.path.join(query_dir, pair_id)
    blast_copy = copy_if_exists(result.get("blast_xml"), f"{pair_prefix}.blast.xml")
    coords_copy = copy_if_exists(result.get("coords"), f"{pair_prefix}.coords")
    snps_copy = copy_if_exists(result.get("snps"), f"{pair_prefix}.snps")
    hl_copy = write_pair_hl(result.get("snps"), f"{pair_prefix}.hl")

    ref_length = read_fasta_length(result.get("fasta"))
    candidates = parse_coords_candidates(result.get("coords"), ref_length)
    stats = parse_variant_stats(result.get("snps"))

    if mode == "sequence":
        query_fasta = genome_files.get("ref_fasta")
        query_gff = genome_files.get("ref_gff")
        source_ref_fasta = result.get("fasta")
        source_ref_gff = None
    else:
        query_fasta = genome_files.get("qry_fasta")
        query_gff = genome_files.get("qry_gff")
        source_ref_fasta = genome_files.get("ref_fasta")
        source_ref_gff = genome_files.get("ref_gff")

    payload = {
        "schema_version": "multi_query_report.v1",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "mode": mode,
        "input": {
            "id": result_id,
            "safe_id": safe_id,
            "location": result.get("location"),
            "sequence_length": len(result.get("sequence", "")) if result.get("sequence") else ref_length,
            "extraction_info": result.get("extraction_info"),
        },
        "parameters": {
            "identity": identity,
            "min_aln_len": min_aln_len,
            "candidate_limit": 3,
            "pairwise_all": False,
        },
        "genomes": {
            "ref": {
                "genome_id": "ref",
                "name": ref_entry_name,
                "role": "ref_source",
                "source_fasta": source_ref_fasta,
                "source_gff": source_ref_gff,
                "artifact_fasta": relpath(ref_fasta_copy, report_dir),
                "artifact_gff": relpath(ref_gff_copy, report_dir),
            },
            "queries": [
                {
                    "genome_id": query_safe,
                    "name": query_entry_name,
                    "role": "query_genome",
                    "source_fasta": query_fasta,
                    "source_gff": query_gff,
                    "has_annotation": bool(query_gff),
                }
            ],
        },
        "pairs": [
            {
                "pair_id": pair_id,
                "type": "ref_query",
                "ref_genome_id": "ref",
                "query_genome_id": query_safe,
                "coords_semantics": {
                    "coords_ref_columns": "query genome BLAST subject",
                    "coords_query_columns": "ref-derived BLAST query sequence",
                    "variants_are_relative_to": "ref_source",
                },
                "artifacts": {
                    "blast_xml": relpath(blast_copy, report_dir),
                    "coords": relpath(coords_copy, report_dir),
                    "snps": relpath(snps_copy, report_dir),
                    "hl": relpath(hl_copy, report_dir),
                    "legacy_report": relpath(legacy_report, report_dir),
                },
                "candidate_ids": [candidate["candidate_id"] for candidate in candidates],
                "stats": {
                    "candidate_count": len(candidates),
                    "total_aln_len": sum(candidate["total_aln_len"] for candidate in candidates),
                    "variants": stats,
                },
            }
        ],
        "candidates": {
            query_safe: candidates,
        },
        "default_selection": {
            "track_order": ["ref", query_safe],
            "selected_candidates": {
                query_safe: candidates[0]["candidate_id"] if candidates else None,
            },
        },
        "statistics": {
            "genome_count": 2,
            "query_count": 1,
            "total_candidate_count": len(candidates),
            "per_query_candidate_counts": {query_safe: len(candidates)},
            "selected_combination": {
                "link_count": len(candidates[0]["blocks"]) if candidates else 0,
                "variants": stats,
            },
        },
        "extra_files": [
            {"label": label, "path": relpath(path, report_dir), "description": desc}
            for label, path, desc in extra_files
            if path
        ],
    }
    return payload


def write_static_report(
    payload: Dict[str, Any],
    output_dir: str,
    legacy_report: Optional[str] = None,
) -> str:
    report_dir = ensure_dir(os.path.join(output_dir, "report"))
    data_path = os.path.join(report_dir, "data.json")
    with open(data_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)

    index_path = os.path.join(report_dir, "index.html")
    if legacy_report and os.path.exists(legacy_report):
        with open(legacy_report, "r", encoding="utf-8", errors="ignore") as handle:
            html = handle.read()
        marker = '<meta name="genescreen-report-schema" content="multi_query_report.v1">'
        if "</head>" in html and marker not in html:
            html = html.replace("</head>", f"  {marker}\n</head>", 1)
        with open(index_path, "w", encoding="utf-8") as handle:
            handle.write(html)
    else:
        with open(index_path, "w", encoding="utf-8") as handle:
            handle.write(
                "<!DOCTYPE html><html lang=\"zh-CN\"><head><meta charset=\"utf-8\">"
                "<title>GeneScreen Report</title></head><body>"
                "<h1>GeneScreen Report</h1><p>Report data is available in data.json.</p>"
                "</body></html>"
            )
    return index_path
