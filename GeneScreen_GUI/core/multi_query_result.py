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


def parse_variants(snps_file: Optional[str]) -> List[Dict[str, Any]]:
    variants: List[Dict[str, Any]] = []
    if not snps_file or not os.path.exists(snps_file):
        return variants
    with open(snps_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
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
            variants.append(
                {
                    "query_genome_pos": query_genome_pos,
                    "ref_source_pos": ref_source_pos,
                    "ref": parts[1],
                    "alt": parts[2],
                    "type": parts[4].upper(),
                    "query_chr": parts[5].split()[0] if len(parts) > 5 else "query",
                    "ref_seq": parts[6].split()[0] if len(parts) > 6 else "ref",
                }
            )
    return variants


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


def build_multi_query_payload(
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
    query_results = result.get("query_results") or []
    query_entries = result.get("queries") or genome_files.get("queries") or []
    result_id = result.get("id", "input")
    safe_id = sanitize_path_segment(result_id, "input")
    ref_entry_name = result_id if mode == "sequence" else (ref_name or "ref")

    ref_dir = ensure_dir(os.path.join(output_dir, "ref"))
    report_dir = ensure_dir(os.path.join(output_dir, "report"))
    ref_fasta_copy = copy_if_exists(result.get("fasta"), os.path.join(ref_dir, f"{safe_id}.fasta"))
    ref_gff_copy = copy_if_exists(result.get("gff"), os.path.join(ref_dir, f"{safe_id}.gff3"))
    ref_length = read_fasta_length(result.get("fasta"))

    tracks = [{
        "track_id": "ref",
        "genome_id": "ref",
        "name": ref_entry_name,
        "role": "ref_source",
        "sequence_id": result_id,
        "length": ref_length,
        "artifact_fasta": relpath(ref_fasta_copy, report_dir),
        "artifact_gff": relpath(ref_gff_copy, report_dir),
    }]
    queries_payload = []
    pairs = []
    candidates_map = {}
    variants_payload = []
    visible_links = []
    query_order = []
    selected_candidates = {}
    best_candidates = {}
    per_query_candidate_counts = {}
    total_candidate_count = 0
    combined_variant_stats = {"SNP": 0, "INS": 0, "DEL": 0, "INDEL": 0, "total": 0}

    def entry_for(index: int, query_result: Dict[str, Any]) -> Dict[str, Any]:
        if index < len(query_entries):
            return query_entries[index]
        return {
            "name": query_result.get("query_name") or f"query_{index + 1}",
            "fasta": query_result.get("query_fasta"),
            "gff": query_result.get("query_gff"),
        }

    for index, query_result in enumerate(query_results):
        entry = entry_for(index, query_result)
        query_name = entry.get("name") or query_result.get("query_name") or f"query_{index + 1}"
        query_safe = sanitize_path_segment(query_name, f"query_{index + 1}")
        pair_id = f"ref__{query_safe}"
        query_dir = ensure_dir(os.path.join(output_dir, "queries", query_safe))
        pair_prefix = os.path.join(query_dir, pair_id)

        blast_copy = copy_if_exists(query_result.get("blast_xml"), f"{pair_prefix}.blast.xml")
        coords_copy = copy_if_exists(query_result.get("coords"), f"{pair_prefix}.coords")
        snps_copy = copy_if_exists(query_result.get("snps"), f"{pair_prefix}.snps")
        hl_copy = write_pair_hl(query_result.get("snps"), f"{pair_prefix}.hl")
        candidates = parse_coords_candidates(query_result.get("coords"), ref_length)
        variants = parse_variants(query_result.get("snps"))
        stats = parse_variant_stats(query_result.get("snps"))
        selected_candidate = candidates[0] if candidates else None

        query_order.append(query_safe)
        selected_candidates[query_safe] = selected_candidate["candidate_id"] if selected_candidate else None
        best_candidates[query_safe] = selected_candidate["candidate_id"] if selected_candidate else None
        candidates_map[query_safe] = candidates
        per_query_candidate_counts[query_safe] = len(candidates)
        total_candidate_count += len(candidates)
        for key in combined_variant_stats:
            combined_variant_stats[key] += stats.get(key, 0)

        queries_payload.append({
            "genome_id": query_safe,
            "name": query_name,
            "role": "query_genome",
            "source_fasta": entry.get("fasta") or query_result.get("query_fasta"),
            "source_gff": entry.get("gff") or query_result.get("query_gff"),
            "has_annotation": bool(entry.get("gff") or query_result.get("query_gff")),
        })
        tracks.append({
            "track_id": query_safe,
            "genome_id": query_safe,
            "name": query_name,
            "role": "query_genome",
            "candidate_count": len(candidates),
            "selected_candidate_id": selected_candidate["candidate_id"] if selected_candidate else None,
            "has_annotation": bool(entry.get("gff") or query_result.get("query_gff")),
        })
        pairs.append({
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
                "legacy_report": relpath(legacy_report, report_dir) if index == 0 else None,
            },
            "candidate_ids": [candidate["candidate_id"] for candidate in candidates],
            "stats": {
                "candidate_count": len(candidates),
                "total_aln_len": sum(candidate["total_aln_len"] for candidate in candidates),
                "variants": stats,
            },
        })
        if selected_candidate:
            for block_index, block in enumerate(selected_candidate.get("blocks", []), start=1):
                visible_links.append({
                    "link_id": f"{pair_id}_link_{block_index}",
                    "pair_id": pair_id,
                    "source_track_id": "ref",
                    "target_track_id": query_safe,
                    "ref_start": block["ref_start"],
                    "ref_end": block["ref_end"],
                    "query_start": block["query_start"],
                    "query_end": block["query_end"],
                    "identity": block["identity"],
                    "strand": block["strand"],
                })
        for variant_index, variant in enumerate(variants, start=1):
            variants_payload.append({
                "variant_id": f"{pair_id}_var_{variant_index}",
                "pair_id": pair_id,
                **variant,
            })

    for pairwise_result in result.get("pairwise_results") or []:
        pair_id = pairwise_result.get("pair_id") or "query__query"
        comparison_id = sanitize_path_segment(
            f"{pairwise_result.get('left_query')}_{pairwise_result.get('left_candidate_id')}__"
            f"{pairwise_result.get('right_query')}_{pairwise_result.get('right_candidate_id')}",
            "pairwise"
        )
        pair_dir = ensure_dir(os.path.join(output_dir, "pairwise", sanitize_path_segment(pair_id, "pairwise")))
        pair_prefix = os.path.join(pair_dir, comparison_id)
        blast_copy = copy_if_exists(pairwise_result.get("blast_xml"), f"{pair_prefix}.blast.xml")
        coords_copy = copy_if_exists(pairwise_result.get("coords"), f"{pair_prefix}.coords")
        snps_copy = copy_if_exists(pairwise_result.get("snps"), f"{pair_prefix}.snps")
        hl_copy = write_pair_hl(pairwise_result.get("snps"), f"{pair_prefix}.hl")
        left_id = sanitize_path_segment(pairwise_result.get("left_query") or "left", "left")
        right_id = sanitize_path_segment(pairwise_result.get("right_query") or "right", "right")
        pairs.append({
            "pair_id": f"{pair_id}::{comparison_id}",
            "type": "query_query",
            "ref_genome_id": left_id,
            "query_genome_id": right_id,
            "left_candidate_id": pairwise_result.get("left_candidate_id"),
            "right_candidate_id": pairwise_result.get("right_candidate_id"),
            "artifacts": {
                "blast_xml": relpath(blast_copy, report_dir),
                "coords": relpath(coords_copy, report_dir),
                "snps": relpath(snps_copy, report_dir),
                "hl": relpath(hl_copy, report_dir),
            },
            "stats": {
                "variants": parse_variant_stats(pairwise_result.get("snps")),
            },
        })

    source_ref_fasta = genome_files.get("ref_fasta")
    source_ref_gff = genome_files.get("ref_gff")
    return {
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
            "candidate_limit": result.get("candidate_limit", 3),
            "pairwise_all": bool(result.get("pairwise_all", False)),
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
            "queries": queries_payload,
        },
        "tracks": tracks,
        "pairs": pairs,
        "candidates": candidates_map,
        "default_selection": {
            "track_order": ["ref"] + query_order,
            "selected_candidates": selected_candidates,
        },
        "overview": {
            "query_order": query_order,
            "best_candidates": best_candidates,
            "candidate_ids": {
                query_id: [candidate["candidate_id"] for candidate in candidates]
                for query_id, candidates in candidates_map.items()
            },
        },
        "detail": {
            "track_order": ["ref"] + query_order,
            "selected_candidates": selected_candidates,
            "visible_pair_ids": [pair["pair_id"] for pair in pairs],
            "visible_links": visible_links,
            "variant_ids": [variant["variant_id"] for variant in variants_payload],
        },
        "variants": variants_payload,
        "statistics": {
            "genome_count": 1 + len(query_order),
            "query_count": len(query_order),
            "total_candidate_count": total_candidate_count,
            "per_query_candidate_counts": per_query_candidate_counts,
            "selected_combination": {
                "link_count": len(visible_links),
                "variants": combined_variant_stats,
            },
        },
        "extra_files": [
            {"label": label, "path": relpath(path, report_dir), "description": desc}
            for label, path, desc in extra_files
            if path
        ],
    }


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
    if len(result.get("query_results") or []) > 1:
        return build_multi_query_payload(
            result=result,
            output_dir=output_dir,
            mode=mode,
            ref_name=ref_name,
            qry_name=qry_name,
            identity=identity,
            min_aln_len=min_aln_len,
            genome_files=genome_files,
            legacy_report=legacy_report,
            extra_files=extra_files,
        )
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
    variants = parse_variants(result.get("snps"))
    selected_candidate = candidates[0] if candidates else None
    visible_links = []
    if selected_candidate:
        for index, block in enumerate(selected_candidate.get("blocks", []), start=1):
            visible_links.append(
                {
                    "link_id": f"{pair_id}_link_{index}",
                    "pair_id": pair_id,
                    "source_track_id": "ref",
                    "target_track_id": query_safe,
                    "ref_start": block["ref_start"],
                    "ref_end": block["ref_end"],
                    "query_start": block["query_start"],
                    "query_end": block["query_end"],
                    "identity": block["identity"],
                    "strand": block["strand"],
                }
            )

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
            "candidate_limit": result.get("candidate_limit", 3),
            "pairwise_all": bool(result.get("pairwise_all", False)),
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
        "tracks": [
            {
                "track_id": "ref",
                "genome_id": "ref",
                "name": ref_entry_name,
                "role": "ref_source",
                "sequence_id": result_id,
                "length": ref_length,
                "artifact_fasta": relpath(ref_fasta_copy, report_dir),
                "artifact_gff": relpath(ref_gff_copy, report_dir),
            },
            {
                "track_id": query_safe,
                "genome_id": query_safe,
                "name": query_entry_name,
                "role": "query_genome",
                "candidate_count": len(candidates),
                "selected_candidate_id": selected_candidate["candidate_id"] if selected_candidate else None,
                "has_annotation": bool(query_gff),
            },
        ],
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
                query_safe: selected_candidate["candidate_id"] if selected_candidate else None,
            },
        },
        "overview": {
            "query_order": [query_safe],
            "best_candidates": {
                query_safe: selected_candidate["candidate_id"] if selected_candidate else None,
            },
            "candidate_ids": {
                query_safe: [candidate["candidate_id"] for candidate in candidates],
            },
        },
        "detail": {
            "track_order": ["ref", query_safe],
            "selected_candidates": {
                query_safe: selected_candidate["candidate_id"] if selected_candidate else None,
            },
            "visible_pair_ids": [pair_id],
            "visible_links": visible_links,
            "variant_ids": [f"var_{index}" for index in range(1, len(variants) + 1)],
        },
        "variants": [
            {"variant_id": f"var_{index}", "pair_id": pair_id, **variant}
            for index, variant in enumerate(variants, start=1)
        ],
        "statistics": {
            "genome_count": 2,
            "query_count": 1,
            "total_candidate_count": len(candidates),
            "per_query_candidate_counts": {query_safe: len(candidates)},
            "selected_combination": {
                "link_count": len(visible_links),
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
    if payload.get("statistics", {}).get("query_count", 0) > 1:
        with open(index_path, "w", encoding="utf-8") as handle:
            handle.write(_render_multi_query_report_html(payload))
    elif legacy_report and os.path.exists(legacy_report):
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


def _render_multi_query_report_html(payload: Dict[str, Any]) -> str:
    embedded = json.dumps(payload, ensure_ascii=False)
    title = payload.get("input", {}).get("id") or "GeneScreen Report"
    return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="genescreen-report-schema" content="multi_query_report.v1">
  <title>{title} - GeneScreen</title>
  <style>
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; font-family: Arial, 'Microsoft YaHei', sans-serif; color: #17202a; background: #f6f8fb; }}
    header {{ padding: 18px 24px; background: #ffffff; border-bottom: 1px solid #d9e0ea; }}
    h1 {{ margin: 0 0 6px; font-size: 22px; font-weight: 700; }}
    .subtle {{ color: #5d6979; font-size: 13px; }}
    main {{ padding: 18px 24px 32px; max-width: 1500px; margin: 0 auto; }}
    .cards {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 12px; margin-bottom: 14px; }}
    .card, .panel {{ background: #fff; border: 1px solid #d9e0ea; border-radius: 6px; }}
    .card {{ padding: 14px 16px; min-height: 96px; }}
    .card h2, .panel h2 {{ margin: 0 0 10px; font-size: 15px; }}
    .metrics {{ display: flex; flex-wrap: wrap; gap: 14px; }}
    .metric strong {{ display: block; font-size: 22px; line-height: 1.1; }}
    .metric span {{ color: #5d6979; font-size: 12px; }}
    .panel {{ padding: 14px 16px; margin-bottom: 14px; }}
    .overview-row, .track-row {{ display: grid; grid-template-columns: 180px 1fr; gap: 12px; align-items: center; min-height: 44px; }}
    .track-label {{ font-weight: 600; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
    .lane {{ position: relative; height: 34px; border-left: 1px solid #c6d0dd; border-right: 1px solid #c6d0dd; background: linear-gradient(#fff, #fff) padding-box; }}
    .candidate {{ position: absolute; top: 9px; height: 16px; min-width: 3px; background: #5b8def; border: 1px solid #2e63c7; cursor: pointer; }}
    .candidate.best {{ background: #27ae60; border-color: #1e874b; }}
    .candidate.selected {{ outline: 2px solid #111827; outline-offset: 1px; }}
    .detail-layout {{ display: grid; grid-template-columns: 260px 1fr; gap: 14px; }}
    .controls {{ border-right: 1px solid #e0e5ec; padding-right: 12px; }}
    .control-row {{ display: grid; grid-template-columns: 1fr 32px 32px; gap: 6px; margin-bottom: 8px; align-items: center; }}
    select, button {{ height: 30px; border: 1px solid #bcc7d5; border-radius: 4px; background: #fff; }}
    button {{ cursor: pointer; }}
    .track-stack {{ position: relative; min-height: 180px; }}
    .track-row {{ border: 1px solid #dde4ee; border-radius: 6px; background: #fbfcfe; padding: 8px; margin-bottom: 8px; cursor: grab; }}
    .track-row.dragging {{ opacity: .5; }}
    svg.links {{ position: absolute; inset: 0; width: 100%; height: 100%; pointer-events: none; overflow: visible; }}
    .empty {{ color: #6b7280; padding: 18px; }}
    @media (max-width: 900px) {{
      main {{ padding: 12px; }}
      .cards, .detail-layout {{ grid-template-columns: 1fr; }}
      .overview-row, .track-row {{ grid-template-columns: 120px 1fr; }}
      .controls {{ border-right: 0; padding-right: 0; }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>{title}</h1>
    <div class="subtle">GeneScreen multi-query static report</div>
  </header>
  <main>
    <section class="cards">
      <div class="card"><h2>输入与候选</h2><div id="inputMetrics" class="metrics"></div></div>
      <div class="card"><h2>当前组合统计</h2><div id="variantMetrics" class="metrics"></div></div>
    </section>
    <section class="panel">
      <h2>Overview</h2>
      <div id="overview"></div>
    </section>
    <section class="panel">
      <h2>Detail</h2>
      <div class="detail-layout">
        <div id="controls" class="controls"></div>
        <div id="detail" class="track-stack"></div>
      </div>
    </section>
  </main>
  <script id="embedded-report-data" type="application/json">{embedded}</script>
  <script>
    let reportData = JSON.parse(document.getElementById('embedded-report-data').textContent);
    fetch('data.json').then(r => r.ok ? r.json() : reportData).then(data => {{ reportData = data; init(); }}).catch(init);

    const state = {{ selected: {{}}, order: [] }};
    function init() {{
      state.selected = Object.assign({{}}, reportData.default_selection.selected_candidates || {{}});
      state.order = (reportData.default_selection.track_order || []).slice();
      renderCards();
      renderOverview();
      renderControls();
      renderDetail();
    }}
    function queryIds() {{ return (reportData.overview && reportData.overview.query_order) || []; }}
    function candidateById(queryId, candidateId) {{
      return (reportData.candidates[queryId] || []).find(c => c.candidate_id === candidateId);
    }}
    function renderCards() {{
      const s = reportData.statistics || {{}};
      const variants = (s.selected_combination && s.selected_combination.variants) || {{}};
      document.getElementById('inputMetrics').innerHTML = [
        metric('Genomes', s.genome_count || 0),
        metric('Queries', s.query_count || 0),
        metric('Candidates', s.total_candidate_count || 0)
      ].join('');
      document.getElementById('variantMetrics').innerHTML = [
        metric('Links', (s.selected_combination && s.selected_combination.link_count) || 0),
        metric('SNP', variants.SNP || 0),
        metric('INS', variants.INS || 0),
        metric('DEL', variants.DEL || 0)
      ].join('');
    }}
    function metric(label, value) {{ return `<div class="metric"><strong>${{value}}</strong><span>${{label}}</span></div>`; }}
    function renderOverview() {{
      const root = document.getElementById('overview');
      root.innerHTML = '';
      for (const qid of queryIds()) {{
        const track = reportData.tracks.find(t => t.track_id === qid) || {{ name: qid }};
        const row = document.createElement('div');
        row.className = 'overview-row';
        row.innerHTML = `<div class="track-label" title="${{track.name}}">${{track.name}}</div><div class="lane"></div>`;
        const lane = row.querySelector('.lane');
        const candidates = reportData.candidates[qid] || [];
        const maxEnd = Math.max(...candidates.map(c => c.query_end || 1), 1);
        for (const candidate of candidates) {{
          const el = document.createElement('div');
          el.className = 'candidate' + (candidate.is_best ? ' best' : '') + (state.selected[qid] === candidate.candidate_id ? ' selected' : '');
          const start = Math.max(0, ((candidate.query_start || 1) / maxEnd) * 100);
          const end = Math.max(start + 0.5, ((candidate.query_end || candidate.query_start || 1) / maxEnd) * 100);
          el.style.left = `${{Math.min(start, 99)}}%`;
          el.style.width = `${{Math.max(0.5, Math.min(end - start, 100 - start))}}%`;
          el.title = `${{candidate.candidate_id}} ${{candidate.query_chr}}:${{candidate.query_start}}-${{candidate.query_end}}`;
          el.onclick = () => {{ state.selected[qid] = candidate.candidate_id; renderOverview(); renderControls(); renderDetail(); }};
          lane.appendChild(el);
        }}
        root.appendChild(row);
      }}
    }}
    function renderControls() {{
      const root = document.getElementById('controls');
      root.innerHTML = '';
      for (const qid of queryIds()) {{
        const track = reportData.tracks.find(t => t.track_id === qid) || {{ name: qid }};
        const row = document.createElement('div');
        row.className = 'control-row';
        const options = (reportData.candidates[qid] || []).map(c => `<option value="${{c.candidate_id}}" ${{state.selected[qid] === c.candidate_id ? 'selected' : ''}}>${{track.name}} · ${{c.candidate_id}}</option>`).join('');
        row.innerHTML = `<select>${{options}}</select><button title="上移">↑</button><button title="下移">↓</button>`;
        row.querySelector('select').onchange = e => {{ state.selected[qid] = e.target.value; renderOverview(); renderDetail(); }};
        row.children[1].onclick = () => moveTrack(qid, -1);
        row.children[2].onclick = () => moveTrack(qid, 1);
        root.appendChild(row);
      }}
    }}
    function moveTrack(trackId, delta) {{
      const index = state.order.indexOf(trackId);
      const target = index + delta;
      if (index <= 0 || target <= 0 || target >= state.order.length) return;
      state.order.splice(index, 1);
      state.order.splice(target, 0, trackId);
      renderDetail();
    }}
    function renderDetail() {{
      const root = document.getElementById('detail');
      root.innerHTML = '';
      if (!state.order.length) {{ root.innerHTML = '<div class="empty">No tracks</div>'; return; }}
      for (const trackId of state.order) {{
        const track = reportData.tracks.find(t => t.track_id === trackId) || {{ name: trackId, role: 'query_genome' }};
        const row = document.createElement('div');
        row.className = 'track-row';
        row.draggable = trackId !== 'ref';
        row.dataset.trackId = trackId;
        const candidate = trackId === 'ref' ? null : candidateById(trackId, state.selected[trackId]);
        const label = candidate ? `${{track.name}} · ${{candidate.candidate_id}} · ${{candidate.query_chr}}:${{candidate.query_start}}-${{candidate.query_end}}` : `${{track.name}} · ref`;
        row.innerHTML = `<div class="track-label" title="${{label}}">${{label}}</div><div class="lane"></div>`;
        const lane = row.querySelector('.lane');
        if (candidate) {{
          const marker = document.createElement('div');
          marker.className = 'candidate selected';
          marker.style.left = '8%';
          marker.style.width = '84%';
          lane.appendChild(marker);
        }}
        wireDrag(row);
        root.appendChild(row);
      }}
    }}
    function wireDrag(row) {{
      row.addEventListener('dragstart', e => {{ row.classList.add('dragging'); e.dataTransfer.setData('text/plain', row.dataset.trackId); }});
      row.addEventListener('dragend', () => row.classList.remove('dragging'));
      row.addEventListener('dragover', e => e.preventDefault());
      row.addEventListener('drop', e => {{
        e.preventDefault();
        const dragged = e.dataTransfer.getData('text/plain');
        const target = row.dataset.trackId;
        if (!dragged || dragged === target || target === 'ref') return;
        state.order = state.order.filter(id => id !== dragged);
        const index = state.order.indexOf(target);
        state.order.splice(index, 0, dragged);
        renderDetail();
      }});
    }}
  </script>
</body>
</html>
"""
