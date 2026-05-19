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
from html import escape
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


def report_candidate_limit(result: Dict[str, Any]) -> Optional[int]:
    if bool(result.get("pairwise_all", False)):
        return None
    return result.get("candidate_limit", 3)


def read_fasta_length(fasta_file: Optional[str]) -> int:
    if not fasta_file or not os.path.exists(fasta_file):
        return 0
    length = 0
    with open(fasta_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(">"):
                length += len(line.strip())
    return length


def parse_coords_candidates(
    coords_file: Optional[str],
    ref_length: int = 0,
    merge_gap: int = 1000,
    min_aln_len: int = 0,
) -> List[Dict[str, Any]]:
    """Build query-genome candidate summaries from the BLAST-derived coords file."""
    if not coords_file or not os.path.exists(coords_file):
        return []

    blocks_by_chr: Dict[str, List[Dict[str, Any]]] = {}
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
            if min_aln_len and int((len1 + len2) / 2) < min_aln_len:
                continue

            query_chr = parts[7].split()[0]
            ref_seq = parts[8].split()[0]
            strand = "+" if e1 >= s1 else "-"
            blocks_by_chr.setdefault(query_chr, []).append(
                {
                    "query_chr": query_chr,
                    "ref_seq": ref_seq,
                    "query_start": min(s1, e1),
                    "query_end": max(s1, e1),
                    "raw_query_start": s1,
                    "raw_query_end": e1,
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

    def build_candidate(group: List[Dict[str, Any]]) -> Dict[str, Any]:
        item = {
            "query_chr": group[0]["query_chr"],
            "ref_seq": group[0]["ref_seq"],
            "blocks": [],
            "total_aln_len": 0,
            "identity_weighted_sum": 0.0,
            "ref_aligned_len": 0,
            "query_start": group[0]["query_start"],
            "query_end": group[0]["query_end"],
            "strand_votes": {"+": 0, "-": 0},
        }
        for block in group:
            block_len = max(int(block["query_aln_len"]), int(block["ref_aln_len"]))
            strand = block.get("strand", "+")
            item["blocks"].append(
                {
                    "query_start": block["raw_query_start"],
                    "query_end": block["raw_query_end"],
                    "ref_start": block["ref_start"],
                    "ref_end": block["ref_end"],
                    "query_aln_len": block["query_aln_len"],
                    "ref_aln_len": block["ref_aln_len"],
                    "identity": block["identity"],
                    "query_chr": block["query_chr"],
                    "ref_seq": block["ref_seq"],
                    "strand": strand,
                }
            )
            item["total_aln_len"] += block_len
            item["identity_weighted_sum"] += float(block["identity"]) * block_len
            item["ref_aligned_len"] += int(block["ref_aln_len"])
            item["query_start"] = min(item["query_start"], block["query_start"])
            item["query_end"] = max(item["query_end"], block["query_end"])
            item["strand_votes"][strand] += block_len
        return item

    candidates: List[Dict[str, Any]] = []
    grouped: List[Dict[str, Any]] = []
    merge_gap = max(0, int(merge_gap or 0))
    for chr_blocks in blocks_by_chr.values():
        sorted_blocks = sorted(chr_blocks, key=lambda b: (b["query_start"], b["query_end"]))
        current: List[Dict[str, Any]] = []
        current_end = 0
        for block in sorted_blocks:
            if not current:
                current = [block]
                current_end = block["query_end"]
                continue
            gap = block["query_start"] - current_end
            if gap <= merge_gap:
                current.append(block)
                current_end = max(current_end, block["query_end"])
            else:
                grouped.append(build_candidate(current))
                current = [block]
                current_end = block["query_end"]
        if current:
            grouped.append(build_candidate(current))

    for item in grouped:
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


def apply_query_flank_metadata(
    candidates: List[Dict[str, Any]],
    upstream: int = 0,
    downstream: int = 0,
) -> List[Dict[str, Any]]:
    if upstream <= 0 and downstream <= 0:
        return candidates
    for candidate in candidates:
        left = min(int(candidate.get("query_start") or 0), int(candidate.get("query_end") or 0))
        right = max(int(candidate.get("query_start") or 0), int(candidate.get("query_end") or 0))
        if left <= 0 or right <= 0:
            continue
        strand = candidate.get("strand", "+")
        if strand == "-":
            region_start = max(1, left - downstream)
            region_end = right + upstream
        else:
            region_start = max(1, left - upstream)
            region_end = right + downstream
        candidate["query_match_start"] = left
        candidate["query_match_end"] = right
        candidate["query_region_start"] = region_start
        candidate["query_region_end"] = region_end
        candidate["query_upstream"] = upstream
        candidate["query_downstream"] = downstream
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
    legacy_report: Optional[Any] = None,
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
    merge_gap = int(result.get("merge_gap") or 1000)

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
    legacy_extra_files = []
    combined_variant_stats = {"SNP": 0, "INS": 0, "DEL": 0, "INDEL": 0, "total": 0}

    def entry_for(index: int, query_result: Dict[str, Any]) -> Dict[str, Any]:
        if index < len(query_entries):
            return query_entries[index]
        return {
            "name": query_result.get("query_name") or f"query_{index + 1}",
            "fasta": query_result.get("query_fasta"),
            "gff": query_result.get("query_gff"),
        }

    def legacy_report_for(index: int, query_safe: str) -> Optional[str]:
        if isinstance(legacy_report, dict):
            return legacy_report.get(query_safe) or legacy_report.get(str(index))
        if isinstance(legacy_report, (list, tuple)):
            return legacy_report[index] if index < len(legacy_report) else None
        return legacy_report if index == 0 else None

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
        candidates = apply_query_flank_metadata(
            parse_coords_candidates(query_result.get("coords"), ref_length, merge_gap, min_aln_len),
            int(result.get("query_upstream") or 0),
            int(result.get("query_downstream") or 0),
        )
        variants = parse_variants(query_result.get("snps"))
        stats = parse_variant_stats(query_result.get("snps"))
        selected_candidate = candidates[0] if candidates else None
        pair_legacy_report = legacy_report_for(index, query_safe)
        if pair_legacy_report:
            legacy_extra_files.append({
                "label": os.path.basename(pair_legacy_report),
                "path": relpath(pair_legacy_report, report_dir),
                "description": f"旧版单报告：{ref_entry_name} vs {query_name}",
            })

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
                "legacy_report": relpath(pair_legacy_report, report_dir),
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
            "candidate_limit": report_candidate_limit(result),
            "pairwise_all": bool(result.get("pairwise_all", False)),
            "merge_gap": merge_gap,
            "query_upstream": result.get("query_upstream", 0),
            "query_downstream": result.get("query_downstream", 0),
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
        "extra_files": legacy_extra_files + [
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
    legacy_report: Optional[Any] = None,
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
    merge_gap = int(result.get("merge_gap") or 1000)
    candidates = apply_query_flank_metadata(
        parse_coords_candidates(result.get("coords"), ref_length, merge_gap, min_aln_len),
        int(result.get("query_upstream") or 0),
        int(result.get("query_downstream") or 0),
    )
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
            "candidate_limit": report_candidate_limit(result),
            "pairwise_all": bool(result.get("pairwise_all", False)),
            "merge_gap": merge_gap,
            "query_upstream": result.get("query_upstream", 0),
            "query_downstream": result.get("query_downstream", 0),
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
    with open(index_path, "w", encoding="utf-8") as handle:
        handle.write(_render_multi_query_report_html(payload))
    return index_path


def _render_multi_query_report_html(payload: Dict[str, Any]) -> str:
    embedded = json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")
    title = str(payload.get("input", {}).get("id") or "GeneScreen Report")
    title_html = escape(title)
    template = r"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="genescreen-report-schema" content="multi_query_report.v1">
  <title>__TITLE__ - GeneScreen Report</title>
  <style>
    * { margin: 0; padding: 0; box-sizing: border-box; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, "Microsoft YaHei", sans-serif;
      line-height: 1.6;
      color: #333;
      background: #f5f5f5;
      padding: 20px;
    }
    .container {
      max-width: 1320px;
      margin: 0 auto;
      background: #fff;
      border-radius: 8px;
      box-shadow: 0 2px 10px rgba(0,0,0,0.1);
      overflow: hidden;
    }
    .header {
      background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
      color: white;
      padding: 30px;
    }
    .header-content {
      display: flex;
      justify-content: space-between;
      gap: 20px;
      align-items: flex-start;
    }
    .header h1 { font-size: 24px; margin-bottom: 10px; }
    .subtitle {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 12px;
      opacity: 0.95;
      font-size: 14px;
    }
    .gen-time { margin-top: 8px; font-size: 13px; opacity: 0.85; }
    .mode-badge {
      display: inline-block;
      padding: 4px 12px;
      border-radius: 20px;
      font-size: 12px;
      font-weight: 600;
      text-transform: uppercase;
      background: rgba(255,255,255,0.22);
      color: white;
    }
    .json-badge {
      display: inline-block;
      color: rgba(255,255,255,0.92);
      border: 1px solid rgba(255,255,255,0.35);
      border-radius: 6px;
      padding: 5px 10px;
      font-size: 12px;
      text-decoration: none;
    }
    .section { padding: 25px 30px; border-bottom: 1px solid #eee; }
    .section:last-child { border-bottom: none; }
    .section h2 {
      font-size: 18px;
      color: #667eea;
      margin-bottom: 15px;
      padding-bottom: 10px;
      border-bottom: 2px solid #667eea;
    }
    .section h3 { font-size: 14px; color: #666; margin: 15px 0 10px; }
    .info-table { width: 100%; border-collapse: collapse; margin-bottom: 15px; }
    .info-table td { padding: 10px 15px; border-bottom: 1px solid #eee; vertical-align: top; }
    .info-table .label { width: 150px; font-weight: 600; color: #666; background: #f9f9f9; }
    .info-table .value { color: #333; }
    .info-table code {
      background: #f0f0f0;
      padding: 2px 6px;
      border-radius: 3px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      word-break: break-all;
    }
    .stats-table { width: 100%; border-collapse: collapse; margin-bottom: 15px; text-align: center; }
    .stats-table th {
      padding: 12px 15px;
      background: #f9f9f9;
      font-weight: 600;
      color: #666;
      border-bottom: 2px solid #eee;
    }
    .stats-table td {
      padding: 15px;
      font-size: 18px;
      font-weight: 500;
      color: #333;
      border-bottom: 1px solid #eee;
    }
    .file-table { width: 100%; border-collapse: collapse; }
    .file-table th, .file-table td { padding: 10px 15px; text-align: left; border-bottom: 1px solid #eee; vertical-align: top; }
    .file-table th { background: #f9f9f9; font-weight: 600; color: #666; }
    .path-code {
      font-size: 11px;
      color: #666;
      word-break: break-all;
      background: #f5f5f5;
      padding: 2px 6px;
      border-radius: 3px;
    }
    .file-na { color: #999; font-size: 12px; }
    .visualization-section { background: #fafafa; }
    .viz-container {
      padding: 20px;
      background: white;
      border: 1px solid #eee;
      border-radius: 8px;
      overflow-x: auto;
    }
    .overview-canvas {
      min-width: 920px;
      display: grid;
      grid-template-columns: 150px minmax(0, 1fr);
      gap: 18px;
      align-items: start;
    }
    .overview-ref {
      position: sticky;
      left: 0;
      min-height: 100%;
      border-left: 3px solid #222;
      border-bottom: 3px solid #222;
      padding: 12px 0 18px 18px;
      font-weight: 600;
      color: #333;
    }
    .overview-ref-title { font-size: 15px; margin-bottom: 8px; }
    .overview-ref-id { font-size: 12px; color: #777; font-weight: 400; word-break: break-all; }
    .overview-body { display: grid; gap: 16px; }
    .chr-group {
      display: grid;
      grid-template-columns: 96px minmax(0, 1fr);
      gap: 16px;
      align-items: center;
      padding: 4px 0;
    }
    .chr-label {
      font-size: 20px;
      color: #444;
      font-weight: 600;
      text-align: right;
      padding-right: 4px;
    }
    .chr-tracks { display: grid; gap: 10px; }
    .overview-track {
      display: grid;
      grid-template-columns: minmax(120px, 180px) minmax(420px, 1fr);
      gap: 12px;
      align-items: center;
    }
    .track-label {
      font-weight: 600;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      color: #333;
    }
    .track-subtext { font-size: 12px; color: #777; margin-top: 2px; font-weight: 400; }
    .lane {
      position: relative;
      height: 36px;
      background: transparent;
    }
    .lane::before {
      content: '';
      position: absolute;
      left: 0;
      right: 0;
      top: 17px;
      height: 3px;
      background: #252525;
      border-radius: 3px;
    }
    .lane-tick {
      position: absolute;
      top: 9px;
      width: 3px;
      height: 19px;
      background: #555;
      border-radius: 3px;
      transform: translateX(-50%);
    }
    .candidate {
      position: absolute;
      top: 12px;
      height: 12px;
      min-width: 3px;
      background: #4f86d9;
      border: 1px solid #2f5f9f;
      border-radius: 999px;
      cursor: pointer;
      transition: transform .12s ease, box-shadow .12s ease;
    }
    .candidate:hover { transform: translateY(-1px); box-shadow: 0 2px 6px rgba(46,99,199,.25); }
    .candidate.best { background: #159447; border-color: #0f6c34; }
    .candidate.selected { outline: 2px solid #111; outline-offset: 3px; }
    .candidate-ref { background: #667eea; border-color: #5364c9; left: 8%; width: 84%; }
    .viz-legend-note {
      margin-top: 15px;
      padding: 10px 15px;
      background: #f5f5f5;
      border-radius: 4px;
      font-size: 13px;
      color: #666;
    }
    .detail-layout {
      display: grid;
      grid-template-columns: 300px minmax(0, 1fr);
      gap: 18px;
      align-items: start;
    }
    .controls {
      background: #fff;
      border: 1px solid #eee;
      border-radius: 8px;
      padding: 14px;
    }
    .control-row {
      display: grid;
      grid-template-columns: 1fr 32px 32px;
      gap: 6px;
      margin-bottom: 10px;
      align-items: center;
    }
    .control-row[draggable="true"] { cursor: grab; }
    .control-row.dragging { opacity: .55; }
    select, button {
      height: 32px;
      border: 1px solid #ddd;
      border-radius: 4px;
      background: #fff;
      color: #333;
    }
    select { width: 100%; padding: 0 8px; }
    button { cursor: pointer; color: #667eea; font-weight: 600; }
    button:hover { background: #f3f0ff; }
    .track-stack { min-width: 760px; }
    .detail-canvas {
      min-width: 820px;
      border: 1px solid #e6e8ef;
      border-radius: 8px;
      background: linear-gradient(180deg, #fff 0%, #fbfcff 100%);
      overflow-x: auto;
    }
    .detail-svg {
      display: block;
      width: 100%;
      min-width: 820px;
      height: auto;
    }
    .detail-track-bg { fill: #fff; }
    .detail-track-bg:nth-of-type(even) { fill: #fafbff; }
    .detail-label { font-size: 13px; font-weight: 600; fill: #333; }
    .detail-subtext { font-size: 11px; fill: #777; }
    .detail-axis { stroke: #252525; stroke-width: 3; stroke-linecap: round; }
    .detail-axis-light { stroke: #a9b1c3; stroke-width: 1; }
    .detail-block { fill: rgba(102,126,234,.18); stroke: rgba(102,126,234,.35); stroke-width: 1; }
    .detail-block.reverse { fill: rgba(231,126,34,.16); stroke: rgba(231,126,34,.35); }
    .detail-match { fill: #4f86d9; stroke: #2f5f9f; stroke-width: 1; }
    .detail-ref-gene { fill: #667eea; stroke: #5364c9; stroke-width: 1; }
    .detail-variant.snp { fill: #f39c12; stroke: #b76d00; }
    .detail-variant.indel { fill: #2d9cdb; stroke: #176a95; }
    .detail-empty {
      padding: 22px;
      border: 1px dashed #d9dce8;
      border-radius: 8px;
      color: #777;
      background: #fff;
    }
    .empty { color: #777; padding: 18px; }
    .copy-toast {
      position: fixed;
      bottom: 20px;
      right: 20px;
      background: #333;
      color: white;
      padding: 10px 20px;
      border-radius: 4px;
      display: none;
      z-index: 1000;
    }
    @media (max-width: 900px) {
      body { padding: 10px; }
      .header-content, .detail-layout { display: block; }
      .section { padding: 20px 16px; }
      .overview-canvas { grid-template-columns: 120px minmax(650px, 1fr); }
      .overview-track { grid-template-columns: 120px minmax(420px, 1fr); }
      .controls { margin-bottom: 14px; }
    }
  </style>
</head>
<body>
<div class="container">
  <div class="header">
    <div class="header-content">
      <div>
        <h1>__TITLE__</h1>
        <div class="subtitle">
          <span class="mode-badge" id="modeBadge">GeneScreen</span>
          <span>1 ref + N query alignment report</span>
        </div>
        <div class="gen-time" id="generatedAt"></div>
      </div>
      <a class="json-badge" href="data.json">data.json</a>
    </div>
  </div>

  <div class="section input-section">
    <h2>初始化</h2>
    <table class="info-table" id="initTable"></table>
  </div>

  <div class="section">
    <h2>输入文件</h2>
    <table class="file-table" id="genomeTable"></table>
  </div>

  <div class="section stats-section">
    <h2>结果统计</h2>
    <table class="stats-table">
      <thead><tr id="statsHead"></tr></thead>
      <tbody><tr id="statsBody"></tr></tbody>
    </table>
  </div>

  <div class="section overview-section visualization-section">
    <h2>宏观比对图</h2>
    <div class="viz-container">
      <div id="overview" class="overview-canvas"></div>
      <div class="viz-legend-note">
        <p>绿色块表示每个查询基因组的最优候选；点击候选块可切换下方局部组合。数据来自 <code>data.json</code>，切换显示不会重新计算。</p>
      </div>
    </div>
  </div>

  <div class="section visualization-section">
    <h2>局部比对图</h2>
    <div class="viz-container">
      <div class="detail-layout">
        <div id="controls" class="controls"></div>
        <div id="detail" class="track-stack"></div>
      </div>
      <div class="viz-legend-note">
        <p>左侧下拉框切换候选组合；上移/下移按钮或拖动轨道可调整显示顺序。</p>
      </div>
    </div>
  </div>

  <div class="section output-section">
    <h2>结果文件</h2>
    <table class="file-table" id="outputTable"></table>
  </div>
</div>
<div id="copyToast" class="copy-toast">已复制到剪贴板</div>

  <script id="embedded-report-data" type="application/json">__EMBEDDED_DATA__</script>
  <script>
    let reportData = JSON.parse(document.getElementById('embedded-report-data').textContent);
    fetch('data.json').then(r => r.ok ? r.json() : reportData).then(data => { reportData = data; init(); }).catch(init);

    const state = { selected: {}, order: [] };
    function esc(value) {
      return String(value ?? '').replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
    }
    function init() {
      state.selected = Object.assign({}, reportData.default_selection.selected_candidates || {});
      state.order = (reportData.default_selection.track_order || []).slice();
      renderHeader();
      renderInit();
      renderGenomes();
      renderStats();
      renderOverview();
      renderControls();
      renderDetail();
      renderOutputs();
    }
    function queryIds() { return (reportData.overview && reportData.overview.query_order) || []; }
    function trackById(trackId) {
      return (reportData.tracks || []).find(t => t.track_id === trackId) || { name: trackId, role: 'query_genome' };
    }
    function candidateById(queryId, candidateId) {
      return (reportData.candidates[queryId] || []).find(c => c.candidate_id === candidateId);
    }
    function numberOr(value, fallback) {
      const n = Number(value);
      return Number.isFinite(n) ? n : fallback;
    }
    function normRange(a, b) {
      const x = numberOr(a, 0);
      const y = numberOr(b, x);
      return [Math.min(x, y), Math.max(x, y)];
    }
    function candidateBounds(candidate) {
      if (!candidate) return { start: 1, end: 2 };
      const points = [];
      if (candidate.query_region_start != null) points.push(Number(candidate.query_region_start));
      if (candidate.query_region_end != null) points.push(Number(candidate.query_region_end));
      if (candidate.query_start != null) points.push(Number(candidate.query_start));
      if (candidate.query_end != null) points.push(Number(candidate.query_end));
      for (const block of candidate.blocks || []) {
        points.push(Number(block.query_start), Number(block.query_end));
      }
      const clean = points.filter(Number.isFinite);
      if (!clean.length) return { start: 1, end: 2 };
      const start = Math.min(...clean);
      const end = Math.max(...clean);
      return end > start ? { start, end } : { start, end: start + 1 };
    }
    function refLength() {
      const input = reportData.input || {};
      if (Number(input.sequence_length) > 0) return Number(input.sequence_length);
      let maxEnd = 1;
      for (const qid of queryIds()) {
        for (const candidate of reportData.candidates[qid] || []) {
          for (const block of candidate.blocks || []) {
            maxEnd = Math.max(maxEnd, Number(block.ref_start) || 1, Number(block.ref_end) || 1);
          }
        }
      }
      return maxEnd;
    }
    function xFor(pos, range, left, width) {
      const span = Math.max(1, range.end - range.start);
      const p = Math.max(range.start, Math.min(range.end, Number(pos) || range.start));
      return left + ((p - range.start) / span) * width;
    }
    function selectedPairForQuery(qid) {
      return (reportData.pairs || []).find(p => p.query_genome_id === qid || p.pair_id === `ref__${qid}`);
    }
    function variantClass(type) {
      return String(type || '').toUpperCase() === 'SNP' ? 'snp' : 'indel';
    }
    function renderHeader() {
      document.getElementById('modeBadge').textContent = reportData.mode || 'GeneScreen';
      document.getElementById('generatedAt').textContent = reportData.generated_at ? `生成时间：${reportData.generated_at}` : '';
    }
    function renderInit() {
      const input = reportData.input || {};
      const params = reportData.parameters || {};
      const rows = [
        ['输入 ID', `<code>${esc(input.id || '-')}</code>`],
        ['模式', esc(reportData.mode || '-')],
        ['序列长度', input.sequence_length ? `${Number(input.sequence_length).toLocaleString()} bp` : '-'],
        ['Identity 阈值', params.identity != null ? `${params.identity}%` : '-'],
        ['最小比对长度', params.min_aln_len != null ? `${params.min_aln_len} bp` : '-'],
        ['Pairwise Top-N', params.pairwise_all ? '全量候选' : (params.candidate_limit ?? '-')],
        ['查询上游延伸', `${params.query_upstream || 0} bp`],
        ['查询下游延伸', `${params.query_downstream || 0} bp`]
      ];
      document.getElementById('initTable').innerHTML = rows.map(([k, v]) => `<tr><td class="label">${k}</td><td class="value">${v}</td></tr>`).join('');
    }
    function renderGenomes() {
      const ref = (reportData.genomes && reportData.genomes.ref) || {};
      const queries = (reportData.genomes && reportData.genomes.queries) || [];
      const rows = [
        '<tr><th>类型</th><th>名称</th><th>FASTA</th><th>注释</th></tr>',
        `<tr><td>参考</td><td>${esc(ref.name || '-')}</td><td><code class="path-code">${esc(ref.source_fasta || ref.artifact_fasta || '-')}</code></td><td><code class="path-code">${esc(ref.source_gff || ref.artifact_gff || '-')}</code></td></tr>`
      ];
      for (const q of queries) {
        rows.push(`<tr><td>查询</td><td>${esc(q.name || q.genome_id)}</td><td><code class="path-code">${esc(q.source_fasta || '-')}</code></td><td>${q.has_annotation ? `<code class="path-code">${esc(q.source_gff || '-')}</code>` : '<span class="file-na">未使用注释</span>'}</td></tr>`);
      }
      document.getElementById('genomeTable').innerHTML = rows.join('');
    }
    function renderStats() {
      const s = reportData.statistics || {};
      const variants = (s.selected_combination && s.selected_combination.variants) || {};
      const indelTotal = variants.INDEL || ((variants.INS || 0) + (variants.DEL || 0));
      const metrics = [
        ['基因组', s.genome_count || 0],
        ['查询基因组', s.query_count || 0],
        ['候选片段', s.total_candidate_count || 0],
        ['可见连接', (s.selected_combination && s.selected_combination.link_count) || 0],
        ['SNP 数量', variants.SNP || 0],
        ['Indel 数量', indelTotal],
        ['INS', variants.INS || 0],
        ['DEL', variants.DEL || 0]
      ];
      document.getElementById('statsHead').innerHTML = metrics.map(m => `<th>${m[0]}</th>`).join('');
      document.getElementById('statsBody').innerHTML = metrics.map(m => `<td>${m[1]}</td>`).join('');
    }
    function renderOverview() {
      const root = document.getElementById('overview');
      const input = reportData.input || {};
      const chrMap = new Map();
      for (const qid of queryIds()) {
        for (const candidate of reportData.candidates[qid] || []) {
          const chr = candidate.query_chr || 'unknown';
          if (!chrMap.has(chr)) chrMap.set(chr, []);
          chrMap.get(chr).push({ qid, candidate, bounds: candidateBounds(candidate) });
        }
      }
      if (!chrMap.size) { root.innerHTML = '<div class="empty">No candidates</div>'; return; }
      const refTitle = esc((reportData.genomes && reportData.genomes.ref && reportData.genomes.ref.name) || 'ref-gene');
      const refId = esc(input.id || input.safe_id || 'A1');
      const body = document.createElement('div');
      body.className = 'overview-body';
      root.innerHTML = `<div class="overview-ref"><div class="overview-ref-title">${refTitle}</div><div class="overview-ref-id">${refId}</div></div>`;
      root.appendChild(body);
      const chrNames = Array.from(chrMap.keys()).sort((a, b) => String(a).localeCompare(String(b), undefined, { numeric: true }));
      for (const chr of chrNames) {
        const entries = chrMap.get(chr);
        const rangeStart = Math.min(...entries.map(e => e.bounds.start));
        const rangeEnd = Math.max(...entries.map(e => e.bounds.end), rangeStart + 1);
        const group = document.createElement('div');
        group.className = 'chr-group';
        group.innerHTML = `<div class="chr-label">chr ${esc(chr)}</div><div class="chr-tracks"></div>`;
        const tracks = group.querySelector('.chr-tracks');
        for (const qid of queryIds()) {
          const candidates = entries.filter(e => e.qid === qid).map(e => e.candidate);
          if (!candidates.length) continue;
          const track = trackById(qid);
          const row = document.createElement('div');
          row.className = 'overview-track';
          row.innerHTML = `<div class="track-label" title="${esc(track.name)}">${esc(track.name)}<div class="track-subtext">${candidates.length} candidates</div></div><div class="lane"></div>`;
          const lane = row.querySelector('.lane');
          for (const candidate of candidates) {
            const bounds = candidateBounds(candidate);
            const left = ((bounds.start - rangeStart) / Math.max(1, rangeEnd - rangeStart)) * 100;
            const width = Math.max(0.8, ((bounds.end - bounds.start) / Math.max(1, rangeEnd - rangeStart)) * 100);
            const tick = document.createElement('div');
            tick.className = 'lane-tick';
            tick.style.left = `${Math.max(0, Math.min(100, left + width / 2))}%`;
            lane.appendChild(tick);
            const el = document.createElement('div');
            el.className = 'candidate' + (candidate.is_best ? ' best' : '') + (state.selected[qid] === candidate.candidate_id ? ' selected' : '');
            const region = candidate.query_region_start ? `${candidate.query_chr}:${candidate.query_region_start}-${candidate.query_region_end}` : `${candidate.query_chr}:${candidate.query_start}-${candidate.query_end}`;
            el.style.left = `${Math.max(0, Math.min(99, left))}%`;
            el.style.width = `${Math.max(0.8, Math.min(width, 100 - left))}%`;
            el.title = `${candidate.candidate_id} ${region} identity=${candidate.identity || '-'} coverage=${candidate.coverage || '-'}`;
            el.onclick = () => { state.selected[qid] = candidate.candidate_id; renderOverview(); renderControls(); renderDetail(); };
            lane.appendChild(el);
          }
          tracks.appendChild(row);
        }
        body.appendChild(group);
      }
    }
    function renderControls() {
      const root = document.getElementById('controls');
      root.innerHTML = '<h3>候选切换</h3>';
      for (const qid of queryIds()) {
        const track = trackById(qid);
        const row = document.createElement('div');
        row.className = 'control-row';
        row.draggable = true;
        row.dataset.trackId = qid;
        const options = (reportData.candidates[qid] || []).map(c => `<option value="${esc(c.candidate_id)}" ${state.selected[qid] === c.candidate_id ? 'selected' : ''}>${esc(track.name)} · ${esc(c.candidate_id)}</option>`).join('');
        row.innerHTML = `<select>${options}</select><button title="上移">↑</button><button title="下移">↓</button>`;
        row.querySelector('select').onchange = e => { state.selected[qid] = e.target.value; renderOverview(); renderDetail(); };
        row.children[1].onclick = () => moveTrack(qid, -1);
        row.children[2].onclick = () => moveTrack(qid, 1);
        wireControlDrag(row);
        root.appendChild(row);
      }
    }
    function moveTrack(trackId, delta) {
      const index = state.order.indexOf(trackId);
      const target = index + delta;
      if (index <= 0 || target <= 0 || target >= state.order.length) return;
      state.order.splice(index, 1);
      state.order.splice(target, 0, trackId);
      renderDetail();
    }
    function wireControlDrag(row) {
      row.addEventListener('dragstart', e => {
        row.classList.add('dragging');
        e.dataTransfer.setData('text/plain', row.dataset.trackId);
      });
      row.addEventListener('dragend', () => row.classList.remove('dragging'));
      row.addEventListener('dragover', e => e.preventDefault());
      row.addEventListener('drop', e => {
        e.preventDefault();
        const dragged = e.dataTransfer.getData('text/plain');
        const target = row.dataset.trackId;
        if (!dragged || dragged === target) return;
        state.order = state.order.filter(id => id !== dragged);
        const index = state.order.indexOf(target);
        state.order.splice(index < 0 ? state.order.length : index, 0, dragged);
        renderControls();
        renderDetail();
      });
    }
    function renderDetail() {
      const root = document.getElementById('detail');
      const order = state.order.length ? state.order.slice() : ['ref'].concat(queryIds());
      const rows = order.map(trackId => {
        const track = trackById(trackId);
        const candidate = trackId === 'ref' ? null : candidateById(trackId, state.selected[trackId]);
        const range = trackId === 'ref' ? { start: 1, end: refLength() } : candidateBounds(candidate);
        return { trackId, track, candidate, range };
      }).filter(row => row.trackId === 'ref' || row.candidate);
      if (!rows.length) { root.innerHTML = '<div class="detail-empty">No tracks</div>'; return; }

      const width = 1040;
      const labelWidth = 250;
      const plotLeft = 285;
      const plotWidth = 690;
      const rowStep = 78;
      const top = 44;
      const height = top + rows.length * rowStep + 34;
      const rowY = row => top + rows.indexOf(row) * rowStep;
      const refRow = rows.find(r => r.trackId === 'ref');
      const fragments = [];
      const ribbons = [];
      const markers = [];

      for (const row of rows) {
        const y = rowY(row);
        const subtitle = row.trackId === 'ref'
          ? `reference · ${row.range.start}-${row.range.end}`
          : `${row.candidate.candidate_id} · ${row.candidate.query_chr}:${row.range.start}-${row.range.end}`;
        fragments.push(`<rect class="detail-track-bg" x="0" y="${y - 30}" width="${width}" height="${rowStep - 8}" rx="8"></rect>`);
        fragments.push(`<text class="detail-label" x="18" y="${y - 5}">${esc(row.track.name || row.trackId)}</text>`);
        fragments.push(`<text class="detail-subtext" x="18" y="${y + 14}">${esc(subtitle)}</text>`);
        fragments.push(`<line class="detail-axis" x1="${plotLeft}" y1="${y}" x2="${plotLeft + plotWidth}" y2="${y}"></line>`);
        fragments.push(`<text class="detail-subtext" x="${plotLeft}" y="${y + 24}">${Number(row.range.start).toLocaleString()}</text>`);
        fragments.push(`<text class="detail-subtext" x="${plotLeft + plotWidth}" y="${y + 24}" text-anchor="end">${Number(row.range.end).toLocaleString()}</text>`);

        if (row.trackId === 'ref') {
          const input = reportData.input || {};
          const info = input.extraction_info || {};
          const geneStart = Number(info.gene_rel_start || 1);
          const geneEnd = Number(info.gene_rel_end || input.sequence_length || row.range.end);
          const x1 = xFor(geneStart, row.range, plotLeft, plotWidth);
          const x2 = xFor(geneEnd, row.range, plotLeft, plotWidth);
          fragments.push(`<rect class="detail-ref-gene" x="${Math.min(x1, x2)}" y="${y - 7}" width="${Math.max(3, Math.abs(x2 - x1))}" height="14" rx="7"></rect>`);
        } else {
          for (const block of row.candidate.blocks || []) {
            const [qs, qe] = normRange(block.query_start, block.query_end);
            const x1 = xFor(qs, row.range, plotLeft, plotWidth);
            const x2 = xFor(qe, row.range, plotLeft, plotWidth);
            fragments.push(`<rect class="detail-match" x="${Math.min(x1, x2)}" y="${y - 6}" width="${Math.max(3, Math.abs(x2 - x1))}" height="12" rx="6"></rect>`);
            if (refRow) {
              const refY = rowY(refRow);
              const rx1 = xFor(block.ref_start, refRow.range, plotLeft, plotWidth);
              const rx2 = xFor(block.ref_end, refRow.range, plotLeft, plotWidth);
              const qx1 = xFor(block.query_start, row.range, plotLeft, plotWidth);
              const qx2 = xFor(block.query_end, row.range, plotLeft, plotWidth);
              const reverse = (Number(block.ref_start) - Number(block.ref_end)) * (Number(block.query_start) - Number(block.query_end)) < 0;
              ribbons.push(`<polygon class="detail-block${reverse ? ' reverse' : ''}" points="${rx1},${refY + 9} ${rx2},${refY + 9} ${qx2},${y - 9} ${qx1},${y - 9}"></polygon>`);
            }
          }
        }
      }

      for (const row of rows) {
        if (row.trackId === 'ref' || !row.candidate || !refRow) continue;
        const pair = selectedPairForQuery(row.trackId);
        const pairId = pair && pair.pair_id;
        const y = rowY(row);
        const refY = rowY(refRow);
        for (const variant of reportData.variants || []) {
          if (pairId && variant.pair_id !== pairId) continue;
          if (variant.query_chr && row.candidate.query_chr && String(variant.query_chr) !== String(row.candidate.query_chr)) continue;
          const vClass = variantClass(variant.type);
          const qPos = Number(variant.query_genome_pos);
          const rPos = Number(variant.ref_source_pos);
          if (Number.isFinite(qPos) && qPos >= row.range.start && qPos <= row.range.end) {
            const x = xFor(qPos, row.range, plotLeft, plotWidth);
            markers.push(`<circle class="detail-variant ${vClass}" cx="${x}" cy="${y - 15}" r="4"><title>${esc(variant.type)} ${esc(row.track.name)}:${qPos}</title></circle>`);
          }
          if (Number.isFinite(rPos) && rPos >= refRow.range.start && rPos <= refRow.range.end) {
            const x = xFor(rPos, refRow.range, plotLeft, plotWidth);
            markers.push(`<circle class="detail-variant ${vClass}" cx="${x}" cy="${refY - 15}" r="4"><title>${esc(variant.type)} ref:${rPos}</title></circle>`);
          }
        }
      }

      root.innerHTML = `<div class="detail-canvas"><svg class="detail-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="multi-track local alignment">${ribbons.join('')}${fragments.join('')}${markers.join('')}</svg></div>`;
    }
    function renderOutputs() {
      const rows = ['<tr><th>文件名</th><th>描述</th><th>路径</th></tr>'];
      rows.push('<tr><td><a href="data.json">data.json</a></td><td>新版报告结构化数据</td><td><code class="path-code">report/data.json</code></td></tr>');
      for (const f of reportData.extra_files || []) {
        rows.push(`<tr><td>${esc(f.label || '-')}</td><td>${esc(f.description || '')}</td><td><code class="path-code">${esc(f.path || '')}</code></td></tr>`);
      }
      document.getElementById('outputTable').innerHTML = rows.join('');
    }
  </script>
</body>
</html>
"""
    return template.replace("__TITLE__", title_html).replace("__EMBEDDED_DATA__", embedded)
