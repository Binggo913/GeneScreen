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
import sys
import traceback
import types
from hashlib import md5
from itertools import permutations, product
from datetime import datetime
from html import escape, unescape
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
        return os.path.relpath(path, start).replace(os.sep, "/")
    except ValueError:
        return str(path).replace(os.sep, "/")


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


VALID_GFF_TYPES = {
    "gene", "mrna", "transcript", "exon", "cds",
    "five_prime_utr", "three_prime_utr", "utr", "5utr", "3utr",
}


def _selection_key(query_order: List[str], selected: Dict[str, str]) -> str:
    return "||".join(f"{qid}={selected.get(qid, '')}" for qid in query_order)


def _order_key(order: List[str]) -> str:
    return "||".join(order)


def _chrom_key(value: Any) -> str:
    text = str(value or "").split()[0]
    text = re.sub(r"^(chr|Chr|CHR)", "", text)
    return text.lstrip("0").lower() or text.lower()


def _chrom_matches(left: Any, right: Any) -> bool:
    return str(left or "") == str(right or "") or _chrom_key(left) == _chrom_key(right)


def _read_fai_lengths(fasta_file: Optional[str]) -> Dict[str, int]:
    if not fasta_file:
        return {}
    candidates = [fasta_file]
    if not str(fasta_file).endswith(".fai"):
        candidates.append(f"{fasta_file}.fai")
    for path in candidates:
        if not path or not os.path.exists(path):
            continue
        lengths: Dict[str, int] = {}
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                try:
                    lengths[parts[0]] = int(parts[1])
                except ValueError:
                    continue
        if lengths:
            return lengths
    return {}


def _lookup_length(lengths: Dict[str, int], name: Any, fallback: int) -> int:
    if not lengths:
        return max(1, int(fallback or 1))
    text = str(name or "")
    candidates = [text, text.split()[0] if text else text]
    for candidate in candidates:
        if candidate in lengths:
            return lengths[candidate]
    target_key = _chrom_key(text)
    for chrom, length in lengths.items():
        if _chrom_key(chrom) == target_key:
            return length
    return max(1, int(fallback or 1))


def _candidate_region(candidate: Dict[str, Any]) -> Tuple[int, int]:
    starts = [
        candidate.get("query_region_start"),
        candidate.get("query_start"),
        candidate.get("query_match_start"),
    ]
    ends = [
        candidate.get("query_region_end"),
        candidate.get("query_end"),
        candidate.get("query_match_end"),
    ]
    points = []
    for value in starts + ends:
        try:
            if value is not None:
                points.append(int(value))
        except (TypeError, ValueError):
            continue
    for block in candidate.get("blocks") or []:
        for key in ("query_start", "query_end"):
            try:
                points.append(int(block.get(key)))
            except (TypeError, ValueError):
                continue
    if not points:
        return 1, 2
    left, right = min(points), max(points)
    return max(1, left), max(max(1, left) + 1, right)


def _candidate_fasta_region(fasta_file: Optional[str]) -> Optional[Tuple[str, int, int]]:
    if not fasta_file or not os.path.exists(fasta_file):
        return None
    try:
        with open(fasta_file, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                match = re.search(r"\bregion=([^:\s]+):([0-9]+)-([0-9]+)", line)
                if not match:
                    return None
                return match.group(1), int(match.group(2)), int(match.group(3))
    except OSError:
        return None
    return None


def _local_to_region(region: Optional[Tuple[str, int, int]], value: Any) -> int:
    try:
        pos = int(value)
    except (TypeError, ValueError):
        return 1
    if not region:
        return max(1, pos)
    return max(1, int(region[1]) + pos - 1)


def _parse_location_text(location: Optional[str]) -> Optional[Tuple[str, int, int]]:
    if not location:
        return None
    match = re.match(r"^([^:]+):([0-9,]+)-([0-9,]+)$", str(location).strip())
    if not match:
        return None
    try:
        return (
            match.group(1),
            int(match.group(2).replace(",", "")),
            int(match.group(3).replace(",", "")),
        )
    except ValueError:
        return None


def _write_gff_records(
    source_gff: Optional[str],
    handle,
    target_seqid: str,
    source_chr: Optional[str] = None,
    region_start: Optional[int] = None,
    region_end: Optional[int] = None,
    relative: bool = False,
) -> int:
    if not source_gff or not os.path.exists(source_gff):
        return 0
    written = 0
    with open(source_gff, "r", encoding="utf-8", errors="ignore") as source:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat_type = parts[2].lower()
            if feat_type not in VALID_GFF_TYPES:
                continue
            try:
                feat_start = int(parts[3])
                feat_end = int(parts[4])
            except ValueError:
                continue
            if source_chr and not _chrom_matches(parts[0], source_chr):
                continue
            if region_start is not None and region_end is not None:
                if feat_end < region_start or feat_start > region_end:
                    continue
                if relative:
                    parts[3] = str(max(1, feat_start - region_start + 1))
                    parts[4] = str(max(1, feat_end - region_start + 1))
            parts[0] = target_seqid
            handle.write("\t".join(parts) + "\n")
            written += 1
    return written


def _run_linkview_for_report(
    input_file: str,
    output_prefix: str,
    k_file: str,
    hl_file: Optional[str],
    gff_file: Optional[str],
    min_identity: float,
    min_aln_len: int,
) -> bool:
    try:
        import argparse
        try:
            from . import LINKVIEW
        except ImportError:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            if module_dir not in sys.path:
                sys.path.insert(0, module_dir)
            if "cairosvg" not in sys.modules:
                cairosvg_stub = types.ModuleType("cairosvg")
                cairosvg_stub.svg2png = lambda *args, **kwargs: None
                sys.modules["cairosvg"] = cairosvg_stub
            import LINKVIEW

        args = argparse.Namespace(
            input=input_file,
            type=2,
            output=output_prefix,
            karyotype=k_file,
            highlight=hl_file if hl_file and os.path.exists(hl_file) else None,
            gff=gff_file if gff_file and os.path.exists(gff_file) else None,
            min_identity=min_identity,
            min_alignment_length=min_aln_len,
            svg_height=400,
            svg_width=1200,
            svg_space=0.2,
            chro_thickness=15,
            label_font_size=18,
            label_angle=0,
            chro_axis=True,
            chro_axis_density=2,
            show_pos_with_label=False,
            bezier=True,
            style="simple",
            svg2png="",
            svg2png_dpi=350,
            no_label=False,
            no_dash=False,
            no_scale=False,
            hl_min1px=False,
            scale=None,
            gap_length=0.2,
            chro_len=None,
            parameter=None,
            max_evalue=1e-5,
            min_bit_score=5000,
        )
        LINKVIEW.args = args
        LINKVIEW.main(args)
        return os.path.exists(f"{output_prefix}.svg")
    except Exception as exc:
        print(f"[WARNING] 新版报告 LINKVIEW 局部图生成失败: {exc}")
        traceback.print_exc()
        return False


def _legacy_legend_svg(svg_width: int, y_top: int) -> str:
    items = [
        ("SNP", "orange", "line", "snp"),
        ("Indel", "blue", "line", "indel"),
        ("5' UTR", "#6B5B7B", "rect", "utr5"),
        ("3' UTR", "#B6AEC9", "rect", "utr3"),
        ("CDS", "#7A7A7A", "rect", "cds"),
    ]
    item_height = 20
    padding_x = 10
    padding_y = 8
    box_width = 120
    box_height = padding_y * 2 + len(items) * item_height
    x = max(20, svg_width - box_width - 30)
    y = y_top
    parts = [
        '<g class="genescreen-legend">',
        f'<rect x="{x}" y="{y}" width="{box_width}" height="{box_height}" fill="white" stroke="#ddd" stroke-width="1" rx="4" opacity="0.95"/>',
    ]
    for index, (label, color, kind, data_type) in enumerate(items):
        yy = y + padding_y + 10 + index * item_height
        icon_x = x + padding_x
        parts.append(
            f'<g class="legend-item" data-type="{data_type}" style="cursor: pointer;">'
            f'<rect x="{x + 2}" y="{yy - 9}" width="{box_width - 4}" height="{item_height - 2}" fill="transparent" class="legend-hitarea"/>'
        )
        if kind == "line":
            parts.append(f'<line class="legend-icon" x1="{icon_x + 6}" y1="{yy - 6}" x2="{icon_x + 6}" y2="{yy + 6}" stroke="{color}" stroke-width="3"/>')
        else:
            parts.append(f'<rect class="legend-icon" x="{icon_x}" y="{yy - 5}" width="12" height="10" rx="2" fill="{color}"/>')
        parts.append(f'<text x="{icon_x + 20}" y="{yy + 4}" font-family="Arial, sans-serif" font-size="10" fill="#666" class="legend-text">{label}</text></g>')
    parts.append("</g>")
    return "".join(parts)


def _process_linkview_svg_file(svg_path: str) -> None:
    try:
        with open(svg_path, "r", encoding="utf-8", errors="ignore") as handle:
            svg = handle.read()
        width_match = re.search(r'width="(\d+)"', svg)
        height_match = re.search(r'height="(\d+)"', svg)
        if not width_match or not height_match:
            return
        width = int(width_match.group(1))
        height = int(height_match.group(1))
        top_margin = 90
        bottom_margin = 70
        svg = re.sub(
            r'<svg\s+width="[^"]*"\s+height="[^"]*"',
            f'<svg width="{width}" height="{height + top_margin + bottom_margin}" viewBox="0 -{top_margin} {width} {height + top_margin + bottom_margin}" preserveAspectRatio="xMidYMid meet"',
            svg,
            count=1,
        )
        if "genescreen-legend" not in svg:
            svg = svg.replace("</svg>", f"{_legacy_legend_svg(width, -top_margin + 35)}</svg>")
        with open(svg_path, "w", encoding="utf-8") as handle:
            handle.write(svg)
    except Exception as exc:
        print(f"[WARNING] 处理 LINKVIEW SVG 失败: {exc}")


def _normalize_track_token(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def _svg_text_records(svg: str) -> List[Tuple[float, str]]:
    records: List[Tuple[float, str]] = []
    for match in re.finditer(r"<text\b([^>]*)>(.*?)</text>", svg, flags=re.IGNORECASE | re.DOTALL):
        attrs = match.group(1) or ""
        body = re.sub(r"<[^>]+>", "", match.group(2) or "")
        text = unescape(body).strip()
        if not text:
            continue
        y_match = re.search(r'\by="([\-0-9.]+)"', attrs)
        try:
            y = float(y_match.group(1)) if y_match else 0.0
        except ValueError:
            y = 0.0
        records.append((y, text))
    return records


def _parse_linkview_svg_order(
    svg_path: str,
    expected_order: List[str],
    tracks: List[Dict[str, Any]],
    track_aliases: Dict[str, str],
) -> List[str]:
    try:
        with open(svg_path, "r", encoding="utf-8", errors="ignore") as handle:
            svg = handle.read()
    except OSError:
        return expected_order

    track_by_id = {str(track.get("track_id")): track for track in tracks}
    candidates = expected_order + [
        str(track.get("track_id")) for track in tracks
        if str(track.get("track_id")) not in expected_order
    ]
    tokens_by_track: Dict[str, List[str]] = {}
    for track_id in candidates:
        track = track_by_id.get(track_id) or {}
        tokens = [
            track_id,
            track_aliases.get(track_id),
            track.get("name"),
            track.get("sequence_id"),
            track.get("genome_id"),
        ]
        tokens_by_track[track_id] = [
            token for token in (_normalize_track_token(item) for item in tokens)
            if token
        ]

    matched: List[Tuple[float, str]] = []
    used = set()
    for y, text in _svg_text_records(svg):
        normalized = _normalize_track_token(text)
        if not normalized:
            continue
        best_track = None
        best_score = 0
        tied = False
        for track_id in candidates:
            if track_id in used:
                continue
            score = 0
            for token in tokens_by_track.get(track_id, []):
                if normalized == token:
                    score = max(score, 5)
                elif len(token) >= 4 and token in normalized:
                    score = max(score, 3)
                elif len(normalized) >= 4 and normalized in token:
                    score = max(score, 2)
            if score > best_score:
                best_track = track_id
                best_score = score
                tied = False
            elif score and score == best_score:
                tied = True
        if best_track and not tied:
            matched.append((y, best_track))
            used.add(best_track)

    inferred = [track_id for _, track_id in sorted(matched, key=lambda item: item[0])]
    return inferred if len(inferred) == len(expected_order) else expected_order


def _write_linkview_relation(
    handle,
    first_start: int,
    first_end: int,
    second_start: int,
    second_end: int,
    first_len: int,
    second_len: int,
    identity: float,
    first_alias: str,
    second_alias: str,
) -> None:
    first_aln_len = max(1, abs(int(first_end) - int(first_start)) + 1)
    second_aln_len = max(1, abs(int(second_end) - int(second_start)) + 1)
    cov_first = first_aln_len / max(1, int(first_len or first_aln_len)) * 100
    cov_second = second_aln_len / max(1, int(second_len or second_aln_len)) * 100
    handle.write(
        f"{int(first_start):>8} {int(first_end):>8}  | {int(second_start):>8} {int(second_end):>8}  | "
        f"{first_aln_len:>8} {second_aln_len:>8}  | {float(identity):>8.2f}  | "
        f"{int(first_len or first_aln_len):>8} {int(second_len or second_aln_len):>8}  | "
        f"{cov_first:>8.2f} {cov_second:>8.2f}  | {first_alias}\t{second_alias}\n"
    )


def _write_hl_marker(handle, alias: str, pos: int, color: str) -> None:
    pos = max(1, int(pos))
    handle.write(f"{alias}\t{pos - 1}\t{pos}\t{color}\n")


def _append_marker(
    markers: List[Dict[str, Any]],
    track_id: str,
    alias: str,
    pos: int,
    color: str,
    variant: Dict[str, Any],
    pair_id: str,
    marker_role: str,
) -> None:
    var_type = str(variant.get("type", "")).upper()
    markers.append({
        "track_id": track_id,
        "alias": alias,
        "position": max(1, int(pos)),
        "color": color,
        "type": var_type,
        "pair_id": pair_id,
        "variant_id": variant.get("variant_id") or f"{pair_id}:{var_type}:{pos}:{marker_role}",
        "marker_role": marker_role,
        "ref": variant.get("ref"),
        "alt": variant.get("alt"),
        "ref_source_pos": variant.get("ref_source_pos"),
        "query_genome_pos": variant.get("query_genome_pos"),
    })


def _pairwise_lookup(result: Dict[str, Any]) -> Dict[Tuple[str, str, str, str], List[Dict[str, Any]]]:
    lookup: Dict[Tuple[str, str, str, str], List[Dict[str, Any]]] = {}
    for pairwise_result in result.get("pairwise_results") or []:
        left_id = sanitize_path_segment(pairwise_result.get("left_query") or "left", "left")
        right_id = sanitize_path_segment(pairwise_result.get("right_query") or "right", "right")
        left_candidate_id = str(pairwise_result.get("left_candidate_id") or "")
        right_candidate_id = str(pairwise_result.get("right_candidate_id") or "")
        lookup.setdefault((left_id, left_candidate_id, right_id, right_candidate_id), []).append(pairwise_result)
    return lookup


def _pairwise_for_selection(
    lookup: Dict[Tuple[str, str, str, str], List[Dict[str, Any]]],
    first_id: str,
    first_candidate_id: str,
    second_id: str,
    second_candidate_id: str,
) -> Tuple[Optional[Dict[str, Any]], bool]:
    direct = lookup.get((first_id, first_candidate_id, second_id, second_candidate_id))
    if direct:
        return direct[0], False
    reverse = lookup.get((second_id, second_candidate_id, first_id, first_candidate_id))
    if reverse:
        return reverse[0], True
    return None, False


def _read_pairwise_blocks(
    pairwise_result: Dict[str, Any],
    left_region: Optional[Tuple[str, int, int]],
    right_region: Optional[Tuple[str, int, int]],
) -> List[Dict[str, Any]]:
    coords_file = pairwise_result.get("coords")
    if not coords_file or not os.path.exists(coords_file):
        return []
    blocks = []
    with open(coords_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("[") or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            try:
                right_start = _local_to_region(right_region, parts[0])
                right_end = _local_to_region(right_region, parts[1])
                left_start = _local_to_region(left_region, parts[2])
                left_end = _local_to_region(left_region, parts[3])
                identity = float(parts[6])
            except ValueError:
                continue
            blocks.append({
                "left_start": left_start,
                "left_end": left_end,
                "right_start": right_start,
                "right_end": right_end,
                "identity": identity,
            })
    return blocks


def _generate_linkview_detail_assets(
    result: Dict[str, Any],
    output_dir: str,
    report_dir: str,
    mode: str,
    ref_length: int,
    genome_files: Dict[str, Any],
    query_order: List[str],
    tracks: List[Dict[str, Any]],
    candidates_map: Dict[str, List[Dict[str, Any]]],
    variants_payload: List[Dict[str, Any]],
    identity: float,
    min_aln_len: int,
) -> Dict[str, Any]:
    if not query_order or not all(candidates_map.get(qid) for qid in query_order):
        return {}

    asset_dir = ensure_dir(os.path.join(report_dir, "linkview_detail"))
    input_id = sanitize_path_segment(result.get("id") or "input", "input")
    ref_alias = input_id
    track_aliases = {"ref": ref_alias}
    used_aliases = {ref_alias}
    track_by_id = {track.get("track_id"): track for track in tracks}
    for qid in query_order:
        base = sanitize_path_segment(track_by_id.get(qid, {}).get("name") or qid, qid)
        alias = base
        suffix = 2
        while alias in used_aliases:
            alias = f"{base}_{suffix}"
            suffix += 1
        track_aliases[qid] = alias
        used_aliases.add(alias)

    query_entries = result.get("queries") or genome_files.get("queries") or []
    query_entry_by_id: Dict[str, Dict[str, Any]] = {}
    for index, qid in enumerate(query_order):
        if index < len(query_entries):
            query_entry_by_id[qid] = query_entries[index]
        else:
            query_entry_by_id[qid] = {}

    ref_gff_source = result.get("gff")
    ref_region = None
    if not ref_gff_source and mode == "location":
        ref_region = _parse_location_text(result.get("location"))
        ref_gff_source = genome_files.get("ref_gff")

    query_lengths: Dict[str, Dict[str, int]] = {}
    for qid in query_order:
        entry = query_entry_by_id.get(qid) or {}
        query_lengths[qid] = _read_fai_lengths(entry.get("fasta") or entry.get("source_fasta"))
    pairwise_lookup = _pairwise_lookup(result)

    default_order = ["ref"] + query_order
    all_orders = list(permutations(default_order))
    candidate_lists = [candidates_map[qid] for qid in query_order]
    svg_map: Dict[str, Dict[str, str]] = {}
    inline_svg_map: Dict[str, Dict[str, str]] = {}
    order_map: Dict[str, Dict[str, List[str]]] = {}
    marker_map: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
    range_map: Dict[str, Dict[str, Dict[str, Any]]] = {}
    views = []

    for combo in product(*candidate_lists):
        selected = {qid: candidate["candidate_id"] for qid, candidate in zip(query_order, combo)}
        selected_candidates = {qid: candidate for qid, candidate in zip(query_order, combo)}
        combo_key = _selection_key(query_order, selected)
        svg_map.setdefault(combo_key, {})
        inline_svg_map.setdefault(combo_key, {})
        order_map.setdefault(combo_key, {})
        marker_map.setdefault(combo_key, {})
        range_map[combo_key] = {
            "ref": {"start": 1, "end": max(1, ref_length), "label": ref_alias}
        }
        for qid, candidate in selected_candidates.items():
            q_left, q_right = _candidate_region(candidate)
            range_map[combo_key][qid] = {
                "start": q_left,
                "end": q_right,
                "chrom": candidate.get("query_chr"),
                "label": track_aliases.get(qid, qid),
            }

        for order_tuple in all_orders:
            order = list(order_tuple)
            order_key = _order_key(order)
            view_hash = md5(f"{combo_key}@@{order_key}".encode("utf-8")).hexdigest()[:16]
            prefix = os.path.join(asset_dir, f"detail_{view_hash}")
            input_file = f"{prefix}.coords"
            k_file = f"{prefix}.k"
            hl_file = f"{prefix}.hl"
            gff_file = f"{prefix}.gff3"

            with open(input_file, "w", encoding="utf-8") as handle:
                handle.write("ref.fasta query.fasta\n")
                handle.write("NUCMER\n\n")
                handle.write("    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n")
                handle.write("=" * 120 + "\n")
                for first_id, second_id in zip(order, order[1:]):
                    if first_id == "ref" or second_id == "ref":
                        qid = second_id if first_id == "ref" else first_id
                        candidate = selected_candidates.get(qid)
                        if not candidate:
                            continue
                        q_chr = candidate.get("query_chr") or "query"
                        q_len = _lookup_length(query_lengths.get(qid, {}), q_chr, _candidate_region(candidate)[1])
                        for block in candidate.get("blocks") or []:
                            q_start = int(block.get("query_start") or 1)
                            q_end = int(block.get("query_end") or q_start)
                            r_start = int(block.get("ref_start") or 1)
                            r_end = int(block.get("ref_end") or r_start)
                            block_identity = float(block.get("identity") or identity)
                            if first_id == "ref":
                                _write_linkview_relation(
                                    handle, r_start, r_end, q_start, q_end,
                                    ref_length, q_len, block_identity, ref_alias, track_aliases[qid]
                                )
                            else:
                                _write_linkview_relation(
                                    handle, q_start, q_end, r_start, r_end,
                                    q_len, ref_length, block_identity, track_aliases[qid], ref_alias
                                )
                        continue

                    first_candidate = selected_candidates.get(first_id)
                    second_candidate = selected_candidates.get(second_id)
                    if not first_candidate or not second_candidate:
                        continue
                    pairwise_result, reversed_pair = _pairwise_for_selection(
                        pairwise_lookup,
                        first_id,
                        str(first_candidate.get("candidate_id") or ""),
                        second_id,
                        str(second_candidate.get("candidate_id") or ""),
                    )
                    if not pairwise_result:
                        continue
                    if reversed_pair:
                        left_id, right_id = second_id, first_id
                        left_candidate, right_candidate = second_candidate, first_candidate
                    else:
                        left_id, right_id = first_id, second_id
                        left_candidate, right_candidate = first_candidate, second_candidate
                    left_region = _candidate_fasta_region(pairwise_result.get("left_candidate_fasta")) or (
                        left_candidate.get("query_chr"), *_candidate_region(left_candidate)
                    )
                    right_region = _candidate_fasta_region(pairwise_result.get("right_candidate_fasta")) or (
                        right_candidate.get("query_chr"), *_candidate_region(right_candidate)
                    )
                    first_range = _candidate_region(first_candidate)
                    second_range = _candidate_region(second_candidate)
                    first_len = max(1, first_range[1] - first_range[0] + 1)
                    second_len = max(1, second_range[1] - second_range[0] + 1)
                    for block in _read_pairwise_blocks(pairwise_result, left_region, right_region):
                        coords_by_track = {
                            left_id: (block["left_start"], block["left_end"]),
                            right_id: (block["right_start"], block["right_end"]),
                        }
                        first_start, first_end = coords_by_track[first_id]
                        second_start, second_end = coords_by_track[second_id]
                        _write_linkview_relation(
                            handle, first_start, first_end, second_start, second_end,
                            first_len, second_len, block.get("identity") or identity,
                            track_aliases[first_id], track_aliases[second_id]
                        )

            with open(k_file, "w", encoding="utf-8") as handle:
                for track_id in order:
                    if track_id == "ref":
                        handle.write(f"{ref_alias}:1:{max(1, ref_length)}\n")
                    else:
                        left, right = _candidate_region(selected_candidates[track_id])
                        handle.write(f"{track_aliases[track_id]}:{left}:{right}\n")

            with open(hl_file, "w", encoding="utf-8") as handle:
                marker_entries: List[Dict[str, Any]] = []
                for first_id, second_id in zip(order, order[1:]):
                    if first_id == "ref" or second_id == "ref":
                        qid = second_id if first_id == "ref" else first_id
                        candidate = selected_candidates.get(qid)
                        if not candidate:
                            continue
                        q_left, q_right = _candidate_region(candidate)
                        pair_id = f"ref__{qid}"
                        for variant in variants_payload:
                            if variant.get("pair_id") != pair_id:
                                continue
                            if variant.get("query_chr") and not _chrom_matches(variant.get("query_chr"), candidate.get("query_chr")):
                                continue
                            color = "orange" if str(variant.get("type", "")).upper() == "SNP" else "blue"
                            try:
                                q_pos = int(variant.get("query_genome_pos"))
                                r_pos = int(variant.get("ref_source_pos"))
                            except (TypeError, ValueError):
                                continue
                            var_type = str(variant.get("type", "")).upper()
                            if var_type == "SNP":
                                if 1 <= r_pos <= ref_length:
                                    _write_hl_marker(handle, ref_alias, r_pos, color)
                                    _append_marker(marker_entries, "ref", ref_alias, r_pos, color, variant, pair_id, "ref")
                                if q_left <= q_pos <= q_right:
                                    _write_hl_marker(handle, track_aliases[qid], q_pos, color)
                                    _append_marker(marker_entries, qid, track_aliases[qid], q_pos, color, variant, pair_id, "query")
                            elif var_type == "INS" and q_left <= q_pos <= q_right:
                                _write_hl_marker(handle, track_aliases[qid], q_pos, color)
                                _append_marker(marker_entries, qid, track_aliases[qid], q_pos, color, variant, pair_id, "query")
                            elif var_type == "DEL" and 1 <= r_pos <= ref_length:
                                _write_hl_marker(handle, ref_alias, r_pos, color)
                                _append_marker(marker_entries, "ref", ref_alias, r_pos, color, variant, pair_id, "ref")
                        continue

                    first_candidate = selected_candidates.get(first_id)
                    second_candidate = selected_candidates.get(second_id)
                    if not first_candidate or not second_candidate:
                        continue
                    pairwise_result, reversed_pair = _pairwise_for_selection(
                        pairwise_lookup,
                        first_id,
                        str(first_candidate.get("candidate_id") or ""),
                        second_id,
                        str(second_candidate.get("candidate_id") or ""),
                    )
                    if not pairwise_result:
                        continue
                    if reversed_pair:
                        left_id, right_id = second_id, first_id
                        left_candidate, right_candidate = second_candidate, first_candidate
                    else:
                        left_id, right_id = first_id, second_id
                        left_candidate, right_candidate = first_candidate, second_candidate
                    left_region = _candidate_fasta_region(pairwise_result.get("left_candidate_fasta")) or (
                        left_candidate.get("query_chr"), *_candidate_region(left_candidate)
                    )
                    right_region = _candidate_fasta_region(pairwise_result.get("right_candidate_fasta")) or (
                        right_candidate.get("query_chr"), *_candidate_region(right_candidate)
                    )
                    left_range = _candidate_region(left_candidate)
                    right_range = _candidate_region(right_candidate)
                    pairwise_pair_id = f"{pairwise_result.get('pair_id') or f'{left_id}__{right_id}'}::{pairwise_result.get('left_candidate_id')}__{pairwise_result.get('right_candidate_id')}"
                    for variant_index, variant in enumerate(parse_variants(pairwise_result.get("snps")), start=1):
                        variant["variant_id"] = f"{pairwise_pair_id}:var_{variant_index}"
                        color = "orange" if str(variant.get("type", "")).upper() == "SNP" else "blue"
                        left_pos = _local_to_region(left_region, variant.get("ref_source_pos"))
                        right_pos = _local_to_region(right_region, variant.get("query_genome_pos"))
                        var_type = str(variant.get("type", "")).upper()
                        if var_type == "SNP":
                            if left_range[0] <= left_pos <= left_range[1]:
                                _write_hl_marker(handle, track_aliases[left_id], left_pos, color)
                                _append_marker(marker_entries, left_id, track_aliases[left_id], left_pos, color, variant, pairwise_pair_id, "left")
                            if right_range[0] <= right_pos <= right_range[1]:
                                _write_hl_marker(handle, track_aliases[right_id], right_pos, color)
                                _append_marker(marker_entries, right_id, track_aliases[right_id], right_pos, color, variant, pairwise_pair_id, "right")
                        elif var_type == "INS" and right_range[0] <= right_pos <= right_range[1]:
                            _write_hl_marker(handle, track_aliases[right_id], right_pos, color)
                            _append_marker(marker_entries, right_id, track_aliases[right_id], right_pos, color, variant, pairwise_pair_id, "right")
                        elif var_type == "DEL" and left_range[0] <= left_pos <= left_range[1]:
                            _write_hl_marker(handle, track_aliases[left_id], left_pos, color)
                            _append_marker(marker_entries, left_id, track_aliases[left_id], left_pos, color, variant, pairwise_pair_id, "left")
                marker_map[combo_key][order_key] = marker_entries

            with open(gff_file, "w", encoding="utf-8") as handle:
                if ref_gff_source:
                    if ref_region:
                        _write_gff_records(
                            ref_gff_source, handle, ref_alias,
                            source_chr=ref_region[0],
                            region_start=ref_region[1],
                            region_end=ref_region[2],
                            relative=True,
                        )
                    else:
                        _write_gff_records(ref_gff_source, handle, ref_alias)
                for qid in query_order:
                    entry = query_entry_by_id.get(qid) or {}
                    candidate = selected_candidates[qid]
                    q_left, q_right = _candidate_region(candidate)
                    _write_gff_records(
                        entry.get("gff") or entry.get("source_gff"),
                        handle,
                        track_aliases[qid],
                        source_chr=candidate.get("query_chr"),
                        region_start=q_left,
                        region_end=q_right,
                        relative=False,
                    )
            if os.path.exists(gff_file) and os.path.getsize(gff_file) == 0:
                os.remove(gff_file)

            if _run_linkview_for_report(input_file, prefix, k_file, hl_file, gff_file, identity, min_aln_len):
                svg_path = f"{prefix}.svg"
                _process_linkview_svg_file(svg_path)
                order_map[combo_key][order_key] = _parse_linkview_svg_order(svg_path, order, tracks, track_aliases)
                svg_map[combo_key][order_key] = relpath(svg_path, report_dir) or ""
                try:
                    with open(svg_path, "r", encoding="utf-8", errors="ignore") as svg_handle:
                        inline_svg_map[combo_key][order_key] = svg_handle.read()
                except OSError:
                    pass

        if svg_map.get(combo_key):
            views.append({
                "selection_key": combo_key,
                "selected_candidates": selected,
                "svg_count": len(svg_map[combo_key]),
            })

    if not views:
        return {}
    return {
        "linkview_svg_map": svg_map,
        "linkview_svg_inline_map": inline_svg_map,
        "linkview_svg_order_map": order_map,
        "linkview_marker_map": marker_map,
        "linkview_views": views,
        "linkview_track_aliases": track_aliases,
        "linkview_track_ranges": range_map,
        "linkview_default_order_key": _order_key(default_order),
        "linkview_layout": {
            "svg_width": 1200,
            "svg_height": 400,
            "top_margin": 90,
        },
    }


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
    detail_assets = _generate_linkview_detail_assets(
        result=result,
        output_dir=output_dir,
        report_dir=report_dir,
        mode=mode,
        ref_length=ref_length,
        genome_files=genome_files,
        query_order=query_order,
        tracks=tracks,
        candidates_map=candidates_map,
        variants_payload=variants_payload,
        identity=identity,
        min_aln_len=min_aln_len,
    )
    detail_payload = {
        "track_order": ["ref"] + query_order,
        "selected_candidates": selected_candidates,
        "visible_pair_ids": [pair["pair_id"] for pair in pairs],
        "visible_links": visible_links,
        "variant_ids": [variant["variant_id"] for variant in variants_payload],
    }
    detail_payload.update(detail_assets)
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
        "detail": detail_payload,
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
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      line-height: 1.6;
      color: #333;
      background: #f5f5f5;
      padding: 20px;
    }
    .container {
      max-width: 1200px;
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
      flex-wrap: nowrap;
    }
    .header-left {
      display: flex;
      flex-direction: column;
      align-items: flex-start;
      min-width: 0;
    }
    .header-right {
      display: flex;
      align-items: flex-start;
      flex: 0 0 auto;
    }
    .header h1 { font-size: 24px; margin-bottom: 10px; }
    .subtitle {
      display: flex;
      flex-wrap: nowrap;
      align-items: center;
      gap: 15px;
      opacity: 0.95;
      font-size: 14px;
      min-width: 0;
    }
    .gen-time { margin-top: 8px; font-size: 13px; opacity: 0.85; }
    .lang-switch {
      display: flex;
      gap: 0;
      border-radius: 6px;
      overflow: hidden;
      border: 1px solid rgba(255,255,255,0.3);
    }
    .lang-btn {
      height: 32px;
      background: transparent;
      border: none;
      color: rgba(255,255,255,0.7);
      padding: 6px 14px;
      font-size: 13px;
      cursor: pointer;
      transition: all 0.2s;
    }
    .lang-btn:hover {
      background: rgba(255,255,255,0.1);
      color: white;
    }
    .lang-btn.active {
      background: rgba(255,255,255,0.2);
      color: white;
      font-weight: 600;
    }
    .mode-badge {
      display: inline-block;
      padding: 4px 12px;
      border-radius: 20px;
      font-size: 12px;
      font-weight: 600;
      text-transform: uppercase;
      flex: 0 0 auto;
    }
    .mode-geneid { background: #e3f2fd; color: #1976d2; }
    .mode-location { background: #f3e5f5; color: #7b1fa2; }
    .mode-sequence { background: #e8f5e9; color: #388e3c; }
    .json-badge {
      color: rgba(255,255,255,0.85);
      font-size: 12px;
      text-decoration: underline;
      white-space: nowrap;
    }
    .json-badge:hover { color: white; }
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
    .info-table-grid { width: 100%; border-collapse: collapse; margin-bottom: 15px; }
    .info-table-grid td { padding: 10px 15px; border-bottom: 1px solid #eee; vertical-align: top; }
    .info-table-grid .label {
      width: 120px;
      font-weight: 600;
      color: #666;
      background: #f9f9f9;
    }
    .info-table-grid .value { color: #333; width: 30%; }
    .info-table-grid code {
      background: #f0f0f0;
      padding: 2px 6px;
      border-radius: 3px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      word-break: break-all;
    }
    .value-list {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
    }
    .value-chip {
      display: inline-block;
      padding: 2px 8px;
      background: #eef1ff;
      color: #4f5fc9;
      border-radius: 999px;
      font-size: 12px;
      line-height: 1.6;
      max-width: 100%;
      overflow-wrap: anywhere;
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
      overflow: visible;
    }
    .overview-canvas {
      width: 100%;
      display: block;
    }
    .overview-body { display: grid; gap: 12px; }
    .chr-group { display: block; padding: 0; }
    .chr-tracks { display: grid; gap: 12px; }
    .overview-track {
      display: grid;
      grid-template-columns: minmax(120px, 190px) minmax(0, 1fr);
      gap: 10px;
      align-items: center;
    }
    .track-label {
      font-family: Arial, sans-serif;
      font-size: 12px;
      font-weight: 600;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      color: #333;
      text-align: right;
    }
    .track-area {
      position: relative;
      min-width: 0;
      height: 20px;
    }
    .track-length {
      position: absolute;
      top: 3px;
      font-family: Arial, sans-serif;
      font-size: 10px;
      color: #666;
      white-space: nowrap;
    }
    .lane {
      position: relative;
      height: 20px;
      background: transparent;
    }
    .lane::before {
      content: '';
      position: absolute;
      left: 0;
      right: 0;
      top: 0;
      height: 20px;
      box-sizing: border-box;
      background: #e0e0e0;
      border: 1px solid #999;
      border-radius: 3px;
    }
    .candidate {
      position: absolute;
      top: 2px;
      height: 16px;
      min-width: 3px;
      background: #4a90d9;
      border: none;
      border-radius: 0;
      cursor: pointer;
      transition: transform .12s ease, box-shadow .12s ease;
    }
    .candidate:hover { background: #2d6cb5; transform: translateY(-1px); }
    .candidate.selected { outline: none; }
    .candidate.selected::after {
      content: '';
      position: absolute;
      left: 50%;
      bottom: -12px;
      transform: translateX(-50%);
      width: 0;
      height: 0;
      border-left: 6px solid transparent;
      border-right: 6px solid transparent;
      border-bottom: 8px solid #e53935;
    }
    .candidate-primary-label {
      position: absolute;
      left: 50%;
      top: -14px;
      transform: translateX(-50%);
      font-family: Arial, sans-serif;
      font-size: 9px;
      color: #1d4ed8;
      font-weight: 700;
      line-height: 1;
      pointer-events: none;
      white-space: nowrap;
    }
    .viz-legend-note {
      margin-top: 15px;
      padding: 10px 15px;
      background: #f5f5f5;
      border-radius: 4px;
      font-size: 13px;
      color: #666;
    }
    .detail-layout { display: block; }
    .controls {
      background: #fff;
      border: 1px solid #eee;
      border-radius: 8px;
      padding: 14px;
    }
    .track-card-list { display: grid; gap: 10px; }
    .control-row.track-card {
      border: 1px solid #e3e6ef;
      border-radius: 8px;
      background: #fff;
      box-shadow: 0 1px 2px rgba(0,0,0,.04);
      overflow: hidden;
    }
    .control-row.track-card.dragging { opacity: .55; }
    .track-card-summary {
      display: grid;
      grid-template-columns: 1fr 12px;
      gap: 8px;
      align-items: center;
      min-height: 36px;
      padding: 8px 10px;
      cursor: grab;
      list-style: none;
      color: #333;
      font-weight: 600;
      font-size: 13px;
    }
    .track-card-summary::-webkit-details-marker { display: none; }
    .track-card-name {
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }
    .track-card-toggle {
      width: 0;
      height: 0;
      border-left: 5px solid transparent;
      border-right: 5px solid transparent;
      border-top: 6px solid #667eea;
      transition: transform .15s ease;
    }
    .track-card:not([open]) .track-card-toggle { transform: rotate(-90deg); }
    .track-card-body {
      display: grid;
      grid-template-columns: 1fr 32px 32px;
      gap: 6px;
      align-items: center;
      padding: 0 10px 10px;
    }
    .track-card-meta {
      min-height: 32px;
      display: flex;
      align-items: center;
      color: #777;
      font-size: 12px;
    }
    select, .track-card-body button {
      height: 32px;
      border: 1px solid #ddd;
      border-radius: 4px;
      background: #fff;
      color: #333;
    }
    select { width: 100%; padding: 0 8px; }
    .track-card-body button { cursor: pointer; color: #667eea; font-weight: 600; }
    .track-card-body button:hover { background: #f3f0ff; }
    .track-stack { min-width: 0; }
    .detail-canvas {
      position: relative;
      border: 1px solid #e6e8ef;
      border-radius: 8px;
      background: #fff;
      overflow: visible;
      padding: 8px 0;
    }
    .detail-linkview-svg {
      min-height: 220px;
      min-width: 0;
      box-sizing: border-box;
    }
    .detail-svg {
      display: block;
      width: 100%;
      height: auto;
    }
    .detail-track-overlays {
      position: absolute;
      inset: 8px 0;
      pointer-events: none;
    }
    .detail-canvas.drag-active .detail-track-overlays {
      pointer-events: auto;
    }
    .detail-linkview-img,
    .detail-linkview-svg svg {
      display: block;
      width: 100%;
      height: auto;
      min-height: 220px;
      object-fit: contain;
    }
    .detail-track-control {
      position: absolute;
      display: inline-flex;
      align-items: center;
      gap: 5px;
      min-height: 26px;
      padding: 2px 4px;
      background: rgba(255,255,255,.94);
      border-radius: 4px;
      cursor: grab;
      user-select: none;
      z-index: 5;
      pointer-events: auto;
    }
    .detail-track-control.dragging { opacity: .55; }
    .detail-track-control.drop-before {
      box-shadow: 0 -3px 0 #667eea;
    }
    .detail-track-control.drop-after {
      box-shadow: 0 3px 0 #667eea;
    }
    .detail-track-control.ref-control {
      background: transparent;
      cursor: grab;
      padding-left: 0;
    }
    .detail-drag-handle {
      width: 12px;
      height: 18px;
      flex: 0 0 12px;
      opacity: .62;
      background-image: radial-gradient(circle, #8b93a7 1.3px, transparent 1.4px);
      background-size: 5px 5px;
      background-position: 0 1px;
    }
    .detail-track-control:hover .detail-drag-handle { opacity: .95; }
    .detail-track-name {
      flex: 1 1 auto;
      min-width: 0;
      max-width: none;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      font-family: Arial, sans-serif;
      font-size: 13px;
      font-weight: 600;
      color: #333;
    }
    .detail-track-select {
      flex: 0 0 62px;
      width: 62px;
      height: 22px;
      padding: 0 16px 0 5px;
      border: 1px solid #d6d9e6;
      border-radius: 4px;
      background: #fff;
      color: #333;
      font-size: 11px;
      cursor: pointer;
    }
    .detail-track-end-label {
      font-family: Arial, sans-serif;
      font-size: 10px;
      fill: #555;
      pointer-events: none;
    }
    .detail-gene-bracket {
      fill: none;
      stroke: #333;
      stroke-width: 2;
      pointer-events: none;
    }
    .detail-gene-label {
      font-family: Arial, sans-serif;
      font-size: 12px;
      fill: #333;
      font-weight: 500;
      pointer-events: none;
    }
    .svg-tooltip {
      position: absolute;
      background: rgba(0, 0, 0, 0.85);
      color: white;
      padding: 8px 12px;
      border-radius: 4px;
      font-size: 12px;
      pointer-events: none;
      opacity: 0;
      transition: opacity 0.2s;
      z-index: 1000;
      max-width: 300px;
      white-space: nowrap;
    }
    .svg-tooltip.visible { opacity: 1; }
    .svg-tooltip .tooltip-title {
      font-weight: 700;
      margin-bottom: 4px;
      color: #ffd700;
    }
    .svg-tooltip .tooltip-row { margin: 2px 0; }
    .tooltip-label { color: #aaa; }
    .tooltip-seq {
      font-family: monospace;
      background: rgba(255,255,255,0.1);
      padding: 1px 4px;
      border-radius: 2px;
    }
    .legend-item { transition: opacity 0.2s; }
    .legend-item:hover { opacity: 0.8; }
    .legend-item.disabled { opacity: 0.4; }
    .legend-item.disabled .legend-text { text-decoration: line-through; }
    .zoom-btn {
      position: absolute;
      right: 10px;
      bottom: 10px;
      width: 36px;
      height: 36px;
      border: none;
      border-radius: 6px;
      background: rgba(102, 126, 234, 0.9);
      color: white;
      cursor: pointer;
      display: flex;
      align-items: center;
      justify-content: center;
      transition: all 0.2s;
      z-index: 100;
    }
    .zoom-btn:hover {
      background: rgba(102, 126, 234, 1);
      transform: scale(1.1);
    }
    .zoom-modal {
      display: none;
      position: fixed;
      inset: 0;
      background: rgba(0, 0, 0, 0.8);
      z-index: 10000;
      overflow: auto;
    }
    .zoom-modal.visible {
      display: flex;
      align-items: center;
      justify-content: center;
    }
    .zoom-modal-content {
      position: relative;
      background: white;
      border-radius: 8px;
      width: 95%;
      max-width: 1400px;
      max-height: 90vh;
      overflow: auto;
      padding: 20px;
    }
    .zoom-modal-close {
      position: absolute;
      top: 10px;
      right: 15px;
      font-size: 28px;
      font-weight: bold;
      color: #666;
      background: none;
      border: none;
      cursor: pointer;
      z-index: 10;
    }
    .zoom-modal-close:hover { color: #333; }
    .zoom-modal-body {
      width: 100%;
      overflow-x: auto;
      position: relative;
    }
    .zoom-modal-body .detail-canvas {
      min-width: 1200px;
      padding: 8px 0;
    }
    .zoom-modal-body svg { min-width: 1200px; }
    .detail-label { font-family: Arial, sans-serif; font-size: 13px; font-weight: 600; fill: #333; }
    .detail-subtext { font-family: Arial, sans-serif; font-size: 11px; fill: #777; }
    .detail-chro { fill: none; stroke: #000; stroke-width: 1.5; }
    .detail-axis { stroke: #000; stroke-width: 1.5; stroke-linecap: round; }
    .detail-axis-light { stroke: #a9b1c3; stroke-width: 1; }
    .detail-block { fill: #b7b7b7; opacity: .5; stroke: none; }
    .detail-block.reverse { fill: #b7b7b7; opacity: .5; }
    .detail-match { fill: #607b8b; opacity: .45; stroke: none; }
    .detail-ref-gene { fill: #667eea; stroke: #5364c9; stroke-width: 1; }
    .detail-variant.snp { fill: orange; stroke: #b76d00; }
    .detail-variant.indel { fill: blue; stroke: #176a95; }
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
      .detail-layout { display: block; }
      .detail-track-name { max-width: 92px; }
      .detail-track-select { width: 62px; }
      .header { padding: 24px 20px; }
      .header-content { gap: 12px; }
      .subtitle { gap: 10px; }
      .section { padding: 20px 16px; }
      .overview-track { grid-template-columns: minmax(104px, 150px) minmax(0, 1fr); }
      .controls { margin-bottom: 14px; }
    }
  </style>
</head>
<body>
<div class="container">
  <div class="header">
    <div class="header-content">
      <div class="header-left">
        <h1 data-zh="GeneScreen 比对分析报告" data-en="GeneScreen Alignment Analysis Report">GeneScreen 比对分析报告</h1>
        <div class="subtitle">
          <span class="mode-badge" id="modeBadge">GeneScreen</span>
          <a class="json-badge" href="data.json" data-zh="data.json" data-en="data.json">data.json</a>
        </div>
        <div class="gen-time" id="generatedAt"></div>
      </div>
      <div class="header-right">
        <div class="lang-switch">
          <button class="lang-btn active" data-lang="zh" onclick="switchLang('zh')">中文</button>
          <button class="lang-btn" data-lang="en" onclick="switchLang('en')">EN</button>
        </div>
      </div>
    </div>
  </div>

  <div class="section input-section">
    <h2 data-zh="初始化" data-en="Initialization">初始化</h2>
    <table class="info-table-grid" id="initTable"></table>
  </div>

  <div class="section">
    <h2 data-zh="输入文件" data-en="Input Files">输入文件</h2>
    <table class="file-table" id="genomeTable"></table>
  </div>

  <div class="section overview-section visualization-section">
    <h2 data-zh="宏观比对图" data-en="Alignment Overview">宏观比对图</h2>
    <div class="viz-container">
      <div id="overview" class="overview-canvas"></div>
      <div class="viz-legend-note">
        <p data-zh="点击比对区域可跳转到对应的局部比对图"
           data-en="Click alignment region to jump to the corresponding local alignment">点击比对区域可跳转到对应的局部比对图</p>
      </div>
    </div>
  </div>

  <div class="section visualization-section">
    <h2 data-zh="局部比对图" data-en="Local Alignment">局部比对图</h2>
    <div class="viz-container">
      <div class="detail-layout">
        <div id="detail" class="track-stack"></div>
      </div>
      <div class="viz-legend-note">
        <p data-zh="轨道名右侧下拉框可切换候选片段；拖动轨道名可调整显示顺序。"
           data-en="Use the dropdown next to a track name to switch candidates; drag track names to reorder tracks.">轨道名右侧下拉框可切换候选片段；拖动轨道名可调整显示顺序。</p>
      </div>
    </div>
  </div>

  <div id="zoom-modal" class="zoom-modal">
    <div class="zoom-modal-content">
      <button class="zoom-modal-close" id="zoom-modal-close">&times;</button>
      <div class="zoom-modal-body" id="zoom-modal-body"></div>
    </div>
  </div>

  <div class="section output-section">
    <h2 data-zh="结果文件" data-en="Output Files">结果文件</h2>
    <table class="file-table" id="outputTable"></table>
  </div>
</div>
<div id="copyToast" class="copy-toast" data-zh="已复制到剪贴板" data-en="Copied to clipboard">已复制到剪贴板</div>

  <script id="embedded-report-data" type="application/json">__EMBEDDED_DATA__</script>
  <script>
    let reportData = JSON.parse(document.getElementById('embedded-report-data').textContent);
    fetch('data.json').then(r => r.ok ? r.json() : reportData).then(data => { reportData = data; init(); }).catch(init);

    const state = { selected: {}, order: [] };
    let detailRenderToken = 0;
    let currentLang = 'zh';
    let resizeAlignTimer = null;
    const i18n = {
      geneId: { zh: '基因ID', en: 'Gene ID' },
      mode: { zh: '模式', en: 'Mode' },
      referenceGenome: { zh: '参考基因组', en: 'Reference Genome' },
      sequenceLength: { zh: '序列长度', en: 'Sequence Length' },
      identityThreshold: { zh: 'Identity 阈值', en: 'Identity Threshold' },
      minAlignmentLength: { zh: '最小比对长度', en: 'Minimum Alignment Length' },
      pairwiseTopN: { zh: 'Pairwise Top-N', en: 'Pairwise Top-N' },
      allCandidates: { zh: '全量候选', en: 'All Candidates' },
      upstream: { zh: '上游延伸', en: 'Upstream Extension' },
      downstream: { zh: '下游延伸', en: 'Downstream Extension' },
      type: { zh: '类型', en: 'Type' },
      name: { zh: '名称', en: 'Name' },
      fasta: { zh: 'FASTA', en: 'FASTA' },
      annotation: { zh: '注释', en: 'Annotation' },
      reference: { zh: '参考', en: 'Reference' },
      query: { zh: '查询', en: 'Query' },
      annotationNotUsed: { zh: '未使用注释', en: 'Annotation not used' },
      genomes: { zh: '基因组', en: 'Genomes' },
      queryGenomes: { zh: '查询基因组', en: 'Query Genomes' },
      candidates: { zh: '候选片段', en: 'Candidates' },
      primaryCandidate: { zh: '主', en: 'Primary' },
      visibleLinks: { zh: '可见连接', en: 'Visible Links' },
      snpCount: { zh: 'SNP 数量', en: 'SNP Count' },
      indelCount: { zh: 'Indel 数量', en: 'Indel Count' },
      noCandidates: { zh: '无候选片段', en: 'No candidates' },
      candidateSwitch: { zh: '候选切换', en: 'Candidate Selection' },
      dragTrack: { zh: '拖动轨道调整顺序', en: 'Drag track to reorder' },
      moveUp: { zh: '上移', en: 'Move up' },
      moveDown: { zh: '下移', en: 'Move down' },
      noTracks: { zh: '无可显示轨道', en: 'No tracks' },
      filename: { zh: '文件名', en: 'File Name' },
      description: { zh: '描述', en: 'Description' },
      path: { zh: '路径', en: 'Path' },
      position: { zh: '位置', en: 'Position' },
      mutation: { zh: '变异', en: 'Mutation' },
      refBase: { zh: '参考', en: 'Ref' },
      qryBase: { zh: '查询', en: 'Qry' },
      insertSeq: { zh: '插入序列', en: 'Inserted Seq' },
      deleteSeq: { zh: '缺失序列', en: 'Deleted Seq' },
      length: { zh: '长度', en: 'Length' },
      insertion: { zh: '插入 (INS)', en: 'Insertion (INS)' },
      deletion: { zh: '缺失 (DEL)', en: 'Deletion (DEL)' },
      utr3: { zh: "3' UTR", en: "3' UTR" },
      utr3Desc: { zh: '三端非翻译区', en: "3' Untranslated Region" },
      utr5: { zh: "5' UTR", en: "5' UTR" },
      utr5Desc: { zh: '五端非翻译区', en: "5' Untranslated Region" },
      cds: { zh: 'CDS', en: 'CDS' },
      cdsDesc: { zh: '编码序列', en: 'Coding Sequence' },
      structuredData: { zh: '新版报告结构化数据', en: 'Structured data for the new report' },
      generatedAt: { zh: '生成时间', en: 'Generated' },
      geneIdMode: { zh: 'Gene ID 模式', en: 'Gene ID Mode' },
      locationMode: { zh: 'Location 模式', en: 'Location Mode' },
      sequenceMode: { zh: 'Sequence 模式', en: 'Sequence Mode' },
      reportTitle: { zh: 'GeneScreen 比对分析报告', en: 'GeneScreen Alignment Analysis Report' }
    };
    function esc(value) {
      return String(value ?? '').replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
    }
    function t(key) {
      return (i18n[key] && i18n[key][currentLang]) || key;
    }
    function modeText(mode) {
      const raw = String(mode || 'GeneScreen');
      const normalized = raw.toLowerCase().replace(/[\s-]+/g, '_');
      if (normalized.includes('gene_id')) return t('geneIdMode');
      if (normalized.includes('location')) return t('locationMode');
      if (normalized.includes('sequence')) return t('sequenceMode');
      return raw;
    }
    function modeClass(mode) {
      const normalized = String(mode || '').toLowerCase().replace(/[\s-]+/g, '_');
      if (normalized.includes('gene_id')) return 'mode-geneid';
      if (normalized.includes('location')) return 'mode-location';
      if (normalized.includes('sequence')) return 'mode-sequence';
      return 'mode-geneid';
    }
    function applyStaticLang() {
      document.documentElement.lang = currentLang === 'en' ? 'en' : 'zh-CN';
      document.querySelectorAll('.lang-btn').forEach(btn => {
        btn.classList.toggle('active', btn.dataset.lang === currentLang);
      });
      document.querySelectorAll('[data-zh][data-en]').forEach(el => {
        el.textContent = el.dataset[currentLang];
      });
    }
    function switchLang(lang) {
      currentLang = lang === 'en' ? 'en' : 'zh';
      renderHeader();
      renderInit();
      renderGenomes();
      renderOverview();
      renderControls();
      renderDetail();
      renderOutputs();
      applyStaticLang();
    }
    function init() {
      state.selected = Object.assign({}, reportData.default_selection.selected_candidates || {});
      state.order = (reportData.default_selection.track_order || []).slice();
      switchLang(currentLang);
      window.addEventListener('resize', () => {
        clearTimeout(resizeAlignTimer);
        resizeAlignTimer = setTimeout(() => {
          const assetOrder = resolvedLinkviewOrder(state.order.length ? state.order.slice() : ['ref'].concat(queryIds()));
          const order = mappedLinkviewOrder(assetOrder);
          document.querySelectorAll('.detail-linkview-svg-root').forEach(svg => {
            ensureDetailGutter(svg, order);
            requestAnimationFrame(() => alignTrackControls(svg, order));
          });
        }, 80);
      });
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
      const modeBadge = document.getElementById('modeBadge');
      modeBadge.textContent = modeText(reportData.mode);
      modeBadge.className = `mode-badge ${modeClass(reportData.mode)}`;
      const generatedText = reportData.generated_at ? `${t('generatedAt')}: ${reportData.generated_at}` : '';
      document.getElementById('generatedAt').textContent = generatedText;
    }
    function renderInit() {
      const input = reportData.input || {};
      const params = reportData.parameters || {};
      const ref = (reportData.genomes && reportData.genomes.ref) || {};
      const queries = (reportData.genomes && reportData.genomes.queries) || [];
      const queryNames = queries.map(q => q.name || q.genome_id).filter(Boolean);
      const queryValue = queryNames.length
        ? `<div class="value-list">${queryNames.map(name => `<span class="value-chip">${esc(name)}</span>`).join('')}</div>`
        : '-';
      const extraction = input.extraction_info || {};
      const upstream = Number(params.query_upstream || extraction.upstream || 0);
      const downstream = Number(params.query_downstream || extraction.downstream || 0);
      const rows = [
        [t('referenceGenome'), esc(ref.name || '-')],
        [t('queryGenomes'), queryValue],
        [t('geneId'), `<code>${esc(input.id || '-')}</code>`],
        [t('identityThreshold'), params.identity != null ? `${params.identity}%` : '-'],
        [t('minAlignmentLength'), params.min_aln_len != null ? `${params.min_aln_len} bp` : '-'],
        [t('pairwiseTopN'), params.pairwise_all ? t('allCandidates') : (params.candidate_limit ?? '-')]
      ];
      if (upstream > 0) rows.push([t('upstream'), `${upstream} bp`]);
      if (downstream > 0) rows.push([t('downstream'), `${downstream} bp`]);
      const htmlRows = [];
      for (let i = 0; i < rows.length; i += 2) {
        const left = rows[i];
        const right = rows[i + 1];
        htmlRows.push(
          `<tr><td class="label">${left[0]}</td><td class="value">${left[1]}</td>` +
          (right ? `<td class="label">${right[0]}</td><td class="value">${right[1]}</td>` : '<td class="label"></td><td class="value"></td>') +
          '</tr>'
        );
      }
      document.getElementById('initTable').innerHTML = htmlRows.join('');
    }
    function renderGenomes() {
      const ref = (reportData.genomes && reportData.genomes.ref) || {};
      const queries = (reportData.genomes && reportData.genomes.queries) || [];
      const rows = [
        `<tr><th>${t('type')}</th><th>${t('name')}</th><th>${t('fasta')}</th><th>${t('annotation')}</th></tr>`,
        `<tr><td>${t('reference')}</td><td>${esc(ref.name || '-')}</td><td><code class="path-code">${esc(ref.source_fasta || ref.artifact_fasta || '-')}</code></td><td><code class="path-code">${esc(ref.source_gff || ref.artifact_gff || '-')}</code></td></tr>`
      ];
      for (const q of queries) {
        rows.push(`<tr><td>${t('query')}</td><td>${esc(q.name || q.genome_id)}</td><td><code class="path-code">${esc(q.source_fasta || '-')}</code></td><td>${q.has_annotation ? `<code class="path-code">${esc(q.source_gff || '-')}</code>` : `<span class="file-na">${t('annotationNotUsed')}</span>`}</td></tr>`);
      }
      document.getElementById('genomeTable').innerHTML = rows.join('');
    }
    function renderOverview() {
      const root = document.getElementById('overview');
      const chrMap = new Map();
      for (const qid of queryIds()) {
        for (const candidate of reportData.candidates[qid] || []) {
          const chr = candidate.query_chr || 'unknown';
          if (!chrMap.has(chr)) chrMap.set(chr, []);
          chrMap.get(chr).push({ qid, candidate, bounds: candidateBounds(candidate) });
        }
      }
      if (!chrMap.size) { root.innerHTML = `<div class="empty">${t('noCandidates')}</div>`; return; }
      const body = document.createElement('div');
      body.className = 'overview-body';
      root.innerHTML = '';
      root.appendChild(body);
      const chrScore = chr => (chrMap.get(chr) || []).reduce((sum, entry) => sum + Number(entry.candidate.total_aln_len || 0), 0);
      const chrNames = Array.from(chrMap.keys()).sort((a, b) => {
        const scoreDelta = chrScore(b) - chrScore(a);
        return scoreDelta || String(a).localeCompare(String(b), undefined, { numeric: true });
      });
      const trackLength = (chr, qid) => {
        const entries = (chrMap.get(chr) || []).filter(e => e.qid === qid);
        return Math.max(...entries.map(e => e.bounds.end), 1);
      };
      const visibleTrackLengths = [];
      for (const chr of chrNames) {
        for (const qid of queryIds()) {
          if ((chrMap.get(chr) || []).some(e => e.qid === qid)) {
            visibleTrackLengths.push(trackLength(chr, qid));
          }
        }
      }
      const maxTrackLength = Math.max(...visibleTrackLengths, 1);
      for (const chr of chrNames) {
        const entries = chrMap.get(chr);
        const group = document.createElement('div');
        group.className = 'chr-group';
        group.innerHTML = `<div class="chr-tracks"></div>`;
        const tracks = group.querySelector('.chr-tracks');
        for (const qid of queryIds()) {
          const candidates = entries.filter(e => e.qid === qid).map(e => e.candidate);
          if (!candidates.length) continue;
          const track = trackById(qid);
          const rowLabel = `${formatChrName(chr)}_${track.name || qid}`;
          const rangeStart = 0;
          const rangeEnd = trackLength(chr, qid);
          const trackScale = Math.max(0.08, rangeEnd / maxTrackLength);
          const lenLabel = formatLength(rangeEnd);
          const row = document.createElement('div');
          row.className = 'overview-track';
          row.innerHTML = `<div class="track-label" title="${esc(rowLabel)}">${esc(rowLabel)}</div><div class="track-area"><div class="lane" style="width:calc((100% - 54px) * ${trackScale})"></div><div class="track-length" style="left:calc((100% - 54px) * ${trackScale} + 5px)">${lenLabel}</div></div>`;
          const lane = row.querySelector('.lane');
          for (const candidate of candidates) {
            const bounds = candidateBounds(candidate);
            const left = ((bounds.start - rangeStart) / Math.max(1, rangeEnd - rangeStart)) * 100;
            const width = Math.max(0.8, ((bounds.end - bounds.start) / Math.max(1, rangeEnd - rangeStart)) * 100);
            const el = document.createElement('div');
            el.className = 'candidate' + (state.selected[qid] === candidate.candidate_id ? ' selected' : '');
            const region = candidate.query_region_start ? `${candidate.query_chr}:${candidate.query_region_start}-${candidate.query_region_end}` : `${candidate.query_chr}:${candidate.query_start}-${candidate.query_end}`;
            el.style.left = `${Math.max(0, Math.min(99, left))}%`;
            el.style.width = `${Math.max(0.8, Math.min(width, 100 - left))}%`;
            el.title = `${candidate.candidate_id} ${region} identity=${candidate.identity || '-'} coverage=${candidate.coverage || '-'}`;
            el.onclick = () => { state.selected[qid] = candidate.candidate_id; renderOverview(); renderControls(); renderDetail(); };
            if (candidate.is_best) {
              const primaryLabel = document.createElement('span');
              primaryLabel.className = 'candidate-primary-label';
              primaryLabel.textContent = t('primaryCandidate');
              el.appendChild(primaryLabel);
            }
            lane.appendChild(el);
          }
          tracks.appendChild(row);
        }
        body.appendChild(group);
      }
    }
    function formatLength(value) {
      const length = Number(value) || 0;
      if (length >= 1e6) return `${(length / 1e6).toFixed(1)}Mb`;
      if (length >= 1e3) return `${(length / 1e3).toFixed(1)}kb`;
      return `${Math.max(0, Math.round(length))}bp`;
    }
    function formatChrName(value) {
      const text = String(value || 'chr');
      return text.toLowerCase().startsWith('chr') ? text : `chr${text}`;
    }
    function renderControls() {
      const root = document.getElementById('controls');
      if (root) root.innerHTML = '';
    }
    function moveTrack(trackId, delta) {
      const index = state.order.indexOf(trackId);
      const target = index + delta;
      if (index < 0 || target < 0 || target >= state.order.length) return;
      state.order.splice(index, 1);
      state.order.splice(target, 0, trackId);
      renderControls();
      renderDetail();
    }
    function isZoomModalOpen() {
      const modal = document.getElementById('zoom-modal');
      return Boolean(modal && modal.classList.contains('visible'));
    }
    function syncOpenZoomModal() {
      if (isZoomModalOpen()) renderZoomModalBody();
    }
    function clearDropHints(root = document) {
      root.querySelectorAll('.detail-track-control.drop-before,.detail-track-control.drop-after').forEach(row => {
        row.classList.remove('drop-before', 'drop-after');
      });
    }
    function applyTrackDrop(dragged, target, placement) {
      if (!dragged || dragged === target) return;
      const nextOrder = state.order.filter(id => id !== dragged);
      const index = nextOrder.indexOf(target);
      if (index < 0) {
        nextOrder.push(dragged);
      } else {
        nextOrder.splice(placement === 'after' ? index + 1 : index, 0, dragged);
      }
      state.order = nextOrder;
      renderDetail();
      syncOpenZoomModal();
    }
    function draggedTrackId(root) {
      const dragging = root.querySelector('.detail-track-control.dragging') || document.querySelector('.detail-track-control.dragging');
      return dragging ? dragging.dataset.trackId : '';
    }
    function nearestTrackControl(root, clientY, draggedId = '') {
      const rows = Array.from(root.querySelectorAll('.detail-track-control'))
        .filter(row => row.dataset.trackId && row.dataset.trackId !== draggedId && !row.classList.contains('dragging'));
      if (!rows.length) return null;
      return rows
        .map(row => {
          const rect = row.getBoundingClientRect();
          const center = rect.top + rect.height / 2;
          return { row, rect, distance: Math.abs(clientY - center), placement: clientY > center ? 'after' : 'before' };
        })
        .sort((a, b) => a.distance - b.distance)[0];
    }
    function markNearestDrop(root, event) {
      const dragged = (event && event.dataTransfer && event.dataTransfer.getData('text/plain')) || draggedTrackId(root);
      const nearest = nearestTrackControl(root, event.clientY, dragged);
      clearDropHints(root);
      if (!nearest) return null;
      nearest.row.classList.add(nearest.placement === 'after' ? 'drop-after' : 'drop-before');
      return nearest;
    }
    function wireDetailTrackControls(root = document) {
      const overlay = root.querySelector('.detail-track-overlays') || root;
      if (!overlay.dataset.dropWired) {
        overlay.dataset.dropWired = '1';
        overlay.addEventListener('dragover', e => {
          e.preventDefault();
          markNearestDrop(root, e);
        });
        overlay.addEventListener('dragleave', e => {
          if (!overlay.contains(e.relatedTarget)) clearDropHints(root);
        });
        overlay.addEventListener('drop', e => {
          e.preventDefault();
          const dragged = e.dataTransfer.getData('text/plain');
          const nearest = markNearestDrop(root, e);
          clearDropHints(root);
          if (nearest) applyTrackDrop(dragged, nearest.row.dataset.trackId, nearest.placement);
        });
      }
      root.querySelectorAll('.detail-track-control').forEach(row => {
        row.draggable = true;
        row.addEventListener('dragstart', e => {
          if (e.target.closest('select')) {
            e.preventDefault();
            return;
          }
          row.classList.add('dragging');
          const canvas = row.closest('.detail-canvas');
          if (canvas) canvas.classList.add('drag-active');
          e.dataTransfer.setData('text/plain', row.dataset.trackId);
        });
        row.addEventListener('dragend', () => {
          row.classList.remove('dragging');
          const canvas = row.closest('.detail-canvas');
          if (canvas) canvas.classList.remove('drag-active');
          clearDropHints(root);
        });
        row.addEventListener('dragover', e => {
          e.preventDefault();
          markNearestDrop(root, e);
        });
        row.addEventListener('dragleave', e => {
          if (!row.contains(e.relatedTarget)) row.classList.remove('drop-before', 'drop-after');
        });
        row.addEventListener('drop', e => {
          e.preventDefault();
          const dragged = e.dataTransfer.getData('text/plain');
          const nearest = markNearestDrop(root, e);
          clearDropHints(root);
          if (nearest) applyTrackDrop(dragged, nearest.row.dataset.trackId, nearest.placement);
        });
        const select = row.querySelector('select');
        if (select) {
          select.addEventListener('mousedown', e => e.stopPropagation());
          select.addEventListener('click', e => e.stopPropagation());
          select.addEventListener('dragstart', e => e.preventDefault());
          select.onchange = e => {
            state.selected[row.dataset.trackId] = e.target.value;
            renderOverview();
            renderDetail();
            syncOpenZoomModal();
          };
        }
      });
    }
    function wireControlDrag(row) {
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
      if (reportData.detail && reportData.detail.linkview_svg_map) {
        renderLinkviewDetail();
        return;
      }
      renderVectorDetail();
    }
    function detailSelectionKey() {
      return queryIds().map(qid => `${qid}=${state.selected[qid] || ''}`).join('||');
    }
    function detailOrderKey(order) {
      return order.join('||');
    }
    function resolvedLinkviewOrder(order) {
      const detail = reportData.detail || {};
      const viewMap = detail.linkview_svg_inline_map || detail.linkview_svg_map || {};
      const selectedViews = viewMap[detailSelectionKey()];
      if (!selectedViews) return order;
      const requestedKey = detailOrderKey(order);
      if (selectedViews[requestedKey]) return order;
      const defaultKey = detail.linkview_default_order_key;
      if (defaultKey && selectedViews[defaultKey]) return defaultKey.split('||').filter(Boolean);
      const firstKey = Object.keys(selectedViews)[0];
      return firstKey ? firstKey.split('||').filter(Boolean) : order;
    }
    function mappedLinkviewOrder(assetOrder) {
      const detail = reportData.detail || {};
      const mappedViews = (detail.linkview_svg_order_map || {})[detailSelectionKey()] || {};
      const mapped = mappedViews[detailOrderKey(assetOrder)];
      return Array.isArray(mapped) && mapped.length ? mapped.slice() : assetOrder;
    }
    function resolveLinkviewSvg(order) {
      const detail = reportData.detail || {};
      const viewMap = detail.linkview_svg_map || {};
      const selectedViews = viewMap[detailSelectionKey()];
      if (!selectedViews) return null;
      const exact = selectedViews[detailOrderKey(order)];
      if (exact) return exact;
      const defaultKey = detail.linkview_default_order_key;
      if (defaultKey && selectedViews[defaultKey]) return selectedViews[defaultKey];
      const firstKey = Object.keys(selectedViews)[0];
      return firstKey ? selectedViews[firstKey] : null;
    }
    function resolveLinkviewInlineSvg(order) {
      const detail = reportData.detail || {};
      const viewMap = detail.linkview_svg_inline_map || {};
      const selectedViews = viewMap[detailSelectionKey()];
      if (!selectedViews) return null;
      const exact = selectedViews[detailOrderKey(order)];
      if (exact) return exact;
      const defaultKey = detail.linkview_default_order_key;
      if (defaultKey && selectedViews[defaultKey]) return selectedViews[defaultKey];
      const firstKey = Object.keys(selectedViews)[0];
      return firstKey ? selectedViews[firstKey] : null;
    }
    function resolveLinkviewMarkers(order) {
      const detail = reportData.detail || {};
      const markerViews = detail.linkview_marker_map || {};
      const selectedMarkers = markerViews[detailSelectionKey()];
      if (!selectedMarkers) return [];
      const exact = selectedMarkers[detailOrderKey(order)];
      if (exact) return exact;
      const defaultKey = detail.linkview_default_order_key;
      if (defaultKey && selectedMarkers[defaultKey]) return selectedMarkers[defaultKey];
      const firstKey = Object.keys(selectedMarkers)[0];
      return firstKey ? selectedMarkers[firstKey] : [];
    }
    function currentTrackRanges(order) {
      const detail = reportData.detail || {};
      const ranges = (detail.linkview_track_ranges || {})[detailSelectionKey()] || {};
      return order.map(trackId => {
        const fallback = trackId === 'ref' ? { start: 1, end: refLength() } : candidateBounds(candidateById(trackId, state.selected[trackId]));
        return { trackId, range: ranges[trackId] || fallback || { start: 1, end: 1 } };
      });
    }
    function normalizeTrackToken(value) {
      return String(value || '').toLowerCase().replace(/[^a-z0-9]+/g, '');
    }
    function matchTokensForTrack(trackId) {
      const detail = reportData.detail || {};
      const aliases = detail.linkview_track_aliases || {};
      const track = trackById(trackId);
      return [trackId, aliases[trackId], track.name, track.sequence_id, track.genome_id]
        .map(normalizeTrackToken)
        .filter(Boolean);
    }
    function inferSvgTrackOrder(svg, fallbackOrder) {
      const labels = Array.from(svg.querySelectorAll('text.label'))
        .map(label => ({ text: normalizeTrackToken(label.textContent), y: Number(label.getAttribute('y') || 0) }))
        .filter(item => item.text)
        .sort((a, b) => a.y - b.y);
      if (!labels.length) return fallbackOrder;
      const allIds = (reportData.tracks || []).map(track => track.track_id);
      const candidates = fallbackOrder.concat(allIds.filter(id => !fallbackOrder.includes(id)));
      const used = new Set();
      const inferred = [];
      for (const label of labels) {
        const matched = candidates.find(trackId => {
          if (used.has(trackId)) return false;
          return matchTokensForTrack(trackId).some(token => token && (label.text.includes(token) || token.includes(label.text)));
        });
        if (matched) {
          inferred.push(matched);
          used.add(matched);
        }
      }
      const chroCount = svg.querySelectorAll('rect.chro').length;
      return inferred.length === chroCount ? inferred : fallbackOrder;
    }
    function makeTrackControls(order, svgHeight, topMargin) {
      const displayHeight = Math.max(1, svgHeight + topMargin);
      return order.map((trackId, index) => {
        const track = trackById(trackId);
        const y = (svgHeight / (order.length + 1)) * (index + 1);
        const topPct = (((y + topMargin) / displayHeight) * 100).toFixed(3);
        const trackTitle = esc(track.name || trackId);
        const handle = `<span class="detail-drag-handle" title="${t('dragTrack')}"></span>`;
        if (trackId === 'ref') {
          return `<div class="detail-track-control ref-control" data-track-id="${esc(trackId)}" title="${trackTitle}" style="left:10px; top:${topPct}%; width:220px; transform:translateY(-50%);">${handle}<span class="detail-track-name">${trackTitle}</span></div>`;
        }
        const options = (reportData.candidates[trackId] || [])
          .map(c => `<option value="${esc(c.candidate_id)}" ${state.selected[trackId] === c.candidate_id ? 'selected' : ''}>${esc(c.candidate_id)}</option>`)
          .join('');
        return `<div class="detail-track-control" data-track-id="${esc(trackId)}" title="${trackTitle}" style="left:10px; top:${topPct}%; width:220px; transform:translateY(-50%);">${handle}<span class="detail-track-name">${trackTitle}</span><select class="detail-track-select" title="${t('candidateSwitch')}">${options}</select></div>`;
      }).join('');
    }
    function syncTrackOverlay(svg, order) {
      const canvas = svg.closest('.detail-canvas');
      const overlay = canvas ? canvas.querySelector('.detail-track-overlays') : null;
      if (!canvas || !overlay) return;
      const currentOrder = Array.from(overlay.querySelectorAll('.detail-track-control')).map(row => row.dataset.trackId).join('||');
      const nextOrder = detailOrderKey(order);
      if (currentOrder === nextOrder) return;
      const detail = reportData.detail || {};
      const layout = detail.linkview_layout || { svg_height: 400, top_margin: 90 };
      overlay.innerHTML = makeTrackControls(order, Number(layout.svg_height || 400), Number(layout.top_margin || 0));
      wireDetailTrackControls(canvas);
    }
    function showSvgTooltip(content, event, tooltip) {
      tooltip.innerHTML = content;
      tooltip.classList.add('visible');
      updateSvgTooltip(event, tooltip);
    }
    function updateSvgTooltip(event, tooltip) {
      const container = tooltip.parentElement;
      const rect = container.getBoundingClientRect();
      let x = event.clientX - rect.left + 15;
      let y = event.clientY - rect.top - 10;
      if (x + 220 > rect.width) x = event.clientX - rect.left - 220;
      tooltip.style.left = `${Math.max(5, x)}px`;
      tooltip.style.top = `${Math.max(5, y)}px`;
    }
    function hideSvgTooltip(tooltip) {
      tooltip.classList.remove('visible');
    }
    function markerRects(svg, type) {
      const color = type === 'snp' ? 'orange' : 'blue';
      return Array.from(svg.querySelectorAll(`rect[fill="${color}"]`));
    }
    function allMarkerRects(svg) {
      return Array.from(svg.querySelectorAll('rect[fill="orange"], rect[fill="blue"]'));
    }
    function centerOfRect(rect) {
      return {
        x: Number(rect.getAttribute('x') || 0) + Number(rect.getAttribute('width') || 0) / 2,
        y: Number(rect.getAttribute('y') || 0) + Number(rect.getAttribute('height') || 0) / 2,
      };
    }
    function removeHoverLines(svg) {
      svg.querySelectorAll('.variant-hover-line').forEach(line => line.remove());
    }
    function drawHoverLines(svg, rect, type) {
      removeHoverLines(svg);
      if (type !== 'snp') return;
      const variantId = rect.dataset.variantId;
      const point = centerOfRect(rect);
      const candidates = markerRects(svg, type).filter(other => other !== rect && (!variantId || other.dataset.variantId === variantId));
      const related = candidates
        .map(other => ({ rect: other, point: centerOfRect(other) }))
        .filter(item => Math.abs(item.point.x - point.x) <= 8 && Math.abs(item.point.y - point.y) > 4);
      const targets = related.length ? related : candidates
        .map(other => ({ rect: other, point: centerOfRect(other), dist: Math.abs(centerOfRect(other).x - point.x) + Math.abs(centerOfRect(other).y - point.y) }))
        .sort((a, b) => a.dist - b.dist)
        .slice(0, 1);
      targets.forEach(item => {
        const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        line.setAttribute('x1', point.x);
        line.setAttribute('y1', point.y);
        line.setAttribute('x2', item.point.x);
        line.setAttribute('y2', item.point.y);
        line.setAttribute('stroke', 'orange');
        line.setAttribute('stroke-width', '2');
        line.setAttribute('stroke-dasharray', '4,2');
        line.setAttribute('pointer-events', 'none');
        line.classList.add('variant-hover-line');
        svg.appendChild(line);
      });
    }
    function variantTooltip(meta) {
      if (!meta) return '';
      const varType = String(meta.type || '').toUpperCase();
      if (varType === 'SNP') {
        return `<div class="tooltip-title">SNP</div><div class="tooltip-row"><span class="tooltip-label">${t('position')}:</span> ${esc(meta.query_genome_pos || meta.position || '')}</div><div class="tooltip-row"><span class="tooltip-label">${t('mutation')}:</span> (${t('refBase')})<span class="tooltip-seq">${esc(meta.ref || '')}</span> → (${t('qryBase')})<span class="tooltip-seq">${esc(meta.alt || '')}</span></div>`;
      }
      if (varType === 'INS') {
        const seq = meta.alt && meta.alt !== '-' ? meta.alt : '';
        return `<div class="tooltip-title">${t('insertion')}</div><div class="tooltip-row"><span class="tooltip-label">${t('position')}:</span> ${esc(meta.ref_source_pos || meta.position || '')}</div><div class="tooltip-row"><span class="tooltip-label">${t('insertSeq')}:</span> <span class="tooltip-seq">${esc(seq)}</span></div><div class="tooltip-row"><span class="tooltip-label">${t('length')}:</span> ${seq.length || 1} bp</div>`;
      }
      const seq = meta.ref && meta.ref !== '-' ? meta.ref : '';
      return `<div class="tooltip-title">${t('deletion')}</div><div class="tooltip-row"><span class="tooltip-label">${t('position')}:</span> ${esc(meta.query_genome_pos || meta.position || '')}</div><div class="tooltip-row"><span class="tooltip-label">${t('deleteSeq')}:</span> <span class="tooltip-seq">${esc(seq)}</span></div><div class="tooltip-row"><span class="tooltip-label">${t('length')}:</span> ${seq.length || 1} bp</div>`;
    }
    function featureTooltip(type, rect, ranges) {
      const bbox = rect.getBBox ? rect.getBBox() : null;
      const width = bbox ? Math.max(1, bbox.width) : 1;
      const chro = rect.closest('svg') ? Array.from(rect.closest('svg').querySelectorAll('rect.chro')).sort((a, b) => Number(a.getAttribute('y') || 0) - Number(b.getAttribute('y') || 0)) : [];
      const trackIndex = chro.findIndex(track => {
        const ty = Number(track.getAttribute('y') || 0);
        const th = Number(track.getAttribute('height') || 0);
        const ry = Number(rect.getAttribute('y') || (bbox ? bbox.y : 0));
        return Math.abs((ty + th / 2) - (ry + th / 2)) < 24;
      });
      const trackRect = trackIndex >= 0 ? chro[trackIndex] : null;
      const row = ranges[trackIndex >= 0 ? trackIndex : 0];
      let start = row && row.range ? Number(row.range.start || 1) : 1;
      let end = row && row.range ? Number(row.range.end || start) : start;
      if (trackRect && bbox) {
        const tx = Number(trackRect.getAttribute('x') || 0);
        const tw = Number(trackRect.getAttribute('width') || 1);
        const span = Math.max(1, end - start + 1);
        const featureStart = Math.round(start + ((bbox.x - tx) / tw) * span);
        const featureEnd = Math.round(start + (((bbox.x + width) - tx) / tw) * span);
        start = Math.max(1, featureStart);
        end = Math.max(start, featureEnd);
      }
      const title = type === 'utr3' ? t('utr3') : type === 'utr5' ? t('utr5') : t('cds');
      const desc = type === 'utr3' ? t('utr3Desc') : type === 'utr5' ? t('utr5Desc') : t('cdsDesc');
      return `<div class="tooltip-title">${title}</div><div class="tooltip-row">${desc}</div><div class="tooltip-row"><span class="tooltip-label">${t('position')}:</span> ${start} - ${end}</div><div class="tooltip-row"><span class="tooltip-label">${t('length')}:</span> ${Math.max(1, end - start + 1)} bp</div>`;
    }
    function svgViewBox(svg) {
      const vb = svg.viewBox && svg.viewBox.baseVal;
      if (vb && vb.width && vb.height) return vb;
      return { x: 0, y: 0, width: Number(svg.getAttribute('width') || 1200), height: Number(svg.getAttribute('height') || 400) };
    }
    function trackPixelLeft(svg, trackRect, overlay) {
      const svgRect = svg.getBoundingClientRect();
      const overlayRect = overlay.getBoundingClientRect();
      const viewBox = svgViewBox(svg);
      const x = Number(trackRect.getAttribute('x') || 0);
      return svgRect.left - overlayRect.left + ((x - viewBox.x) / viewBox.width) * svgRect.width;
    }
    function ensureDetailGutter(svg, order) {
      const canvas = svg.closest('.detail-canvas');
      const host = canvas ? canvas.querySelector('.detail-linkview-svg') : null;
      const overlay = canvas ? canvas.querySelector('.detail-track-overlays') : null;
      if (!canvas || !host || !overlay) return;
      host.style.paddingLeft = '0px';
      host.style.paddingRight = '0px';
      const chros = Array.from(svg.querySelectorAll('rect.chro')).sort((a, b) => Number(a.getAttribute('y') || 0) - Number(b.getAttribute('y') || 0));
      if (!chros.length) return;
      const controls = Array.from(overlay.querySelectorAll('.detail-track-control'));
      const requiredWidth = Math.max(0, ...controls.map(control => control.getBoundingClientRect().width || 220)) + 24;
      const firstTrackLeft = Math.min(...chros.map(rect => trackPixelLeft(svg, rect, overlay)));
      const deficit = Math.ceil(requiredWidth - firstTrackLeft);
      if (deficit > 0) host.style.paddingLeft = `${deficit}px`;
    }
    function alignTrackControls(svg, order) {
      const canvas = svg.closest('.detail-canvas');
      const overlay = canvas ? canvas.querySelector('.detail-track-overlays') : null;
      if (!canvas || !overlay) return;
      const chros = Array.from(svg.querySelectorAll('rect.chro')).sort((a, b) => Number(a.getAttribute('y') || 0) - Number(b.getAttribute('y') || 0));
      if (!chros.length) return;
      const svgRect = svg.getBoundingClientRect();
      const overlayRect = overlay.getBoundingClientRect();
      const viewBox = svgViewBox(svg);
      chros.forEach((rect, index) => {
        const trackId = order[index];
        const control = Array.from(overlay.querySelectorAll('.detail-track-control')).find(item => item.dataset.trackId === trackId);
        if (!control) return;
        const x = Number(rect.getAttribute('x') || 0);
        const y = Number(rect.getAttribute('y') || 0);
        const height = Number(rect.getAttribute('height') || 0);
        const controlWidth = control.getBoundingClientRect().width || 220;
        const xPx = svgRect.left - overlayRect.left + ((x - viewBox.x) / viewBox.width) * svgRect.width;
        const yPx = svgRect.top - overlayRect.top + (((y + height / 2) - viewBox.y) / viewBox.height) * svgRect.height;
        const left = Math.max(10, Math.min(overlayRect.width - controlWidth - 8, xPx - controlWidth - 14));
        control.style.left = `${left}px`;
        control.style.top = `${Math.max(14, yPx)}px`;
        control.style.transform = 'translateY(-50%)';
      });
    }
    function nudgeScaleBar(svg) {
      svg.querySelectorAll('.detail-scale-shift').forEach(el => el.classList.remove('detail-scale-shift'));
      const texts = Array.from(svg.querySelectorAll('text'));
      const scaleTexts = texts.filter(text => {
        const label = (text.textContent || '').trim();
        const fontSize = Number(String(text.getAttribute('font-size') || '').replace(/[^\d.]/g, '')) || 0;
        return /^\d[\d,.\s]*bp$/i.test(label) && fontSize >= 14;
      });
      scaleTexts.forEach(text => {
        const x = Number(text.getAttribute('x') || 0);
        const y = Number(text.getAttribute('y') || 0);
        const shift = 22;
        text.setAttribute('transform', `${text.getAttribute('transform') || ''} translate(0 ${shift})`.trim());
        text.classList.add('detail-scale-shift');
        svg.querySelectorAll('line,path').forEach(el => {
          const y1 = Number(el.getAttribute('y1') || el.getAttribute('y') || NaN);
          const y2 = Number(el.getAttribute('y2') || el.getAttribute('y') || NaN);
          const x1 = Number(el.getAttribute('x1') || el.getAttribute('x') || NaN);
          const x2 = Number(el.getAttribute('x2') || el.getAttribute('x') || NaN);
          if (Number.isFinite(y1) && Number.isFinite(y2) && Math.abs(Math.max(y1, y2) - y) < 36 && Number.isFinite(x1) && Number.isFinite(x2) && Math.abs(((x1 + x2) / 2) - x) < 180) {
            el.setAttribute('transform', `${el.getAttribute('transform') || ''} translate(0 ${shift})`.trim());
            el.classList.add('detail-scale-shift');
          }
        });
      });
    }
    function addGeneBracket(svg, order, ranges) {
      svg.querySelectorAll('.detail-gene-bracket,.detail-gene-label').forEach(el => el.remove());
      const refIndex = order.indexOf('ref');
      if (refIndex < 0) return;
      const input = reportData.input || {};
      const info = input.extraction_info || {};
      const geneStart = Number(info.gene_rel_start || 0);
      const geneEnd = Number(info.gene_rel_end || 0);
      if (!geneStart || !geneEnd || !input.id) return;
      const chros = Array.from(svg.querySelectorAll('rect.chro')).sort((a, b) => Number(a.getAttribute('y') || 0) - Number(b.getAttribute('y') || 0));
      const rect = chros[refIndex];
      const row = ranges[refIndex];
      if (!rect || !row || !row.range) return;
      const start = Number(row.range.start || 1);
      const end = Number(row.range.end || start);
      const span = Math.max(1, end - start + 1);
      const x = Number(rect.getAttribute('x') || 0);
      const y = Number(rect.getAttribute('y') || 0);
      const width = Number(rect.getAttribute('width') || 0);
      const height = Number(rect.getAttribute('height') || 0);
      const x1 = x + ((Math.max(start, geneStart) - start) / span) * width;
      const x2 = x + ((Math.min(end, geneEnd) - start + 1) / span) * width;
      if (!Number.isFinite(x1) || !Number.isFinite(x2) || Math.abs(x2 - x1) < 2) return;
      const left = Math.min(x1, x2);
      const right = Math.max(x1, x2);
      const bracketY = y - 18;
      const bracketHeight = 14;
      const curl = 10;
      const mid = (left + right) / 2;
      const d = (right - left) > curl * 4
        ? `M ${left},${bracketY} Q ${left},${bracketY - bracketHeight} ${left + curl},${bracketY - bracketHeight} L ${mid - curl},${bracketY - bracketHeight} Q ${mid},${bracketY - bracketHeight} ${mid},${bracketY - bracketHeight - 6} Q ${mid},${bracketY - bracketHeight} ${mid + curl},${bracketY - bracketHeight} L ${right - curl},${bracketY - bracketHeight} Q ${right},${bracketY - bracketHeight} ${right},${bracketY}`
        : `M ${left},${bracketY} Q ${mid},${bracketY - bracketHeight - 6} ${right},${bracketY}`;
      const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      path.setAttribute('d', d);
      path.setAttribute('class', 'detail-gene-bracket');
      svg.appendChild(path);
      const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      label.setAttribute('x', mid);
      label.setAttribute('y', bracketY - bracketHeight - 10);
      label.setAttribute('text-anchor', 'middle');
      label.setAttribute('class', 'detail-gene-label');
      label.textContent = input.id;
      svg.appendChild(label);
    }
    function enhanceLinkviewSvg(svg, assetOrder, tooltip) {
      if (!svg) return;
      const mappedOrder = mappedLinkviewOrder(assetOrder);
      const effectiveOrder = inferSvgTrackOrder(svg, mappedOrder);
      syncTrackOverlay(svg, effectiveOrder);
      svg.classList.add('detail-linkview-svg-root');
      svg.setAttribute('width', '100%');
      svg.style.overflow = 'visible';
      svg.querySelectorAll('text.label').forEach(label => label.remove());
      svg.querySelectorAll('.detail-track-end-label').forEach(label => label.remove());
      const ranges = currentTrackRanges(effectiveOrder);
      const chros = Array.from(svg.querySelectorAll('rect.chro')).sort((a, b) => Number(a.getAttribute('y') || 0) - Number(b.getAttribute('y') || 0));
      chros.forEach((rect, index) => {
        const row = ranges[index];
        if (!row) return;
        const x = Number(rect.getAttribute('x') || 0);
        const y = Number(rect.getAttribute('y') || 0);
        const width = Number(rect.getAttribute('width') || 0);
        const height = Number(rect.getAttribute('height') || 0);
        const labelY = y + height + 17;
        [['start', x + 2, 'start'], ['end', x + width - 2, 'end']].forEach(([kind, labelX, anchor]) => {
          const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
          text.setAttribute('x', labelX);
          text.setAttribute('y', labelY);
          text.setAttribute('text-anchor', anchor);
          text.setAttribute('class', 'detail-track-end-label');
          text.textContent = Number(row.range[kind] || 0).toLocaleString();
          svg.appendChild(text);
        });
      });
      addGeneBracket(svg, effectiveOrder, ranges);
      nudgeScaleBar(svg);
      requestAnimationFrame(() => {
        ensureDetailGutter(svg, effectiveOrder);
        requestAnimationFrame(() => alignTrackControls(svg, effectiveOrder));
      });
      const markers = resolveLinkviewMarkers(assetOrder);
      allMarkerRects(svg).forEach((rect, index) => {
        const meta = markers[index] || {};
        rect.dataset.variantId = meta.variant_id || '';
        rect.dataset.variantType = String(meta.type || '').toLowerCase();
        rect.dataset.markerIndex = String(index);
      });
      [['.UTR3', 'utr3'], ['.UTR5', 'utr5'], ['.exon', 'cds']].forEach(([selector, type]) => {
        svg.querySelectorAll(selector).forEach(el => {
          el.style.cursor = 'pointer';
          el.addEventListener('mouseenter', event => showSvgTooltip(featureTooltip(type, el, ranges), event, tooltip));
          el.addEventListener('mousemove', event => updateSvgTooltip(event, tooltip));
          el.addEventListener('mouseleave', () => hideSvgTooltip(tooltip));
        });
      });
      svg.querySelectorAll('.legend-item').forEach(item => {
        item.addEventListener('click', () => {
          const type = item.getAttribute('data-type');
          item.classList.toggle('disabled');
          const visible = !item.classList.contains('disabled');
          let elements = [];
          if (type === 'snp') elements = markerRects(svg, 'snp');
          if (type === 'indel') elements = markerRects(svg, 'indel');
          if (type === 'utr5') elements = Array.from(svg.querySelectorAll('.UTR5'));
          if (type === 'utr3') elements = Array.from(svg.querySelectorAll('.UTR3'));
          if (type === 'cds') elements = Array.from(svg.querySelectorAll('.exon'));
          elements.forEach(el => { el.style.display = visible ? '' : 'none'; });
        });
      });
      [['snp', 'SNP'], ['indel', 'Indel']].forEach(([type, label]) => {
        markerRects(svg, type).forEach(rect => {
          rect.style.cursor = 'pointer';
          rect.addEventListener('mouseenter', event => {
            drawHoverLines(svg, rect, type);
            const meta = markers[Number(rect.dataset.markerIndex || -1)];
            showSvgTooltip(variantTooltip(meta) || `<div class="tooltip-title">${label}</div>`, event, tooltip);
          });
          rect.addEventListener('mousemove', event => updateSvgTooltip(event, tooltip));
          rect.addEventListener('mouseleave', () => {
            hideSvgTooltip(tooltip);
            removeHoverLines(svg);
          });
        });
      });
    }
    function renderZoomModalBody() {
      const modal = document.getElementById('zoom-modal');
      const body = document.getElementById('zoom-modal-body');
      if (!modal || !body) return false;
      const requestedOrder = state.order.length ? state.order.slice() : ['ref'].concat(queryIds());
      const assetOrder = resolvedLinkviewOrder(requestedOrder);
      const displayOrder = mappedLinkviewOrder(assetOrder);
      const detail = reportData.detail || {};
      const layout = detail.linkview_layout || { svg_height: 400, top_margin: 90 };
      const svgHeight = Number(layout.svg_height || 400);
      const topMargin = Number(layout.top_margin || 0);
      const inlineSvg = resolveLinkviewInlineSvg(assetOrder);
      const currentSvg = document.querySelector('#detail .detail-linkview-svg svg');
      if (!inlineSvg && !currentSvg) return false;
      body.innerHTML = `<div class="detail-canvas zoom-detail-canvas"><div id="modal-tooltip" class="svg-tooltip"></div><div class="detail-linkview-svg"></div><div class="detail-track-overlays">${makeTrackControls(displayOrder, svgHeight, topMargin)}</div></div>`;
      const host = body.querySelector('.detail-linkview-svg');
      if (inlineSvg) {
        host.innerHTML = inlineSvg;
      } else {
        const clone = currentSvg.cloneNode(true);
        clone.id = 'modal-svg';
        host.appendChild(clone);
      }
      const svg = host.querySelector('svg');
      if (!svg) return false;
      svg.style.width = '100%';
      svg.style.minWidth = '1200px';
      enhanceLinkviewSvg(svg, assetOrder, body.querySelector('#modal-tooltip'));
      wireDetailTrackControls(body);
      return true;
    }
    function openZoomModal() {
      const modal = document.getElementById('zoom-modal');
      if (!modal) return;
      modal.classList.add('visible');
      document.body.style.overflow = 'hidden';
      if (!renderZoomModalBody()) closeZoomModal();
    }
    function closeZoomModal() {
      const modal = document.getElementById('zoom-modal');
      if (!modal) return;
      modal.classList.remove('visible');
      document.body.style.overflow = '';
    }
    function renderLinkviewDetail() {
      const root = document.getElementById('detail');
      const requestedOrder = state.order.length ? state.order.slice() : ['ref'].concat(queryIds());
      const assetOrder = resolvedLinkviewOrder(requestedOrder);
      const displayOrder = mappedLinkviewOrder(assetOrder);
      const svgPath = resolveLinkviewSvg(assetOrder);
      if (!svgPath) {
        renderVectorDetail();
        return;
      }
      const detail = reportData.detail || {};
      const layout = detail.linkview_layout || { svg_height: 400, top_margin: 90 };
      const svgHeight = Number(layout.svg_height || 400);
      const topMargin = Number(layout.top_margin || 0);
      const token = ++detailRenderToken;
      root.innerHTML = `<div class="detail-canvas"><div id="svg-tooltip" class="svg-tooltip"></div><div class="detail-linkview-svg">Loading...</div><div class="detail-track-overlays">${makeTrackControls(displayOrder, svgHeight, topMargin)}</div><button id="zoom-btn" class="zoom-btn" title="Zoom"><svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="11" cy="11" r="8"></circle><path d="M21 21l-4.35-4.35"></path><path d="M11 8v6M8 11h6"></path></svg></button></div>`;
      wireDetailTrackControls();
      const zoomBtn = document.getElementById('zoom-btn');
      if (zoomBtn) zoomBtn.onclick = openZoomModal;
      const modal = document.getElementById('zoom-modal');
      const closeBtn = document.getElementById('zoom-modal-close');
      if (closeBtn) closeBtn.onclick = closeZoomModal;
      if (modal) modal.onclick = event => { if (event.target === modal) closeZoomModal(); };
      const inlineSvg = resolveLinkviewInlineSvg(assetOrder);
      if (inlineSvg) {
        const host = root.querySelector('.detail-linkview-svg');
        host.innerHTML = inlineSvg;
        enhanceLinkviewSvg(host.querySelector('svg'), assetOrder, root.querySelector('#svg-tooltip'));
        return;
      }
      fetch(svgPath)
        .then(response => response.ok ? response.text() : Promise.reject(new Error('SVG fetch failed')))
        .then(svgText => {
          if (token !== detailRenderToken) return;
          const host = root.querySelector('.detail-linkview-svg');
          host.innerHTML = svgText;
          enhanceLinkviewSvg(host.querySelector('svg'), assetOrder, root.querySelector('#svg-tooltip'));
        })
        .catch(() => {
          if (token !== detailRenderToken) return;
          const host = root.querySelector('.detail-linkview-svg');
          host.innerHTML = `<img class="detail-linkview-img" src="${esc(svgPath)}" alt="LINKVIEW local alignment">`;
        });
    }
    function renderVectorDetail() {
      const root = document.getElementById('detail');
      const order = state.order.length ? state.order.slice() : ['ref'].concat(queryIds());
      const rows = order.map(trackId => {
        const track = trackById(trackId);
        const candidate = trackId === 'ref' ? null : candidateById(trackId, state.selected[trackId]);
        const range = trackId === 'ref' ? { start: 1, end: refLength() } : candidateBounds(candidate);
        return { trackId, track, candidate, range };
      }).filter(row => row.trackId === 'ref' || row.candidate);
      if (!rows.length) { root.innerHTML = `<div class="detail-empty">${t('noTracks')}</div>`; return; }

      const width = 1200;
      const plotLeft = 170;
      const plotWidth = 980;
      const rowStep = 72;
      const top = 52;
      const height = top + rows.length * rowStep + 36;
      const rowY = row => top + rows.indexOf(row) * rowStep;
      const refRow = rows.find(r => r.trackId === 'ref');
      const fragments = [];
      const ribbons = [];
      const markers = [];
      const overlays = [];

      for (const row of rows) {
        const y = rowY(row);
        const subtitle = row.trackId === 'ref'
          ? `${t('reference')} · ${row.range.start}-${row.range.end}`
          : `${row.candidate.candidate_id} · ${row.candidate.query_chr}:${row.range.start}-${row.range.end}`;
        const overlayTop = ((y / height) * 100).toFixed(3);
        const overlayWidth = `calc(${((plotLeft / width) * 100).toFixed(3)}% - 16px)`;
        const trackTitle = esc(row.track.name || row.trackId);
        const handle = `<span class="detail-drag-handle" title="${t('dragTrack')}"></span>`;
        if (row.trackId === 'ref') {
          overlays.push(`<div class="detail-track-control" data-track-id="${esc(row.trackId)}" title="${trackTitle}" style="left:8px; top:${overlayTop}%; width:${overlayWidth}; transform:translateY(-50%);">${handle}<span class="detail-track-name">${trackTitle}</span></div>`);
        } else {
          const options = (reportData.candidates[row.trackId] || [])
            .map(c => `<option value="${esc(c.candidate_id)}" ${state.selected[row.trackId] === c.candidate_id ? 'selected' : ''}>${esc(c.candidate_id)}</option>`)
            .join('');
          overlays.push(`<div class="detail-track-control" data-track-id="${esc(row.trackId)}" title="${trackTitle}" style="left:8px; top:${overlayTop}%; width:${overlayWidth}; transform:translateY(-50%);">${handle}<span class="detail-track-name">${trackTitle}</span><select class="detail-track-select" title="${t('candidateSwitch')}">${options}</select></div>`);
        }
        fragments.push(`<text class="detail-subtext" x="${plotLeft}" y="${y + 27}">${esc(subtitle)}</text>`);
        fragments.push(`<rect class="detail-chro" x="${plotLeft}" y="${y - 8}" width="${plotWidth}" height="16" rx="0"></rect>`);
        fragments.push(`<text class="detail-subtext" x="${plotLeft}" y="${y + 16}">${Number(row.range.start).toLocaleString()}</text>`);
        fragments.push(`<text class="detail-subtext" x="${plotLeft + plotWidth}" y="${y + 16}" text-anchor="end">${Number(row.range.end).toLocaleString()}</text>`);

        if (row.trackId === 'ref') {
          const input = reportData.input || {};
          const info = input.extraction_info || {};
          const geneStart = Number(info.gene_rel_start || 1);
          const geneEnd = Number(info.gene_rel_end || input.sequence_length || row.range.end);
          const x1 = xFor(geneStart, row.range, plotLeft, plotWidth);
          const x2 = xFor(geneEnd, row.range, plotLeft, plotWidth);
          fragments.push(`<rect class="detail-ref-gene" x="${Math.min(x1, x2)}" y="${y - 8}" width="${Math.max(3, Math.abs(x2 - x1))}" height="16" rx="3"></rect>`);
        } else {
          for (const block of row.candidate.blocks || []) {
            const [qs, qe] = normRange(block.query_start, block.query_end);
            const x1 = xFor(qs, row.range, plotLeft, plotWidth);
            const x2 = xFor(qe, row.range, plotLeft, plotWidth);
            fragments.push(`<rect class="detail-match" x="${Math.min(x1, x2)}" y="${y - 8}" width="${Math.max(3, Math.abs(x2 - x1))}" height="16" rx="3"></rect>`);
            if (refRow) {
              const refY = rowY(refRow);
              const rx1 = xFor(block.ref_start, refRow.range, plotLeft, plotWidth);
              const rx2 = xFor(block.ref_end, refRow.range, plotLeft, plotWidth);
              const qx1 = xFor(block.query_start, row.range, plotLeft, plotWidth);
              const qx2 = xFor(block.query_end, row.range, plotLeft, plotWidth);
              const reverse = (Number(block.ref_start) - Number(block.ref_end)) * (Number(block.query_start) - Number(block.query_end)) < 0;
              ribbons.push(`<polygon class="detail-block${reverse ? ' reverse' : ''}" points="${rx1},${refY + 10} ${rx2},${refY + 10} ${qx2},${y - 10} ${qx1},${y - 10}"></polygon>`);
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
            markers.push(`<rect class="detail-variant ${vClass}" x="${x - 2}" y="${y - 24}" width="4" height="18"><title>${esc(variant.type)} ${esc(row.track.name)}:${qPos}</title></rect>`);
          }
          if (Number.isFinite(rPos) && rPos >= refRow.range.start && rPos <= refRow.range.end) {
            const x = xFor(rPos, refRow.range, plotLeft, plotWidth);
            markers.push(`<rect class="detail-variant ${vClass}" x="${x - 2}" y="${refY - 24}" width="4" height="18"><title>${esc(variant.type)} ref:${rPos}</title></rect>`);
          }
        }
      }

      root.innerHTML = `<div class="detail-canvas"><svg class="detail-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="multi-track local alignment">${ribbons.join('')}${fragments.join('')}${markers.join('')}</svg><div class="detail-track-overlays">${overlays.join('')}</div></div>`;
      wireDetailTrackControls();
    }
    function renderOutputs() {
      const rows = [`<tr><th>${t('filename')}</th><th>${t('description')}</th><th>${t('path')}</th></tr>`];
      rows.push(`<tr><td><a href="data.json">data.json</a></td><td>${t('structuredData')}</td><td><code class="path-code">report/data.json</code></td></tr>`);
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
