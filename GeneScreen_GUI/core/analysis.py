#!/usr/bin/env python3
"""
GeneScreen 1.0 - 序列分析模块

整合自 CLI 版本:
- SequenceExtractor: 使用 pyfaidx 从参考基因组提取序列
- BlastAligner: 使用 BLAST+ 进行序列比对和变异检测

改动（相比 CLI 版本）：
- 移除 CLI 相关代码
- 添加类型注解
- 使用 Database 记录分析历史
"""

import os
import re
import hashlib
import shutil
import tempfile
from itertools import combinations, product
from pathlib import Path
from typing import Optional, Dict, List, Any, Tuple

from pyfaidx import Fasta
from Bio.Blast import NCBIXML

from utils import run_subprocess

try:
    from .multi_query_result import parse_coords_candidates
except ImportError:
    from core.multi_query_result import parse_coords_candidates


def run_cmd(cmd: str, error_msg: str = "命令执行失败") -> bool:
    """执行 shell 命令（Windows 下隐藏黑窗口）"""
    print(f"[CMD] {cmd}")
    result = run_subprocess(cmd, shell=True, capture_output=False)
    if result.returncode != 0:
        print(f"[ERROR] {error_msg}")
        return False
    return True


def _blast_safe_fasta_path(fasta_path: str) -> str:
    """Return a whitespace-free FASTA path for BLAST tools that misparse spaces."""
    if not re.search(r"\s", fasta_path):
        return fasta_path
    source = Path(fasta_path)
    stat = source.stat()
    key_src = f"{source.resolve()}|{stat.st_mtime_ns}|{stat.st_size}"
    digest = hashlib.sha1(key_src.encode("utf-8")).hexdigest()[:16]
    safe_dir = Path(tempfile.gettempdir()) / "genescreen_blast_inputs" / digest
    safe_dir.mkdir(parents=True, exist_ok=True)
    suffix = "".join(source.suffixes) or ".fasta"
    safe_path = safe_dir / f"input{suffix}"
    if not safe_path.exists() or safe_path.stat().st_size != stat.st_size:
        shutil.copyfile(source, safe_path)
    return str(safe_path)


def ensure_blast_db(fasta_path: str) -> Optional[str]:
    """确保 BLAST 数据库存在，返回可传给 blastn -db 的数据库前缀。"""
    db_prefix = _blast_safe_fasta_path(fasta_path)
    db_file = f"{db_prefix}.nin"
    if not os.path.exists(db_file):
        print(f"[INFO] 创建 BLAST 数据库: {fasta_path}")
        cmd = f'makeblastdb -in "{db_prefix}" -dbtype nucl -out "{db_prefix}"'
        if not run_cmd(cmd, "创建 BLAST 数据库失败"):
            return None
    return db_prefix


def sanitize_path_segment(text: Any) -> str:
    """清理用于输出子目录名的输入 ID。"""
    return re.sub(r'[<>:"/\\|?*]', "_", str(text)).strip()


def _normalize_query_entries(
    query_genome: Optional[str],
    qry_name: str = "",
    qry_gff: Optional[str] = None,
    query_genomes: Optional[List[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    entries: List[Dict[str, Any]] = []
    if query_genomes:
        for index, genome in enumerate(query_genomes, start=1):
            fasta = genome.get("fasta_path") or genome.get("fasta") or genome.get("path")
            if not fasta:
                continue
            name = genome.get("name") or genome.get("display_name") or f"query_{index}"
            entries.append({
                "name": name,
                "fasta": fasta,
                "gff": genome.get("annotation_path") or genome.get("gff"),
            })
    elif query_genome:
        entries.append({
            "name": qry_name or "query",
            "fasta": query_genome,
            "gff": qry_gff,
        })
    return entries


def _extract_candidate_fasta(
    query_entry: Dict[str, Any],
    candidate: Dict[str, Any],
    output_dir: str,
    upstream: int = 0,
    downstream: int = 0,
) -> Optional[str]:
    chrom = candidate.get("query_chr")
    start = int(candidate.get("query_start") or 0)
    end = int(candidate.get("query_end") or 0)
    fasta_path = query_entry.get("fasta")
    if not chrom or not fasta_path or start <= 0 or end <= 0:
        return None
    safe_query = sanitize_path_segment(query_entry.get("name") or "query") or "query"
    candidate_id = sanitize_path_segment(candidate.get("candidate_id") or "candidate") or "candidate"
    out_dir = os.path.join(output_dir, "queries", safe_query, "candidates")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{candidate_id}.fasta")
    left, right = min(start, end), max(start, end)
    strand = candidate.get("strand", "+")
    if strand == "-":
        extract_left = max(1, left - downstream)
        extract_right = right + upstream
    else:
        extract_left = max(1, left - upstream)
        extract_right = right + downstream
    try:
        fasta = Fasta(fasta_path)
        chrom_len = len(fasta[chrom])
        extract_right = min(chrom_len, extract_right)
        seq = fasta[chrom][extract_left - 1:extract_right].seq
        fasta.close()
    except Exception as exc:
        print(f"[WARNING] 无法提取候选序列 {query_entry.get('name')}:{chrom}:{extract_left}-{extract_right}: {exc}")
        return None
    with open(out_path, "w", encoding="utf-8") as handle:
        seq_id = sanitize_path_segment(candidate_id) or "candidate"
        handle.write(
            f">{seq_id} query={safe_query} region={chrom}:{extract_left}-{extract_right} "
            f"match={left}-{right} strand={strand} upstream={upstream} downstream={downstream}\n"
        )
        handle.write(f"{seq}\n")
    return out_path


def _precompute_pairwise(
    query_results: List[Dict[str, Any]],
    query_entries: List[Dict[str, Any]],
    aligner: "BlastAligner",
    output_dir: str,
    identity: float,
    candidate_limit: Optional[int] = 3,
    pairwise_all: bool = False,
    merge_gap: int = 1000,
    min_aln_len: int = 0,
    upstream: int = 0,
    downstream: int = 0,
) -> List[Dict[str, Any]]:
    bundles = []
    for index, query_result in enumerate(query_results):
        if index >= len(query_entries):
            continue
        entry = query_entries[index]
        candidates = parse_coords_candidates(
            query_result.get("coords"),
            merge_gap=merge_gap,
            min_aln_len=min_aln_len,
        )
        if not pairwise_all:
            candidates = candidates[:candidate_limit or 3]
        prepared = []
        for candidate in candidates:
            fasta = _extract_candidate_fasta(entry, candidate, output_dir, upstream, downstream)
            if fasta:
                prepared.append({"candidate": candidate, "fasta": fasta})
        bundles.append({"entry": entry, "safe": sanitize_path_segment(entry.get("name") or "query") or "query", "candidates": prepared})

    pairwise_results = []
    original_output_dir = aligner.output_dir
    for left, right in combinations(bundles, 2):
        if not left["candidates"] or not right["candidates"]:
            continue
        pair_safe = f"{left['safe']}__{right['safe']}"
        pair_dir = os.path.join(output_dir, "pairwise", pair_safe)
        os.makedirs(pair_dir, exist_ok=True)
        aligner.output_dir = pair_dir
        for left_candidate, right_candidate in product(left["candidates"], right["candidates"]):
            left_id = sanitize_path_segment(left_candidate["candidate"].get("candidate_id") or "candidate")
            right_id = sanitize_path_segment(right_candidate["candidate"].get("candidate_id") or "candidate")
            prefix = f"{left['safe']}_{left_id}__{right['safe']}_{right_id}"
            aligned = aligner.align(right_candidate["fasta"], left_candidate["fasta"], prefix, identity)
            if not aligned:
                continue
            aligned.update({
                "pair_id": pair_safe,
                "left_query": left["entry"].get("name"),
                "right_query": right["entry"].get("name"),
                "left_candidate_id": left_candidate["candidate"].get("candidate_id"),
                "right_candidate_id": right_candidate["candidate"].get("candidate_id"),
                "left_candidate_fasta": left_candidate["fasta"],
                "right_candidate_fasta": right_candidate["fasta"],
                "pair_dir": pair_dir,
            })
            pairwise_results.append(aligned)
    aligner.output_dir = original_output_dir
    return pairwise_results


class SequenceExtractor:
    """使用 pyfaidx 从参考基因组提取序列"""

    def __init__(self, genome: str, annotation: Optional[str], output_dir: str):
        self.genome = genome
        self.annotation = annotation
        self.output_dir = output_dir
        self._fasta: Optional[Fasta] = None

    @property
    def fasta(self) -> Fasta:
        """懒加载 FASTA 文件"""
        if self._fasta is None:
            self._fasta = Fasta(self.genome)
        return self._fasta

    def extract_by_gene_id(
        self,
        gene_id: str,
        upstream: int = 0,
        downstream: int = 0
    ) -> Tuple[Optional[str], Optional[str], Optional[Dict[str, Any]]]:
        """通过 Gene ID 提取序列，支持模糊匹配和上下游延伸"""
        print(f"[INFO] 提取基因 {gene_id} 的序列...")
        if upstream > 0 or downstream > 0:
            print(f"[INFO] 上游延伸: {upstream} bp, 下游延伸: {downstream} bp")

        original_gff = os.path.join(self.output_dir, f"{gene_id}.original.gff3")
        linkview_gff = os.path.join(self.output_dir, f"{gene_id}.linkview.gff3")

        if not self.annotation or not os.path.exists(self.annotation):
            print(f"[WARNING] 参考基因组无注释文件，Gene ID 模式需要 GFF 注释来定位基因")
            return None, None, None

        # 使用 ripgrep 过滤 GFF（模糊匹配）
        cmd = f'rg "{gene_id}" "{self.annotation}" > "{original_gff}"'
        if not run_cmd(cmd, f"提取 {gene_id} 注释失败"):
            self._grep_fallback(gene_id, self.annotation, original_gff)

        # 检查是否找到匹配
        if not os.path.exists(original_gff) or os.path.getsize(original_gff) == 0:
            print(f"[INFO] 尝试模糊匹配 {gene_id}...")
            self._fuzzy_grep(gene_id, self.annotation, original_gff)

        gene_info = self._parse_gff_info(original_gff, gene_id)
        if not gene_info:
            print(f"[ERROR] 未找到 {gene_id} 的基因记录")
            return None, None, None

        gene_start = gene_info["start"]
        gene_end = gene_info["end"]
        strand = gene_info["strand"]
        chrom = gene_info["chrom"]

        if strand == "+":
            theory_start = gene_start - upstream
            theory_end = gene_end + downstream
        else:
            theory_start = gene_start - downstream
            theory_end = gene_end + upstream

        chrom_normalized = self._normalize_chrom(chrom)
        chrom_length = len(self.fasta[chrom_normalized])
        actual_start = max(1, theory_start)
        actual_end = min(chrom_length, theory_end)

        if actual_start != theory_start or actual_end != theory_end:
            print(f"[INFO] 边界裁剪: [{theory_start}, {theory_end}] -> [{actual_start}, {actual_end}]")

        gene_rel_start = gene_start - actual_start + 1
        gene_rel_end = gene_end - actual_start + 1

        self._generate_linkview_gff(gene_id, chrom_normalized, actual_start, actual_end, linkview_gff)
        fasta_file = self._extract_fasta_region(gene_id, chrom_normalized, actual_start, actual_end)

        extraction_info = {
            "gene_start": gene_start,
            "gene_end": gene_end,
            "gene_rel_start": gene_rel_start,
            "gene_rel_end": gene_rel_end,
            "strand": strand,
            "actual_start": actual_start,
            "actual_end": actual_end,
            "upstream": upstream,
            "downstream": downstream,
            "chrom": chrom,
        }

        return fasta_file, linkview_gff, extraction_info

    def _parse_gff_info(self, original_gff: str, gene_id: str) -> Optional[Dict[str, Any]]:
        """解析 GFF 获取基因信息（chrom, start, end, strand）"""
        gene_info = None
        mrna_records = []

        def id_matches(attributes: str, target_id: str) -> bool:
            match_keys = ("ID", "Parent", "Name", "gene_id", "locus_tag")
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    key = key.strip()
                    value = value.strip()
                    if key in match_keys:
                        pattern = r'(^|[,;._:])' + re.escape(target_id) + r'($|[,;._:])'
                        if re.search(pattern, value):
                            return True
            return False

        with open(original_gff, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue

                feature_type = fields[2]
                attributes = fields[8]

                if not id_matches(attributes, gene_id):
                    continue

                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6] if fields[6] in ["+", "-"] else "+"

                if feature_type == "gene":
                    gene_info = {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                    }
                    break
                if feature_type == "mRNA":
                    mrna_records.append({
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "range": end - start,
                    })

        if not gene_info and mrna_records:
            best_mrna = max(mrna_records, key=lambda x: x["range"])
            gene_info = {
                "chrom": best_mrna["chrom"],
                "start": best_mrna["start"],
                "end": best_mrna["end"],
                "strand": best_mrna["strand"],
            }
            print(f"[INFO] 未找到 gene 记录，使用 mRNA 记录作为兜底")

        return gene_info

    def _generate_linkview_gff(
        self,
        gene_id: str,
        chrom: str,
        actual_start: int,
        actual_end: int,
        linkview_gff: str
    ):
        """从全量 GFF 按坐标区间提取 feature，生成相对坐标的 linkview.gff3"""
        region_len = actual_end - actual_start + 1
        allowed_types = {
            "gene", "mrna", "transcript", "exon", "cds",
            "five_prime_utr", "three_prime_utr", "intron"
        }

        def chrom_matches(feat_chrom: str, target_chrom: str) -> bool:
            if feat_chrom == target_chrom:
                return True
            feat_num = re.search(r'(\d+)', feat_chrom)
            target_num = re.search(r'(\d+)', target_chrom)
            if feat_num and target_num:
                return int(feat_num.group(1)) == int(target_num.group(1))
            return False

        def feature_belongs_to_gene(attributes: str, target_gene_id: str) -> bool:
            id_match = re.search(r'ID=([^;]+)', attributes)
            if id_match and target_gene_id in id_match.group(1):
                return True
            parent_match = re.search(r'Parent=([^;]+)', attributes)
            if parent_match and target_gene_id in parent_match.group(1):
                return True
            gene_id_match = re.search(r'gene_id=([^;]+)', attributes)
            if gene_id_match and target_gene_id in gene_id_match.group(1):
                return True
            return False

        with open(self.annotation, "r", encoding="utf-8", errors="ignore") as f_in, \
             open(linkview_gff, "w") as f_out:
            for line in f_in:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue

                feat_chrom = fields[0]
                if not chrom_matches(feat_chrom, chrom):
                    continue

                feat_type = fields[2]
                if feat_type.lower() not in allowed_types:
                    continue

                attributes = fields[8]
                if not feature_belongs_to_gene(attributes, gene_id):
                    continue

                feat_start = int(fields[3])
                feat_end = int(fields[4])
                if feat_end < actual_start or feat_start > actual_end:
                    continue

                clipped_start = max(feat_start, actual_start)
                clipped_end = min(feat_end, actual_end)
                rel_start = clipped_start - actual_start + 1
                rel_end = clipped_end - actual_start + 1
                rel_start = max(1, min(rel_start, region_len))
                rel_end = max(1, min(rel_end, region_len))

                f_out.write(
                    f"{gene_id}\t{fields[1]}\t{fields[2]}\t{rel_start}\t{rel_end}\t"
                    f"{fields[5]}\t{fields[6]}\t{fields[7]}\t{fields[8]}\n"
                )

    def _extract_fasta_region(self, gene_id: str, chrom: str, start: int, end: int) -> Optional[str]:
        """提取指定区域的序列"""
        fasta_file = os.path.join(self.output_dir, f"{gene_id}.fasta")
        try:
            seq = self.fasta[chrom][start - 1 : end]
            with open(fasta_file, "w") as f:
                f.write(f">{gene_id}\n{str(seq)}\n")
            print(f"[INFO] 已生成序列文件: {fasta_file} (长度: {len(str(seq))} bp)")
            return fasta_file
        except Exception as e:
            print(f"[ERROR] 提取序列失败: {e}")
            return None

    def _normalize_chrom(self, chrom: str) -> str:
        """标准化染色体名称，自动匹配基因组中的实际名称"""
        available_chroms = list(self.fasta.keys())
        
        if chrom in available_chroms:
            return chrom
        
        num_match = re.search(r'(\d+)', chrom)
        if not num_match:
            return chrom
        
        num = num_match.group(1)
        num_int = int(num)
        
        candidates = [
            str(num_int), num,
            f"Chr{num_int}", f"chr{num_int}",
            f"Chr{num.zfill(2)}", f"chr{num.zfill(2)}",
            f"Chr{num}", f"chr{num}",
        ]
        
        for candidate in candidates:
            if candidate in available_chroms:
                if candidate != chrom:
                    print(f"[INFO] 染色体名称转换: {chrom} -> {candidate}")
                return candidate
        
        return chrom

    def extract_by_location(
        self, chrom: str, start: int, end: int, name: Optional[str] = None
    ) -> Optional[str]:
        """通过染色体坐标提取序列"""
        loc_name = name or f"{chrom}_{start}_{end}"
        print(f"[INFO] 提取区域 {chrom}:{start}-{end} 的序列...")

        chrom = self._normalize_chrom(chrom)
        fasta_file = os.path.join(self.output_dir, f"{loc_name}.fasta")

        try:
            seq = self.fasta[chrom][start - 1 : end]
            with open(fasta_file, "w") as f:
                f.write(f">{loc_name}\n{str(seq)}\n")
            print(f"[INFO] 已生成序列文件: {fasta_file}")
            return fasta_file
        except Exception as e:
            print(f"[ERROR] 提取序列失败: {e}")
            return None

    def _grep_fallback(self, pattern: str, input_file: str, output_file: str):
        """纯 Python grep 回退"""
        with open(input_file, "r", encoding="utf-8", errors="ignore") as f_in:
            with open(output_file, "w") as f_out:
                for line in f_in:
                    if pattern in line:
                        f_out.write(line)

    def _fuzzy_grep(self, gene_id: str, input_file: str, output_file: str):
        """模糊匹配 Gene ID"""
        core_match = re.search(r'(Os\d+g\d+)', gene_id, re.IGNORECASE)
        if not core_match:
            return
        
        core_id = core_match.group(1)
        print(f"[INFO] 使用核心 ID 搜索: {core_id}")
        
        with open(input_file, "r", encoding="utf-8", errors="ignore") as f_in:
            with open(output_file, "w") as f_out:
                for line in f_in:
                    if core_id.lower() in line.lower():
                        f_out.write(line)

    def _parse_gff(
        self, original_gff: str, linkview_gff: str, gene_id: str
    ) -> Optional[Dict[str, Any]]:
        """解析 GFF 文件"""
        start_offset = None
        gene_info = {}

        with open(original_gff, "r") as f_in, open(linkview_gff, "w") as f_out:
            for line in f_in:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue

                if start_offset is None:
                    start_offset = int(fields[3])

                rel_start = int(fields[3]) - start_offset + 1
                rel_end = int(fields[4]) - start_offset + 1
                f_out.write(
                    f"{gene_id}\t{fields[1]}\t{fields[2]}\t{rel_start}\t{rel_end}\t"
                    f"{fields[5]}\t{fields[6]}\t{fields[7]}\t{fields[8]}\n"
                )

                if fields[2] == "gene":
                    gene_info = {
                        "chrom": fields[0],
                        "start": int(fields[3]),
                        "end": int(fields[4]),
                    }

        if not gene_info:
            print(f"[ERROR] 未找到 {gene_id} 的 gene 记录")
            return None
        return gene_info

    def _extract_fasta(self, gene_id: str, gene_info: Dict[str, Any]) -> Optional[str]:
        """使用 pyfaidx 提取序列"""
        fasta_file = os.path.join(self.output_dir, f"{gene_id}.fasta")

        try:
            chrom = gene_info["chrom"]
            start = gene_info["start"]
            end = gene_info["end"]
            seq = self.fasta[chrom][start - 1 : end]

            with open(fasta_file, "w") as f:
                f.write(f">{gene_id}\n{str(seq)}\n")

            print(f"[INFO] 已生成序列文件: {fasta_file}")
            return fasta_file
        except Exception as e:
            print(f"[ERROR] 提取序列失败: {e}")
            return None


class BlastAligner:
    """使用 BLAST+ 进行序列比对"""

    def __init__(self, output_dir: str):
        self.output_dir = output_dir

    def align(
        self,
        reference: str,
        query: str,
        prefix_name: str,
        identity_threshold: float = 90
    ) -> Optional[Dict[str, str]]:
        """执行 BLAST 比对"""
        # 清理 prefix_name 中的非法路径字符
        import re
        safe_prefix_name = re.sub(r'[<>:"/\\|?*]', "_", prefix_name)
        prefix = os.path.join(self.output_dir, safe_prefix_name)

        db_prefix = ensure_blast_db(reference)
        if not db_prefix:
            return None

        xml_file = f"{prefix}.blast.xml"
        print(f"[INFO] 运行 BLAST 比对...")
        cmd = (
            f'blastn -query "{query}" -db "{db_prefix}" '
            f'-out "{xml_file}" -outfmt 5 '
            f"-perc_identity {identity_threshold} "
            f"-num_threads 2"
        )
        if not run_cmd(cmd, "BLAST 比对失败"):
            return None

        coords_file = f"{prefix}.coords"
        snps_file = f"{prefix}.snps"
        self._parse_blast_xml(xml_file, coords_file, snps_file)

        print(f"[INFO] 比对完成，结果文件: {coords_file}, {snps_file}")
        return {
            "blast_xml": xml_file,
            "coords": coords_file,
            "snps": snps_file,
            "prefix": prefix,
        }

    def _parse_blast_xml(self, xml_file: str, coords_file: str, snps_file: str):
        """解析 BLAST XML 输出"""
        with open(xml_file) as f:
            blast_records = NCBIXML.parse(f)

            with open(coords_file, "w") as f_coords, open(snps_file, "w") as f_snps:
                f_coords.write(
                    "[S1]\t[E1]\t[S2]\t[E2]\t[LEN1]\t[LEN2]\t[%IDY]\t[REF]\t[QUERY]\n"
                )
                f_snps.write("[QRY_POS]\t[REF]\t[ALT]\t[REF_POS]\t[TYPE]\t[QRY_NAME]\t[REF_NAME]\n")

                for record in blast_records:
                    query_name = record.query
                    for alignment in record.alignments:
                        ref_name = alignment.hit_def
                        for hsp in alignment.hsps:
                            identity = (hsp.identities / hsp.align_length) * 100
                            f_coords.write(
                                f"{hsp.sbjct_start}\t{hsp.sbjct_end}\t"
                                f"{hsp.query_start}\t{hsp.query_end}\t"
                                f"{abs(hsp.sbjct_end - hsp.sbjct_start) + 1}\t"
                                f"{abs(hsp.query_end - hsp.query_start) + 1}\t"
                                f"{identity:.2f}\t{ref_name}\t{query_name}\n"
                            )

                            variants = self._extract_variants(hsp, ref_name, query_name)
                            for var in variants:
                                f_snps.write(
                                    f"{var['ref_pos']}\t{var['ref']}\t{var['alt']}\t"
                                    f"{var['query_pos']}\t{var['type']}\t"
                                    f"{ref_name}\t{query_name}\n"
                                )

    def _extract_variants(
        self, hsp, ref_name: str, query_name: str
    ) -> List[Dict[str, Any]]:
        """从 HSP 中提取变异"""
        variants = []
        ref_pos = hsp.sbjct_start
        query_pos = hsp.query_start
        ref_step = 1 if hsp.sbjct_end >= hsp.sbjct_start else -1
        query_step = 1 if hsp.query_end >= hsp.query_start else -1

        query_seq = hsp.query
        ref_seq = hsp.sbjct

        i = 0
        while i < len(query_seq):
            q_base = query_seq[i]
            r_base = ref_seq[i]

            if q_base == "-":
                # BLAST query is the ref-derived sequence; a gap there means
                # the query genome has an insertion relative to that ref.
                ins_seq = ""
                while i < len(query_seq) and query_seq[i] == "-":
                    ins_seq += ref_seq[i]
                    i += 1
                variants.append({
                    "ref_pos": ref_pos,
                    "query_pos": query_pos,
                    "ref": "-",
                    "alt": ins_seq,
                    "type": "INS",
                })
                ref_pos += ref_step * len(ins_seq)
                continue

            elif r_base == "-":
                # BLAST subject is the query-genome hit; a gap there means
                # the query genome has a deletion relative to the ref.
                del_seq = ""
                while i < len(ref_seq) and ref_seq[i] == "-":
                    del_seq += query_seq[i]
                    i += 1
                variants.append({
                    "ref_pos": ref_pos,
                    "query_pos": query_pos,
                    "ref": del_seq,
                    "alt": "-",
                    "type": "DEL",
                })
                query_pos += query_step * len(del_seq)
                continue

            elif q_base != r_base:
                variants.append({
                    "ref_pos": ref_pos,
                    "query_pos": query_pos,
                    "ref": q_base,
                    "alt": r_base,
                    "type": "SNP",
                })

            ref_pos += ref_step
            query_pos += query_step
            i += 1

        return variants



# ======================= 分析处理器 =======================

class GeneIDProcessor:
    """Gene ID 模式处理器"""

    def __init__(
        self,
        ref_genome: str,
        ref_annotation: str,
        query_genome: str,
        output_dir: str,
        ref_name: str = "",
        qry_name: str = "",
        identity: float = 90,
        qry_gff: Optional[str] = None,
        upstream: int = 0,
        downstream: int = 0,
        min_aln_len: int = 100,
        merge_gap: int = 1000,
        query_genomes: Optional[List[Dict[str, Any]]] = None,
        candidate_limit: Optional[int] = 3,
        pairwise_all: bool = False,
    ):
        self.extractor = SequenceExtractor(ref_genome, ref_annotation, output_dir)
        self.aligner = BlastAligner(output_dir)
        self.query_genome = query_genome
        self.query_entries = _normalize_query_entries(query_genome, qry_name, qry_gff, query_genomes)
        self.output_dir = output_dir
        self.ref_name = ref_name
        self.qry_name = self.query_entries[0]["name"] if self.query_entries else qry_name
        self.identity = identity
        self.upstream = upstream
        self.downstream = downstream
        self.min_aln_len = min_aln_len
        self.merge_gap = merge_gap
        self.candidate_limit = candidate_limit
        self.pairwise_all = pairwise_all
        self.genome_files = {
            'ref_fasta': ref_genome,
            'ref_gff': ref_annotation,
            'qry_fasta': self.query_entries[0]["fasta"] if self.query_entries else query_genome,
            'qry_gff': self.query_entries[0]["gff"] if self.query_entries else qry_gff,
            'queries': self.query_entries,
        }

    def _align_queries(self, fasta_file: str, item_id: str) -> Optional[Dict[str, Any]]:
        query_results = []
        original_output_dir = self.aligner.output_dir
        for entry in self.query_entries:
            safe_query = sanitize_path_segment(entry["name"]) or "query"
            query_dir = os.path.join(self.output_dir, "queries", safe_query)
            os.makedirs(query_dir, exist_ok=True)
            self.aligner.output_dir = query_dir
            pair_prefix = f"ref__{safe_query}"
            aligned = self.aligner.align(entry["fasta"], fasta_file, pair_prefix, self.identity)
            if not aligned:
                continue
            aligned["query_name"] = entry["name"]
            aligned["query_fasta"] = entry["fasta"]
            aligned["query_gff"] = entry.get("gff")
            aligned["query_dir"] = query_dir
            query_results.append(aligned)
        self.aligner.output_dir = original_output_dir
        if not query_results:
            return None
        result = dict(query_results[0])
        result["query_results"] = query_results
        result["queries"] = self.query_entries
        result["pairwise_results"] = _precompute_pairwise(
            query_results, self.query_entries, self.aligner, self.output_dir,
            self.identity, self.candidate_limit, self.pairwise_all, self.merge_gap, self.min_aln_len
        )
        result["candidate_limit"] = self.candidate_limit
        result["pairwise_all"] = self.pairwise_all
        result["merge_gap"] = self.merge_gap
        return result

    def process(self, gene_id: str) -> Optional[Dict[str, Any]]:
        """处理单个基因"""
        print(f"\n{'='*50}")
        print(f"[INFO] 处理基因: {gene_id}")
        print(f"{'='*50}")

        fasta_file, gff_file, extraction_info = self.extractor.extract_by_gene_id(
            gene_id, self.upstream, self.downstream
        )
        if not fasta_file:
            return None

        result = self._align_queries(fasta_file, gene_id)
        if not result:
            return None

        result["fasta"] = fasta_file
        result["gff"] = gff_file
        result["id"] = gene_id
        result["output_dir"] = self.output_dir
        result["mode"] = "gene_id"
        result["ref_name"] = self.ref_name
        result["qry_name"] = self.qry_name
        result["identity"] = self.identity
        result["min_aln_len"] = self.min_aln_len
        result["genome_files"] = self.genome_files
        if extraction_info:
            result["extraction_info"] = extraction_info

        return result

    def set_output_dir(self, output_dir: str) -> None:
        self.output_dir = output_dir
        self.extractor.output_dir = output_dir
        self.aligner.output_dir = output_dir


class LocationProcessor:
    """Location 模式处理器"""

    def __init__(
        self,
        ref_genome: str,
        query_genome: str,
        output_dir: str,
        ref_name: str = "",
        qry_name: str = "",
        identity: float = 90,
        ref_gff: Optional[str] = None,
        qry_gff: Optional[str] = None,
        min_aln_len: int = 100,
        merge_gap: int = 1000,
        query_genomes: Optional[List[Dict[str, Any]]] = None,
        candidate_limit: Optional[int] = 3,
        pairwise_all: bool = False,
    ):
        self.extractor = SequenceExtractor(ref_genome, None, output_dir)
        self.aligner = BlastAligner(output_dir)
        self.ref_genome = ref_genome
        self.query_genome = query_genome
        self.query_entries = _normalize_query_entries(query_genome, qry_name, qry_gff, query_genomes)
        self.output_dir = output_dir
        self.ref_name = ref_name
        self.qry_name = self.query_entries[0]["name"] if self.query_entries else qry_name
        self.identity = identity
        self.min_aln_len = min_aln_len
        self.merge_gap = merge_gap
        self.candidate_limit = candidate_limit
        self.pairwise_all = pairwise_all
        self.genome_files = {
            'ref_fasta': ref_genome,
            'qry_fasta': self.query_entries[0]["fasta"] if self.query_entries else query_genome,
            'ref_gff': ref_gff,
            'qry_gff': self.query_entries[0]["gff"] if self.query_entries else qry_gff,
            'queries': self.query_entries,
        }

    def _align_queries(self, fasta_file: str, item_id: str) -> Optional[Dict[str, Any]]:
        query_results = []
        original_output_dir = self.aligner.output_dir
        for entry in self.query_entries:
            safe_query = sanitize_path_segment(entry["name"]) or "query"
            query_dir = os.path.join(self.output_dir, "queries", safe_query)
            os.makedirs(query_dir, exist_ok=True)
            self.aligner.output_dir = query_dir
            pair_prefix = f"ref__{safe_query}"
            aligned = self.aligner.align(entry["fasta"], fasta_file, pair_prefix, self.identity)
            if not aligned:
                continue
            aligned["query_name"] = entry["name"]
            aligned["query_fasta"] = entry["fasta"]
            aligned["query_gff"] = entry.get("gff")
            aligned["query_dir"] = query_dir
            query_results.append(aligned)
        self.aligner.output_dir = original_output_dir
        if not query_results:
            return None
        result = dict(query_results[0])
        result["query_results"] = query_results
        result["queries"] = self.query_entries
        result["pairwise_results"] = _precompute_pairwise(
            query_results, self.query_entries, self.aligner, self.output_dir,
            self.identity, self.candidate_limit, self.pairwise_all, self.merge_gap, self.min_aln_len
        )
        result["candidate_limit"] = self.candidate_limit
        result["pairwise_all"] = self.pairwise_all
        result["merge_gap"] = self.merge_gap
        return result

    def process(
        self,
        chrom: str,
        start: int,
        end: int,
        name: Optional[str] = None
    ) -> Optional[Dict[str, Any]]:
        """处理单个区域"""
        loc_name = name or f"{chrom}_{start}_{end}"
        print(f"\n{'='*50}")
        print(f"[INFO] 处理区域: {chrom}:{start}-{end}")
        print(f"{'='*50}")

        fasta_file = self.extractor.extract_by_location(chrom, start, end, loc_name)
        if not fasta_file:
            return None

        result = self._align_queries(fasta_file, loc_name)
        if not result:
            return None

        result["fasta"] = fasta_file
        result["id"] = loc_name
        result["location"] = f"{chrom}:{start}-{end}"
        result["output_dir"] = self.output_dir
        result["mode"] = "location"
        result["ref_name"] = self.ref_name
        result["qry_name"] = self.qry_name
        result["identity"] = self.identity
        result["min_aln_len"] = self.min_aln_len
        result["genome_files"] = self.genome_files

        return result

    def set_output_dir(self, output_dir: str) -> None:
        self.output_dir = output_dir
        self.extractor.output_dir = output_dir
        self.aligner.output_dir = output_dir


class SequenceProcessor:
    """Sequence 模式处理器"""

    def __init__(
        self,
        ref_genome: str,
        output_dir: str,
        ref_name: str = "",
        identity: float = 90,
        ref_gff: Optional[str] = None,
        min_aln_len: int = 100,
        merge_gap: int = 1000,
        query_genomes: Optional[List[Dict[str, Any]]] = None,
        upstream: int = 0,
        downstream: int = 0,
        candidate_limit: Optional[int] = 3,
        pairwise_all: bool = False,
    ):
        self.aligner = BlastAligner(output_dir)
        self.ref_genome = ref_genome
        self.query_entries = _normalize_query_entries(
            ref_genome, ref_name, ref_gff, query_genomes
        )
        self.output_dir = output_dir
        self.ref_name = self.query_entries[0]["name"] if self.query_entries else ref_name
        self.identity = identity
        self.min_aln_len = min_aln_len
        self.merge_gap = merge_gap
        self.upstream = upstream
        self.downstream = downstream
        self.candidate_limit = candidate_limit
        self.pairwise_all = pairwise_all
        self.genome_files = {
            'ref_fasta': ref_genome,
            'ref_gff': ref_gff,
            'queries': self.query_entries,
        }

    def _align_queries(self, fasta_file: str, item_id: str) -> Optional[Dict[str, Any]]:
        query_results = []
        original_output_dir = self.aligner.output_dir
        for entry in self.query_entries:
            safe_query = sanitize_path_segment(entry["name"]) or "query"
            query_dir = os.path.join(self.output_dir, "queries", safe_query)
            os.makedirs(query_dir, exist_ok=True)
            self.aligner.output_dir = query_dir
            pair_prefix = f"ref__{safe_query}"
            aligned = self.aligner.align(entry["fasta"], fasta_file, pair_prefix, self.identity)
            if not aligned:
                continue
            aligned["query_name"] = entry["name"]
            aligned["query_fasta"] = entry["fasta"]
            aligned["query_gff"] = entry.get("gff")
            aligned["query_dir"] = query_dir
            query_results.append(aligned)
        self.aligner.output_dir = original_output_dir
        if not query_results:
            return None
        result = dict(query_results[0])
        result["query_results"] = query_results
        result["queries"] = self.query_entries
        result["pairwise_results"] = _precompute_pairwise(
            query_results, self.query_entries, self.aligner, self.output_dir,
            self.identity, self.candidate_limit, self.pairwise_all,
            self.merge_gap, self.min_aln_len, self.upstream, self.downstream
        )
        result["candidate_limit"] = self.candidate_limit
        result["pairwise_all"] = self.pairwise_all
        result["merge_gap"] = self.merge_gap
        result["query_upstream"] = self.upstream
        result["query_downstream"] = self.downstream
        return result

    def set_output_dir(self, output_dir: str):
        """设置输出目录（同时更新 aligner）"""
        self.output_dir = output_dir
        self.aligner = BlastAligner(output_dir)

    def process(self, seq_file: str) -> Optional[Dict[str, Any]]:
        """处理序列文件"""
        print(f"\n{'='*50}")
        print(f"[INFO] 处理序列文件: {seq_file}")
        print(f"{'='*50}")

        seq_id, processed_file = self._prepare_fasta(seq_file)

        sequence = ""
        with open(processed_file, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    sequence += line.strip()

        result = self._align_queries(processed_file, seq_id)
        if not result:
            return None

        result["fasta"] = processed_file
        result["id"] = seq_id
        result["sequence"] = sequence
        result["mode"] = "sequence"
        result["ref_name"] = self.ref_name
        result["identity"] = self.identity
        result["min_aln_len"] = self.min_aln_len
        result["query_upstream"] = self.upstream
        result["query_downstream"] = self.downstream
        result["genome_files"] = self.genome_files

        return result

    def process_sequence_text(self, sequence: str, seq_id: str = "query_seq") -> Optional[Dict[str, Any]]:
        """处理序列文本（GUI 用）"""
        # 清理 seq_id 中的非法路径字符用于文件名和 FASTA header
        import re
        safe_seq_id = re.sub(r'[<>:"/\\|?*]', "_", seq_id)
        fasta_file = os.path.join(self.output_dir, f"{safe_seq_id}.fasta")
        with open(fasta_file, "w") as f:
            f.write(f">{safe_seq_id}\n{sequence}\n")
        
        result = self._align_queries(fasta_file, safe_seq_id)
        if not result:
            return None

        result["fasta"] = fasta_file
        result["id"] = safe_seq_id  # 使用清理后的 ID
        result["sequence"] = sequence
        result["mode"] = "sequence"
        result["ref_name"] = self.ref_name
        result["identity"] = self.identity
        result["min_aln_len"] = self.min_aln_len
        result["query_upstream"] = self.upstream
        result["query_downstream"] = self.downstream
        result["genome_files"] = self.genome_files

        return result

    def _prepare_fasta(self, seq_file: str) -> Tuple[str, str]:
        """确保 FASTA 文件有 ID 行"""
        with open(seq_file, "r") as f:
            lines = f.readlines()

        has_id = any(line.startswith(">") for line in lines)

        if has_id:
            seq_id = None
            for line in lines:
                if line.startswith(">"):
                    seq_id = line.strip().lstrip(">").split()[0]
                    break
            out_path = os.path.join(self.output_dir, f"{seq_id}.fasta")
            with open(out_path, "w") as f:
                f.writelines(lines)
            return seq_id, out_path
        else:
            seq_id = "query_seq"
            out_path = os.path.join(self.output_dir, f"{seq_id}.fasta")
            with open(out_path, "w") as f:
                f.write(f">{seq_id}\n")
                for line in lines:
                    if line.strip():
                        f.write(line)
            return seq_id, out_path
