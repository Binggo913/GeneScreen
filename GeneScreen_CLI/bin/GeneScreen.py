#!/usr/bin/env python3
"""
GeneScreen 1.0 - 基因组比对与变异分析工具

使用 BLAST+ + pyfaidx + ripgrep，实现跨平台兼容 (Windows/macOS/Linux)

支持三种输入模式：
1. Gene ID 模式: 从参考基因组提取指定基因，与目标基因组比对
2. Location 模式: 从参考基因组提取指定区域，与目标基因组比对
3. Sequence 模式: 用户提供序列，与参考基因组比对
"""

import os
import subprocess
import argparse
import re
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path

from pyfaidx import Fasta
from Bio.Blast import NCBIXML

from GeneScreenVisualizer import GeneIDVisualizer, LocationVisualizer, SequenceVisualizer
from GenomeManager import GenomeManager
from multi_query_result import parse_coords_candidates


# 全局基因组管理器
genome_manager = GenomeManager()


@dataclass
class GenomeEntry:
    name: str
    fasta: str
    annotation: str = None
    source: str = "preloaded"


def _parse_ref_entry(values):
    if not values:
        raise ValueError("-ref 需要指定参考基因组")
    if len(values) == 1:
        return {"name": values[0], "fasta": None, "annotation": None, "source": "preloaded"}
    if len(values) == 3:
        return {"name": values[0], "fasta": values[1], "annotation": values[2], "source": "explicit"}
    raise ValueError("-ref 支持 1 个值（入库名称）或 3 个值（名称 FASTA GFF）")


def _parse_query_entries(entries):
    parsed = []
    for values in entries or []:
        if len(values) == 1:
            parsed.append({"name": values[0], "fasta": None, "annotation": None, "source": "preloaded"})
        elif len(values) in (2, 3):
            parsed.append({
                "name": values[0],
                "fasta": values[1],
                "annotation": values[2] if len(values) == 3 else None,
                "source": "explicit",
            })
        else:
            raise ValueError("-qry 支持 1 个值（入库名称）、2 个值（名称 FASTA）或 3 个值（名称 FASTA GFF）")
    return parsed


def _resolve_genome_entry(entry, annotation_source=None, require_annotation=False, role="基因组"):
    if entry["source"] == "explicit":
        fasta = entry["fasta"]
        annotation = entry.get("annotation")
    else:
        fasta, annotation = genome_manager.get(entry["name"], annotation_source=annotation_source)
    if not fasta:
        raise ValueError(f"{role} '{entry['name']}' 不可用")
    if require_annotation and not annotation:
        raise ValueError(f"{role} '{entry['name']}' 需要注释文件")
    return GenomeEntry(entry["name"], fasta, annotation, entry["source"])


# ======================= 工具函数 =======================
def run_cmd(cmd, error_msg="命令执行失败"):
    """执行 shell 命令"""
    print(f"[CMD] {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"[ERROR] {error_msg}")
        return False
    return True


def ensure_blast_db(fasta_path):
    """确保 BLAST 数据库存在"""
    db_file = f"{fasta_path}.nin"
    if not os.path.exists(db_file):
        print(f"[INFO] 创建 BLAST 数据库: {fasta_path}")
        cmd = f'makeblastdb -in "{fasta_path}" -dbtype nucl'
        if not run_cmd(cmd, "创建 BLAST 数据库失败"):
            return False
    return True


def sanitize_path_segment(text):
    """清理用于输出子目录名的输入 ID。"""
    return re.sub(r'[<>:"/\\|?*]', "_", str(text)).strip()


def normalize_query_entries(query_genome=None, qry_name="", qry_gff=None, query_entries=None):
    normalized = []
    if query_entries:
        for index, entry in enumerate(query_entries, start=1):
            if isinstance(entry, GenomeEntry):
                normalized.append({
                    "name": entry.name,
                    "fasta": entry.fasta,
                    "gff": entry.annotation,
                })
            else:
                normalized.append({
                    "name": entry.get("name") or f"query_{index}",
                    "fasta": entry.get("fasta") or entry.get("fasta_path"),
                    "gff": entry.get("annotation") or entry.get("annotation_path") or entry.get("gff"),
                })
    elif query_genome:
        normalized.append({"name": qry_name or "query", "fasta": query_genome, "gff": qry_gff})
    return [entry for entry in normalized if entry.get("fasta")]


def extract_candidate_fasta(query_entry, candidate, output_dir, upstream=0, downstream=0):
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
        handle.write(
            f">{safe_query}_{candidate_id}|{chrom}:{extract_left}-{extract_right}|"
            f"match={left}-{right}|strand={strand}|upstream={upstream}|downstream={downstream}\n"
        )
        handle.write(f"{seq}\n")
    return out_path


def precompute_pairwise(
    query_results,
    query_entries,
    aligner,
    output_dir,
    identity,
    candidate_limit=3,
    pairwise_all=False,
    upstream=0,
    downstream=0,
):
    bundles = []
    for index, query_result in enumerate(query_results):
        if index >= len(query_entries):
            continue
        entry = query_entries[index]
        candidates = parse_coords_candidates(query_result.get("coords"))
        if not pairwise_all:
            candidates = candidates[:candidate_limit]
        prepared = []
        for candidate in candidates:
            fasta = extract_candidate_fasta(entry, candidate, output_dir, upstream, downstream)
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


# ======================= 序列提取模块 =======================
class SequenceExtractor:
    """使用 pyfaidx 从参考基因组提取序列"""

    def __init__(self, genome, annotation, output_dir):
        self.genome = genome
        self.annotation = annotation
        self.output_dir = output_dir
        self._fasta = None

    @property
    def fasta(self):
        """懒加载 FASTA 文件"""
        if self._fasta is None:
            self._fasta = Fasta(self.genome)
        return self._fasta

    def extract_by_gene_id(self, gene_id, upstream=0, downstream=0):
        """通过 Gene ID 提取序列，支持模糊匹配和上下游延伸
        
        参数:
            gene_id: 基因 ID
            upstream: 上游延伸长度 (bp)
            downstream: 下游延伸长度 (bp)
        
        返回:
            (fasta_file, linkview_gff, extraction_info) 或 (None, None, None)
            extraction_info 包含: gene_rel_start, gene_rel_end, strand, actual_start, actual_end
        """
        print(f"[INFO] 提取基因 {gene_id} 的序列...")
        if upstream > 0 or downstream > 0:
            print(f"[INFO] 上游延伸: {upstream} bp, 下游延伸: {downstream} bp")

        original_gff = os.path.join(self.output_dir, f"{gene_id}.original.gff3")
        linkview_gff = os.path.join(self.output_dir, f"{gene_id}.linkview.gff3")

        # 检查注释文件是否存在
        if not self.annotation or not os.path.exists(self.annotation):
            print(f"[WARNING] 参考基因组无注释文件，Gene ID 模式需要 GFF 注释来定位基因")
            print(f"[WARNING] 请使用 Location 模式 (-loc) 或为参考基因组添加注释文件")
            return None, None, None

        # 使用 ripgrep 过滤 GFF（模糊匹配）
        cmd = f'rg "{gene_id}" "{self.annotation}" > "{original_gff}"'
        if not run_cmd(cmd, f"提取 {gene_id} 注释失败"):
            # 回退到纯 Python
            self._grep_fallback(gene_id, self.annotation, original_gff)

        # 检查是否找到匹配
        if not os.path.exists(original_gff) or os.path.getsize(original_gff) == 0:
            # 尝试更宽松的匹配（去掉前缀/后缀）
            print(f"[INFO] 尝试模糊匹配 {gene_id}...")
            self._fuzzy_grep(gene_id, self.annotation, original_gff)

        # 阶段1: 解析 GFF 获取基因信息（带 ID 匹配校验）
        gene_info = self._parse_gff_info(original_gff, gene_id)
        if not gene_info:
            print(f"[ERROR] 未找到 {gene_id} 的基因记录")
            return None, None, None

        # 计算理论提取区域（考虑链方向）
        gene_start = gene_info["start"]
        gene_end = gene_info["end"]
        strand = gene_info["strand"]
        chrom = gene_info["chrom"]

        if strand == "+":
            # 正链：上游在左侧，下游在右侧
            theory_start = gene_start - upstream
            theory_end = gene_end + downstream
        else:
            # 负链：上游在右侧，下游在左侧
            theory_start = gene_start - downstream
            theory_end = gene_end + upstream

        # 边界防御：截断到染色体边界
        chrom_normalized = self._normalize_chrom(chrom)
        chrom_length = len(self.fasta[chrom_normalized])
        
        actual_start = max(1, theory_start)
        actual_end = min(chrom_length, theory_end)

        if actual_start != theory_start or actual_end != theory_end:
            print(f"[INFO] 边界裁剪: [{theory_start}, {theory_end}] -> [{actual_start}, {actual_end}]")

        # 计算基因相对位置（基于 actual_start）
        gene_rel_start = gene_start - actual_start + 1
        gene_rel_end = gene_end - actual_start + 1

        # 阶段2: 从全量 GFF 按坐标区间提取 feature，生成 linkview.gff3
        self._generate_linkview_gff(gene_id, chrom_normalized, actual_start, actual_end, linkview_gff)

        # 使用 pyfaidx 提取序列
        fasta_file = self._extract_fasta_region(gene_id, chrom_normalized, actual_start, actual_end)
        
        # 构建提取信息
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

    def _normalize_chrom(self, chrom):
        """标准化染色体名称，自动匹配基因组中的实际名称
        
        支持的输入格式：Chr1, chr1, Chr01, chr01, 1
        """
        # 获取基因组中的染色体列表
        available_chroms = list(self.fasta.keys())
        
        # 如果完全匹配，直接返回
        if chrom in available_chroms:
            return chrom
        
        # 提取数字部分
        import re
        num_match = re.search(r'(\d+)', chrom)
        if not num_match:
            return chrom  # 无法解析，返回原值
        
        num = num_match.group(1)
        num_int = int(num)
        
        # 尝试各种格式
        candidates = [
            str(num_int),           # 1
            num,                    # 01
            f"Chr{num_int}",        # Chr1
            f"chr{num_int}",        # chr1
            f"Chr{num.zfill(2)}",   # Chr01
            f"chr{num.zfill(2)}",   # chr01
            f"Chr{num}",            # Chr1 或 Chr01
            f"chr{num}",            # chr1 或 chr01
        ]
        
        for candidate in candidates:
            if candidate in available_chroms:
                if candidate != chrom:
                    print(f"[INFO] 染色体名称转换: {chrom} -> {candidate}")
                return candidate
        
        # 都没匹配到，返回原值
        return chrom

    def extract_by_location(self, chrom, start, end, name=None):
        """通过染色体坐标提取序列"""
        loc_name = name or f"{chrom}_{start}_{end}"
        print(f"[INFO] 提取区域 {chrom}:{start}-{end} 的序列...")

        # 标准化染色体名称
        chrom = self._normalize_chrom(chrom)

        fasta_file = os.path.join(self.output_dir, f"{loc_name}.fasta")

        # 使用 pyfaidx 提取序列
        try:
            seq = self.fasta[chrom][start - 1 : end]  # pyfaidx 是 0-based
            with open(fasta_file, "w") as f:
                f.write(f">{loc_name}\n{str(seq)}\n")
            print(f"[INFO] 已生成序列文件: {fasta_file}")
            return fasta_file
        except Exception as e:
            print(f"[ERROR] 提取序列失败: {e}")
            return None

    def _grep_fallback(self, pattern, input_file, output_file):
        """纯 Python grep 回退"""
        with open(input_file, "r", encoding="utf-8", errors="ignore") as f_in, open(output_file, "w") as f_out:
            for line in f_in:
                if pattern in line:
                    f_out.write(line)

    def _fuzzy_grep(self, gene_id, input_file, output_file):
        """模糊匹配 Gene ID
        
        支持的匹配模式：
        - 完整 ID: AGIS_Os01g000010
        - 简化 ID: Os01g000010
        - RAP-DB 格式: Os01g0100100
        """
        import re
        
        # 提取核心 ID 部分（如 Os01g000010 或 Os01g0100100）
        # 匹配模式：Os + 染色体号 + g + 数字
        core_match = re.search(r'(Os\d+g\d+)', gene_id, re.IGNORECASE)
        if not core_match:
            return
        
        core_id = core_match.group(1)
        print(f"[INFO] 使用核心 ID 搜索: {core_id}")
        
        with open(input_file, "r", encoding="utf-8", errors="ignore") as f_in, open(output_file, "w") as f_out:
            for line in f_in:
                # 不区分大小写匹配
                if core_id.lower() in line.lower():
                    f_out.write(line)

    def _parse_gff_info(self, original_gff, gene_id):
        """阶段1: 解析 GFF 文件，只获取基因信息 (chrom, start, end, strand)
        
        优先从 gene 记录获取，若无则从 mRNA 记录中选取坐标范围最大的
        只匹配 ID 或 Parent 精确包含 gene_id 的记录（按分隔符边界匹配）
        """
        import re
        gene_info = None
        mRNA_records = []

        def id_matches(attributes, target_id):
            """检查 GFF attributes 是否精确匹配目标 gene_id
            
            检查字段：ID, Parent, Name, gene_id, locus_tag
            使用分隔符边界匹配，避免 Os01g000010 误匹配 Os01g0000100
            分隔符包括：逗号、分号、点号、下划线、冒号、行首、行尾
            """
            match_keys = ("ID", "Parent", "Name", "gene_id", "locus_tag")
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    key = key.strip()
                    value = value.strip()
                    if key in match_keys:
                        # 使用正则边界匹配：target_id 前后必须是分隔符或字符串边界
                        # 支持 ID=gene:LOC_Os01g01010 格式
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
                
                # 只处理 ID 匹配的记录
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
                    break  # 找到匹配的 gene 记录，直接返回
                elif feature_type == "mRNA":
                    mRNA_records.append({
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "range": end - start,
                    })

        # 兜底：若无 gene 记录，从 mRNA 中选取坐标范围最大的
        if not gene_info and mRNA_records:
            best_mRNA = max(mRNA_records, key=lambda x: x["range"])
            gene_info = {
                "chrom": best_mRNA["chrom"],
                "start": best_mRNA["start"],
                "end": best_mRNA["end"],
                "strand": best_mRNA["strand"],
            }
            print(f"[INFO] 未找到 gene 记录，使用 mRNA 记录作为兜底")

        return gene_info

    def _generate_linkview_gff(self, gene_id, chrom, actual_start, actual_end, linkview_gff):
        """阶段2: 从全量 GFF 按坐标区间提取 feature，生成相对坐标的 linkview.gff3
        
        使用坐标范围过滤（从全量 annotation 提取）：
        - 只保留同一染色体上与 [actual_start, actual_end] 有重叠的 feature
        - 只保留与目标 gene_id 相关的 feature（通过 ID/Parent 属性匹配）
        - 跨界 feature 会被裁剪到区间范围内
        - 过滤掉 region/chromosome 等无 Parent 的顶层记录，避免 LINKVIEW 报错
        """
        import re
        region_len = actual_end - actual_start + 1
        
        # 允许的 feature 类型（LINKVIEW 支持的基因结构类型，小写比对）
        allowed_types = {"gene", "mrna", "transcript", "exon", "cds", 
                         "five_prime_utr", "three_prime_utr", "intron"}

        def chrom_matches(feat_chrom, target_chrom):
            """检查染色体名称是否匹配（支持不同命名格式）
            
            支持：Chr01 vs 1, chr1 vs Chr01, Chr1 vs chr01 等
            """
            if feat_chrom == target_chrom:
                return True
            
            # 提取数字部分比较
            feat_num = re.search(r'(\d+)', feat_chrom)
            target_num = re.search(r'(\d+)', target_chrom)
            if feat_num and target_num:
                return int(feat_num.group(1)) == int(target_num.group(1))
            return False

        def feature_belongs_to_gene(attributes, target_gene_id):
            """检查 feature 是否属于目标基因
            
            通过 ID 或 Parent 属性中是否包含 gene_id 来判断
            """
            # 检查 ID 属性（gene 记录）
            id_match = re.search(r'ID=([^;]+)', attributes)
            if id_match:
                id_value = id_match.group(1)
                if target_gene_id in id_value:
                    return True
            
            # 检查 Parent 属性（mRNA、exon、CDS、UTR 等）
            parent_match = re.search(r'Parent=([^;]+)', attributes)
            if parent_match:
                parent_value = parent_match.group(1)
                if target_gene_id in parent_value:
                    return True
            
            # 检查 gene_id 属性
            gene_id_match = re.search(r'gene_id=([^;]+)', attributes)
            if gene_id_match:
                gene_id_value = gene_id_match.group(1)
                if target_gene_id in gene_id_value:
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

                # 染色体过滤（支持不同命名格式）
                feat_chrom = fields[0]
                if not chrom_matches(feat_chrom, chrom):
                    continue

                # feature 类型过滤：只保留基因结构相关类型（小写比对）
                feat_type = fields[2]
                if feat_type.lower() not in allowed_types:
                    continue

                # Gene ID 过滤：只保留与目标基因相关的 feature
                attributes = fields[8]
                if not feature_belongs_to_gene(attributes, gene_id):
                    continue

                # 坐标范围过滤：只保留与 [actual_start, actual_end] 有重叠的记录
                feat_start = int(fields[3])
                feat_end = int(fields[4])
                if feat_end < actual_start or feat_start > actual_end:
                    continue  # 无重叠，跳过

                # 裁剪：将跨界 feature 裁剪到 [actual_start, actual_end] 范围内
                clipped_start = max(feat_start, actual_start)
                clipped_end = min(feat_end, actual_end)

                # 计算相对坐标（基于 actual_start，1-based）
                rel_start = clipped_start - actual_start + 1
                rel_end = clipped_end - actual_start + 1

                # 确保相对坐标在有效范围内 [1, region_len]
                rel_start = max(1, min(rel_start, region_len))
                rel_end = max(1, min(rel_end, region_len))

                f_out.write(
                    f"{gene_id}\t{fields[1]}\t{fields[2]}\t{rel_start}\t{rel_end}\t"
                    f"{fields[5]}\t{fields[6]}\t{fields[7]}\t{fields[8]}\n"
                )

    def _extract_fasta_region(self, gene_id, chrom, start, end):
        """使用 pyfaidx 提取指定区域的序列"""
        fasta_file = os.path.join(self.output_dir, f"{gene_id}.fasta")

        try:
            # pyfaidx 使用 0-based 半开区间，GFF 使用 1-based 闭区间
            seq = self.fasta[chrom][start - 1 : end]

            with open(fasta_file, "w") as f:
                f.write(f">{gene_id}\n{str(seq)}\n")

            print(f"[INFO] 已生成序列文件: {fasta_file} (长度: {len(str(seq))} bp)")
            return fasta_file
        except Exception as e:
            print(f"[ERROR] 提取序列失败: {e}")
            return None

    def _parse_gff(self, original_gff, linkview_gff, gene_id):
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

    def _extract_fasta(self, gene_id, gene_info):
        """使用 pyfaidx 提取序列"""
        fasta_file = os.path.join(self.output_dir, f"{gene_id}.fasta")

        try:
            chrom = gene_info["chrom"]
            start = gene_info["start"]
            end = gene_info["end"]
            seq = self.fasta[chrom][start - 1 : end]  # pyfaidx 是 0-based

            with open(fasta_file, "w") as f:
                f.write(f">{gene_id}\n{str(seq)}\n")

            print(f"[INFO] 已生成序列文件: {fasta_file}")
            return fasta_file
        except Exception as e:
            print(f"[ERROR] 提取序列失败: {e}")
            return None


# ======================= BLAST 比对模块 =======================
class BlastAligner:
    """使用 BLAST+ 进行序列比对"""

    def __init__(self, output_dir):
        self.output_dir = output_dir

    def align(self, reference, query, prefix_name, identity_threshold=90):
        """
        执行 BLAST 比对

        参数:
        - reference: 参考基因组（作为数据库）
        - query: 查询序列
        - prefix_name: 输出前缀
        - identity_threshold: 最小 identity 阈值
        """
        prefix = os.path.join(self.output_dir, prefix_name)

        # 确保 BLAST 数据库存在
        if not ensure_blast_db(reference):
            return None

        # 运行 blastn
        xml_file = f"{prefix}.blast.xml"
        print(f"[INFO] 运行 BLAST 比对...")
        cmd = (
            f'blastn -query "{query}" -db "{reference}" '
            f'-out "{xml_file}" -outfmt 5 '
            f"-perc_identity {identity_threshold} "
            f"-num_threads 2"
        )
        if not run_cmd(cmd, "BLAST 比对失败"):
            return None

        # 解析结果
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

    def _parse_blast_xml(self, xml_file, coords_file, snps_file):
        """解析 BLAST XML 输出"""
        with open(xml_file) as f:
            blast_records = NCBIXML.parse(f)

            with open(coords_file, "w") as f_coords, open(snps_file, "w") as f_snps:
                # 写入 header
                f_coords.write(
                    "[S1]\t[E1]\t[S2]\t[E2]\t[LEN1]\t[LEN2]\t[%IDY]\t[REF]\t[QUERY]\n"
                )
                f_snps.write("[QRY_POS]\t[REF]\t[ALT]\t[REF_POS]\t[TYPE]\t[QRY_NAME]\t[REF_NAME]\n")

                for record in blast_records:
                    query_name = record.query
                    for alignment in record.alignments:
                        ref_name = alignment.hit_def
                        for hsp in alignment.hsps:
                            # 写入 coords
                            identity = (hsp.identities / hsp.align_length) * 100
                            f_coords.write(
                                f"{hsp.sbjct_start}\t{hsp.sbjct_end}\t"
                                f"{hsp.query_start}\t{hsp.query_end}\t"
                                f"{abs(hsp.sbjct_end - hsp.sbjct_start) + 1}\t"
                                f"{abs(hsp.query_end - hsp.query_start) + 1}\t"
                                f"{identity:.2f}\t{ref_name}\t{query_name}\n"
                            )

                            # 提取变异
                            variants = self._extract_variants(hsp, ref_name, query_name)
                            for var in variants:
                                f_snps.write(
                                    f"{var['ref_pos']}\t{var['ref']}\t{var['alt']}\t"
                                    f"{var['query_pos']}\t{var['type']}\t"
                                    f"{ref_name}\t{query_name}\n"
                                )

    def _extract_variants(self, hsp, ref_name, query_name):
        """从 HSP 中提取变异"""
        variants = []
        ref_pos = hsp.sbjct_start
        query_pos = hsp.query_start
        ref_step = 1 if hsp.sbjct_end >= hsp.sbjct_start else -1
        query_step = 1 if hsp.query_end >= hsp.query_start else -1

        query_seq = hsp.query
        ref_seq = hsp.sbjct
        midline = hsp.match

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
                variants.append(
                    {
                        "ref_pos": ref_pos,
                        "query_pos": query_pos,
                        "ref": "-",
                        "alt": ins_seq,
                        "type": "INS",
                    }
                )
                ref_pos += ref_step * len(ins_seq)
                continue

            elif r_base == "-":
                # BLAST subject is the query-genome hit; a gap there means
                # the query genome has a deletion relative to the ref.
                del_seq = ""
                while i < len(ref_seq) and ref_seq[i] == "-":
                    del_seq += query_seq[i]
                    i += 1
                variants.append(
                    {
                        "ref_pos": ref_pos,
                        "query_pos": query_pos,
                        "ref": del_seq,
                        "alt": "-",
                        "type": "DEL",
                    }
                )
                query_pos += query_step * len(del_seq)
                continue

            elif q_base != r_base:
                # SNP
                variants.append(
                    {
                        "ref_pos": ref_pos,
                        "query_pos": query_pos,
                        "ref": q_base,
                        "alt": r_base,
                        "type": "SNP",
                    }
                )

            ref_pos += ref_step
            query_pos += query_step
            i += 1

        return variants


# ======================= 输入处理器 =======================
class GeneIDProcessor:
    """Gene ID 模式处理器"""

    def __init__(self, ref_genome, ref_annotation, query_genome, output_dir, ref_name="", qry_name="", identity=90, qry_gff=None, upstream=0, downstream=0, min_aln_len=100, merge_gap=1000, query_entries=None, candidate_limit=3, pairwise_all=False):
        self.extractor = SequenceExtractor(ref_genome, ref_annotation, output_dir)
        self.aligner = BlastAligner(output_dir)
        self.query_genome = query_genome
        self.query_entries = normalize_query_entries(query_genome, qry_name, qry_gff, query_entries)
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
        # 保存基因组文件路径
        self.genome_files = {
            'ref_fasta': ref_genome,
            'ref_gff': ref_annotation,
            'qry_fasta': self.query_entries[0]["fasta"] if self.query_entries else query_genome,
            'qry_gff': self.query_entries[0]["gff"] if self.query_entries else qry_gff,
            'queries': self.query_entries,
        }
        # 查询基因组索引文件（用于宏观概览图）
        self.qry_index_file = f"{query_genome}.fai" if query_genome else None

    def _align_queries(self, fasta_file, item_id):
        query_results = []
        original_output_dir = self.aligner.output_dir
        for entry in self.query_entries:
            safe_query = sanitize_path_segment(entry["name"]) or "query"
            query_dir = os.path.join(self.output_dir, "queries", safe_query)
            os.makedirs(query_dir, exist_ok=True)
            self.aligner.output_dir = query_dir
            aligned = self.aligner.align(entry["fasta"], fasta_file, f"ref__{safe_query}", self.identity)
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
        result["pairwise_results"] = precompute_pairwise(
            query_results, self.query_entries, self.aligner, self.output_dir,
            self.identity, self.candidate_limit, self.pairwise_all
        )
        result["candidate_limit"] = self.candidate_limit
        result["pairwise_all"] = self.pairwise_all
        return result

    def process(self, gene_id):
        print(f"\n{'='*50}")
        print(f"[INFO] 处理基因: {gene_id}")
        if self.upstream > 0 or self.downstream > 0:
            print(f"[INFO] 上游: {self.upstream} bp, 下游: {self.downstream} bp")
        print(f"{'='*50}")

        # 提取基因序列（带上下游延伸）
        fasta_file, gff_file, extraction_info = self.extractor.extract_by_gene_id(
            gene_id, self.upstream, self.downstream
        )
        if not fasta_file:
            return None

        # 比对
        result = self._align_queries(fasta_file, gene_id)
        if not result:
            return None

        result["fasta"] = fasta_file
        result["gff"] = gff_file
        result["id"] = gene_id

        # 添加提取信息到结果
        if extraction_info:
            result["extraction_info"] = extraction_info

        # 可视化（传递 min_aln_len, merge_gap, ref_index_file）
        visualizer = GeneIDVisualizer(
            self.output_dir,
            self.min_aln_len,
            self.merge_gap,
            self.qry_index_file,
            self.identity,
        )
        visualizer.visualize(result, self.ref_name, self.qry_name, self.identity, genome_files=self.genome_files)

        return result

    def set_output_dir(self, output_dir):
        self.output_dir = output_dir
        self.extractor.output_dir = output_dir
        self.aligner.output_dir = output_dir


class LocationProcessor:
    """Location 模式处理器"""

    def __init__(self, ref_genome, query_genome, output_dir, ref_name="", qry_name="", identity=90, ref_gff=None, qry_gff=None, min_aln_len=100, merge_gap=1000, query_entries=None, candidate_limit=3, pairwise_all=False):
        self.extractor = SequenceExtractor(ref_genome, None, output_dir)
        self.aligner = BlastAligner(output_dir)
        self.ref_genome = ref_genome
        self.query_genome = query_genome
        self.query_entries = normalize_query_entries(query_genome, qry_name, qry_gff, query_entries)
        self.output_dir = output_dir
        self.ref_name = ref_name
        self.qry_name = self.query_entries[0]["name"] if self.query_entries else qry_name
        self.identity = identity
        self.min_aln_len = min_aln_len
        self.merge_gap = merge_gap
        self.candidate_limit = candidate_limit
        self.pairwise_all = pairwise_all
        # 保存基因组文件路径
        self.genome_files = {
            'ref_fasta': ref_genome,
            'qry_fasta': self.query_entries[0]["fasta"] if self.query_entries else query_genome,
            'ref_gff': ref_gff,  # 参考基因组 GFF（用于提取区域内的基因注释）
            'qry_gff': self.query_entries[0]["gff"] if self.query_entries else qry_gff,
            'queries': self.query_entries,
        }
        # 查询基因组索引文件（用于宏观概览图）
        self.qry_index_file = f"{query_genome}.fai" if query_genome else None

    def _align_queries(self, fasta_file, item_id):
        query_results = []
        original_output_dir = self.aligner.output_dir
        for entry in self.query_entries:
            safe_query = sanitize_path_segment(entry["name"]) or "query"
            query_dir = os.path.join(self.output_dir, "queries", safe_query)
            os.makedirs(query_dir, exist_ok=True)
            self.aligner.output_dir = query_dir
            aligned = self.aligner.align(entry["fasta"], fasta_file, f"ref__{safe_query}", self.identity)
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
        result["pairwise_results"] = precompute_pairwise(
            query_results, self.query_entries, self.aligner, self.output_dir,
            self.identity, self.candidate_limit, self.pairwise_all
        )
        result["candidate_limit"] = self.candidate_limit
        result["pairwise_all"] = self.pairwise_all
        return result

    def process(self, chrom, start, end, name=None):
        loc_name = name or f"{chrom}_{start}_{end}"
        print(f"\n{'='*50}")
        print(f"[INFO] 处理区域: {chrom}:{start}-{end}")
        print(f"{'='*50}")

        # 提取序列
        fasta_file = self.extractor.extract_by_location(chrom, start, end, loc_name)
        if not fasta_file:
            return None

        # 比对
        result = self._align_queries(fasta_file, loc_name)
        if not result:
            return None

        result["fasta"] = fasta_file
        result["id"] = loc_name
        result["location"] = f"{chrom}:{start}-{end}"  # 保存原始位置信息

        # 可视化（传递 min_aln_len, merge_gap, ref_index_file）
        visualizer = LocationVisualizer(
            self.output_dir,
            self.min_aln_len,
            self.merge_gap,
            self.qry_index_file,
            self.identity,
        )
        visualizer.visualize(result, self.ref_name, self.qry_name, self.identity, genome_files=self.genome_files)

        return result

    def set_output_dir(self, output_dir):
        self.output_dir = output_dir
        self.extractor.output_dir = output_dir
        self.aligner.output_dir = output_dir


class SequenceProcessor:
    """Sequence 模式处理器"""

    def __init__(
        self,
        ref_genome,
        output_dir,
        ref_name="",
        identity=90,
        ref_gff=None,
        min_aln_len=100,
        merge_gap=1000,
        query_entries=None,
        upstream=0,
        downstream=0,
        candidate_limit=3,
        pairwise_all=False,
    ):
        self.aligner = BlastAligner(output_dir)
        self.ref_genome = ref_genome
        self.query_entries = normalize_query_entries(ref_genome, ref_name, ref_gff, query_entries)
        self.output_dir = output_dir
        self.ref_name = self.query_entries[0]["name"] if self.query_entries else ref_name
        self.identity = identity
        self.min_aln_len = min_aln_len
        self.merge_gap = merge_gap
        self.upstream = upstream
        self.downstream = downstream
        self.candidate_limit = candidate_limit
        self.pairwise_all = pairwise_all
        # 保存基因组文件路径
        self.genome_files = {
            'ref_fasta': ref_genome,
            'ref_gff': ref_gff,
            'queries': self.query_entries,
        }
        # 参考基因组索引文件（用于宏观概览图）
        self.ref_index_file = f"{ref_genome}.fai" if ref_genome else None

    def _align_queries(self, fasta_file, item_id):
        query_results = []
        original_output_dir = self.aligner.output_dir
        for entry in self.query_entries:
            safe_query = sanitize_path_segment(entry["name"]) or "query"
            query_dir = os.path.join(self.output_dir, "queries", safe_query)
            os.makedirs(query_dir, exist_ok=True)
            self.aligner.output_dir = query_dir
            aligned = self.aligner.align(entry["fasta"], fasta_file, f"ref__{safe_query}", self.identity)
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
        result["pairwise_results"] = precompute_pairwise(
            query_results, self.query_entries, self.aligner, self.output_dir,
            self.identity, self.candidate_limit, self.pairwise_all,
            self.upstream, self.downstream
        )
        result["candidate_limit"] = self.candidate_limit
        result["pairwise_all"] = self.pairwise_all
        result["query_upstream"] = self.upstream
        result["query_downstream"] = self.downstream
        return result

    def process(self, seq_file):
        """处理序列文件（支持多序列 FASTA）"""
        print(f"\n{'='*50}")
        print(f"[INFO] 处理序列文件: {seq_file}")
        print(f"{'='*50}")

        # 解析多序列 FASTA
        sequences = self._parse_multi_fasta(seq_file)
        if not sequences:
            print(f"[ERROR] 未找到有效序列")
            return None
        
        print(f"[INFO] 共解析到 {len(sequences)} 个序列")
        
        results = []
        base_output_dir = self.output_dir
        multi_mode = len(sequences) > 1
        for i, seq in enumerate(sequences):
            seq_id = seq["seq_id"]
            sequence = seq["sequence"]
            print(f"\n[INFO] 处理序列 {seq_id} ({i+1}/{len(sequences)})")

            if multi_mode:
                item_output_dir = os.path.join(base_output_dir, sanitize_path_segment(seq_id) or "sequence")
                os.makedirs(item_output_dir, exist_ok=True)
                self.set_output_dir(item_output_dir)

            result = self._process_single(seq_id, sequence)
            if result:
                results.append(result)

        if multi_mode:
            self.set_output_dir(base_output_dir)
        
        print(f"\n[INFO] 完成 {len(results)}/{len(sequences)} 个序列")
        return results[-1] if results else None

    def _process_single(self, seq_id, sequence):
        """处理单个序列"""
        # 清理 seq_id 中的非法路径字符
        import re
        safe_seq_id = re.sub(r'[<>:"/\\|?*]', "_", seq_id)
        
        # 写入临时 FASTA 文件
        fasta_file = os.path.join(self.output_dir, f"{safe_seq_id}.fasta")
        with open(fasta_file, "w") as f:
            f.write(f">{safe_seq_id}\n{sequence}\n")

        # 比对
        result = self._align_queries(fasta_file, safe_seq_id)
        if not result:
            return None

        result["fasta"] = fasta_file
        result["id"] = safe_seq_id
        result["sequence"] = sequence
        result["query_upstream"] = self.upstream
        result["query_downstream"] = self.downstream

        # 可视化
        visualizer = SequenceVisualizer(
            self.output_dir,
            self.min_aln_len,
            self.merge_gap,
            self.ref_index_file,
            self.identity,
        )
        visualizer.visualize(result, self.ref_name, None, self.identity, genome_files=self.genome_files)

        return result

    def _parse_multi_fasta(self, seq_file):
        """
        解析多序列 FASTA 文件
        返回: [{"seq_id": "xxx", "sequence": "ATGC..."}, ...]
        """
        sequences = []
        current_id = None
        current_seq = []
        auto_idx = 0
        
        with open(seq_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    # 保存上一个序列
                    if current_seq:
                        seq_str = ''.join(current_seq)
                        if seq_str:
                            if current_id is None:
                                auto_idx += 1
                                current_id = f"query_seq_{auto_idx}"
                            sequences.append({"seq_id": current_id, "sequence": seq_str})
                    # 开始新序列
                    header = line[1:].strip()
                    current_id = header.split()[0] if header else None
                    current_seq = []
                else:
                    current_seq.append(line)
        
        # 保存最后一个序列
        if current_seq:
            seq_str = ''.join(current_seq)
            if seq_str:
                if current_id is None:
                    auto_idx += 1
                    current_id = f"query_seq_{auto_idx}"
                sequences.append({"seq_id": current_id, "sequence": seq_str})
        
        return sequences

    def set_output_dir(self, output_dir):
        self.output_dir = output_dir
        self.aligner.output_dir = output_dir


# ======================= 主程序 =======================
def main():
    parser = argparse.ArgumentParser(
        description="GeneScreen 1.0 - 基因组比对与变异分析工具 (跨平台版)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # Gene ID 模式（单个）
  python GeneScreen.py -ref Nippon -qry ZS97 -gid LOC_Os06g10990 -o output/

  # Gene ID 模式（多个）
  python GeneScreen.py -ref Nippon -qry ZS97 -gid LOC_Os06g10990 LOC_Os06g10991 -o output/

  # Gene ID 模式（从文件读取）
  python GeneScreen.py -ref Nippon -qry ZS97 -gid gene_list.txt -o output/

  # Gene ID 模式 + 上下游延伸（上游 2kb，下游 1kb）
  python GeneScreen.py -ref Nippon -qry ZS97 -gid LOC_Os06g10990 -u 2000 -d 1000 -o output/

  # Gene ID 模式 + 多查询基因组（已入库 / 混合显式路径）
  python GeneScreen.py -ref Nippon -qry ZS97 -qry MH63 -gid LOC_Os06g10990 -o output/
  python GeneScreen.py -ref Nippon -qry q1 q1.fa q1.gff3 -qry q2 q2.fa -gid LOC_Os06g10990 -o output/

  # Location 模式（单个）
  python GeneScreen.py -ref Nippon -qry ZS97 -loc Chr1:1000-2000 -o output/

  # Location 模式（多个）
  python GeneScreen.py -ref Nippon -qry ZS97 -loc Chr1:1000-2000 Chr2:3000-4000 -o output/

  # Location 模式（从文件读取）
  python GeneScreen.py -ref Nippon -qry ZS97 -loc positions.txt -o output/

  # Sequence 模式（支持多序列 FASTA，查询基因组可多选）
  python GeneScreen.py -qry ZS97 -seq query.fasta -o output/
  python GeneScreen.py -qry ZS97 -qry MH63 -seq query.fasta -u 2000 -d 1000 -o output/
        """,
    )

    # 基因组参数
    parser.add_argument("-ref", nargs='+', help="参考基因组条目：名称，或 名称 FASTA GFF；Gene ID/Location 必填")
    parser.add_argument("-qry", nargs='+', action='append', help="目标基因组条目，可重复：名称，或 名称 FASTA [GFF]")
    parser.add_argument("-ra", help="参考基因组注释文件（覆盖默认）")
    parser.add_argument("-ra-source", help="参考基因组注释版本（如 igv, ensembl_plants）")

    # 输入模式（互斥）
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-gid", nargs='+', help="Gene ID（多个值或文件路径）")
    input_group.add_argument("-loc", nargs='+', help="位置（多个区域字符串或文件路径，格式: Chr1:1000-2000）")
    input_group.add_argument("-seq", help="序列文件（FASTA 格式，支持多序列）")

    # 比对参数
    parser.add_argument(
        "--identity",
        type=float,
        default=90,
        help="最小 identity 阈值 (默认: 90)",
    )

    # 可视化参数
    parser.add_argument(
        "--min-aln-len",
        type=int,
        default=100,
        help="最小比对长度阈值 bp，过滤 max(LEN1, LEN2) < 此值的比对块 (默认: 100)",
    )
    parser.add_argument(
        "--merge-gap",
        type=int,
        default=1000,
        help="合并间距阈值 bp，间距小于此值的比对块合并为一组 (默认: 1000)",
    )
    parser.add_argument(
        "--candidate-limit",
        type=int,
        default=3,
        help="每个 query 用于 query-query 预计算的候选数量 (默认: 3)",
    )
    parser.add_argument(
        "--pairwise-all",
        action="store_true",
        help="query-query 预计算使用所有通过阈值的候选，覆盖 --candidate-limit",
    )

    # Gene ID 模式专用参数：上下游延伸
    parser.add_argument(
        "-u", "--upstream",
        type=int,
        default=0,
        help="上游延伸长度 bp（Gene ID 为参考基因，Sequence 为查询命中片段，默认: 0）",
    )
    parser.add_argument(
        "-d", "--downstream",
        type=int,
        default=0,
        help="下游延伸长度 bp（Gene ID 为参考基因，Sequence 为查询命中片段，默认: 0）",
    )

    # 输出
    parser.add_argument("-o", "--output", required=True, help="输出目录")

    args = parser.parse_args()

    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)

    # 校验 -u/-d 参数非负
    if args.upstream < 0:
        raise ValueError(f"-u/--upstream 参数不能为负数，当前值: {args.upstream}")
    if args.downstream < 0:
        raise ValueError(f"-d/--downstream 参数不能为负数，当前值: {args.downstream}")
    if args.candidate_limit < 1:
        raise ValueError(f"--candidate-limit 参数必须大于 0，当前值: {args.candidate_limit}")

    # 获取基因组路径（支持注释版本选择）
    ra_source = getattr(args, 'ra_source', None)
    query_entries = _parse_query_entries(args.qry)
    ref_entry = None
    ref_genome = None
    ref_annotation = None
    ref_name = ""
    if args.ref:
        ref_entry_raw = _parse_ref_entry(args.ref)
        ref_entry = _resolve_genome_entry(ref_entry_raw, annotation_source=ra_source, role="参考基因组")
        ref_genome = ref_entry.fasta
        ref_annotation = ref_entry.annotation
        ref_name = ref_entry.name
        if args.ra:
            ref_annotation = args.ra
    elif args.gid or args.loc:
        raise ValueError("Gene ID/Location 模式需要指定 -ref 参考基因组")

    # 根据模式处理
    if args.seq:
        sequence_target = ref_entry
        resolved_queries = None
        if query_entries:
            resolved_queries = [
                _resolve_genome_entry(entry, role="目标基因组")
                for entry in query_entries
            ]
            sequence_target = resolved_queries[0]
        if not sequence_target:
            raise ValueError("Sequence 模式需要指定 -qry 查询基因组")
        processor = SequenceProcessor(
            sequence_target.fasta, args.output, sequence_target.name, args.identity,
            ref_gff=sequence_target.annotation,
            min_aln_len=args.min_aln_len, merge_gap=args.merge_gap,
            query_entries=resolved_queries,
            upstream=args.upstream,
            downstream=args.downstream,
            candidate_limit=args.candidate_limit,
            pairwise_all=args.pairwise_all
        )
        processor.process(args.seq)

    elif args.gid:
        if not query_entries:
            raise ValueError("Gene ID 模式需要指定 -qry 目标基因组")
        if not ref_annotation:
            raise ValueError("Gene ID 模式需要参考基因组注释文件")

        resolved_queries = [
            _resolve_genome_entry(entry, role="目标基因组")
            for entry in query_entries
        ]
        query_entry = resolved_queries[0]
        qry_genome = query_entry.fasta
        qry_annotation = query_entry.annotation
        qry_name = query_entry.name

        processor = GeneIDProcessor(
            ref_genome, ref_annotation, qry_genome, args.output, ref_name, qry_name, args.identity,
            qry_gff=qry_annotation, upstream=args.upstream, downstream=args.downstream,
            min_aln_len=args.min_aln_len, merge_gap=args.merge_gap,
            query_entries=resolved_queries,
            candidate_limit=args.candidate_limit,
            pairwise_all=args.pairwise_all
        )

        # 判断是文件还是 Gene ID 列表
        gene_ids = []
        if len(args.gid) == 1 and os.path.exists(args.gid[0]):
            # 从文件读取（跳过注释行和空行）
            with open(args.gid[0], encoding="utf-8-sig") as f:
                gene_ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        else:
            # 直接使用命令行参数
            gene_ids = args.gid

        multi_mode = len(gene_ids) > 1
        for gene_id in gene_ids:
            if multi_mode:
                item_output_dir = os.path.join(args.output, sanitize_path_segment(gene_id) or "gene")
                os.makedirs(item_output_dir, exist_ok=True)
                processor.set_output_dir(item_output_dir)
            processor.process(gene_id)

    elif args.loc:
        if not query_entries:
            raise ValueError("Location 模式需要指定 -qry 目标基因组")

        # 提示 -u/-d 在 Location 模式下无效
        if args.upstream > 0 or args.downstream > 0:
            print(f"[WARNING] -u/-d 参数仅在 Gene ID 模式下有效，当前 Location 模式将忽略这些参数")

        resolved_queries = [
            _resolve_genome_entry(entry, role="目标基因组")
            for entry in query_entries
        ]
        query_entry = resolved_queries[0]
        qry_genome = query_entry.fasta
        qry_annotation = query_entry.annotation
        qry_name = query_entry.name

        processor = LocationProcessor(
            ref_genome, qry_genome, args.output, ref_name, qry_name, args.identity,
            ref_gff=ref_annotation, qry_gff=qry_annotation,
            min_aln_len=args.min_aln_len, merge_gap=args.merge_gap,
            query_entries=resolved_queries,
            candidate_limit=args.candidate_limit,
            pairwise_all=args.pairwise_all
        )

        # 判断是文件还是区域字符串列表
        locations = []
        
        if len(args.loc) == 1 and os.path.exists(args.loc[0]):
            # 从文件读取
            with open(args.loc[0], encoding="utf-8-sig") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 3:
                        locations.append({
                            'chrom': parts[0],
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'name': parts[3] if len(parts) > 3 else None
                        })
        else:
            # 解析区域字符串列表
            for loc_str in args.loc:
                match = re.match(r"(\w+):(\d+)-(\d+)", loc_str)
                if match:
                    locations.append({
                        'chrom': match.group(1),
                        'start': int(match.group(2)),
                        'end': int(match.group(3)),
                        'name': None
                    })
                else:
                    print(f"[WARNING] 无效的位置格式: {loc_str}，跳过")
        
        if not locations:
            raise ValueError("未找到有效的位置输入")
        
        multi_mode = len(locations) > 1
        for loc in locations:
            if multi_mode:
                loc_name = loc['name'] or f"{loc['chrom']}_{loc['start']}_{loc['end']}"
                item_output_dir = os.path.join(args.output, sanitize_path_segment(loc_name) or "location")
                os.makedirs(item_output_dir, exist_ok=True)
                processor.set_output_dir(item_output_dir)
            processor.process(loc['chrom'], loc['start'], loc['end'], loc['name'])

    print(f"\n[INFO] 处理完成，结果保存在: {args.output}")


if __name__ == "__main__":
    main()
