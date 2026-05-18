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
from pathlib import Path

from pyfaidx import Fasta
from Bio.Blast import NCBIXML

from GeneScreenVisualizer import GeneIDVisualizer, LocationVisualizer, SequenceVisualizer
from GenomeManager import GenomeManager


# 全局基因组管理器
genome_manager = GenomeManager()


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

    def __init__(self, ref_genome, ref_annotation, query_genome, output_dir, ref_name="", qry_name="", identity=90, qry_gff=None, upstream=0, downstream=0, min_aln_len=100, merge_gap=1000):
        self.extractor = SequenceExtractor(ref_genome, ref_annotation, output_dir)
        self.aligner = BlastAligner(output_dir)
        self.query_genome = query_genome
        self.output_dir = output_dir
        self.ref_name = ref_name
        self.qry_name = qry_name
        self.identity = identity
        self.upstream = upstream
        self.downstream = downstream
        self.min_aln_len = min_aln_len
        self.merge_gap = merge_gap
        # 保存基因组文件路径
        self.genome_files = {
            'ref_fasta': ref_genome,
            'ref_gff': ref_annotation,
            'qry_fasta': query_genome,
            'qry_gff': qry_gff  # 查询基因组 GFF（用于提取比对区域的注释）
        }
        # 查询基因组索引文件（用于宏观概览图）
        self.qry_index_file = f"{query_genome}.fai" if query_genome else None

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
        result = self.aligner.align(self.query_genome, fasta_file, gene_id)
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


class LocationProcessor:
    """Location 模式处理器"""

    def __init__(self, ref_genome, query_genome, output_dir, ref_name="", qry_name="", identity=90, ref_gff=None, qry_gff=None, min_aln_len=100, merge_gap=1000):
        self.extractor = SequenceExtractor(ref_genome, None, output_dir)
        self.aligner = BlastAligner(output_dir)
        self.ref_genome = ref_genome
        self.query_genome = query_genome
        self.output_dir = output_dir
        self.ref_name = ref_name
        self.qry_name = qry_name
        self.identity = identity
        self.min_aln_len = min_aln_len
        self.merge_gap = merge_gap
        # 保存基因组文件路径
        self.genome_files = {
            'ref_fasta': ref_genome,
            'qry_fasta': query_genome,
            'ref_gff': ref_gff,  # 参考基因组 GFF（用于提取区域内的基因注释）
            'qry_gff': qry_gff   # 查询基因组 GFF（用于提取比对区域的注释）
        }
        # 查询基因组索引文件（用于宏观概览图）
        self.qry_index_file = f"{query_genome}.fai" if query_genome else None

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
        result = self.aligner.align(self.query_genome, fasta_file, loc_name)
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


class SequenceProcessor:
    """Sequence 模式处理器"""

    def __init__(self, ref_genome, output_dir, ref_name="", identity=90, ref_gff=None, min_aln_len=100, merge_gap=1000):
        self.aligner = BlastAligner(output_dir)
        self.ref_genome = ref_genome
        self.output_dir = output_dir
        self.ref_name = ref_name
        self.identity = identity
        self.min_aln_len = min_aln_len
        self.merge_gap = merge_gap
        # 保存基因组文件路径
        self.genome_files = {
            'ref_fasta': ref_genome,
            'ref_gff': ref_gff  # 参考基因组 GFF（用于提取比对区域的注释）
        }
        # 参考基因组索引文件（用于宏观概览图）
        self.ref_index_file = f"{ref_genome}.fai" if ref_genome else None

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
        for i, seq in enumerate(sequences):
            seq_id = seq["seq_id"]
            sequence = seq["sequence"]
            print(f"\n[INFO] 处理序列 {seq_id} ({i+1}/{len(sequences)})")
            
            result = self._process_single(seq_id, sequence)
            if result:
                results.append(result)
        
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
        result = self.aligner.align(self.ref_genome, fasta_file, safe_seq_id)
        if not result:
            return None

        result["fasta"] = fasta_file
        result["id"] = safe_seq_id
        result["sequence"] = sequence

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

  # Location 模式（单个）
  python GeneScreen.py -ref Nippon -qry ZS97 -loc Chr1:1000-2000 -o output/

  # Location 模式（多个）
  python GeneScreen.py -ref Nippon -qry ZS97 -loc Chr1:1000-2000 Chr2:3000-4000 -o output/

  # Location 模式（从文件读取）
  python GeneScreen.py -ref Nippon -qry ZS97 -loc positions.txt -o output/

  # Sequence 模式（支持多序列 FASTA）
  python GeneScreen.py -ref Nippon -seq query.fasta -o output/
        """,
    )

    # 基因组参数
    parser.add_argument("-ref", required=True, help="参考基因组（内置名称或路径）")
    parser.add_argument("-qry", help="目标基因组（Gene ID/Location 模式必需）")
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

    # Gene ID 模式专用参数：上下游延伸
    parser.add_argument(
        "-u", "--upstream",
        type=int,
        default=0,
        help="上游延伸长度 bp（仅 Gene ID 模式有效，默认: 0）",
    )
    parser.add_argument(
        "-d", "--downstream",
        type=int,
        default=0,
        help="下游延伸长度 bp（仅 Gene ID 模式有效，默认: 0）",
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

    # 获取基因组路径（支持注释版本选择）
    ra_source = getattr(args, 'ra_source', None)
    ref_genome, ref_annotation = genome_manager.get(args.ref, annotation_source=ra_source)
    if not ref_genome:
        print(f"\n基因组 '{args.ref}' 不可用。")
        print("可以使用以下命令管理基因组:")
        print(f"  python GenomeManager.py search {args.ref}")
        print(f"  python GenomeManager.py download <id>")
        print(f"  python GenomeManager.py add <name> <fasta>")
        return

    if args.ra:
        ref_annotation = args.ra

    # 根据模式处理
    if args.seq:
        # 提示 -u/-d 在 Sequence 模式下无效
        if args.upstream > 0 or args.downstream > 0:
            print(f"[WARNING] -u/-d 参数仅在 Gene ID 模式下有效，当前 Sequence 模式将忽略这些参数")
        processor = SequenceProcessor(ref_genome, args.output, args.ref, args.identity, ref_gff=ref_annotation, 
                                       min_aln_len=args.min_aln_len, merge_gap=args.merge_gap)
        processor.process(args.seq)

    elif args.gid:
        if not args.qry:
            raise ValueError("Gene ID 模式需要指定 -qry 目标基因组")
        if not ref_annotation:
            raise ValueError("Gene ID 模式需要参考基因组注释文件")

        qry_genome, qry_annotation = genome_manager.get(args.qry)
        if not qry_genome:
            print(f"[ERROR] 目标基因组 '{args.qry}' 不可用")
            return

        processor = GeneIDProcessor(ref_genome, ref_annotation, qry_genome, args.output, args.ref, args.qry, args.identity, 
                                     qry_gff=qry_annotation, upstream=args.upstream, downstream=args.downstream,
                                     min_aln_len=args.min_aln_len, merge_gap=args.merge_gap)

        # 判断是文件还是 Gene ID 列表
        gene_ids = []
        if len(args.gid) == 1 and os.path.exists(args.gid[0]):
            # 从文件读取（跳过注释行和空行）
            with open(args.gid[0], encoding="utf-8-sig") as f:
                gene_ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        else:
            # 直接使用命令行参数
            gene_ids = args.gid

        for gene_id in gene_ids:
            processor.process(gene_id)

    elif args.loc:
        if not args.qry:
            raise ValueError("Location 模式需要指定 -qry 目标基因组")

        # 提示 -u/-d 在 Location 模式下无效
        if args.upstream > 0 or args.downstream > 0:
            print(f"[WARNING] -u/-d 参数仅在 Gene ID 模式下有效，当前 Location 模式将忽略这些参数")

        qry_genome, qry_annotation = genome_manager.get(args.qry)
        if not qry_genome:
            print(f"[ERROR] 目标基因组 '{args.qry}' 不可用")
            return

        processor = LocationProcessor(ref_genome, qry_genome, args.output, args.ref, args.qry, args.identity, 
                                        ref_gff=ref_annotation, qry_gff=qry_annotation,
                                        min_aln_len=args.min_aln_len, merge_gap=args.merge_gap)

        # 判断是文件还是区域字符串列表
        import re
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
        
        for loc in locations:
            processor.process(loc['chrom'], loc['start'], loc['end'], loc['name'])

    print(f"\n[INFO] 处理完成，结果保存在: {args.output}")


if __name__ == "__main__":
    main()
