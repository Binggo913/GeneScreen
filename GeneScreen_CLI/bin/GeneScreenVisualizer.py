#!/usr/bin/env python3
"""
GeneScreen 1.0 可视化与报告模块

为三种输入模式提供独立的可视化处理逻辑，生成 HTML 报告
支持 LINKVIEW 可视化工具生成比对图
"""

import os
import json
import datetime
import subprocess
import base64
import re
from dataclasses import dataclass, field
from typing import List, Optional, Dict

from multi_query_result import build_single_query_payload, write_static_report


# ============================================================
# 数据结构定义
# ============================================================

@dataclass
class AlignmentBlock:
    """
    单个比对块
    
    对应 coords 文件中的一行记录
    字段来源: [S1] [E1] [S2] [E2] [LEN1] [LEN2] [%IDY] [REF] [QUERY]
    """
    ref_start: int      # S1: 参考序列起始位置
    ref_end: int        # E1: 参考序列结束位置
    qry_start: int      # S2: 查询序列起始位置
    qry_end: int        # E2: 查询序列结束位置
    ref_aln_len: int    # LEN1: 参考序列比对长度
    qry_aln_len: int    # LEN2: 查询序列比对长度
    identity: float     # %IDY: 比对一致性百分比
    ref_name: str       # REF: 参考序列名（染色体名）
    qry_name: str       # QUERY: 查询序列名


@dataclass
class AlignmentGroup:
    """
    比对组（多个相近的比对块）
    
    将同一染色体上间距小于 merge_gap 的比对块合并为一组
    """
    chr_name: str                           # 染色体名
    blocks: List[AlignmentBlock] = field(default_factory=list)  # 该组包含的比对块
    region_start: int = 0                   # 组的起始位置（最小 ref_start）
    region_end: int = 0                     # 组的结束位置（最大 ref_end）
    total_aln_len: int = 0                  # 总比对长度（sum of LEN1）


@dataclass
class VisualizationResult:
    """
    可视化结果
    
    包含宏观概览图和多个详细比对图的路径
    """
    overview_path: Optional[str] = None     # 宏观概览图路径
    detail_paths: List[str] = field(default_factory=list)       # 详细图路径列表
    detail_groups: List[AlignmentGroup] = field(default_factory=list)  # 对应的比对组信息


class LinkviewVisualizer:
    """
    LINKVIEW 可视化生成器
    使用 LINKVIEW.py 工具生成比对可视化图
    """
    
    def __init__(self, output_dir: str, min_aln_len: int = 100,
                 merge_gap: int = 1000, ref_index_file: Optional[str] = None,
                 min_identity: float = 90.0):
        """
        初始化可视化器
        
        Args:
            output_dir: 输出目录
            min_aln_len: 最小比对长度阈值，过滤 max(LEN1, LEN2) < 此值的比对块 (默认 100bp)
            merge_gap: 合并间距阈值，间距小于此值的比对块合并为一组 (默认 1000bp)
            ref_index_file: 参考基因组索引文件路径 (.fai)，用于读取染色体长度
            min_identity: 最小一致性阈值，过滤 %IDY < 此值的比对块 (默认 90)
        """
        self.output_dir = output_dir
        self.min_aln_len = min_aln_len
        self.merge_gap = merge_gap
        self.ref_index_file = ref_index_file
        self.min_identity = min_identity
        # 颜色配置
        self.colors = {
            'snp': 'orange',    # SNP 颜色（橙色）
            'indel': 'blue',    # Indel 颜色（蓝色）
        }
    
    def _get_chromosome_lengths(self) -> Dict[str, int]:
        """
        从 ref_index_file (.fai) 读取染色体长度
        
        .fai 格式: name\tlength\toffset\tlinebases\tlinewidth
        
        Returns:
            Dict[str, int]: {chr_name: length, ...}
        """
        chr_lengths: Dict[str, int] = {}
        
        if not self.ref_index_file or not os.path.exists(self.ref_index_file):
            return chr_lengths
        
        try:
            with open(self.ref_index_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        chr_name = parts[0]
                        try:
                            length = int(parts[1])
                            chr_lengths[chr_name] = length
                        except ValueError:
                            continue
        except Exception as e:
            print(f"[ERROR] 读取染色体长度失败: {e}")
        
        return chr_lengths
    
    def _parse_coords_file(self, coords_file: str) -> Dict[str, List[AlignmentBlock]]:
        """
        解析 coords 文件，返回按染色体分组的比对块
        
        coords 文件格式（TSV，首行表头以 [S1] 开头）:
        [S1]  [E1]  [S2]  [E2]  [LEN1]  [LEN2]  [%IDY]  [REF]  [QUERY]
        
        过滤规则: int((LEN1 + LEN2) / 2) > self.min_aln_len 且 %IDY > self.min_identity
        
        Returns:
            Dict[str, List[AlignmentBlock]]: {chr_name: [AlignmentBlock, ...], ...}
        """
        chr_alignments: Dict[str, List[AlignmentBlock]] = {}
        
        if not coords_file or not os.path.exists(coords_file):
            return chr_alignments
        
        try:
            with open(coords_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    # 跳过空行、表头行（以 [ 开头）、注释行（以 # 开头）
                    if not line or line.startswith('[') or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 9:
                        continue
                    
                    try:
                        ref_start = int(parts[0])
                        ref_end = int(parts[1])
                        qry_start = int(parts[2])
                        qry_end = int(parts[3])
                        ref_aln_len = int(parts[4]) if parts[4] else abs(ref_end - ref_start) + 1
                        qry_aln_len = int(parts[5]) if parts[5] else abs(qry_end - qry_start) + 1
                        identity = float(parts[6]) if parts[6] else 0.0
                        
                        # 过滤短比对和低一致性比对（与 LINKVIEW 规则一致）
                        alignment_length = int((ref_aln_len + qry_aln_len) / 2)
                        if alignment_length <= self.min_aln_len:
                            continue
                        if identity <= self.min_identity:
                            continue
                        
                        # REF 字段可能包含额外信息，使用 _extract_chr_name 统一处理
                        ref_name = self._extract_chr_name(parts[7]) if parts[7] else "chr"
                        qry_name = parts[8] if parts[8] else "query"
                        
                        block = AlignmentBlock(
                            ref_start=ref_start,
                            ref_end=ref_end,
                            qry_start=qry_start,
                            qry_end=qry_end,
                            ref_aln_len=ref_aln_len,
                            qry_aln_len=qry_aln_len,
                            identity=identity,
                            ref_name=ref_name,
                            qry_name=qry_name
                        )
                        
                        if ref_name not in chr_alignments:
                            chr_alignments[ref_name] = []
                        chr_alignments[ref_name].append(block)
                        
                    except (ValueError, IndexError) as e:
                        # 跳过解析失败的行
                        continue
                        
        except Exception as e:
            print(f"[ERROR] 解析 coords 文件失败: {e}")
        
        return chr_alignments
    
    def _cluster_alignments(self, blocks: List[AlignmentBlock], chr_name: str) -> List[AlignmentGroup]:
        """
        将同一染色体上的比对块按间距分组
        
        算法：
        1. 按 ref_start 升序排序
        2. 遍历比对块，如果与前一个块的间距 < merge_gap，合并到同一组
        3. 否则创建新组
        4. 返回按 total_aln_len 降序排列的组列表
        
        Args:
            blocks: 同一染色体上的比对块列表
            chr_name: 染色体名
            
        Returns:
            List[AlignmentGroup]: 按 total_aln_len 降序排列的比对组列表
        """
        if not blocks:
            return []
        
        # 按 ref_start 升序排序（处理反向比对，取 min(ref_start, ref_end)）
        sorted_blocks = sorted(blocks, key=lambda b: min(b.ref_start, b.ref_end))
        
        groups: List[AlignmentGroup] = []
        current_group: Optional[AlignmentGroup] = None
        
        for block in sorted_blocks:
            block_start = min(block.ref_start, block.ref_end)
            block_end = max(block.ref_start, block.ref_end)
            
            if current_group is None:
                # 创建第一个组
                current_group = AlignmentGroup(
                    chr_name=chr_name,
                    blocks=[block],
                    region_start=block_start,
                    region_end=block_end,
                    total_aln_len=block.ref_aln_len
                )
            else:
                # 计算与当前组的间距
                gap = block_start - current_group.region_end
                
                if gap < self.merge_gap:
                    # 间距小于阈值，合并到当前组
                    current_group.blocks.append(block)
                    current_group.region_end = max(current_group.region_end, block_end)
                    current_group.total_aln_len += block.ref_aln_len
                else:
                    # 间距大于阈值，保存当前组，创建新组
                    groups.append(current_group)
                    current_group = AlignmentGroup(
                        chr_name=chr_name,
                        blocks=[block],
                        region_start=block_start,
                        region_end=block_end,
                        total_aln_len=block.ref_aln_len
                    )
        
        # 保存最后一个组
        if current_group is not None:
            groups.append(current_group)
        
        # 按 total_aln_len 降序排列
        groups.sort(key=lambda g: g.total_aln_len, reverse=True)
        
        return groups
    
    def generate_k_file(self, gene_id, fasta_file, coords_file, 
                        alignment_group: Optional[AlignmentGroup] = None,
                        output_suffix: str = ""):
        """
        生成 .k 文件（定义显示区域）
        
        .k 文件格式：每行一个区域，格式为 name:start:end
        - 第一行：基因序列（相对坐标）
        - 第二行：主比对染色体区域
        
        Args:
            gene_id: 基因 ID
            fasta_file: 基因序列 FASTA 文件
            coords_file: 比对坐标文件
            alignment_group: 指定的比对组（可选），如果指定则仅输出该组区域
            output_suffix: 输出文件后缀（可选），用于区分不同组的文件
        
        逻辑：
        - 如果指定 alignment_group，使用该组的区域
        - 否则找到总比对长度最长的染色体，合并该染色体上的所有比对块范围
        """
        suffix = f".{output_suffix}" if output_suffix else ""
        k_file = os.path.join(self.output_dir, f"{gene_id}{suffix}.k")
        
        try:
            # 获取基因序列长度
            gene_len = 0
            with open(fasta_file, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        gene_len += len(line.strip())
            
            # 从 coords 文件解析比对区域，按染色体分组
            chr_alignments = {}  # chr -> [(start, end, aln_len), ...]
            
            if coords_file and os.path.exists(coords_file):
                with open(coords_file, 'r') as cf:
                    for line in cf:
                        line = line.strip()
                        if not line or line.startswith('[') or line.startswith('#'):
                            continue
                        parts = line.split('\t')
                        if len(parts) >= 8:
                            try:
                                ref_start = int(parts[0])
                                ref_end = int(parts[1])
                                ref_aln_len = int(parts[4]) if len(parts) > 4 else abs(ref_end - ref_start) + 1
                                qry_aln_len = int(parts[5]) if len(parts) > 5 else ref_aln_len
                                
                                # 过滤短比对：使用 max(LEN1, LEN2) 与 self.min_aln_len 比较
                                if max(ref_aln_len, qry_aln_len) < self.min_aln_len:
                                    continue
                                    
                                ref_name = parts[7].split()[0] if parts[7] else "chr"
                                chr_name = self._extract_chr_name(ref_name)
                                
                                if chr_name not in chr_alignments:
                                    chr_alignments[chr_name] = []
                                chr_alignments[chr_name].append((min(ref_start, ref_end), max(ref_start, ref_end), ref_aln_len))
                            except (ValueError, IndexError):
                                continue
            
            # 找到总比对长度最长的染色体
            main_chr = None
            max_total_len = 0
            for chr_name, alns in chr_alignments.items():
                total_len = sum(a[2] for a in alns)
                if total_len > max_total_len:
                    max_total_len = total_len
                    main_chr = chr_name
            
            with open(k_file, 'w') as kf:
                # 第一行：基因序列
                kf.write(f"{gene_id}:1:{gene_len}\n")
                
                # 第二行：主比对染色体区域
                if alignment_group:
                    # 使用指定的比对组区域
                    kf.write(f"{alignment_group.chr_name}:{alignment_group.region_start}:{alignment_group.region_end}\n")
                elif main_chr and chr_alignments[main_chr]:
                    # 合并该染色体上的所有比对块
                    alns = chr_alignments[main_chr]
                    min_start = min(a[0] for a in alns)
                    max_end = max(a[1] for a in alns)
                    kf.write(f"{main_chr}:{min_start}:{max_end}\n")
            
            print(f"[INFO] 已生成 .k 文件: {k_file}")
            return k_file
        except Exception as e:
            print(f"[ERROR] 生成 .k 文件失败: {e}")
            return None
    
    def generate_hl_file(self, gene_id, snps_file,
                         alignment_group: Optional[AlignmentGroup] = None,
                         output_suffix: str = ""):
        """
        生成 .hl 文件（高亮标记）
        
        .hl 文件格式：BED 格式 + 颜色
        chr  start  end  color
        
        Args:
            gene_id: 基因 ID
            snps_file: SNP/Indel 文件
            alignment_group: 指定的比对组（可选），如果指定则仅输出该组区域内的变异
            output_suffix: 输出文件后缀（可选），用于区分不同组的文件
        
        标记逻辑：
        - SNP：两个轨道都标注（两边都有碱基）
        - INS（查询相对 ref 插入）：只在查询基因组轨道标注
        - DEL（查询相对 ref 缺失）：只在 ref 来源序列轨道标注
        """
        suffix = f".{output_suffix}" if output_suffix else ""
        hl_file = os.path.join(self.output_dir, f"{gene_id}{suffix}.hl")
        
        try:
            with open(hl_file, 'w') as hf:
                if snps_file and os.path.exists(snps_file):
                    with open(snps_file, 'r') as sf:
                        for line in sf:
                            line = line.strip()
                            if not line or line.startswith('[') or line.startswith('#'):
                                continue
                            parts = line.split('\t')
                            if len(parts) >= 5:
                                try:
                                    ref_pos = int(parts[0])
                                    qry_pos = int(parts[3])
                                    var_type = parts[4].upper() if len(parts) > 4 else 'SNP'
                                    
                                    # 获取染色体/序列名
                                    ref_name = parts[5].split()[0] if len(parts) > 5 else "ref"
                                    qry_name = parts[6] if len(parts) > 6 else gene_id
                                    
                                    chr_name = self._extract_chr_name(ref_name)
                                    
                                    # 如果指定了 alignment_group，过滤区域外的变异
                                    if alignment_group:
                                        # 检查染色体是否匹配
                                        if chr_name != alignment_group.chr_name:
                                            continue
                                        # 检查位置是否在区域内
                                        if ref_pos < alignment_group.region_start or ref_pos > alignment_group.region_end:
                                            continue
                                    
                                    # 根据类型选择颜色
                                    if var_type == 'SNP':
                                        color = self.colors['snp']
                                    else:  # INS, DEL
                                        color = self.colors['indel']
                                    
                                    # 根据变异类型决定标注位置
                                    if var_type == 'SNP':
                                        # SNP：两个轨道都标注
                                        hf.write(f"{qry_name}\t{qry_pos-1}\t{qry_pos}\t{color}\n")
                                        hf.write(f"{chr_name}\t{ref_pos-1}\t{ref_pos}\t{color}\n")
                                    elif var_type == 'INS':
                                        # INS：查询基因组相对 ref 多出的碱基，标在查询基因组轨道
                                        hf.write(f"{chr_name}\t{ref_pos-1}\t{ref_pos}\t{color}\n")
                                    elif var_type == 'DEL':
                                        # DEL：查询基因组相对 ref 缺失的碱基，标在 ref 来源序列轨道
                                        hf.write(f"{qry_name}\t{qry_pos-1}\t{qry_pos}\t{color}\n")
                                except (ValueError, IndexError):
                                    continue
            
            print(f"[INFO] 已生成 .hl 文件: {hl_file}")
            return hl_file
        except Exception as e:
            print(f"[ERROR] 生成 .hl 文件失败: {e}")
            return None
    
    def _extract_chr_name(self, ref_name):
        """从 ref_name 中提取染色体名称
        
        LINKVIEW 解析 coords 文件时只取空格前的第一个字段，
        所以这里也只返回第一个字段以保持一致
        """
        if not ref_name:
            return "chr"
        # 只取第一个空格前的部分，与 LINKVIEW 解析逻辑一致
        return ref_name.split()[0] if ref_name else "chr"
    
    def _generate_linkview_input(self, gene_id, coords_file, fasta_file):
        """
        生成 LINKVIEW nucmer coords 格式的输入文件 (-t 2)
        
        格式: S1 E1 | S2 E2 | LEN1 LEN2 | %IDY | LEN_R LEN_Q | COV_R COV_Q | TAGS
        只输出主染色体上的比对（总比对长度最长的染色体）
        """
        linkview_input = os.path.join(self.output_dir, f"{gene_id}.linkview.coords")
        
        # 获取查询序列长度
        qry_len = 0
        with open(fasta_file, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    qry_len += len(line.strip())
        
        try:
            # 按染色体分组
            chr_alignments = {}  # chr -> [(ref_start, ref_end, qry_start, qry_end, aln_len, identity, qry_name, ref_full_name), ...]
            
            with open(coords_file, 'r') as cf:
                for line in cf:
                    line = line.strip()
                    if not line or line.startswith('[') or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 9:
                        try:
                            ref_start = int(parts[0])
                            ref_end = int(parts[1])
                            qry_start = int(parts[2])
                            qry_end = int(parts[3])
                            ref_aln_len = int(parts[4]) if len(parts) > 4 else abs(ref_end - ref_start) + 1
                            qry_aln_len = int(parts[5]) if len(parts) > 5 else abs(qry_end - qry_start) + 1
                            identity = float(parts[6]) if len(parts) > 6 else 0.0
                            
                            # 过滤短比对：与 LINKVIEW 的长度规则一致
                            alignment_length = int((ref_aln_len + qry_aln_len) / 2)
                            if alignment_length <= self.min_aln_len:
                                continue
                            if identity <= self.min_identity:
                                continue
                            if identity < self.min_identity:
                                continue
                            
                            ref_full_name = parts[7] if len(parts) > 7 else "chr"
                            ref_name = self._extract_chr_name(ref_full_name)
                            qry_name = parts[8]
                            
                            if ref_name not in chr_alignments:
                                chr_alignments[ref_name] = []
                            chr_alignments[ref_name].append((ref_start, ref_end, qry_start, qry_end, 
                                                           ref_aln_len, qry_aln_len, identity, qry_name, ref_full_name))
                        except (ValueError, IndexError):
                            continue
            
            # 找到总比对长度最长的染色体
            main_chr = None
            max_total_len = 0
            for chr_name, alns in chr_alignments.items():
                total_len = sum(a[4] for a in alns)
                if total_len > max_total_len:
                    max_total_len = total_len
                    main_chr = chr_name

            chr_lengths = self._get_chromosome_lengths()
            if not chr_lengths:
                print(f"[WARNING] 未读取到染色体长度(.fai)，跳过 LINKVIEW 输入: {linkview_input}")
                return None
            
            # 写入 nucmer coords 格式
            line_count = 0
            with open(linkview_input, 'w') as lf:
                # 写入头部（LINKVIEW 需要跳过前5行）
                lf.write("ref.fasta query.fasta\n")
                lf.write("NUCMER\n")
                lf.write("\n")
                lf.write("    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n")
                lf.write("=" * 120 + "\n")
                
                if main_chr and chr_alignments[main_chr]:
                    ref_len = chr_lengths.get(main_chr)
                    if not ref_len:
                        print(f"[WARNING] 未找到染色体长度: {main_chr}")
                        return None
                    for ref_start, ref_end, qry_start, qry_end, ref_aln_len, qry_aln_len, identity, qry_name, ref_full_name in chr_alignments[main_chr]:
                        cov_r = ref_aln_len / ref_len * 100
                        cov_q = qry_aln_len / qry_len * 100
                        
                        lf.write(f"{ref_start:>8} {ref_end:>8}  | {qry_start:>8} {qry_end:>8}  | {ref_aln_len:>8} {qry_aln_len:>8}  | {identity:>8.2f}  | {ref_len:>8} {qry_len:>8}  | {cov_r:>8.2f} {cov_q:>8.2f}  | {main_chr}\t{qry_name}\n")
                        line_count += 1

            if line_count == 0:
                print(f"[WARNING] 未生成有效的 LINKVIEW 比对行: {linkview_input}")
                return None
            
            return linkview_input
        except Exception as e:
            print(f"[ERROR] 生成 LINKVIEW 输入文件失败: {e}")
            return None
    
    def run_linkview(self, gene_id, coords_file, gff_file, k_file, hl_file, output_prefix, fasta_file):
        """
        运行 LINKVIEW.py 生成可视化图
        
        返回生成的图片路径（SVG 或 PNG）
        """
        output_svg = f"{output_prefix}.svg"
        output_png = f"{output_prefix}.png"
        
        # 生成 LINKVIEW 专用输入文件 (nucmer coords 格式)
        linkview_input = self._generate_linkview_input(gene_id, coords_file, fasta_file)
        if not linkview_input:
            return None
        
        # 获取同目录下的 LINKVIEW.py 路径
        script_dir = os.path.dirname(os.path.abspath(__file__))
        linkview_path = os.path.join(script_dir, "LINKVIEW.py")
        
        # 构建 LINKVIEW 命令 (使用 -t 2 nucmer coords 格式)
        cmd = [
            "python", linkview_path,
            "-t", "2",
            linkview_input,
            "-k", k_file,
            "-hl", hl_file,
            "-o", output_prefix,
            "--min_identity", str(self.min_identity),
            "--min_alignment_length", str(self.min_aln_len),
            "--svg_height", "400",
            "--chro_axis",
            "--bezier",
            "-s",
            "--style", "simple"
        ]
        
        # 只有存在 GFF 文件时才添加 -g 参数
        if gff_file and os.path.exists(gff_file):
            cmd.extend(["-g", gff_file])
        
        cmd_str = " ".join(cmd)
        print(f"[CMD] {cmd_str}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"[WARNING] LINKVIEW 运行失败: {result.stderr}")
                return None
            
            # 检查输出文件
            if os.path.exists(output_svg):
                print(f"[INFO] LINKVIEW 生成 SVG: {output_svg}")
                return output_svg
            elif os.path.exists(output_png):
                print(f"[INFO] LINKVIEW 生成 PNG: {output_png}")
                return output_png
            else:
                print(f"[WARNING] LINKVIEW 未生成图片文件")
                return None
        except FileNotFoundError:
            print(f"[WARNING] LINKVIEW.py 未找到: {linkview_path}")
            return None
        except Exception as e:
            print(f"[WARNING] LINKVIEW 运行异常: {e}")
            return None
    
    def generate_visualization(self, gene_id, fasta_file, coords_file, snps_file, gff_file) -> VisualizationResult:
        """
        完整的可视化流程
        
        统一走多图流程，1个组就生成1个图，多个组就生成多个图
        
        Returns:
            VisualizationResult: 包含所有生成的图片路径和对应的比对组信息
        """
        result = VisualizationResult()
        
        # 1. 解析 coords 文件
        chr_alignments = self._parse_coords_file(coords_file)
        
        if not chr_alignments:
            print(f"[WARNING] 没有找到有效的比对数据: {coords_file}")
            return result
        
        # 2. 对每个染色体进行聚类分组
        all_groups: List[AlignmentGroup] = []
        for chr_name, blocks in chr_alignments.items():
            groups = self._cluster_alignments(blocks, chr_name)
            all_groups.extend(groups)
        
        if not all_groups:
            print(f"[WARNING] 没有生成有效的比对组")
            return result
        
        # 3. 按 total_aln_len 降序排列
        all_groups.sort(key=lambda g: g.total_aln_len, reverse=True)
        
        # 4. 为每个组分配索引（同染色体内按 region_start 升序，从 1 开始）
        chr_group_indices: Dict[str, List[AlignmentGroup]] = {}
        for group in all_groups:
            if group.chr_name not in chr_group_indices:
                chr_group_indices[group.chr_name] = []
            chr_group_indices[group.chr_name].append(group)
        
        # 用 id(group) 作为 key，避免 AlignmentGroup 不可 hash 的问题
        group_index_map: Dict[int, int] = {}
        for chr_name, groups in chr_group_indices.items():
            sorted_groups = sorted(groups, key=lambda g: g.region_start)
            for idx, group in enumerate(sorted_groups, start=1):  # 从 1 开始
                group_index_map[id(group)] = idx
        
        # 5. 为每个组生成可视化图
        for group in all_groups:
            group_index = group_index_map[id(group)]
            image_path = self._generate_group_visualization(
                gene_id, fasta_file, coords_file, snps_file, gff_file,
                group, group_index
            )
            if image_path:
                result.detail_paths.append(image_path)
                result.detail_groups.append(group)
        
        # 6. 生成宏观概览图
        if result.detail_groups:
            overview_path = self.generate_overview_chart(gene_id, result.detail_groups)
            if overview_path:
                result.overview_path = overview_path
        
        return result
    
    def generate_overview_chart(self, gene_id: str, groups: List[AlignmentGroup]) -> Optional[str]:
        """
        生成宏观概览图
        
        显示所有染色体轨道，标注比对区域位置
        
        Args:
            gene_id: 基因 ID
            groups: 所有比对组列表
        
        Returns:
            生成的 SVG 文件路径，失败返回 None
        
        布局规则：
        - 轨道按总比对长度降序排列
        - 最长染色体铺满宽度，其它按比例缩放
        - 比对片段最小宽度 3px
        """
        if not groups:
            return None
        
        # 创建输出目录
        gene_dir = os.path.join(self.output_dir, gene_id)
        overview_dir = os.path.join(gene_dir, "overview")
        os.makedirs(overview_dir, exist_ok=True)
        
        output_path = os.path.join(overview_dir, f"{gene_id}.overview.svg")
        
        # 获取染色体长度
        chr_lengths = self._get_chromosome_lengths()
        if not chr_lengths:
            print(f"[WARNING] 未能读取染色体长度，将使用比对区域估算")
        
        # 按染色体汇总比对信息（仅使用 .fai 提供的长度）
        chr_info: Dict[str, dict] = {}  # {chr_name: {length, total_aln_len, groups}}
        for group in groups:
            chr_name = group.chr_name
            if chr_name not in chr_lengths:
                continue
            if chr_name not in chr_info:
                chr_info[chr_name] = {
                    'length': chr_lengths[chr_name],
                    'total_aln_len': 0,
                    'groups': []
                }
            chr_info[chr_name]['total_aln_len'] += group.total_aln_len
            chr_info[chr_name]['groups'].append(group)
        
        # 按总比对长度降序排列染色体
        sorted_chrs = sorted(chr_info.keys(), key=lambda c: chr_info[c]['total_aln_len'], reverse=True)
        
        # SVG 参数 - 使用百分比布局，让 SVG 自适应容器宽度
        # viewBox 使用 1000 作为基准宽度，实际渲染时会按容器宽度缩放
        svg_width = 1000
        
        # 在绘图区域内分配空间（使用百分比概念）
        label_width = 80  # 染色体名称标签宽度
        len_label_width = 60  # 长度标签宽度
        padding_x = 20  # 左右边距
        draw_width = svg_width - label_width - len_label_width - padding_x * 2  # 轨道可用宽度
        track_x = padding_x + label_width  # 轨道起始 X 坐标
        
        padding_top = 15
        padding_bottom = 15
        track_height = 20
        track_gap = 20
        
        # 找到最长染色体用于缩放
        max_chr_len = max(chr_info[c]['length'] for c in sorted_chrs) if sorted_chrs else 1
        
        # 计算 SVG 高度
        svg_height = padding_top + padding_bottom + len(sorted_chrs) * (track_height + track_gap) - track_gap + track_height
        
        # 生成 SVG - 使用 width="100%" 让 SVG 填满容器宽度
        svg_parts = []
        svg_parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="100%" viewBox="0 0 {svg_width} {svg_height}" preserveAspectRatio="xMidYMid meet">')
        svg_parts.append('<style>')
        svg_parts.append('  .track { fill: #e0e0e0; stroke: #999; stroke-width: 1; }')
        svg_parts.append('  .alignment { fill: #4a90d9; stroke: none; cursor: pointer; }')
        svg_parts.append('  .alignment:hover { fill: #2d6cb5; }')
        svg_parts.append('  .alignment-primary { fill: #2563eb; stroke: #1d4ed8; stroke-width: 2; cursor: pointer; }')
        svg_parts.append('  .alignment-primary:hover { fill: #1d4ed8; }')
        svg_parts.append('  .label { font-family: Arial, sans-serif; font-size: 12px; fill: #333; }')
        svg_parts.append('  .primary-label { font-family: Arial, sans-serif; font-size: 9px; fill: #1d4ed8; font-weight: bold; }')
        svg_parts.append('</style>')
        
        # 找出主比对区域（总比对长度最长的 group）
        primary_group = groups[0] if groups else None
        
        # 绘制每个染色体轨道
        y = padding_top
        for chr_name in sorted_chrs:
            info = chr_info[chr_name]
            chr_len = info['length']
            
            # 计算轨道宽度（最长染色体占满 draw_width，其他按比例缩放）
            track_width = draw_width * chr_len / max_chr_len
            
            # 染色体标签
            svg_parts.append(f'<text x="{track_x - 10}" y="{y + track_height/2 + 4}" text-anchor="end" class="label">{chr_name}</text>')
            
            # 染色体轨道背景
            svg_parts.append(f'<rect x="{track_x}" y="{y}" width="{track_width}" height="{track_height}" class="track" rx="3"/>')
            
            # 染色体长度标签（轨道右侧）
            len_label = f"{chr_len/1e6:.1f}Mb" if chr_len >= 1e6 else f"{chr_len/1e3:.1f}kb"
            svg_parts.append(f'<text x="{track_x + track_width + 5}" y="{y + track_height/2 + 4}" class="label" font-size="10" fill="#666">{len_label}</text>')
            
            # 绘制比对区域
            for group in info['groups']:
                # 计算位置
                x_start = track_x + (group.region_start / chr_len) * track_width
                x_end = track_x + (group.region_end / chr_len) * track_width
                width = max(x_end - x_start, 3)  # 最小宽度 3px
                
                # 判断是否为主比对区域
                is_primary = (primary_group is not None and 
                             group.chr_name == primary_group.chr_name and 
                             group.region_start == primary_group.region_start and 
                             group.region_end == primary_group.region_end)
                
                # 比对区域矩形（主比对使用特殊样式）
                rect_class = 'alignment-primary' if is_primary else 'alignment'
                svg_parts.append(f'<rect x="{x_start}" y="{y + 2}" width="{width}" height="{track_height - 4}" class="{rect_class}" '
                               f'data-chr="{chr_name}" data-start="{group.region_start}" data-end="{group.region_end}" '
                               f'data-len="{group.total_aln_len}" data-primary="{str(is_primary).lower()}"/>')
                
                # 为主比对添加"主/Primary"标签（支持国际化）
                if is_primary:
                    label_x = x_start + width / 2
                    label_y = y - 2
                    svg_parts.append(f'<text x="{label_x}" y="{label_y}" text-anchor="middle" class="primary-label" data-zh="主" data-en="Primary">主</text>')
            
            y += track_height + track_gap
        
        svg_parts.append('</svg>')
        
        # 写入文件
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write('\n'.join(svg_parts))
            print(f"[INFO] 已生成宏观比对图: {output_path}")
            return output_path
        except Exception as e:
            print(f"[ERROR] 生成宏观比对图失败: {e}")
            return None

    def _generate_group_visualization(self, gene_id: str, fasta_file: str, coords_file: str, 
                                       snps_file: str, gff_file: str, 
                                       group: AlignmentGroup, group_index: int) -> Optional[str]:
        """
        为单个比对组生成可视化图
        
        目录结构：
        output_dir/{gene_id}/
        ├── report/              # HTML 报告（任务 6 实现）
        ├── overview/            # 宏观概览图（任务 5 实现）
        ├── detail/{chr_name}/   # 局部图、.k、.hl
        └── linkview/            # 中间文件
        
        Args:
            gene_id: 基因 ID
            fasta_file: 基因序列 FASTA 文件
            coords_file: 比对坐标文件
            snps_file: SNP/Indel 文件
            gff_file: GFF 注释文件
            group: 比对组
            group_index: 组索引（同染色体内按 region_start 升序的索引，从 1 开始）
        
        Returns:
            生成的图片路径，失败返回 None
        
        文件命名: gene_id.{group_index}.svg
        """
        # 创建目录结构
        gene_dir = os.path.join(self.output_dir, gene_id)
        detail_dir = os.path.join(gene_dir, "detail", group.chr_name)
        linkview_dir = os.path.join(gene_dir, "linkview")
        report_dir = os.path.join(gene_dir, "report")
        overview_dir = os.path.join(gene_dir, "overview")
        
        os.makedirs(detail_dir, exist_ok=True)
        os.makedirs(linkview_dir, exist_ok=True)
        os.makedirs(report_dir, exist_ok=True)
        os.makedirs(overview_dir, exist_ok=True)
        
        # 文件名基础：gene_id.group_index
        file_base = f"{gene_id}.{group_index}"
        
        # 1. 生成该组的 .k 文件（存放在 detail 目录）
        k_file = self._generate_k_file_for_group(
            gene_id, fasta_file, group, 
            os.path.join(detail_dir, f"{file_base}.k")
        )
        if not k_file:
            return None
        
        # 2. 生成该组的 .hl 文件（存放在 detail 目录）
        hl_file = self._generate_hl_file_for_group(
            gene_id, snps_file, group,
            os.path.join(detail_dir, f"{file_base}.hl")
        )
        if not hl_file:
            return None
        
        # 3. 生成 LINKVIEW 输入文件（存放在 linkview 目录）
        linkview_input = self._generate_linkview_input_for_group(
            gene_id, coords_file, fasta_file, group,
            os.path.join(linkview_dir, f"{file_base}.linkview.coords")
        )
        if not linkview_input:
            return None
        
        # 4. 运行 LINKVIEW（输出到 detail 目录）
        output_prefix = os.path.join(detail_dir, file_base)
        return self._run_linkview_with_input(
            linkview_input, gff_file, k_file, hl_file, output_prefix
        )
    
    def _generate_k_file_for_group(self, gene_id: str, fasta_file: str, 
                                    group: AlignmentGroup, output_path: str) -> Optional[str]:
        """为指定比对组生成 .k 文件"""
        try:
            # 获取基因序列长度
            gene_len = 0
            with open(fasta_file, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        gene_len += len(line.strip())
            
            with open(output_path, 'w') as kf:
                kf.write(f"{gene_id}:1:{gene_len}\n")
                kf.write(f"{group.chr_name}:{group.region_start}:{group.region_end}\n")
            
            print(f"[INFO] 已生成 .k 文件: {output_path}")
            return output_path
        except Exception as e:
            print(f"[ERROR] 生成 .k 文件失败: {e}")
            return None
    
    def _generate_hl_file_for_group(self, gene_id: str, snps_file: str,
                                     group: AlignmentGroup, output_path: str) -> Optional[str]:
        """为指定比对组生成 .hl 文件（仅包含区域内的变异）"""
        try:
            with open(output_path, 'w') as hf:
                if snps_file and os.path.exists(snps_file):
                    with open(snps_file, 'r') as sf:
                        for line in sf:
                            line = line.strip()
                            if not line or line.startswith('[') or line.startswith('#'):
                                continue
                            parts = line.split('\t')
                            if len(parts) >= 5:
                                try:
                                    ref_pos = int(parts[0])
                                    qry_pos = int(parts[3])
                                    var_type = parts[4].upper() if len(parts) > 4 else 'SNP'
                                    
                                    ref_name = parts[5].split()[0] if len(parts) > 5 else "ref"
                                    qry_name = parts[6] if len(parts) > 6 else gene_id
                                    chr_name = self._extract_chr_name(ref_name)
                                    
                                    # 过滤：只保留该组区域内的变异
                                    if chr_name != group.chr_name:
                                        continue
                                    if ref_pos < group.region_start or ref_pos > group.region_end:
                                        continue
                                    
                                    color = self.colors['snp'] if var_type == 'SNP' else self.colors['indel']
                                    
                                    if var_type == 'SNP':
                                        hf.write(f"{qry_name}\t{qry_pos-1}\t{qry_pos}\t{color}\n")
                                        hf.write(f"{chr_name}\t{ref_pos-1}\t{ref_pos}\t{color}\n")
                                    elif var_type == 'INS':
                                        hf.write(f"{chr_name}\t{ref_pos-1}\t{ref_pos}\t{color}\n")
                                    elif var_type == 'DEL':
                                        hf.write(f"{qry_name}\t{qry_pos-1}\t{qry_pos}\t{color}\n")
                                except (ValueError, IndexError):
                                    continue
            
            print(f"[INFO] 已生成 .hl 文件: {output_path}")
            return output_path
        except Exception as e:
            print(f"[ERROR] 生成 .hl 文件失败: {e}")
            return None
    
    def _generate_linkview_input_for_group(self, gene_id: str, coords_file: str, 
                                            fasta_file: str, group: AlignmentGroup,
                                            output_path: str) -> Optional[str]:
        """为指定比对组生成 LINKVIEW 输入文件"""
        try:
            # 获取查询序列长度
            qry_len = 0
            with open(fasta_file, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        qry_len += len(line.strip())
            
            line_count = 0
            chr_lengths = self._get_chromosome_lengths()
            ref_len = chr_lengths.get(group.chr_name) if chr_lengths else None
            if not ref_len:
                print(f"[WARNING] 未找到染色体长度: {group.chr_name}")
                return None

            with open(output_path, 'w') as lf:
                # 写入头部
                lf.write("ref.fasta query.fasta\n")
                lf.write("NUCMER\n")
                lf.write("\n")
                lf.write("    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n")
                lf.write("=" * 120 + "\n")
                
                # 只输出该组的比对块
                for block in group.blocks:
                    cov_r = block.ref_aln_len / ref_len * 100
                    cov_q = block.qry_aln_len / qry_len * 100
                    
                    lf.write(f"{block.ref_start:>8} {block.ref_end:>8}  | {block.qry_start:>8} {block.qry_end:>8}  | "
                            f"{block.ref_aln_len:>8} {block.qry_aln_len:>8}  | {block.identity:>8.2f}  | "
                            f"{ref_len:>8} {qry_len:>8}  | {cov_r:>8.2f} {cov_q:>8.2f}  | "
                            f"{group.chr_name}\t{block.qry_name}\n")
                    line_count += 1

            if line_count == 0:
                print(f"[WARNING] 未生成有效的 LINKVIEW 比对行: {output_path}")
                return None
            
            print(f"[INFO] 已生成 LINKVIEW 输入文件: {output_path}")
            return output_path
        except Exception as e:
            print(f"[ERROR] 生成 LINKVIEW 输入文件失败: {e}")
            return None
    
    def _run_linkview_with_input(self, linkview_input: str, gff_file: str,
                                  k_file: str, hl_file: str, output_prefix: str) -> Optional[str]:
        """使用指定的输入文件运行 LINKVIEW"""
        output_svg = f"{output_prefix}.svg"
        output_png = f"{output_prefix}.png"
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        linkview_path = os.path.join(script_dir, "LINKVIEW.py")
        
        cmd = [
            "python", linkview_path,
            "-t", "2",
            linkview_input,
            "-k", k_file,
            "-hl", hl_file,
            "-o", output_prefix,
            "--min_identity", str(self.min_identity),
            "--min_alignment_length", str(self.min_aln_len),
            "--svg_height", "400",
            "--chro_axis",
            "--bezier",
            "-s",
            "--style", "simple"
        ]
        
        if gff_file and os.path.exists(gff_file):
            cmd.extend(["-g", gff_file])
        
        cmd_str = " ".join(cmd)
        print(f"[CMD] {cmd_str}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"[WARNING] LINKVIEW 运行失败: {result.stderr}")
                return None
            
            if os.path.exists(output_svg):
                print(f"[INFO] LINKVIEW 生成 SVG: {output_svg}")
                return output_svg
            elif os.path.exists(output_png):
                print(f"[INFO] LINKVIEW 生成 PNG: {output_png}")
                return output_png
            else:
                print(f"[WARNING] LINKVIEW 未生成图片文件")
                return None
        except FileNotFoundError:
            print(f"[WARNING] LINKVIEW.py 未找到: {linkview_path}")
            return None
        except Exception as e:
            print(f"[WARNING] LINKVIEW 运行异常: {e}")
            return None

class BaseVisualizer:
    """可视化基类"""

    def __init__(self, output_dir: str, min_aln_len: int = 100,
                 merge_gap: int = 1000, ref_index_file: Optional[str] = None,
                 min_identity: float = 90.0):
        """
        初始化可视化基类
        
        Args:
            output_dir: 输出目录
            min_aln_len: 最小比对长度阈值 (默认 100bp)
            merge_gap: 合并间距阈值 (默认 1000bp)
            ref_index_file: 参考基因组索引文件路径 (.fai)
            min_identity: 最小一致性阈值 (默认 90)
        """
        self.output_dir = output_dir
        self.min_aln_len = min_aln_len
        self.merge_gap = merge_gap
        self.ref_index_file = ref_index_file
        self.min_identity = min_identity
        self.linkview_viz = LinkviewVisualizer(
            output_dir,
            min_aln_len,
            merge_gap,
            ref_index_file,
            min_identity,
        )

    def visualize(self, result):
        """可视化入口，子类实现具体逻辑"""
        raise NotImplementedError

    def _legacy_query_entry(self, result, genome_files, index, query_result):
        query_entries = result.get("queries") or (genome_files or {}).get("queries") or []
        if index < len(query_entries):
            return query_entries[index] or {}
        return {
            "name": query_result.get("query_name") or f"query_{index + 1}",
            "fasta": query_result.get("query_fasta"),
            "gff": query_result.get("query_gff"),
        }

    def _legacy_query_safe(self, query_name, fallback):
        safe = re.sub(r'[<>:"/\\|?*\s]+', "_", str(query_name or fallback)).strip("._")
        return safe or fallback

    def _legacy_report_suffix_for_result(self, result, genome_files):
        query_results = result.get("query_results") or []
        if not query_results:
            return None
        query_entry = self._legacy_query_entry(result, genome_files, 0, query_results[0])
        query_name = query_entry.get("name") or query_results[0].get("query_name") or "query_1"
        return self._legacy_query_safe(query_name, "query_1")

    def _legacy_report_path(self, report_stem, report_suffix=None):
        if not report_suffix:
            return os.path.join(self.output_dir, f"{report_stem}.report.html")
        report_dir = os.path.join(self.output_dir, "report", "single_reports", report_suffix)
        os.makedirs(report_dir, exist_ok=True)
        return os.path.join(report_dir, f"{report_stem}.report.html")

    def _single_query_result(self, result, query_result, query_entry):
        single_result = dict(result)
        single_result.update(query_result)
        single_result["query_results"] = [query_result]
        single_result["queries"] = [query_entry]
        single_result["pairwise_results"] = []
        return single_result

    def _single_query_genome_files(self, mode, genome_files, query_entry, query_result):
        files = dict(genome_files or {})
        query_fasta = query_entry.get("fasta") or query_result.get("query_fasta")
        query_gff = query_entry.get("gff") or query_result.get("query_gff")
        if mode == "sequence":
            files["ref_fasta"] = query_fasta
            files["ref_gff"] = query_gff
        else:
            files["qry_fasta"] = query_fasta
            files["qry_gff"] = query_gff
        files["queries"] = [query_entry]
        return files

    def _build_legacy_report_map(
        self,
        result,
        mode,
        first_report_file,
        ref_genome="",
        qry_genome="",
        identity=90,
        input_source=None,
        genome_files=None,
    ):
        query_results = result.get("query_results") or []
        if len(query_results) <= 1:
            return first_report_file

        legacy_reports = {}
        for index, query_result in enumerate(query_results):
            query_entry = self._legacy_query_entry(result, genome_files, index, query_result)
            query_name = query_entry.get("name") or query_result.get("query_name") or f"query_{index + 1}"
            query_safe = self._legacy_query_safe(query_name, f"query_{index + 1}")
            if index == 0:
                legacy_reports[query_safe] = first_report_file
                continue

            single_result = self._single_query_result(result, query_result, query_entry)
            single_files = self._single_query_genome_files(mode, genome_files, query_entry, query_result)
            if mode == "sequence":
                legacy_reports[query_safe] = self.visualize(
                    single_result,
                    query_name,
                    None,
                    identity,
                    input_source,
                    single_files,
                    legacy_only=True,
                    report_suffix=query_safe,
                )
            else:
                legacy_reports[query_safe] = self.visualize(
                    single_result,
                    ref_genome,
                    query_name,
                    identity,
                    input_source,
                    single_files,
                    legacy_only=True,
                    report_suffix=query_safe,
                )
        return legacy_reports

    def _has_multiple_query_reports(self, result):
        """Return True only when the current input has multiple query genomes."""
        return len(result.get("query_results") or []) > 1

    def _finalize_single_query_report(
        self,
        result,
        mode,
        ref_name="",
        qry_name="",
        identity=90,
        genome_files=None,
        legacy_report=None,
        extra_files=None,
    ):
        """Write the shared N=1 multi-query payload and static report entry."""
        payload = build_single_query_payload(
            result=result,
            output_dir=self.output_dir,
            mode=mode,
            ref_name=ref_name,
            qry_name=qry_name,
            identity=identity,
            min_aln_len=self.min_aln_len,
            genome_files=genome_files,
            legacy_report=legacy_report,
            extra_files=extra_files,
        )
        report_index = write_static_report(payload, self.output_dir, legacy_report)
        result["multi_query_report"] = payload
        result["report_data"] = os.path.join(self.output_dir, "report", "data.json")
        result["report_index"] = report_index
        return report_index

    def _escape_path_for_js(self, path):
        """转义路径用于 JavaScript，兼容 Windows 和 Linux"""
        if not path:
            return ''
        # Windows 路径需要双重转义反斜杠
        return path.replace('\\', '\\\\').replace("'", "\\'")

    def _generate_legend_svg(self, svg_width, y_top):
        """
        生成图例 SVG 代码（右上角，支持点击切换显示/隐藏，带勾选框）
        
        参数:
        - svg_width: SVG 总宽度
        - y_top: 顶部 Y 坐标
        
        颜色与 LINKVIEW 生成的 SVG 保持一致：
        - SNP: orange (橙色)
        - Indel: blue (蓝色)
        - 5' UTR: #6B5B7B (深紫灰色)
        - 3' UTR: #B6AEC9 (浅紫灰色)
        - CDS: #7A7A7A (灰色)
        """
        # (label, color, shape, data-type)
        items = [
            ('SNP', 'orange', 'line', 'snp'),
            ('Indel', 'blue', 'line', 'indel'),
            ("5' UTR", '#6B5B7B', 'rect', 'utr5'),
            ("3' UTR", '#B6AEC9', 'rect', 'utr3'),
            ('CDS', '#7A7A7A', 'rect', 'cds'),
        ]
        
        # 图例框尺寸
        item_height = 20
        padding_x = 10
        padding_y = 8
        icon_width = 12
        text_offset = 20
        max_text_width = 60
        box_width = padding_x * 2 + icon_width + text_offset + max_text_width
        box_height = padding_y * 2 + len(items) * item_height
        
        # 图例框位置（右上角，距离右边 30px，与轨道保持距离）
        box_x = svg_width - box_width - 30
        box_y = y_top
        
        parts = []
        
        # 绘制图例背景框
        parts.append(f'''
            <rect x="{box_x}" y="{box_y}" width="{box_width}" height="{box_height}" 
                  fill="white" stroke="#ddd" stroke-width="1" rx="4" opacity="0.95"/>
        ''')
        
        # 绘制图例项（每项是一个可点击的组）
        current_y = box_y + padding_y + 10
        for label, color, shape, data_type in items:
            icon_x = box_x + padding_x
            
            # 包装成可点击的组
            parts.append(f'''
                <g class="legend-item" data-type="{data_type}" style="cursor: pointer;">
                    <rect x="{box_x + 2}" y="{current_y - 9}" width="{box_width - 4}" height="{item_height - 2}" 
                          fill="transparent" class="legend-hitarea"/>
            ''')
            
            if shape == 'line':
                # 竖线图标
                parts.append(f'''
                    <line class="legend-icon" x1="{icon_x + 6}" y1="{current_y - 6}" x2="{icon_x + 6}" y2="{current_y + 6}" 
                          stroke="{color}" stroke-width="3"/>
                ''')
            else:
                # 矩形图标
                parts.append(f'''
                    <rect class="legend-icon" x="{icon_x}" y="{current_y - 5}" width="{icon_width}" height="10" 
                          fill="{color}" rx="2"/>
                ''')
            
            # 文本标签
            parts.append(f'''
                    <text x="{icon_x + text_offset}" y="{current_y + 4}" 
                          font-family="Arial, sans-serif" font-size="10" fill="#666" class="legend-text">{label}</text>
                </g>
            ''')
            
            current_y += item_height
        
        return '\n'.join(parts)

    def _generate_gene_bracket_svg(self, svg_width, svg_height, gene_rel_start, gene_rel_end, gene_id, seq_len=None):
        """
        生成基因区域大括号 SVG 代码（用于 LINKVIEW 生成的 SVG）
        
        在上方轨道下边缘绘制水平大括号，标注基因区域范围
        
        参数:
        - svg_width: SVG 总宽度
        - svg_height: SVG 总高度
        - gene_rel_start: 基因相对起始位置 (1-based)
        - gene_rel_end: 基因相对结束位置 (1-based)
        - gene_id: 基因 ID（用于标签显示）
        - seq_len: 序列总长度（用于计算 X 坐标）
        
        返回: SVG 字符串
        """
        # LINKVIEW 生成的 SVG 布局参数（从实际 SVG 分析得出）
        # 上方轨道：x 从 120 到 1080，y 从 133.33 到 148.33
        padding = 120  # LINKVIEW 默认 padding
        draw_width = svg_width - 2 * padding  # 960
        
        # 上方轨道下边缘 Y 坐标
        track_bottom_y = 148.33
        
        # 使用传入的 seq_len，如果没有则用 gene_rel_end 作为近似
        if seq_len is None:
            seq_len = gene_rel_end
        
        # 如果 gene_rel_start == 1 且 gene_rel_end == seq_len，说明没有上下游，不绘制
        # 但这个判断应该在调用前完成，这里假设已经判断过需要绘制
        
        # 计算基因区域在 SVG 中的 X 坐标
        # 上方轨道显示范围是 [1, seq_len]
        bracket_x_start = padding + (gene_rel_start - 1) / seq_len * draw_width
        bracket_x_end = padding + gene_rel_end / seq_len * draw_width
        
        # 大括号参数
        bracket_y = track_bottom_y + 8  # 大括号顶部距离轨道下边缘
        bracket_height = 15
        curl_width = 10
        label_offset = 20
        color = "#333333"
        
        mid_x = (bracket_x_start + bracket_x_end) / 2
        bracket_width = bracket_x_end - bracket_x_start
        
        parts = []
        
        # 使用 SVG path 绘制水平大括号 ⎵
        if bracket_width > curl_width * 4:
            # 正常大括号
            path = f'''M {bracket_x_start},{bracket_y}
                       Q {bracket_x_start},{bracket_y + bracket_height} {bracket_x_start + curl_width},{bracket_y + bracket_height}
                       L {mid_x - curl_width},{bracket_y + bracket_height}
                       Q {mid_x},{bracket_y + bracket_height} {mid_x},{bracket_y + bracket_height + 6}
                       Q {mid_x},{bracket_y + bracket_height} {mid_x + curl_width},{bracket_y + bracket_height}
                       L {bracket_x_end - curl_width},{bracket_y + bracket_height}
                       Q {bracket_x_end},{bracket_y + bracket_height} {bracket_x_end},{bracket_y}'''
        else:
            # 窄区域：简化为弧形
            path = f'''M {bracket_x_start},{bracket_y}
                       Q {mid_x},{bracket_y + bracket_height + 6} {bracket_x_end},{bracket_y}'''
        
        parts.append(f'''
            <path d="{path}" fill="none" stroke="{color}" stroke-width="2"/>
        ''')
        
        # 基因 ID 标签
        label_y = bracket_y + bracket_height + label_offset
        parts.append(f'''
            <text x="{mid_x}" y="{label_y}" text-anchor="middle" 
                  font-family="Arial, sans-serif" font-size="12" fill="{color}" font-weight="500">
                {gene_id}
            </text>
        ''')
        
        return '\n'.join(parts)

    def _generate_linkview_visualization_section(self, gene_id, fasta_file, coords_file, snps_file, gff_file, qry_gff_file=None, extraction_info=None):
        """
        使用 LINKVIEW 生成比对可视化模块
        
        双区域布局：
        - 上方区域: 宏观概览图，默认固定显示
        - 下方区域: 局部比对图，默认显示总比对长度最长的组
        
        参数:
        - gff_file: 上方轨道的 GFF 文件（相对坐标）
        - qry_gff_file: 下方轨道的 GFF 文件（绝对坐标），可选
        - extraction_info: 提取信息（可选，包含 gene_rel_start, gene_rel_end 用于绘制大括号）
        """
        # 合并 GFF 文件用于 LINKVIEW 渲染
        merged_gff = self._merge_gff_files(gff_file, qry_gff_file, gene_id) if qry_gff_file else gff_file
        
        # 生成可视化（返回 VisualizationResult）
        viz_result = self.linkview_viz.generate_visualization(
            gene_id, fasta_file, coords_file, snps_file, merged_gff
        )
        
        if not viz_result.detail_paths:
            print(f"[WARNING] LINKVIEW 生成失败")
            return self._generate_visualization_error_section("LINKVIEW 可视化生成失败，请检查比对结果")
        
        # 计算有效索引（对应实际存在的 SVG 文件）
        valid_indices = []
        for i, path in enumerate(viz_result.detail_paths):
            if path.endswith('.svg') and os.path.exists(path):
                valid_indices.append(i)
        
        # 生成双区域布局
        html_parts = []
        
        # 1. 宏观概览图区域 + 区域统计卡片
        if viz_result.overview_path and os.path.exists(viz_result.overview_path):
            html_parts.append(self._embed_overview_section(viz_result, gene_id, valid_indices, snps_file, fasta_file))
        
        # 2. 局部比对图区域
        html_parts.append(self._embed_multi_image_section(
            viz_result, gene_id, snps_file,
            gff_file, qry_gff_file, fasta_file, extraction_info
        ))
        
        return '\n'.join(html_parts)
    
    def _embed_overview_section(self, viz_result: VisualizationResult, gene_id: str, 
                                 valid_indices: Optional[List[int]] = None,
                                 snps_file: str = None, fasta_file: str = None) -> str:
        """
        嵌入宏观概览图区域 + 区域统计卡片
        
        - 显示所有染色体轨道和比对区域
        - 点击比对片段可切换到对应的局部图
        - 在宏观图下方显示区域统计卡片（随选中区域动态更新）
        
        Args:
            viz_result: 可视化结果
            gene_id: 基因 ID
            valid_indices: 有效的组索引列表（对应实际存在的 SVG 文件）
            snps_file: SNP/Indel 文件路径（用于统计）
            fasta_file: FASTA 文件路径（用于计算基因长度）
        """
        if not viz_result.overview_path or not os.path.exists(viz_result.overview_path):
            return ""
        
        # 读取 SVG 内容
        with open(viz_result.overview_path, 'r', encoding='utf-8') as f:
            svg_content = f.read()
        
        # 计算基因长度
        gene_len = 0
        if fasta_file and os.path.exists(fasta_file):
            with open(fasta_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if not line.startswith('>'):
                        gene_len += len(line.strip())
        
        # 解析 SNP/Indel 数据，按区域统计
        all_variants = []
        if snps_file and os.path.exists(snps_file):
            with open(snps_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('[') or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 6:
                        try:
                            ref_pos = int(parts[0])
                            var_type = parts[4].upper() if len(parts) > 4 else 'SNP'
                            chr_name = self.linkview_viz._extract_chr_name(parts[5]) if len(parts) > 5 else ""
                            all_variants.append({
                                'ref_pos': ref_pos,
                                'type': var_type,
                                'chr_name': chr_name
                            })
                        except (ValueError, IndexError):
                            continue
        
        # 为每个区域计算统计数据
        valid_set = set(valid_indices) if valid_indices else set(range(len(viz_result.detail_groups)))
        region_stats = {}
        for i, group in enumerate(viz_result.detail_groups):
            if i not in valid_set:
                continue
            # 统计该区域内的 SNP 和 Indel
            snp_count = 0
            indel_count = 0
            for v in all_variants:
                if v['chr_name'] == group.chr_name and group.region_start <= v['ref_pos'] <= group.region_end:
                    if v['type'] == 'SNP':
                        snp_count += 1
                    else:
                        indel_count += 1
            
            region_stats[i] = {
                'chr_name': group.chr_name,
                'region_start': group.region_start,
                'region_end': group.region_end,
                'aln_len': group.total_aln_len,
                'snp_count': snp_count,
                'indel_count': indel_count
            }
        
        # 构建 group 索引映射
        group_map = {}
        for i, group in enumerate(viz_result.detail_groups):
            if i in valid_set:
                key = f"{group.chr_name}:{group.region_start}:{group.region_end}"
                group_map[key] = i
        
        # 默认显示第一个有效区域的统计
        default_idx = valid_indices[0] if valid_indices else 0
        default_stats = region_stats.get(default_idx, {})
        # 如果没有有效统计数据，使用空值兜底
        if not default_stats and viz_result.detail_groups:
            # 使用第一个 group 的信息兜底
            first_group = viz_result.detail_groups[0]
            default_stats = {
                'chr_name': first_group.chr_name,
                'region_start': first_group.region_start,
                'region_end': first_group.region_end,
                'aln_len': first_group.total_aln_len,
                'snp_count': 0,
                'indel_count': 0
            }
        default_region_name = f"{default_stats.get('chr_name', '')}:{default_stats.get('region_start', 0):,}-{default_stats.get('region_end', 0):,}"
        
        # 注入点击事件脚本（三联动：更新统计卡片 + 切换局部图 + 滚动）
        click_script = f'''
        <script>
        var groupMap = {json.dumps(group_map)};
        var regionStats = {json.dumps(region_stats)};
        var geneLen = {gene_len};
        
        window.updateRegionStats = function(idx) {{
            idx = parseInt(idx, 10);
            if (isNaN(idx)) return;
            var stats = regionStats[idx];
            if (!stats) return;
            
            var regionName = stats.chr_name + ':' + stats.region_start.toLocaleString() + '-' + stats.region_end.toLocaleString();
            var lang = document.documentElement.lang || 'zh';
            
            // 更新标题
            var titleEl = document.getElementById('region-stats-title');
            if (titleEl) {{
                titleEl.textContent = regionName + (lang === 'en' ? ' - Statistics' : '-结果统计');
                titleEl.setAttribute('data-zh', regionName + '-结果统计');
                titleEl.setAttribute('data-en', regionName + ' - Statistics');
            }}
            
            // 更新统计数据
            document.getElementById('stats-gene-len').textContent = geneLen.toLocaleString() + ' bp';
            document.getElementById('stats-aln-len').textContent = stats.aln_len.toLocaleString() + ' bp';
            document.getElementById('stats-snp-count').textContent = stats.snp_count;
            document.getElementById('stats-indel-count').textContent = stats.indel_count;
        }}
        
        document.querySelectorAll('.alignment, .alignment-primary').forEach(function(el) {{
            el.addEventListener('click', function() {{
                var chr = this.getAttribute('data-chr');
                var start = this.getAttribute('data-start');
                var end = this.getAttribute('data-end');
                var key = chr + ':' + start + ':' + end;
                if (groupMap[key] !== undefined) {{
                    var idx = groupMap[key];
                    
                    // 1. 更新统计卡片
                    updateRegionStats(idx);
                    
                    // 2. 切换局部比对图
                    var selector = document.getElementById('group-selector');
                    if (selector) {{
                        selector.value = idx;
                        switchVisualization(idx);
                    }}
                    
                    // 3. 不自动滚动页面
                }}
            }});
        }});
        </script>
        '''
        
        html = f'''
        <div class="section overview-section">
            <h2 data-zh="宏观比对图" data-en="Alignment Overview">宏观比对图</h2>
            <div class="overview-container" style="width:100%; overflow-x: auto;">
                {svg_content}
            </div>
            <div class="overview-note" style="font-size: 12px; color: #666; margin-top: 5px;">
                <p data-zh="点击比对区域可跳转到对应的局部比对图" data-en="Click alignment region to jump to the corresponding local alignment">点击比对区域可跳转到对应的局部比对图</p>
            </div>
        </div>
        
        <div class="section stats-section" id="region-stats-section">
            <h2 id="region-stats-title" data-zh="{default_region_name}-结果统计" data-en="{default_region_name} - Statistics">{default_region_name}-结果统计</h2>
            <table class="stats-table">
                <tr>
                    <th data-zh="基因长度" data-en="Gene Length">基因长度</th>
                    <th data-zh="比对长度" data-en="Alignment Length">比对长度</th>
                    <th data-zh="SNP 数量" data-en="SNP Count">SNP 数量</th>
                    <th data-zh="Indel 数量" data-en="Indel Count">Indel 数量</th>
                </tr>
                <tr>
                    <td id="stats-gene-len">{gene_len:,} bp</td>
                    <td id="stats-aln-len">{default_stats.get('aln_len', 0):,} bp</td>
                    <td id="stats-snp-count">{default_stats.get('snp_count', 0)}</td>
                    <td id="stats-indel-count">{default_stats.get('indel_count', 0)}</td>
                </tr>
            </table>
        </div>
        {click_script}
        '''
        
        return html
    
    def _embed_multi_image_section(self, viz_result: VisualizationResult, gene_id: str, 
                                    snps_file: str, gff_file: str, qry_gff_file: str,
                                    fasta_file: str, extraction_info: dict) -> str:
        """
        嵌入多图切换界面（局部比对图区域）
        
        - 默认显示第一个图（总比对长度最长的组）
        - 右上角下拉框切换不同组的图
        - 显示当前区域信息
        """
        # 读取所有 SVG 内容
        svg_contents = []
        for i, (path, group) in enumerate(zip(viz_result.detail_paths, viz_result.detail_groups)):
            if path.endswith('.svg') and os.path.exists(path):
                with open(path, 'r', encoding='utf-8') as f:
                    svg_content = f.read()
                # 处理 SVG（添加 viewBox、图例等）
                svg_content = self._process_svg_for_embed(svg_content, gene_id, extraction_info)
                svg_contents.append({
                    'content': svg_content,
                    'group': group,
                    'index': i
                })
        
        if not svg_contents:
            return self._generate_visualization_error_section("没有有效的可视化图片")
        
        # 默认显示第一个有效的 SVG（总比对长度最长的组）
        default_idx = svg_contents[0]['index']
        
        # 生成下拉选项
        options_html = ""
        for item in svg_contents:
            group = item['group']
            idx = item['index']
            label = f"{group.chr_name}:{group.region_start:,}-{group.region_end:,} ({group.total_aln_len:,}bp)"
            selected = "selected" if idx == default_idx else ""
            options_html += f'<option value="{idx}" {selected}>{label}</option>\n'
        
        # 生成 SVG 容器（默认只显示第一个有效的）
        svg_containers = ""
        for item in svg_contents:
            idx = item['index']
            group = item['group']
            display = "block" if idx == default_idx else "none"
            region_info_zh = f"{group.chr_name}:{group.region_start:,}-{group.region_end:,} (比对长度: {group.total_aln_len:,}bp)"
            region_info_en = f"{group.chr_name}:{group.region_start:,}-{group.region_end:,} (Alignment Length: {group.total_aln_len:,}bp)"
            svg_containers += f'''
                <div class="svg-container" id="svg-container-{idx}" style="display: {display};" data-region-zh="{region_info_zh}" data-region-en="{region_info_en}">
                    <div class="region-info" style="font-size: 12px; color: #666; margin-bottom: 5px;" data-zh="{region_info_zh}" data-en="{region_info_en}">{region_info_zh}</div>
                    {item['content']}
                </div>
            '''
        
        # 解析变异数据用于 tooltip
        variants = self._parse_snps_for_tooltip(snps_file) if snps_file else []
        ref_features = self._parse_gff_for_tooltip(gff_file, fasta_file) if gff_file else {'utr3': [], 'utr5': [], 'cds': []}
        qry_features = self._parse_gff_for_tooltip(qry_gff_file, None) if qry_gff_file else {'utr3': [], 'utr5': [], 'cds': []}
        
        for feat in ref_features['utr3'] + ref_features['utr5'] + ref_features['cds']:
            feat['track'] = 'top'
        for feat in qry_features['utr3'] + qry_features['utr5'] + qry_features['cds']:
            feat['track'] = 'bottom'
        
        features = {'top': ref_features, 'bottom': qry_features}
        variants_json = json.dumps(variants, ensure_ascii=False)
        features_json = json.dumps(features, ensure_ascii=False)
        
        # 国际化文本
        region_count_zh = f"共 {len(svg_contents)} 个比对区域，按总比对长度降序排列"
        region_count_en = f"{len(svg_contents)} alignment region(s), sorted by total alignment length"
        
        html = f"""
        <div class="section visualization-section">
            <h2 data-zh="局部比对图" data-en="Local Alignment">局部比对图</h2>
            <div class="viz-controls" style="margin-bottom: 10px; display: flex; justify-content: space-between; align-items: center;">
                <span style="font-size: 12px; color: #666;" data-zh="{region_count_zh}" data-en="{region_count_en}">{region_count_zh}</span>
                <div>
                    <label for="group-selector" data-zh="选择比对区域：" data-en="Select region:">选择比对区域：</label>
                    <select id="group-selector" onchange="switchVisualization(this.value)">
                        {options_html}
                    </select>
                </div>
            </div>
            <div class="viz-wrapper">
                <div class="viz-container" id="viz-container" style="width:100%; position: relative; overflow: hidden;">
                    <div id="svg-tooltip" class="svg-tooltip"></div>
                    {svg_containers}
                    <button id="zoom-btn" class="zoom-btn" title="Zoom">
                        <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <circle cx="11" cy="11" r="8"></circle>
                            <path d="M21 21l-4.35-4.35"></path>
                            <path d="M11 8v6M8 11h6"></path>
                        </svg>
                    </button>
                </div>
            </div>
            <div class="viz-legend-note">
                <p data-zh="上方轨道为基因序列，下方轨道为查询基因组比对区域。橙色标记表示 SNP，蓝色标记表示 Indel。鼠标悬停可查看详情。" 
                   data-en="Upper track shows gene sequence, lower track shows aligned region in query genome. Orange marks indicate SNPs, blue marks indicate Indels. Hover for details.">上方轨道为基因序列，下方轨道为查询基因组比对区域。橙色标记表示 SNP，蓝色标记表示 Indel。鼠标悬停可查看详情。</p>
            </div>
        </div>
        
        <!-- 模态框 -->
        <div id="zoom-modal" class="zoom-modal">
            <div class="zoom-modal-content">
                <button class="zoom-modal-close" id="zoom-modal-close">&times;</button>
                <div class="zoom-modal-body" id="zoom-modal-body"></div>
            </div>
        </div>
        
        <style>
            .viz-wrapper {{
                position: relative;
            }}
            .svg-tooltip {{
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
            }}
            .svg-tooltip.visible {{
                opacity: 1;
            }}
            .svg-tooltip .tooltip-title {{
                font-weight: bold;
                margin-bottom: 4px;
                color: #ffd700;
            }}
            .svg-tooltip .tooltip-row {{
                margin: 2px 0;
            }}
            .svg-tooltip .tooltip-label {{
                color: #aaa;
            }}
            .svg-tooltip .tooltip-seq {{
                font-family: monospace;
                background: rgba(255,255,255,0.1);
                padding: 1px 4px;
                border-radius: 2px;
            }}
            /* 放大按钮 */
            .zoom-btn {{
                position: absolute;
                bottom: 10px;
                right: 10px;
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
            }}
            .zoom-btn:hover {{
                background: rgba(102, 126, 234, 1);
                transform: scale(1.1);
            }}
            /* 模态框 */
            .zoom-modal {{
                display: none;
                position: fixed;
                top: 0;
                left: 0;
                width: 100%;
                height: 100%;
                background: rgba(0, 0, 0, 0.8);
                z-index: 10000;
                overflow: auto;
            }}
            .zoom-modal.visible {{
                display: flex;
                align-items: center;
                justify-content: center;
            }}
            .zoom-modal-content {{
                position: relative;
                background: white;
                border-radius: 8px;
                width: 95%;
                max-width: 1400px;
                max-height: 90vh;
                overflow: auto;
                padding: 20px;
            }}
            .zoom-modal-close {{
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
            }}
            .zoom-modal-close:hover {{
                color: #333;
            }}
            .zoom-modal-body {{
                width: 100%;
                overflow-x: auto;
                position: relative;
            }}
            .zoom-modal-body svg {{
                min-width: 1200px;
            }}
            .legend-item {{
                transition: opacity 0.2s;
            }}
            .legend-item:hover {{
                opacity: 0.8;
            }}
            .legend-item.disabled {{
                opacity: 0.4;
            }}
            .legend-item.disabled .legend-text {{
                text-decoration: line-through;
            }}
        </style>
        
        <script>
        function switchVisualization(index) {{
            var indexStr = String(index);
            // 隐藏所有 SVG 容器
            document.querySelectorAll('.svg-container').forEach(function(el) {{
                el.style.display = 'none';
            }});
            // 显示选中的容器
            var container = document.getElementById('svg-container-' + indexStr);
            if (container) {{
                container.style.display = 'block';
                // 重新初始化当前容器的交互
                initSvgInteraction(container.querySelector('svg'));
            }}
            // 同步更新统计卡片（如果 updateRegionStats 函数存在）
            if (typeof updateRegionStats === 'function') {{
                updateRegionStats(indexStr);
            }}
        }}
        
        (function() {{
            const variants = {variants_json};
            const features = {features_json};
            const tooltip = document.getElementById('svg-tooltip');
            
            // 显示/隐藏状态
            const visibility = {{
                snp: true,
                indel: true,
                utr5: true,
                utr3: true,
                cds: true
            }};
            
            // 国际化文本
            const i18n = {{
                zh: {{
                    position: '位置',
                    mutation: '变异',
                    refBase: '参考',
                    qryBase: '查询',
                    insertSeq: '插入序列',
                    deleteSeq: '缺失序列',
                    length: '长度',
                    sequence: '序列',
                    insertion: '插入 (INS)',
                    deletion: '缺失 (DEL)',
                    utr3: "3' UTR",
                    utr3Desc: '三端非翻译区',
                    utr5: "5' UTR",
                    utr5Desc: '五端非翻译区',
                    cds: 'CDS',
                    cdsDesc: '编码序列'
                }},
                en: {{
                    position: 'Position',
                    mutation: 'Mutation',
                    refBase: 'Ref',
                    qryBase: 'Qry',
                    insertSeq: 'Inserted Seq',
                    deleteSeq: 'Deleted Seq',
                    length: 'Length',
                    sequence: 'Sequence',
                    insertion: 'Insertion (INS)',
                    deletion: 'Deletion (DEL)',
                    utr3: "3' UTR",
                    utr3Desc: "3' Untranslated Region",
                    utr5: "5' UTR",
                    utr5Desc: "5' Untranslated Region",
                    cds: 'CDS',
                    cdsDesc: 'Coding Sequence'
                }}
            }};
            
            function getLang() {{
                return typeof currentLang !== 'undefined' ? currentLang : 'zh';
            }}
            
            function t(key) {{
                const lang = getLang();
                return i18n[lang] && i18n[lang][key] ? i18n[lang][key] : i18n['zh'][key];
            }}
            
            function showTooltip(content, e, tooltipEl) {{
                const tt = tooltipEl || tooltip;
                tt.innerHTML = content;
                tt.classList.add('visible');
                updateTooltipPosition(e, tt);
            }}
            
            function updateTooltipPosition(e, tooltipEl) {{
                const tt = tooltipEl || tooltip;
                const container = tt.parentElement;
                const containerRect = container.getBoundingClientRect();
                let x = e.clientX - containerRect.left + 15;
                let y = e.clientY - containerRect.top - 10;
                if (x + 200 > containerRect.width) {{
                    x = e.clientX - containerRect.left - 200;
                }}
                tt.style.left = x + 'px';
                tt.style.top = y + 'px';
            }}
            
            function hideTooltip(tooltipEl) {{
                const tt = tooltipEl || tooltip;
                tt.classList.remove('visible');
            }}
            
            // 初始化 SVG 交互
            window.initSvgInteraction = function(svg, tooltipEl) {{
                if (!svg) return;
                const tt = tooltipEl || tooltip;
                
                // 获取 SVG 的高度中点
                const svgRect = svg.getBBox ? svg.getBBox() : {{ height: 400 }};
                const midY = svgRect.height / 2;
                
                function getTrack(el) {{
                    const bbox = el.getBBox ? el.getBBox() : null;
                    if (!bbox) return 'top';
                    return bbox.y < midY ? 'top' : 'bottom';
                }}
                
                function getTrackFeatures(track, type) {{
                    if (features[track] && features[track][type]) {{
                        return features[track][type];
                    }}
                    return [];
                }}
                
                // 按轨道分组统计索引
                let topUtr3Idx = 0, bottomUtr3Idx = 0;
                let topUtr5Idx = 0, bottomUtr5Idx = 0;
                let topCdsIdx = 0, bottomCdsIdx = 0;
                
                // UTR3 交互
                svg.querySelectorAll('.UTR3').forEach((el) => {{
                    el.style.cursor = 'pointer';
                    const track = getTrack(el);
                    const trackFeatures = getTrackFeatures(track, 'utr3');
                    const idx = track === 'top' ? topUtr3Idx++ : bottomUtr3Idx++;
                    
                    el.addEventListener('mouseenter', function(e) {{
                        const feat = trackFeatures[idx] || {{}};
                        let content = `<div class="tooltip-title">${{t('utr3')}}</div>`;
                        content += `<div class="tooltip-row">${{t('utr3Desc')}}</div>`;
                        if (feat.start) {{
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{feat.start}} - ${{feat.end}}</div>`;
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{feat.length}} bp</div>`;
                            if (feat.seq) {{
                                content += `<div class="tooltip-row"><span class="tooltip-label">${{t('sequence')}}:</span> <span class="tooltip-seq">${{feat.seq}}</span></div>`;
                            }}
                        }}
                        showTooltip(content, e, tt);
                    }});
                    el.addEventListener('mousemove', (e) => updateTooltipPosition(e, tt));
                    el.addEventListener('mouseleave', () => hideTooltip(tt));
                }});
                
                // UTR5 交互
                svg.querySelectorAll('.UTR5').forEach((el) => {{
                    el.style.cursor = 'pointer';
                    const track = getTrack(el);
                    const trackFeatures = getTrackFeatures(track, 'utr5');
                    const idx = track === 'top' ? topUtr5Idx++ : bottomUtr5Idx++;
                    
                    el.addEventListener('mouseenter', function(e) {{
                        const feat = trackFeatures[idx] || {{}};
                        let content = `<div class="tooltip-title">${{t('utr5')}}</div>`;
                        content += `<div class="tooltip-row">${{t('utr5Desc')}}</div>`;
                        if (feat.start) {{
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{feat.start}} - ${{feat.end}}</div>`;
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{feat.length}} bp</div>`;
                            if (feat.seq) {{
                                content += `<div class="tooltip-row"><span class="tooltip-label">${{t('sequence')}}:</span> <span class="tooltip-seq">${{feat.seq}}</span></div>`;
                            }}
                        }}
                        showTooltip(content, e, tt);
                    }});
                    el.addEventListener('mousemove', (e) => updateTooltipPosition(e, tt));
                    el.addEventListener('mouseleave', () => hideTooltip(tt));
                }});
                
                // CDS/exon 交互
                svg.querySelectorAll('.exon').forEach((el) => {{
                    el.style.cursor = 'pointer';
                    const track = getTrack(el);
                    const trackFeatures = getTrackFeatures(track, 'cds');
                    const idx = track === 'top' ? topCdsIdx++ : bottomCdsIdx++;
                    
                    el.addEventListener('mouseenter', function(e) {{
                        const feat = trackFeatures[idx] || {{}};
                        let content = `<div class="tooltip-title">${{t('cds')}}</div>`;
                        content += `<div class="tooltip-row">${{t('cdsDesc')}}</div>`;
                        if (feat.start) {{
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{feat.start}} - ${{feat.end}}</div>`;
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{feat.length}} bp</div>`;
                            if (feat.seq) {{
                                content += `<div class="tooltip-row"><span class="tooltip-label">${{t('sequence')}}:</span> <span class="tooltip-seq">${{feat.seq}}</span></div>`;
                            }}
                        }}
                        showTooltip(content, e, tt);
                    }});
                    el.addEventListener('mousemove', (e) => updateTooltipPosition(e, tt));
                    el.addEventListener('mouseleave', () => hideTooltip(tt));
                }});
                
                // SNP/Indel 交互
                if (variants.length > 0) {{
                    const orangeRects = Array.from(svg.querySelectorAll('rect[fill="orange"]'));
                    const blueRects = Array.from(svg.querySelectorAll('rect[fill="blue"]'));
                    const snps = variants.filter(v => v.type === 'SNP');
                    const indels = variants.filter(v => v.type !== 'SNP');
                    
                    // 创建连线元素（用于 SNP 悬停时显示）
                    let snpConnectLine = null;
                    
                    // 绘制 SNP 连线的函数
                    function drawSnpConnectLine(rect1, rect2) {{
                        if (snpConnectLine) {{
                            snpConnectLine.remove();
                        }}
                        
                        // 获取两个 rect 的中心坐标
                        const x1 = parseFloat(rect1.getAttribute('x')) + parseFloat(rect1.getAttribute('width')) / 2;
                        const y1 = parseFloat(rect1.getAttribute('y')) + parseFloat(rect1.getAttribute('height')) / 2;
                        const x2 = parseFloat(rect2.getAttribute('x')) + parseFloat(rect2.getAttribute('width')) / 2;
                        const y2 = parseFloat(rect2.getAttribute('y')) + parseFloat(rect2.getAttribute('height')) / 2;
                        
                        // 创建连线
                        snpConnectLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                        snpConnectLine.setAttribute('x1', x1);
                        snpConnectLine.setAttribute('y1', y1);
                        snpConnectLine.setAttribute('x2', x2);
                        snpConnectLine.setAttribute('y2', y2);
                        snpConnectLine.setAttribute('stroke', 'orange');
                        snpConnectLine.setAttribute('stroke-width', '2');
                        snpConnectLine.setAttribute('stroke-dasharray', '4,2');
                        snpConnectLine.setAttribute('pointer-events', 'none');
                        snpConnectLine.classList.add('snp-connect-line');
                        svg.appendChild(snpConnectLine);
                    }}
                    
                    // 移除连线的函数
                    function removeSnpConnectLine() {{
                        if (snpConnectLine) {{
                            snpConnectLine.remove();
                            snpConnectLine = null;
                        }}
                    }}
                    
                    orangeRects.forEach((rect, idx) => {{
                        rect.style.cursor = 'pointer';
                        rect.classList.add('variant-snp');
                        
                        rect.addEventListener('mouseenter', function(e) {{
                            const variant = snps[Math.floor(idx / 2)];
                            if (variant) {{
                                const content = `
                                    <div class="tooltip-title">SNP</div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{variant.qry_pos}}</div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('mutation')}}:</span> (${{t('refBase')}})<span class="tooltip-seq">${{variant.ref_base}}</span> → (${{t('qryBase')}})<span class="tooltip-seq">${{variant.alt_base}}</span></div>
                                `;
                                showTooltip(content, e, tt);
                                
                                // 绘制连线：找到配对的 rect（上下轨道）
                                const pairIdx = (idx % 2 === 0) ? idx + 1 : idx - 1;
                                if (pairIdx >= 0 && pairIdx < orangeRects.length) {{
                                    drawSnpConnectLine(rect, orangeRects[pairIdx]);
                                }}
                            }}
                        }});
                        rect.addEventListener('mousemove', (e) => updateTooltipPosition(e, tt));
                        rect.addEventListener('mouseleave', () => {{
                            hideTooltip(tt);
                            removeSnpConnectLine();
                        }});
                    }});
                    
                    blueRects.forEach((rect, idx) => {{
                        rect.style.cursor = 'pointer';
                        rect.classList.add('variant-indel');
                        
                        rect.addEventListener('mouseenter', function(e) {{
                            const variant = indels[idx];
                            if (variant) {{
                                let content = '';
                                if (variant.type === 'INS') {{
                                    content = `
                                        <div class="tooltip-title">${{t('insertion')}}</div>
                                        <div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{variant.ref_pos}}</div>
                                        <div class="tooltip-row"><span class="tooltip-label">${{t('insertSeq')}}:</span> <span class="tooltip-seq">${{variant.seq}}</span></div>
                                        <div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{variant.length}} bp</div>
                                    `;
                                }} else {{
                                    content = `
                                        <div class="tooltip-title">${{t('deletion')}}</div>
                                        <div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{variant.qry_pos}}</div>
                                        <div class="tooltip-row"><span class="tooltip-label">${{t('deleteSeq')}}:</span> <span class="tooltip-seq">${{variant.seq}}</span></div>
                                        <div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{variant.length}} bp</div>
                                    `;
                                }}
                                showTooltip(content, e, tt);
                            }}
                        }});
                        rect.addEventListener('mousemove', (e) => updateTooltipPosition(e, tt));
                        rect.addEventListener('mouseleave', () => hideTooltip(tt));
                    }});
                }}
                
                // 图例点击切换（带勾选框状态）
                svg.querySelectorAll('.legend-item').forEach(item => {{
                    item.addEventListener('click', function() {{
                        const type = item.getAttribute('data-type');
                        visibility[type] = !visibility[type];
                        item.classList.toggle('disabled', !visibility[type]);
                        
                        // 切换勾选框状态
                        const checkmark = item.querySelector('.legend-checkmark');
                        if (checkmark) {{
                            checkmark.style.display = visibility[type] ? '' : 'none';
                        }}
                        
                        let elements = [];
                        switch(type) {{
                            case 'snp':
                                elements = svg.querySelectorAll('rect[fill="orange"]');
                                break;
                            case 'indel':
                                elements = svg.querySelectorAll('rect[fill="blue"]');
                                break;
                            case 'utr5':
                                elements = svg.querySelectorAll('.UTR5');
                                break;
                            case 'utr3':
                                elements = svg.querySelectorAll('.UTR3');
                                break;
                            case 'cds':
                                elements = svg.querySelectorAll('.exon');
                                break;
                        }}
                        
                        elements.forEach(el => {{
                            el.style.display = visibility[type] ? '' : 'none';
                        }});
                    }});
                }});
            }};
            
            // 初始化默认显示的 SVG
            const defaultContainer = document.querySelector('.svg-container[style*="display: block"]');
            if (defaultContainer) {{
                initSvgInteraction(defaultContainer.querySelector('svg'));
            }}
            
            // 放大按钮和模态框
            const zoomBtn = document.getElementById('zoom-btn');
            const zoomModal = document.getElementById('zoom-modal');
            const zoomModalClose = document.getElementById('zoom-modal-close');
            const zoomModalBody = document.getElementById('zoom-modal-body');
            
            if (zoomBtn && zoomModal) {{
                zoomBtn.addEventListener('click', function() {{
                    const activeContainer = document.querySelector('.svg-container[style*="display: block"]');
                    const svg = activeContainer ? activeContainer.querySelector('svg') : null;
                    if (svg) {{
                        const svgClone = svg.cloneNode(true);
                        svgClone.id = 'modal-svg';
                        svgClone.style.width = '100%';
                        svgClone.style.minWidth = '1200px';
                        zoomModalBody.innerHTML = '<div id="modal-tooltip" class="svg-tooltip"></div>';
                        zoomModalBody.appendChild(svgClone);
                        zoomModal.classList.add('visible');
                        document.body.style.overflow = 'hidden';
                        
                        // 为模态框中的 SVG 添加交互
                        const modalTooltip = document.getElementById('modal-tooltip');
                        initSvgInteraction(svgClone, modalTooltip);
                    }}
                }});
                
                zoomModalClose.addEventListener('click', function() {{
                    zoomModal.classList.remove('visible');
                    document.body.style.overflow = '';
                }});
                
                zoomModal.addEventListener('click', function(e) {{
                    if (e.target === zoomModal) {{
                        zoomModal.classList.remove('visible');
                        document.body.style.overflow = '';
                    }}
                }});
            }}
        }})();
        </script>
        """
        
        return html
    
    def _process_svg_for_embed(self, svg_content: str, gene_id: str, extraction_info: dict) -> str:
        """处理 SVG 内容用于嵌入（添加 viewBox、图例等）"""
        import re
        width_match = re.search(r'width="(\d+)"', svg_content)
        height_match = re.search(r'height="(\d+)"', svg_content)
        svg_width = 1200
        svg_height = 400
        if width_match and height_match:
            svg_width, svg_height = int(width_match.group(1)), int(height_match.group(1))
            # 给图注预留顶部空间，避免覆盖绘图区
            legend_top_margin = 90
            svg_content = re.sub(
                r'<svg\s+width="[^"]*"\s+height="[^"]*"',
                f'<svg viewBox="0 -{legend_top_margin} {svg_width} {svg_height + legend_top_margin}" width="100%" preserveAspectRatio="xMidYMid meet" style="overflow:visible"',
                svg_content
            )
            # 添加图例（放在顶部空白区）
            legend_svg = self._generate_legend_svg(svg_width, -legend_top_margin + 35)
            svg_content = svg_content.replace('</svg>', f'{legend_svg}</svg>')
            
            # 如果有 extraction_info 且有上下游延伸，注入基因区域大括号
            if extraction_info and (extraction_info.get("upstream", 0) > 0 or extraction_info.get("downstream", 0) > 0):
                seq_len = extraction_info.get("actual_end", 1) - extraction_info.get("actual_start", 1) + 1
                bracket_svg = self._generate_gene_bracket_svg(
                    svg_width, svg_height,
                    extraction_info.get("gene_rel_start", 1),
                    extraction_info.get("gene_rel_end", 1),
                    gene_id,
                    seq_len
                )
                svg_content = svg_content.replace('</svg>', f'{bracket_svg}</svg>')
        
        return svg_content
    
    def _generate_visualization_error_section(self, message):
        """生成可视化错误提示"""
        return f"""
        <div class="section visualization-section">
            <h2>比对可视化</h2>
            <div class="viz-container" style="padding: 40px; text-align: center; color: #999;">
                <p>{message}</p>
            </div>
        </div>
        """
    
    def _parse_snps_for_tooltip(self, snps_file):
        """
        解析 snps 文件，生成 tooltip 数据
        
        返回: list of dict，每个 dict 包含变异信息
        """
        variants = []
        if not snps_file or not os.path.exists(snps_file):
            return variants
        
        with open(snps_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('[') or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 5:
                    try:
                        ref_pos = int(parts[0])
                        ref_base = parts[1]
                        alt_base = parts[2]
                        qry_pos = int(parts[3])
                        var_type = parts[4].upper()
                        
                        variant = {
                            'ref_pos': ref_pos,
                            'qry_pos': qry_pos,
                            'ref_base': ref_base,
                            'alt_base': alt_base,
                            'type': var_type
                        }
                        
                        # 计算 Indel 长度
                        if var_type == 'INS':
                            variant['length'] = len(alt_base)
                            variant['seq'] = alt_base
                        elif var_type == 'DEL':
                            variant['length'] = len(ref_base)
                            variant['seq'] = ref_base
                        
                        variants.append(variant)
                    except (ValueError, IndexError):
                        continue
        return variants

    def _parse_gff_for_tooltip(self, gff_file, fasta_file=None):
        """
        解析 GFF 文件，生成 UTR/CDS 的 tooltip 数据
        
        返回: dict，包含 utr3, utr5, cds 列表
        """
        features = {
            'utr3': [],
            'utr5': [],
            'cds': []
        }
        
        if not gff_file or not os.path.exists(gff_file):
            return features
        
        # 读取序列用于提取特征序列
        sequence = ""
        if fasta_file and os.path.exists(fasta_file):
            with open(fasta_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if not line.startswith('>'):
                        sequence += line.strip()
        
        with open(gff_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 8:
                    try:
                        feat_type = parts[2].lower()
                        start = int(parts[3])
                        end = int(parts[4])
                        strand = parts[6] if len(parts) > 6 else '+'
                        length = end - start + 1
                        
                        # 提取序列片段（最多显示 20bp）
                        seq = ""
                        if sequence and start > 0 and end <= len(sequence):
                            seq = sequence[start-1:end]
                            if len(seq) > 20:
                                seq = seq[:10] + "..." + seq[-7:]
                        
                        feature = {
                            'start': start,
                            'end': end,
                            'length': length,
                            'strand': strand,
                            'seq': seq
                        }
                        
                        if 'three_prime_utr' in feat_type or '3_utr' in feat_type or feat_type == 'three_prime_utr':
                            features['utr3'].append(feature)
                        elif 'five_prime_utr' in feat_type or '5_utr' in feat_type or feat_type == 'five_prime_utr':
                            features['utr5'].append(feature)
                        elif feat_type == 'cds':
                            features['cds'].append(feature)
                    except (ValueError, IndexError):
                        continue
        
        return features

    def _embed_image_section(self, image_path, gene_id, snps_file=None, gff_file=None, qry_gff_file=None, fasta_file=None, extraction_info=None):
        """
        将图片嵌入 HTML
        SVG 直接嵌入代码（支持交互），PNG 用 Base64
        在 SVG 右上角添加图例，添加鼠标悬停 tooltip
        
        参数:
        - gff_file: 上方轨道的 GFF 文件（相对坐标）
        - qry_gff_file: 下方轨道的 GFF 文件（绝对坐标），可选
        - extraction_info: 提取信息（可选，包含 gene_rel_start, gene_rel_end 用于绘制大括号）
        """
        # 解析变异数据用于 tooltip
        variants = self._parse_snps_for_tooltip(snps_file) if snps_file else []
        
        # 解析 GFF 数据用于 UTR/CDS tooltip
        # 分别解析上方轨道（ref）和下方轨道（qry）的 GFF，并添加轨道标识
        ref_features = self._parse_gff_for_tooltip(gff_file, fasta_file) if gff_file else {'utr3': [], 'utr5': [], 'cds': []}
        qry_features = self._parse_gff_for_tooltip(qry_gff_file, None) if qry_gff_file else {'utr3': [], 'utr5': [], 'cds': []}
        
        # 为每个特征添加轨道标识
        for feat in ref_features['utr3'] + ref_features['utr5'] + ref_features['cds']:
            feat['track'] = 'top'
        for feat in qry_features['utr3'] + qry_features['utr5'] + qry_features['cds']:
            feat['track'] = 'bottom'
        
        # 分开存储，不合并（JavaScript 会根据元素位置匹配）
        features = {
            'top': ref_features,
            'bottom': qry_features
        }
        
        html = """
        <div class="section visualization-section">
            <h2>比对可视化</h2>
            <div class="viz-wrapper">
                <div class="viz-container" id="viz-container" style="width:100%; position: relative; overflow: hidden;">
                    <div id="svg-tooltip" class="svg-tooltip"></div>
        """
        
        try:
            if image_path.endswith('.svg'):
                # SVG 直接嵌入代码，支持后续添加交互
                with open(image_path, 'r', encoding='utf-8') as f:
                    svg_content = f.read()
                # 提取原始宽高，添加 viewBox 保持比例
                import re
                width_match = re.search(r'width="(\d+)"', svg_content)
                height_match = re.search(r'height="(\d+)"', svg_content)
                svg_width = 1200
                svg_height = 400
                if width_match and height_match:
                    svg_width, svg_height = int(width_match.group(1)), int(height_match.group(1))
                    # 给图注预留顶部空间，避免覆盖绘图区
                    legend_top_margin = 90
                    svg_content = re.sub(
                        r'<svg\s+width="[^"]*"\s+height="[^"]*"',
                        f'<svg id="alignment-svg" viewBox="0 -{legend_top_margin} {svg_width} {svg_height + legend_top_margin}" width="100%" preserveAspectRatio="xMidYMid meet" style="overflow:visible"',
                        svg_content
                    )
                    # 在 </svg> 前注入图例（放在顶部空白区）
                    legend_svg = self._generate_legend_svg(svg_width, -legend_top_margin + 35)
                    svg_content = svg_content.replace('</svg>', f'{legend_svg}</svg>')
                    
                    # 如果有 extraction_info 且有上下游延伸，注入基因区域大括号
                    if extraction_info and (extraction_info.get("upstream", 0) > 0 or extraction_info.get("downstream", 0) > 0):
                        # 计算序列总长度
                        seq_len = extraction_info.get("actual_end", 1) - extraction_info.get("actual_start", 1) + 1
                        bracket_svg = self._generate_gene_bracket_svg(
                            svg_width, svg_height,
                            extraction_info.get("gene_rel_start", 1),
                            extraction_info.get("gene_rel_end", 1),
                            gene_id,
                            seq_len
                        )
                        svg_content = svg_content.replace('</svg>', f'{bracket_svg}</svg>')
                
                html += svg_content
                
                # 添加放大按钮
                html += """
                    <button id="zoom-btn" class="zoom-btn" title="放大查看">
                        <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <circle cx="11" cy="11" r="8"></circle>
                            <path d="M21 21l-4.35-4.35"></path>
                            <path d="M11 8v6M8 11h6"></path>
                        </svg>
                    </button>
                </div>
                """
            else:
                # PNG 用 Base64 嵌入
                with open(image_path, 'rb') as f:
                    img_data = base64.b64encode(f.read()).decode()
                html += f'<img src="data:image/png;base64,{img_data}" width="100%">'
                html += "</div>"
        except Exception as e:
            print(f"[WARNING] 嵌入图片失败: {e}")
            html += f'<p>图片加载失败: {image_path}</p>'
            html += "</div>"
        
        html += """
            </div>
            <div class="viz-legend-note">
                <p data-zh="上方轨道为基因序列，下方轨道为查询基因组比对区域。橙色标记表示 SNP，蓝色标记表示 Indel（插入/缺失）。鼠标悬停可查看详情。" 
                   data-en="Upper track shows gene sequence, lower track shows aligned region in query genome. Orange marks indicate SNPs, blue marks indicate Indels. Hover for details.">
                   上方轨道为基因序列，下方轨道为查询基因组比对区域。橙色标记表示 SNP，蓝色标记表示 Indel（插入/缺失）。鼠标悬停可查看详情。
                </p>
            </div>
        </div>
        
        <!-- 模态框 -->
        <div id="zoom-modal" class="zoom-modal">
            <div class="zoom-modal-content">
                <button class="zoom-modal-close" id="zoom-modal-close">&times;</button>
                <div class="zoom-modal-body" id="zoom-modal-body"></div>
            </div>
        </div>
        """
        
        # 添加 tooltip 样式和交互脚本
        variants_json = json.dumps(variants, ensure_ascii=False)
        features_json = json.dumps(features, ensure_ascii=False)
        
        html += f"""
        <style>
            .viz-wrapper {{
                position: relative;
            }}
            .svg-tooltip {{
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
            }}
            .svg-tooltip.visible {{
                opacity: 1;
            }}
            .svg-tooltip .tooltip-title {{
                font-weight: bold;
                margin-bottom: 4px;
                color: #ffd700;
            }}
            .svg-tooltip .tooltip-row {{
                margin: 2px 0;
            }}
            .svg-tooltip .tooltip-label {{
                color: #aaa;
            }}
            .svg-tooltip .tooltip-seq {{
                font-family: monospace;
                background: rgba(255,255,255,0.1);
                padding: 1px 4px;
                border-radius: 2px;
            }}
            /* 放大按钮 */
            .zoom-btn {{
                position: absolute;
                bottom: 10px;
                right: 10px;
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
            }}
            .zoom-btn:hover {{
                background: rgba(102, 126, 234, 1);
                transform: scale(1.1);
            }}
            /* 模态框 */
            .zoom-modal {{
                display: none;
                position: fixed;
                top: 0;
                left: 0;
                width: 100%;
                height: 100%;
                background: rgba(0, 0, 0, 0.8);
                z-index: 10000;
                overflow: auto;
            }}
            .zoom-modal.visible {{
                display: flex;
                align-items: center;
                justify-content: center;
            }}
            .zoom-modal-content {{
                position: relative;
                background: white;
                border-radius: 8px;
                width: 95%;
                max-width: 1400px;
                max-height: 90vh;
                overflow: auto;
                padding: 20px;
            }}
            .zoom-modal-close {{
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
            }}
            .zoom-modal-close:hover {{
                color: #333;
            }}
            .zoom-modal-body {{
                width: 100%;
                overflow-x: auto;
                position: relative;
            }}
            .zoom-modal-body svg {{
                min-width: 1200px;
            }}
            .legend-item {{
                transition: opacity 0.2s;
            }}
            .legend-item:hover {{
                opacity: 0.8;
            }}
            .legend-item.disabled {{
                opacity: 0.4;
            }}
            .legend-item.disabled .legend-text {{
                text-decoration: line-through;
            }}
        </style>
        <script>
        (function() {{
            const variants = {variants_json};
            const features = {features_json};
            const tooltip = document.getElementById('svg-tooltip');
            const svg = document.getElementById('alignment-svg');
            
            if (!svg || !tooltip) return;
            
            // 显示/隐藏状态
            const visibility = {{
                snp: true,
                indel: true,
                utr5: true,
                utr3: true,
                cds: true
            }};
            
            // 国际化文本
            const i18n = {{
                zh: {{
                    position: '位置',
                    mutation: '变异',
                    refBase: '参考',
                    qryBase: '查询',
                    insertSeq: '插入序列',
                    deleteSeq: '缺失序列',
                    length: '长度',
                    sequence: '序列',
                    insertion: '插入 (INS)',
                    deletion: '缺失 (DEL)',
                    utr3: "3' UTR",
                    utr3Desc: '三端非翻译区',
                    utr5: "5' UTR",
                    utr5Desc: '五端非翻译区',
                    cds: 'CDS',
                    cdsDesc: '编码序列'
                }},
                en: {{
                    position: 'Position',
                    mutation: 'Mutation',
                    refBase: 'Ref',
                    qryBase: 'Qry',
                    insertSeq: 'Inserted Seq',
                    deleteSeq: 'Deleted Seq',
                    length: 'Length',
                    sequence: 'Sequence',
                    insertion: 'Insertion (INS)',
                    deletion: 'Deletion (DEL)',
                    utr3: "3' UTR",
                    utr3Desc: "3' Untranslated Region",
                    utr5: "5' UTR",
                    utr5Desc: "5' Untranslated Region",
                    cds: 'CDS',
                    cdsDesc: 'Coding Sequence'
                }}
            }};
            
            // 获取当前语言
            function getLang() {{
                return typeof currentLang !== 'undefined' ? currentLang : 'zh';
            }}
            
            function t(key) {{
                const lang = getLang();
                return i18n[lang] && i18n[lang][key] ? i18n[lang][key] : i18n['zh'][key];
            }}
            
            // 辅助函数：显示 tooltip
            function showTooltip(content, e) {{
                tooltip.innerHTML = content;
                tooltip.classList.add('visible');
                updateTooltipPosition(e);
            }}
            
            function updateTooltipPosition(e) {{
                const container = tooltip.parentElement;
                const containerRect = container.getBoundingClientRect();
                let x = e.clientX - containerRect.left + 15;
                let y = e.clientY - containerRect.top - 10;
                
                // 防止超出右边界
                if (x + 200 > containerRect.width) {{
                    x = e.clientX - containerRect.left - 200;
                }}
                
                tooltip.style.left = x + 'px';
                tooltip.style.top = y + 'px';
            }}
            
            function hideTooltip() {{
                tooltip.classList.remove('visible');
            }}
            
            // 扩大点击区域的辅助函数
            function addHitArea(element, padding) {{
                const bbox = element.getBBox ? element.getBBox() : null;
                if (!bbox) return;
                
                // 创建透明的扩大点击区域
                const hitArea = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
                hitArea.setAttribute('x', bbox.x - padding);
                hitArea.setAttribute('y', bbox.y - padding);
                hitArea.setAttribute('width', bbox.width + padding * 2);
                hitArea.setAttribute('height', bbox.height + padding * 2);
                hitArea.setAttribute('fill', 'transparent');
                hitArea.style.cursor = 'pointer';
                
                element.parentNode.insertBefore(hitArea, element);
                return hitArea;
            }}
            
            // ========== UTR 交互 ==========
            // 获取 SVG 的高度中点，用于判断元素属于上方还是下方轨道
            const svgRect = svg.getBBox ? svg.getBBox() : {{ height: 400 }};
            const midY = svgRect.height / 2;
            
            // 辅助函数：根据元素 Y 坐标判断属于哪个轨道
            function getTrack(el) {{
                const bbox = el.getBBox ? el.getBBox() : null;
                if (!bbox) return 'top';
                return bbox.y < midY ? 'top' : 'bottom';
            }}
            
            // 辅助函数：获取对应轨道的特征数据
            function getTrackFeatures(track, type) {{
                if (features[track] && features[track][type]) {{
                    return features[track][type];
                }}
                return [];
            }}
            
            const utr3Elements = svg.querySelectorAll('.UTR3');
            const utr5Elements = svg.querySelectorAll('.UTR5');
            
            // 按轨道分组统计索引
            let topUtr3Idx = 0, bottomUtr3Idx = 0;
            let topUtr5Idx = 0, bottomUtr5Idx = 0;
            let topCdsIdx = 0, bottomCdsIdx = 0;
            
            utr3Elements.forEach((el) => {{
                el.style.cursor = 'pointer';
                const track = getTrack(el);
                const trackFeatures = getTrackFeatures(track, 'utr3');
                const idx = track === 'top' ? topUtr3Idx++ : bottomUtr3Idx++;
                
                el.addEventListener('mouseenter', function(e) {{
                    const feat = trackFeatures[idx] || {{}};
                    let content = `<div class="tooltip-title">${{t('utr3')}}</div>`;
                    content += `<div class="tooltip-row">${{t('utr3Desc')}}</div>`;
                    if (feat.start) {{
                        content += `<div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{feat.start}} - ${{feat.end}}</div>`;
                        content += `<div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{feat.length}} bp</div>`;
                        if (feat.seq) {{
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('sequence')}}:</span> <span class="tooltip-seq">${{feat.seq}}</span></div>`;
                        }}
                    }}
                    showTooltip(content, e);
                }});
                el.addEventListener('mousemove', updateTooltipPosition);
                el.addEventListener('mouseleave', hideTooltip);
            }});
            
            utr5Elements.forEach((el) => {{
                el.style.cursor = 'pointer';
                const track = getTrack(el);
                const trackFeatures = getTrackFeatures(track, 'utr5');
                const idx = track === 'top' ? topUtr5Idx++ : bottomUtr5Idx++;
                
                el.addEventListener('mouseenter', function(e) {{
                    const feat = trackFeatures[idx] || {{}};
                    let content = `<div class="tooltip-title">${{t('utr5')}}</div>`;
                    content += `<div class="tooltip-row">${{t('utr5Desc')}}</div>`;
                    if (feat.start) {{
                        content += `<div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{feat.start}} - ${{feat.end}}</div>`;
                        content += `<div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{feat.length}} bp</div>`;
                        if (feat.seq) {{
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('sequence')}}:</span> <span class="tooltip-seq">${{feat.seq}}</span></div>`;
                        }}
                    }}
                    showTooltip(content, e);
                }});
                el.addEventListener('mousemove', updateTooltipPosition);
                el.addEventListener('mouseleave', hideTooltip);
            }});
            
            // ========== CDS/exon 交互 ==========
            const exonElements = svg.querySelectorAll('.exon');
            
            exonElements.forEach((el) => {{
                el.style.cursor = 'pointer';
                const track = getTrack(el);
                const trackFeatures = getTrackFeatures(track, 'cds');
                const idx = track === 'top' ? topCdsIdx++ : bottomCdsIdx++;
                
                el.addEventListener('mouseenter', function(e) {{
                    const feat = trackFeatures[idx] || {{}};
                    let content = `<div class="tooltip-title">${{t('cds')}}</div>`;
                    content += `<div class="tooltip-row">${{t('cdsDesc')}}</div>`;
                    if (feat.start) {{
                        content += `<div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{feat.start}} - ${{feat.end}}</div>`;
                        content += `<div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{feat.length}} bp</div>`;
                        if (feat.seq) {{
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('sequence')}}:</span> <span class="tooltip-seq">${{feat.seq}}</span></div>`;
                        }}
                    }}
                    showTooltip(content, e);
                }});
                el.addEventListener('mousemove', updateTooltipPosition);
                el.addEventListener('mouseleave', hideTooltip);
            }});
            
            // ========== SNP/Indel 交互 ==========
            if (variants.length > 0) {{
                const orangeRects = Array.from(svg.querySelectorAll('rect[fill="orange"]'));
                const blueRects = Array.from(svg.querySelectorAll('rect[fill="blue"]'));
                const snps = variants.filter(v => v.type === 'SNP');
                const indels = variants.filter(v => v.type !== 'SNP');
                
                // 创建连线元素（用于 SNP 悬停时显示）
                let snpConnectLine = null;
                
                // 绘制 SNP 连线的函数
                function drawSnpConnectLine(rect1, rect2) {{
                    if (snpConnectLine) {{
                        snpConnectLine.remove();
                    }}
                    
                    // 获取两个 rect 的中心坐标
                    const x1 = parseFloat(rect1.getAttribute('x')) + parseFloat(rect1.getAttribute('width')) / 2;
                    const y1 = parseFloat(rect1.getAttribute('y')) + parseFloat(rect1.getAttribute('height')) / 2;
                    const x2 = parseFloat(rect2.getAttribute('x')) + parseFloat(rect2.getAttribute('width')) / 2;
                    const y2 = parseFloat(rect2.getAttribute('y')) + parseFloat(rect2.getAttribute('height')) / 2;
                    
                    // 创建连线
                    snpConnectLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                    snpConnectLine.setAttribute('x1', x1);
                    snpConnectLine.setAttribute('y1', y1);
                    snpConnectLine.setAttribute('x2', x2);
                    snpConnectLine.setAttribute('y2', y2);
                    snpConnectLine.setAttribute('stroke', 'orange');
                    snpConnectLine.setAttribute('stroke-width', '2');
                    snpConnectLine.setAttribute('stroke-dasharray', '4,2');
                    snpConnectLine.setAttribute('pointer-events', 'none');
                    snpConnectLine.classList.add('snp-connect-line');
                    svg.appendChild(snpConnectLine);
                }}
                
                // 移除连线的函数
                function removeSnpConnectLine() {{
                    if (snpConnectLine) {{
                        snpConnectLine.remove();
                        snpConnectLine = null;
                    }}
                }}
                
                // SNP 交互
                orangeRects.forEach((rect, idx) => {{
                    rect.style.cursor = 'pointer';
                    rect.classList.add('variant-snp');
                    
                    rect.addEventListener('mouseenter', function(e) {{
                        const variant = snps[Math.floor(idx / 2)]; // 每个 SNP 有两个标记
                        if (variant) {{
                            const content = `
                                <div class="tooltip-title">SNP</div>
                                <div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{variant.qry_pos}}</div>
                                <div class="tooltip-row"><span class="tooltip-label">${{t('mutation')}}:</span> (${{t('refBase')}})<span class="tooltip-seq">${{variant.ref_base}}</span> → (${{t('qryBase')}})<span class="tooltip-seq">${{variant.alt_base}}</span></div>
                            `;
                            showTooltip(content, e);
                            
                            // 绘制连线：找到配对的 rect（上下轨道）
                            const pairIdx = (idx % 2 === 0) ? idx + 1 : idx - 1;
                            if (pairIdx >= 0 && pairIdx < orangeRects.length) {{
                                drawSnpConnectLine(rect, orangeRects[pairIdx]);
                            }}
                        }}
                    }});
                    rect.addEventListener('mousemove', updateTooltipPosition);
                    rect.addEventListener('mouseleave', () => {{
                        hideTooltip();
                        removeSnpConnectLine();
                    }});
                }});
                
                // Indel 交互
                blueRects.forEach((rect, idx) => {{
                    rect.style.cursor = 'pointer';
                    rect.classList.add('variant-indel');
                    
                    rect.addEventListener('mouseenter', function(e) {{
                        const variant = indels[idx];
                        if (variant) {{
                            let content = '';
                            if (variant.type === 'INS') {{
                                content = `
                                    <div class="tooltip-title">${{t('insertion')}}</div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{variant.ref_pos}}</div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('insertSeq')}}:</span> <span class="tooltip-seq">${{variant.seq}}</span></div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{variant.length}} bp</div>
                                `;
                            }} else {{
                                content = `
                                    <div class="tooltip-title">${{t('deletion')}}</div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{variant.qry_pos}}</div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('deleteSeq')}}:</span> <span class="tooltip-seq">${{variant.seq}}</span></div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{variant.length}} bp</div>
                                `;
                            }}
                            showTooltip(content, e);
                        }}
                    }});
                    rect.addEventListener('mousemove', updateTooltipPosition);
                    rect.addEventListener('mouseleave', hideTooltip);
                }});
            }}
            
            // ========== 图例点击切换（带勾选框状态） ==========
            const legendItems = svg.querySelectorAll('.legend-item');
            
            legendItems.forEach(item => {{
                item.addEventListener('click', function() {{
                    const type = item.getAttribute('data-type');
                    visibility[type] = !visibility[type];
                    
                    // 切换图例样式
                    item.classList.toggle('disabled', !visibility[type]);
                    
                    // 切换勾选框状态
                    const checkmark = item.querySelector('.legend-checkmark');
                    if (checkmark) {{
                        checkmark.style.display = visibility[type] ? '' : 'none';
                    }}
                    
                    // 切换对应元素的显示/隐藏
                    let elements = [];
                    switch(type) {{
                        case 'snp':
                            elements = svg.querySelectorAll('rect[fill="orange"]');
                            break;
                        case 'indel':
                            elements = svg.querySelectorAll('rect[fill="blue"]');
                            break;
                        case 'utr5':
                            elements = svg.querySelectorAll('.UTR5');
                            break;
                        case 'utr3':
                            elements = svg.querySelectorAll('.UTR3');
                            break;
                        case 'cds':
                            elements = svg.querySelectorAll('.exon');
                            break;
                    }}
                    
                    elements.forEach(el => {{
                        el.style.display = visibility[type] ? '' : 'none';
                    }});
                }});
            }});
            
            // ========== 放大按钮和模态框 ==========
            const zoomBtn = document.getElementById('zoom-btn');
            const zoomModal = document.getElementById('zoom-modal');
            const zoomModalClose = document.getElementById('zoom-modal-close');
            const zoomModalBody = document.getElementById('zoom-modal-body');
            
            if (zoomBtn && zoomModal && svg) {{
                zoomBtn.addEventListener('click', function() {{
                    // 克隆 SVG 到模态框
                    const svgClone = svg.cloneNode(true);
                    svgClone.id = 'modal-svg';
                    svgClone.style.width = '100%';
                    svgClone.style.minWidth = '1200px';
                    zoomModalBody.innerHTML = '<div id="modal-tooltip" class="svg-tooltip"></div>';
                    zoomModalBody.appendChild(svgClone);
                    zoomModal.classList.add('visible');
                    document.body.style.overflow = 'hidden';
                    
                    // 为模态框中的 SVG 添加交互
                    setupModalInteraction(svgClone);
                }});
                
                zoomModalClose.addEventListener('click', function() {{
                    zoomModal.classList.remove('visible');
                    document.body.style.overflow = '';
                }});
                
                zoomModal.addEventListener('click', function(e) {{
                    if (e.target === zoomModal) {{
                        zoomModal.classList.remove('visible');
                        document.body.style.overflow = '';
                    }}
                }});
            }}
            
            // 为模态框 SVG 设置交互
            function setupModalInteraction(modalSvg) {{
                const modalTooltip = document.getElementById('modal-tooltip');
                const modalVisibility = {{ snp: true, indel: true, utr5: true, utr3: true, cds: true }};
                
                function showModalTooltip(content, e) {{
                    modalTooltip.innerHTML = content;
                    modalTooltip.classList.add('visible');
                    updateModalTooltipPosition(e);
                }}
                
                function updateModalTooltipPosition(e) {{
                    const container = zoomModalBody;
                    const containerRect = container.getBoundingClientRect();
                    let x = e.clientX - containerRect.left + 15;
                    let y = e.clientY - containerRect.top - 10;
                    if (x + 200 > containerRect.width) x = e.clientX - containerRect.left - 200;
                    modalTooltip.style.left = x + 'px';
                    modalTooltip.style.top = y + 'px';
                }}
                
                function hideModalTooltip() {{
                    modalTooltip.classList.remove('visible');
                }}
                
                // 获取模态框 SVG 的高度中点
                const modalSvgRect = modalSvg.getBBox ? modalSvg.getBBox() : {{ height: 400 }};
                const modalMidY = modalSvgRect.height / 2;
                
                // 辅助函数：根据元素 Y 坐标判断属于哪个轨道
                function getModalTrack(el) {{
                    const bbox = el.getBBox ? el.getBBox() : null;
                    if (!bbox) return 'top';
                    return bbox.y < modalMidY ? 'top' : 'bottom';
                }}
                
                // 按轨道分组统计索引
                let modalTopUtr3Idx = 0, modalBottomUtr3Idx = 0;
                let modalTopUtr5Idx = 0, modalBottomUtr5Idx = 0;
                let modalTopCdsIdx = 0, modalBottomCdsIdx = 0;
                
                // UTR 交互
                modalSvg.querySelectorAll('.UTR3').forEach((el) => {{
                    el.style.cursor = 'pointer';
                    const track = getModalTrack(el);
                    const trackFeatures = features[track] ? features[track].utr3 || [] : [];
                    const idx = track === 'top' ? modalTopUtr3Idx++ : modalBottomUtr3Idx++;
                    
                    el.addEventListener('mouseenter', function(e) {{
                        const feat = trackFeatures[idx] || {{}};
                        let content = `<div class="tooltip-title">${{t('utr3')}}</div><div class="tooltip-row">${{t('utr3Desc')}}</div>`;
                        if (feat.start) {{
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{feat.start}} - ${{feat.end}}</div>`;
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{feat.length}} bp</div>`;
                            if (feat.seq) content += `<div class="tooltip-row"><span class="tooltip-label">${{t('sequence')}}:</span> <span class="tooltip-seq">${{feat.seq}}</span></div>`;
                        }}
                        showModalTooltip(content, e);
                    }});
                    el.addEventListener('mousemove', updateModalTooltipPosition);
                    el.addEventListener('mouseleave', hideModalTooltip);
                }});
                
                modalSvg.querySelectorAll('.UTR5').forEach((el) => {{
                    el.style.cursor = 'pointer';
                    const track = getModalTrack(el);
                    const trackFeatures = features[track] ? features[track].utr5 || [] : [];
                    const idx = track === 'top' ? modalTopUtr5Idx++ : modalBottomUtr5Idx++;
                    
                    el.addEventListener('mouseenter', function(e) {{
                        const feat = trackFeatures[idx] || {{}};
                        let content = `<div class="tooltip-title">${{t('utr5')}}</div><div class="tooltip-row">${{t('utr5Desc')}}</div>`;
                        if (feat.start) {{
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{feat.start}} - ${{feat.end}}</div>`;
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{feat.length}} bp</div>`;
                            if (feat.seq) content += `<div class="tooltip-row"><span class="tooltip-label">${{t('sequence')}}:</span> <span class="tooltip-seq">${{feat.seq}}</span></div>`;
                        }}
                        showModalTooltip(content, e);
                    }});
                    el.addEventListener('mousemove', updateModalTooltipPosition);
                    el.addEventListener('mouseleave', hideModalTooltip);
                }});
                
                // CDS 交互
                modalSvg.querySelectorAll('.exon').forEach((el) => {{
                    el.style.cursor = 'pointer';
                    const track = getModalTrack(el);
                    const trackFeatures = features[track] ? features[track].cds || [] : [];
                    const idx = track === 'top' ? modalTopCdsIdx++ : modalBottomCdsIdx++;
                    
                    el.addEventListener('mouseenter', function(e) {{
                        const feat = trackFeatures[idx] || {{}};
                        let content = `<div class="tooltip-title">${{t('cds')}}</div><div class="tooltip-row">${{t('cdsDesc')}}</div>`;
                        if (feat.start) {{
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{feat.start}} - ${{feat.end}}</div>`;
                            content += `<div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{feat.length}} bp</div>`;
                            if (feat.seq) content += `<div class="tooltip-row"><span class="tooltip-label">${{t('sequence')}}:</span> <span class="tooltip-seq">${{feat.seq}}</span></div>`;
                        }}
                        showModalTooltip(content, e);
                    }});
                    el.addEventListener('mousemove', updateModalTooltipPosition);
                    el.addEventListener('mouseleave', hideModalTooltip);
                }});
                
                // SNP 交互（带连线）
                const modalOrangeRects = Array.from(modalSvg.querySelectorAll('rect[fill="orange"]'));
                const snps = variants.filter(v => v.type === 'SNP');
                
                // 创建连线元素（用于 SNP 悬停时显示）
                let modalSnpConnectLine = null;
                
                // 绘制 SNP 连线的函数
                function drawModalSnpConnectLine(rect1, rect2) {{
                    if (modalSnpConnectLine) {{
                        modalSnpConnectLine.remove();
                    }}
                    
                    const x1 = parseFloat(rect1.getAttribute('x')) + parseFloat(rect1.getAttribute('width')) / 2;
                    const y1 = parseFloat(rect1.getAttribute('y')) + parseFloat(rect1.getAttribute('height')) / 2;
                    const x2 = parseFloat(rect2.getAttribute('x')) + parseFloat(rect2.getAttribute('width')) / 2;
                    const y2 = parseFloat(rect2.getAttribute('y')) + parseFloat(rect2.getAttribute('height')) / 2;
                    
                    modalSnpConnectLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                    modalSnpConnectLine.setAttribute('x1', x1);
                    modalSnpConnectLine.setAttribute('y1', y1);
                    modalSnpConnectLine.setAttribute('x2', x2);
                    modalSnpConnectLine.setAttribute('y2', y2);
                    modalSnpConnectLine.setAttribute('stroke', 'orange');
                    modalSnpConnectLine.setAttribute('stroke-width', '2');
                    modalSnpConnectLine.setAttribute('stroke-dasharray', '4,2');
                    modalSnpConnectLine.setAttribute('pointer-events', 'none');
                    modalSnpConnectLine.classList.add('snp-connect-line');
                    modalSvg.appendChild(modalSnpConnectLine);
                }}
                
                function removeModalSnpConnectLine() {{
                    if (modalSnpConnectLine) {{
                        modalSnpConnectLine.remove();
                        modalSnpConnectLine = null;
                    }}
                }}
                
                modalOrangeRects.forEach((rect, idx) => {{
                    rect.style.cursor = 'pointer';
                    rect.addEventListener('mouseenter', function(e) {{
                        const variant = snps[Math.floor(idx / 2)];
                        if (variant) {{
                            const content = `
                                <div class="tooltip-title">SNP</div>
                                <div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{variant.qry_pos}}</div>
                                <div class="tooltip-row"><span class="tooltip-label">${{t('mutation')}}:</span> (${{t('refBase')}})<span class="tooltip-seq">${{variant.ref_base}}</span> → (${{t('qryBase')}})<span class="tooltip-seq">${{variant.alt_base}}</span></div>
                            `;
                            showModalTooltip(content, e);
                            
                            // 绘制连线
                            const pairIdx = (idx % 2 === 0) ? idx + 1 : idx - 1;
                            if (pairIdx >= 0 && pairIdx < modalOrangeRects.length) {{
                                drawModalSnpConnectLine(rect, modalOrangeRects[pairIdx]);
                            }}
                        }}
                    }});
                    rect.addEventListener('mousemove', updateModalTooltipPosition);
                    rect.addEventListener('mouseleave', () => {{
                        hideModalTooltip();
                        removeModalSnpConnectLine();
                    }});
                }});
                
                // Indel 交互
                const modalBlueRects = modalSvg.querySelectorAll('rect[fill="blue"]');
                const indels = variants.filter(v => v.type !== 'SNP');
                modalBlueRects.forEach((rect, idx) => {{
                    rect.style.cursor = 'pointer';
                    rect.addEventListener('mouseenter', function(e) {{
                        const variant = indels[idx];
                        if (variant) {{
                            let content = '';
                            if (variant.type === 'INS') {{
                                content = `
                                    <div class="tooltip-title">${{t('insertion')}}</div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{variant.ref_pos}}</div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('insertSeq')}}:</span> <span class="tooltip-seq">${{variant.seq}}</span></div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{variant.length}} bp</div>
                                `;
                            }} else {{
                                content = `
                                    <div class="tooltip-title">${{t('deletion')}}</div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('position')}}:</span> ${{variant.qry_pos}}</div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('deleteSeq')}}:</span> <span class="tooltip-seq">${{variant.seq}}</span></div>
                                    <div class="tooltip-row"><span class="tooltip-label">${{t('length')}}:</span> ${{variant.length}} bp</div>
                                `;
                            }}
                            showModalTooltip(content, e);
                        }}
                    }});
                    rect.addEventListener('mousemove', updateModalTooltipPosition);
                    rect.addEventListener('mouseleave', hideModalTooltip);
                }});
                
                // 图例点击切换（带勾选框状态）
                modalSvg.querySelectorAll('.legend-item').forEach(item => {{
                    item.addEventListener('click', function() {{
                        const type = item.getAttribute('data-type');
                        modalVisibility[type] = !modalVisibility[type];
                        item.classList.toggle('disabled', !modalVisibility[type]);
                        
                        // 切换勾选框状态
                        const checkmark = item.querySelector('.legend-checkmark');
                        if (checkmark) {{
                            checkmark.style.display = modalVisibility[type] ? '' : 'none';
                        }}
                        
                        let elements = [];
                        switch(type) {{
                            case 'snp': elements = modalSvg.querySelectorAll('rect[fill="orange"]'); break;
                            case 'indel': elements = modalSvg.querySelectorAll('rect[fill="blue"]'); break;
                            case 'utr5': elements = modalSvg.querySelectorAll('.UTR5'); break;
                            case 'utr3': elements = modalSvg.querySelectorAll('.UTR3'); break;
                            case 'cds': elements = modalSvg.querySelectorAll('.exon'); break;
                        }}
                        elements.forEach(el => {{ el.style.display = modalVisibility[type] ? '' : 'none'; }});
                    }});
                }});
            }}
        }})();
        </script>
        """
        
        return html

    def _generate_mode_section(self, mode_type, ref_genome, qry_genome, params, genome_files=None):
        """
        生成模块0：初始化信息（两行两列布局）
        """
        html = """
        <div class="section mode-section">
            <h2>初始化</h2>
            <table class="info-table-grid">
        """
        
        # 收集所有要显示的项目
        items = []
        if ref_genome:
            items.append(('参考基因组', ref_genome))
        if qry_genome:
            items.append(('查询基因组', qry_genome))
        for key, value in params.items():
            items.append((key, value))
        
        # 两列布局，每行两个元素
        for i in range(0, len(items), 2):
            html += '<tr>'
            # 第一个元素
            html += f'<td class="label">{items[i][0]}</td><td class="value">{items[i][1]}</td>'
            # 第二个元素（如果存在）
            if i + 1 < len(items):
                html += f'<td class="label">{items[i+1][0]}</td><td class="value">{items[i+1][1]}</td>'
            else:
                html += '<td class="label"></td><td class="value"></td>'
            html += '</tr>'
        
        html += """
            </table>
        </div>
        """
        return html

    def _generate_stats_section(self, fasta_file, coords_file, snps_file):
        """
        生成结果统计模块（支持国际化）
        
        参数:
        - fasta_file: 基因序列文件路径
        - coords_file: 比对坐标文件路径
        - snps_file: SNP/Indel 文件路径
        
        返回: HTML 字符串
        """
        # 计算基因长度
        gene_len = 0
        if fasta_file and os.path.exists(fasta_file):
            with open(fasta_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if not line.startswith('>'):
                        gene_len += len(line.strip())
        
        # 计算比对长度（从 coords 文件）
        align_len = 0
        if coords_file and os.path.exists(coords_file):
            with open(coords_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('[') or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 6:
                        try:
                            # 使用查询序列的比对长度
                            qry_aln_len = int(parts[5])
                            align_len += qry_aln_len
                        except (ValueError, IndexError):
                            continue
        
        # 统计 SNP 和 Indel 数量
        snp_count = 0
        indel_count = 0
        if snps_file and os.path.exists(snps_file):
            with open(snps_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('[') or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        var_type = parts[4].upper() if len(parts) > 4 else 'SNP'
                        if var_type == 'SNP':
                            snp_count += 1
                        else:  # INS, DEL
                            indel_count += 1
        
        # 生成 HTML（支持国际化）
        html = """
        <div class="section stats-section">
            <h2 data-zh="结果统计" data-en="Statistics">结果统计</h2>
            <table class="stats-table">
                <tr>
                    <th data-zh="基因长度" data-en="Gene Length">基因长度</th>
                    <th data-zh="比对长度" data-en="Alignment Length">比对长度</th>
                    <th data-zh="SNP 数量" data-en="SNP Count">SNP 数量</th>
                    <th data-zh="Indel 数量" data-en="Indel Count">Indel 数量</th>
                </tr>
                <tr>
        """
        html += f'<td>{gene_len:,} bp</td>'
        html += f'<td>{align_len:,} bp</td>'
        html += f'<td>{snp_count}</td>'
        html += f'<td>{indel_count}</td>'
        html += """
                </tr>
            </table>
        </div>
        """
        return html

    def _extract_region_gff(self, ref_gff, region_id, chr_name, start, end, output_dir):
        """
        从参考基因组 GFF 中提取指定区域的注释，并转换为相对坐标
        
        参数:
        - ref_gff: 参考基因组 GFF 文件路径
        - region_id: 区域 ID（用于输出文件命名）
        - chr_name: 染色体名称
        - start: 区域起始位置（1-based）
        - end: 区域终止位置（1-based）
        - output_dir: 输出目录
        
        返回:
        - 生成的相对坐标 GFF 文件路径，如果没有注释则返回 None
        """
        if not ref_gff or not os.path.exists(ref_gff):
            return None
        
        output_gff = os.path.join(output_dir, f"{region_id}.linkview.gff3")
        features_found = False
        
        # 只保留基因结构相关的特征类型（LINKVIEW 需要有 Parent 属性的特征）
        valid_types = {'gene', 'mrna', 'transcript', 'exon', 'cds', 
                      'five_prime_utr', 'three_prime_utr', 'utr', 
                      '5utr', '3utr', 'five_prime_UTR', 'three_prime_UTR'}
        
        try:
            with open(ref_gff, 'r', encoding='utf-8') as f_in, \
                 open(output_gff, 'w', encoding='utf-8') as f_out:
                for line in f_in:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue
                    
                    feat_type = parts[2].lower()
                    
                    # 过滤掉不需要的特征类型（如 chromosome, region 等）
                    if feat_type not in valid_types and parts[2] not in valid_types:
                        continue
                    
                    feat_chr = parts[0]
                    feat_start = int(parts[3])
                    feat_end = int(parts[4])
                    
                    # 检查染色体名称是否匹配（支持多种格式）
                    chr_match = False
                    if feat_chr == chr_name:
                        chr_match = True
                    elif feat_chr.replace('Chr', '').replace('chr', '') == chr_name.replace('Chr', '').replace('chr', ''):
                        chr_match = True
                    
                    if not chr_match:
                        continue
                    
                    # 检查是否在区域范围内（有重叠即可）
                    if feat_end < start or feat_start > end:
                        continue
                    
                    # 转换为相对坐标
                    rel_start = max(1, feat_start - start + 1)
                    rel_end = feat_end - start + 1
                    
                    # 写入相对坐标的 GFF
                    f_out.write(f"{region_id}\t{parts[1]}\t{parts[2]}\t{rel_start}\t{rel_end}\t{parts[5]}\t{parts[6]}\t{parts[7]}\t{parts[8]}\n")
                    features_found = True
            
            if features_found:
                print(f"[INFO] 已从 GFF 提取区域注释: {output_gff}")
                return output_gff
            else:
                # 没有找到注释，删除空文件
                os.remove(output_gff)
                print(f"[INFO] 区域 {chr_name}:{start}-{end} 没有找到基因注释")
                return None
                
        except Exception as e:
            print(f"[WARNING] 提取区域 GFF 失败: {e}")
            return None

    def _extract_aligned_region_gff(self, gff_file, coords_file, region_id, output_dir, track_name="qry"):
        """
        从 coords 文件中提取比对区域，然后从 GFF 中提取该区域的注释
        
        用于：
        - GeneID 模式：从 qry_gff 提取查询基因组比对区域的注释
        - Location 模式：从 qry_gff 提取查询基因组比对区域的注释
        - Sequence 模式：从 ref_gff 提取参考基因组比对区域的注释
        
        参数:
        - gff_file: GFF 文件路径
        - coords_file: 比对坐标文件路径
        - region_id: 区域 ID（用于输出文件命名）
        - output_dir: 输出目录
        - track_name: 轨道名称（用于 GFF 输出的 seqid）
        
        返回:
        - 生成的相对坐标 GFF 文件路径，如果没有注释则返回 None
        """
        if not gff_file or not os.path.exists(gff_file):
            return None
        if not coords_file or not os.path.exists(coords_file):
            return None
        
        # 从 coords 文件中提取主比对区域（总比对长度最长的染色体）
        chr_alignments = {}  # chr -> [(start, end, aln_len), ...]
        
        try:
            with open(coords_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('[') or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 8:
                        try:
                            ref_start = int(parts[0])
                            ref_end = int(parts[1])
                            ref_aln_len = int(parts[4]) if len(parts) > 4 else abs(ref_end - ref_start) + 1
                            qry_aln_len = int(parts[5]) if len(parts) > 5 else ref_aln_len
                            
                            # 过滤短比对：使用 max(LEN1, LEN2) 与 self.min_aln_len 比较
                            if max(ref_aln_len, qry_aln_len) < self.min_aln_len:
                                continue
                            
                            # 提取染色体名称（ref_name 在第8列）
                            ref_name = parts[7].split()[0] if parts[7] else "chr"
                            
                            if ref_name not in chr_alignments:
                                chr_alignments[ref_name] = []
                            chr_alignments[ref_name].append((min(ref_start, ref_end), max(ref_start, ref_end), ref_aln_len))
                        except (ValueError, IndexError):
                            continue
            
            # 找到总比对长度最长的染色体
            main_chr = None
            max_total_len = 0
            for chr_name, alns in chr_alignments.items():
                total_len = sum(a[2] for a in alns)
                if total_len > max_total_len:
                    max_total_len = total_len
                    main_chr = chr_name
            
            if not main_chr:
                print(f"[INFO] 未找到有效的比对区域")
                return None
            
            # 计算比对区域范围
            alns = chr_alignments[main_chr]
            region_start = min(a[0] for a in alns)
            region_end = max(a[1] for a in alns)
            
            print(f"[INFO] 比对区域: {main_chr}:{region_start}-{region_end}")
            
            # 从 GFF 中提取该区域的注释
            output_gff = os.path.join(output_dir, f"{region_id}.{track_name}.linkview.gff3")
            features_found = False
            
            # 只保留基因结构相关的特征类型（LINKVIEW 需要有 Parent 属性的特征）
            valid_types = {'gene', 'mrna', 'transcript', 'exon', 'cds', 
                          'five_prime_utr', 'three_prime_utr', 'utr', 
                          '5utr', '3utr', 'five_prime_UTR', 'three_prime_UTR'}
            
            with open(gff_file, 'r', encoding='utf-8') as f_in, \
                 open(output_gff, 'w', encoding='utf-8') as f_out:
                for line in f_in:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue
                    
                    feat_type = parts[2].lower()
                    
                    # 过滤掉不需要的特征类型（如 chromosome, region 等）
                    if feat_type not in valid_types and parts[2] not in valid_types:
                        continue
                    
                    feat_chr = parts[0]
                    feat_start = int(parts[3])
                    feat_end = int(parts[4])
                    
                    # 检查染色体名称是否匹配（支持多种格式）
                    chr_match = False
                    # 标准化染色体名称进行比较
                    feat_chr_norm = feat_chr.replace('Chr', '').replace('chr', '').lstrip('0')
                    main_chr_norm = main_chr.replace('Chr', '').replace('chr', '').lstrip('0')
                    # 处理类似 "1 dna:chromosome..." 的格式
                    if ' ' in main_chr_norm:
                        main_chr_norm = main_chr_norm.split()[0]
                    
                    if feat_chr == main_chr:
                        chr_match = True
                    elif feat_chr_norm == main_chr_norm:
                        chr_match = True
                    
                    if not chr_match:
                        continue
                    
                    # 检查是否在区域范围内（有重叠即可）
                    if feat_end < region_start or feat_start > region_end:
                        continue
                    
                    # 保持绝对坐标（不转换为相对坐标）
                    # LINKVIEW 的 get_gene_structure 函数使用 .k 文件中的绝对坐标范围
                    # 来检查 GFF 中的基因坐标是否有重叠，所以 GFF 必须使用绝对坐标
                    # seqid 使用 main_chr，与 .k 文件中的 name 部分匹配
                    f_out.write(f"{main_chr}\t{parts[1]}\t{parts[2]}\t{feat_start}\t{feat_end}\t{parts[5]}\t{parts[6]}\t{parts[7]}\t{parts[8]}\n")
                    features_found = True
            
            if features_found:
                print(f"[INFO] 已从 GFF 提取比对区域注释: {output_gff}")
                return output_gff
            else:
                # 没有找到注释，删除空文件
                os.remove(output_gff)
                print(f"[INFO] 比对区域 {main_chr}:{region_start}-{region_end} 没有找到基因注释")
                return None
                
        except Exception as e:
            print(f"[WARNING] 提取比对区域 GFF 失败: {e}")
            return None

    def _merge_gff_files(self, gff1, gff2, region_id):
        """
        合并两个 GFF 文件
        
        参数:
        - gff1: 第一个 GFF 文件路径（通常是参考基因组/上方轨道的注释）
        - gff2: 第二个 GFF 文件路径（通常是查询基因组/下方轨道的注释）
        - region_id: 区域 ID（用于输出文件命名）
        
        返回:
        - 合并后的 GFF 文件路径，如果都没有则返回 None
        """
        # 检查文件是否存在（处理 None 和空字符串）
        gff1_exists = gff1 and os.path.exists(gff1)
        gff2_exists = gff2 and os.path.exists(gff2)
        
        if not gff1_exists and not gff2_exists:
            return None
        if not gff1_exists:
            return gff2
        if not gff2_exists:
            return gff1
        
        output_gff = os.path.join(self.output_dir, f"{region_id}.merged.linkview.gff3")
        
        try:
            with open(output_gff, 'w', encoding='utf-8') as f_out:
                # 写入第一个 GFF 的内容
                with open(gff1, 'r', encoding='utf-8') as f_in:
                    for line in f_in:
                        if not line.startswith('#'):
                            f_out.write(line)
                
                # 写入第二个 GFF 的内容
                with open(gff2, 'r', encoding='utf-8') as f_in:
                    for line in f_in:
                        if not line.startswith('#'):
                            f_out.write(line)
            
            print(f"[INFO] 已合并 GFF 文件: {output_gff}")
            return output_gff
        except Exception as e:
            print(f"[WARNING] 合并 GFF 文件失败: {e}")
            return gff1  # 返回第一个文件作为降级

    def _parse_location(self, location):
        """
        解析位置字符串
        
        支持格式:
        - Chr01:1000-5000
        - Chr01:1000:5000
        - 1:1000-5000
        
        返回: (chr_name, start, end) 或 (None, None, None)
        """
        import re
        # 尝试匹配 chr:start-end 或 chr:start:end
        match = re.match(r'^([^:]+):(\d+)[-:](\d+)$', location)
        if match:
            return match.group(1), int(match.group(2)), int(match.group(3))
        return None, None, None

    def _generate_header(self, mode_type, mode_class):
        """
        生成报告标题栏
        
        参数:
        - mode_type: 模式类型显示名称
        - mode_class: CSS 类名
        """
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        # 模式类型的英文映射
        mode_en_map = {
            "Gene ID 模式": "Gene ID Mode",
            "Location 模式": "Location Mode",
            "Sequence 模式": "Sequence Mode"
        }
        mode_type_en = mode_en_map.get(mode_type, mode_type)
        
        return f"""
        <div class="header">
            <div class="header-content">
                <div class="header-left">
                    <h1 data-zh="GeneScreen 比对分析报告" data-en="GeneScreen Alignment Analysis Report">GeneScreen 比对分析报告</h1>
                    <div class="subtitle">
                        <span class="mode-badge {mode_class}" data-zh="{mode_type}" data-en="{mode_type_en}">{mode_type}</span>
                    </div>
                    <div class="gen-time" data-zh="生成时间: {now}" data-en="Generated: {now}">生成时间: {now}</div>
                </div>
                <div class="header-right">
                    <div class="lang-switch">
                        <button class="lang-btn active" data-lang="zh" onclick="switchLang('zh')">中文</button>
                        <button class="lang-btn" data-lang="en" onclick="switchLang('en')">EN</button>
                    </div>
                </div>
            </div>
        </div>
        """

    def _generate_input_files_section(self, ref_genome, qry_genome, genome_files):
        """
        生成输入文件模块（基因组文件）
        
        参数:
        - ref_genome: 参考基因组名称
        - qry_genome: 查询基因组名称
        - genome_files: 基因组文件路径字典
        """
        if not genome_files:
            return ""
        
        html = """
        <div class="section input-files-section">
            <h2>输入文件</h2>
            <table class="file-table">
                <tr><th>名称</th><th>描述</th><th>路径</th></tr>
        """
        
        # 参考基因组 FASTA
        ref_fasta = genome_files.get('ref_fasta', '')
        if ref_fasta:
            abs_path = os.path.abspath(ref_fasta)
            js_path = self._escape_path_for_js(abs_path)
            exists = os.path.exists(ref_fasta)
            name = f"参考基因组 ({ref_genome})" if ref_genome else "参考基因组"
            if exists:
                html += f'''<tr>
                    <td>{name}</td>
                    <td>FASTA 序列文件</td>
                    <td class="path-cell"><code class="path-code path-clickable" onclick="copyToClipboard('{js_path}')" title="点击复制路径">{abs_path}</code></td>
                </tr>'''
            else:
                html += f'''<tr>
                    <td>{name}</td>
                    <td>FASTA 序列文件</td>
                    <td><span class="file-na">文件不存在: {abs_path}</span></td>
                </tr>'''
        
        # 参考基因组 GFF
        ref_gff = genome_files.get('ref_gff', '')
        if ref_gff:
            abs_path = os.path.abspath(ref_gff)
            js_path = self._escape_path_for_js(abs_path)
            exists = os.path.exists(ref_gff)
            name = f"参考基因组注释 ({ref_genome})" if ref_genome else "参考基因组注释"
            if exists:
                html += f'''<tr>
                    <td>{name}</td>
                    <td>GFF/GFF3 注释文件</td>
                    <td class="path-cell"><code class="path-code path-clickable" onclick="copyToClipboard('{js_path}')" title="点击复制路径">{abs_path}</code></td>
                </tr>'''
            else:
                html += f'''<tr>
                    <td>{name}</td>
                    <td>GFF/GFF3 注释文件</td>
                    <td><span class="file-na">文件不存在: {abs_path}</span></td>
                </tr>'''
        
        # 查询基因组 FASTA
        qry_fasta = genome_files.get('qry_fasta', '')
        if qry_fasta:
            abs_path = os.path.abspath(qry_fasta)
            js_path = self._escape_path_for_js(abs_path)
            exists = os.path.exists(qry_fasta)
            name = f"查询基因组 ({qry_genome})" if qry_genome else "查询基因组"
            if exists:
                html += f'''<tr>
                    <td>{name}</td>
                    <td>FASTA 序列文件</td>
                    <td class="path-cell"><code class="path-code path-clickable" onclick="copyToClipboard('{js_path}')" title="点击复制路径">{abs_path}</code></td>
                </tr>'''
            else:
                html += f'''<tr>
                    <td>{name}</td>
                    <td>FASTA 序列文件</td>
                    <td><span class="file-na">文件不存在: {abs_path}</span></td>
                </tr>'''
        
        # 查询基因组 GFF
        qry_gff = genome_files.get('qry_gff', '')
        if qry_gff:
            abs_path = os.path.abspath(qry_gff)
            js_path = self._escape_path_for_js(abs_path)
            exists = os.path.exists(qry_gff)
            name = f"查询基因组注释 ({qry_genome})" if qry_genome else "查询基因组注释"
            if exists:
                html += f'''<tr>
                    <td>{name}</td>
                    <td>GFF/GFF3 注释文件</td>
                    <td class="path-cell"><code class="path-code path-clickable" onclick="copyToClipboard('{js_path}')" title="点击复制路径">{abs_path}</code></td>
                </tr>'''
            else:
                html += f'''<tr>
                    <td>{name}</td>
                    <td>GFF/GFF3 注释文件</td>
                    <td><span class="file-na">文件不存在: {abs_path}</span></td>
                </tr>'''
        
        html += "</table></div>"
        return html

    def _generate_input_files_section_with_extra(self, ref_genome, qry_genome, genome_files, extra_files=None):
        """
        生成输入文件模块（基因组文件 + 额外输入文件）
        
        参数:
        - ref_genome: 参考基因组名称
        - qry_genome: 查询基因组名称
        - genome_files: 基因组文件路径字典
        - extra_files: 额外的输入文件列表 [(名称, 描述, 路径), ...]
        """
        has_genome = genome_files and any(genome_files.values())
        has_extra = extra_files and len(extra_files) > 0
        
        if not has_genome and not has_extra:
            return ""
        
        html = """
        <div class="section input-files-section">
            <h2>输入文件</h2>
            <table class="file-table">
                <tr><th>名称</th><th>描述</th><th>路径</th></tr>
        """
        
        # 基因组文件
        if genome_files:
            # 参考基因组 FASTA
            ref_fasta = genome_files.get('ref_fasta', '')
            if ref_fasta:
                abs_path = os.path.abspath(ref_fasta)
                js_path = self._escape_path_for_js(abs_path)
                exists = os.path.exists(ref_fasta)
                name = f"参考基因组 ({ref_genome})" if ref_genome else "参考基因组"
                if exists:
                    html += f'''<tr>
                        <td>{name}</td>
                        <td>FASTA 序列文件</td>
                        <td class="path-cell"><code class="path-code path-clickable" onclick="copyToClipboard('{js_path}')" title="点击复制路径">{abs_path}</code></td>
                    </tr>'''
                else:
                    html += f'''<tr>
                        <td>{name}</td>
                        <td>FASTA 序列文件</td>
                        <td><span class="file-na">文件不存在: {abs_path}</span></td>
                    </tr>'''
            
            # 参考基因组 GFF
            ref_gff = genome_files.get('ref_gff', '')
            if ref_gff:
                abs_path = os.path.abspath(ref_gff)
                js_path = self._escape_path_for_js(abs_path)
                exists = os.path.exists(ref_gff)
                name = f"参考基因组注释 ({ref_genome})" if ref_genome else "参考基因组注释"
                if exists:
                    html += f'''<tr>
                        <td>{name}</td>
                        <td>GFF/GFF3 注释文件</td>
                        <td class="path-cell"><code class="path-code path-clickable" onclick="copyToClipboard('{js_path}')" title="点击复制路径">{abs_path}</code></td>
                    </tr>'''
                else:
                    html += f'''<tr>
                        <td>{name}</td>
                        <td>GFF/GFF3 注释文件</td>
                        <td><span class="file-na">文件不存在: {abs_path}</span></td>
                    </tr>'''
            
            # 查询基因组 FASTA
            qry_fasta = genome_files.get('qry_fasta', '')
            if qry_fasta:
                abs_path = os.path.abspath(qry_fasta)
                js_path = self._escape_path_for_js(abs_path)
                exists = os.path.exists(qry_fasta)
                name = f"查询基因组 ({qry_genome})" if qry_genome else "查询基因组"
                if exists:
                    html += f'''<tr>
                        <td>{name}</td>
                        <td>FASTA 序列文件</td>
                        <td class="path-cell"><code class="path-code path-clickable" onclick="copyToClipboard('{js_path}')" title="点击复制路径">{abs_path}</code></td>
                    </tr>'''
                else:
                    html += f'''<tr>
                        <td>{name}</td>
                        <td>FASTA 序列文件</td>
                        <td><span class="file-na">文件不存在: {abs_path}</span></td>
                    </tr>'''
            
            # 查询基因组 GFF
            qry_gff = genome_files.get('qry_gff', '')
            if qry_gff:
                abs_path = os.path.abspath(qry_gff)
                js_path = self._escape_path_for_js(abs_path)
                exists = os.path.exists(qry_gff)
                name = f"查询基因组注释 ({qry_genome})" if qry_genome else "查询基因组注释"
                if exists:
                    html += f'''<tr>
                        <td>{name}</td>
                        <td>GFF/GFF3 注释文件</td>
                        <td class="path-cell"><code class="path-code path-clickable" onclick="copyToClipboard('{js_path}')" title="点击复制路径">{abs_path}</code></td>
                    </tr>'''
                else:
                    html += f'''<tr>
                        <td>{name}</td>
                        <td>GFF/GFF3 注释文件</td>
                        <td><span class="file-na">文件不存在: {abs_path}</span></td>
                    </tr>'''
        
        # 额外的输入文件
        if extra_files:
            for name, desc, path in extra_files:
                if path:
                    abs_path = os.path.abspath(path)
                    js_path = self._escape_path_for_js(abs_path)
                    exists = os.path.exists(path)
                    if exists:
                        html += f'''<tr>
                            <td>{name}</td>
                            <td>{desc}</td>
                            <td class="path-cell"><code class="path-code path-clickable" onclick="copyToClipboard('{js_path}')" title="点击复制路径">{abs_path}</code></td>
                        </tr>'''
                    else:
                        html += f'''<tr>
                            <td>{name}</td>
                            <td>{desc}</td>
                            <td><span class="file-na">文件不存在: {abs_path}</span></td>
                        </tr>'''
        
        html += "</table></div>"
        return html

    def _generate_input_section(self, input_type, input_value, input_source=None, genome_files=None):
        """
        生成模块1：输入信息
        
        参数:
        - input_type: 输入类型 (Gene ID / Location / Sequence)
        - input_value: 用户实际输入的内容
        - input_source: 输入来源信息 {
            'type': 'direct' | 'file',  # 直接输入或文件输入
            'file_path': str,           # 如果是文件，文件的绝对路径
            'file_name': str            # 文件名
          }
        - genome_files: 使用的基因组/注释文件 {
            'ref_fasta': str,      # 参考基因组 FASTA 路径
            'ref_gff': str,        # 参考基因组 GFF 路径
            'qry_fasta': str,      # 查询基因组 FASTA 路径（可选）
          }
        """
        html = """
        <div class="section input-section">
            <h2>输入信息</h2>
            <table class="info-table">
        """
        html += f'<tr><td class="label">输入类型</td><td class="value">{input_type}</td></tr>'
        
        # 显示用户输入
        if input_source and input_source.get('type') == 'file':
            # 文件输入模式
            file_path = input_source.get('file_path', '')
            file_name = input_source.get('file_name', os.path.basename(file_path))
            dir_path = os.path.dirname(os.path.abspath(file_path)) if file_path else ''
            
            html += f'<tr><td class="label">输入方式</td><td class="value">文件输入</td></tr>'
            html += f'<tr><td class="label">输入文件</td><td class="value">'
            html += f'<code>{file_name}</code>'
            if dir_path:
                html += f' <button class="open-folder-btn" onclick="copyPath(\'{dir_path.replace(chr(92), chr(92)+chr(92))}\')" title="复制文件夹路径">📁 复制路径</button>'
            html += '</td></tr>'
            html += f'<tr><td class="label">文件路径</td><td class="value"><code class="path-code">{file_path}</code></td></tr>'
            
            # 显示文件内容预览
            if input_value:
                preview = input_value if len(input_value) <= 500 else input_value[:500] + "..."
                html += f'<tr><td class="label">文件内容</td><td class="value"><pre class="content-preview">{preview}</pre></td></tr>'
        else:
            # 直接输入模式
            html += f'<tr><td class="label">输入方式</td><td class="value">直接输入</td></tr>'
            html += f'<tr><td class="label">输入值</td><td class="value"><code>{input_value}</code></td></tr>'
        
        html += "</table>"
        
        # 使用的基因组/注释文件
        if genome_files:
            html += """
            <h3 style="margin-top: 20px; font-size: 14px; color: #666;">使用的基因组文件</h3>
            <table class="file-table">
                <tr><th>类型</th><th>文件路径</th><th>操作</th></tr>
            """
            
            file_types = [
                ('ref_fasta', '参考基因组 FASTA'),
                ('ref_gff', '参考基因组注释 GFF'),
                ('qry_fasta', '查询基因组 FASTA'),
            ]
            
            for key, label in file_types:
                path = genome_files.get(key, '')
                if path:
                    abs_path = os.path.abspath(path) if path else ""
                    dir_path = os.path.dirname(abs_path) if abs_path else ""
                    file_exists = os.path.exists(path) if path else False
                    
                    html += f'<tr><td>{label}</td>'
                    html += f'<td class="path-cell"><code class="path-code" title="{abs_path}">{self._truncate_path(abs_path, 50)}</code></td>'
                    if file_exists and dir_path:
                        html += f'''<td><button class="open-folder-btn" onclick="copyPath('{dir_path.replace(chr(92), chr(92)+chr(92))}')" title="复制文件夹路径">📁 复制路径</button></td>'''
                    else:
                        html += '<td><span class="file-na">文件不存在</span></td>'
                    html += '</tr>'
            
            html += "</table>"
        
        html += "</div>"
        return html

    def _generate_output_section(self, output_files):
        """
        生成结果文件模块
        
        参数:
        - output_files: 生成的文件列表 [(文件名, 文件路径, 描述), ...]
        """
        if not output_files:
            return ""
        
        html = """
        <div class="section output-section">
            <h2>结果文件</h2>
            <table class="file-table">
                <tr><th>文件名</th><th>描述</th><th>路径</th></tr>
        """
        
        for name, path, desc in output_files:
            abs_path = os.path.abspath(path) if path else ""
            rel_path = os.path.basename(path) if path else ""
            file_exists = os.path.exists(path) if path else False
            js_abs_path = self._escape_path_for_js(abs_path)
            
            if file_exists:
                html += f'''<tr>
                    <td><a href="{rel_path}" class="file-name-link" title="点击查看文件">{name}</a></td>
                    <td>{desc}</td>
                    <td class="path-cell"><code class="path-code path-clickable" onclick="copyToClipboard('{js_abs_path}')" title="点击复制路径">{abs_path}</code></td>
                </tr>'''
            else:
                html += f'''<tr>
                    <td><code class="file-name-disabled">{name}</code></td>
                    <td>{desc}</td>
                    <td><span class="file-na">未生成</span></td>
                </tr>'''
        
        html += "</table></div>"
        return html

    def _truncate_path(self, path, max_len=40):
        """截断过长的路径，保留开头和结尾"""
        if len(path) <= max_len:
            return path
        return path[:15] + "..." + path[-(max_len-18):]

    def _get_html_header(self, title):
        """生成 HTML 头部"""
        return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title} - GeneScreen Report</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.6; 
            color: #333; 
            background: #f5f5f5;
            padding: 20px;
        }}
        .container {{ 
            max-width: 1200px; 
            margin: 0 auto; 
            background: #fff;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
        }}
        .header-content {{
            display: flex;
            justify-content: space-between;
            align-items: flex-start;
        }}
        .header-left {{
            display: flex;
            flex-direction: column;
            align-items: flex-start;
        }}
        .header-right {{
            display: flex;
            align-items: flex-start;
        }}
        .header h1 {{ font-size: 24px; margin-bottom: 10px; }}
        .header .subtitle {{ 
            display: flex;
            align-items: center;
            gap: 15px;
            opacity: 0.95; 
            font-size: 14px; 
        }}
        .header .gen-time {{
            margin-top: 8px;
            font-size: 13px;
            opacity: 0.85;
            text-align: left;
        }}
        .lang-switch {{
            display: flex;
            gap: 0;
            border-radius: 6px;
            overflow: hidden;
            border: 1px solid rgba(255,255,255,0.3);
        }}
        .lang-btn {{
            background: transparent;
            border: none;
            color: rgba(255,255,255,0.7);
            padding: 6px 14px;
            font-size: 13px;
            cursor: pointer;
            transition: all 0.2s;
        }}
        .lang-btn:hover {{
            background: rgba(255,255,255,0.1);
            color: white;
        }}
        .lang-btn.active {{
            background: rgba(255,255,255,0.2);
            color: white;
            font-weight: 600;
        }}
        .section {{ padding: 25px 30px; border-bottom: 1px solid #eee; }}
        .section:last-child {{ border-bottom: none; }}
        .section h2 {{ 
            font-size: 18px; 
            color: #667eea; 
            margin-bottom: 15px;
            padding-bottom: 10px;
            border-bottom: 2px solid #667eea;
        }}
        .section h3 {{ font-size: 14px; color: #666; margin-bottom: 10px; }}
        .info-table {{ width: 100%; border-collapse: collapse; margin-bottom: 15px; }}
        .info-table td {{ padding: 10px 15px; border-bottom: 1px solid #eee; }}
        .info-table .label {{ 
            width: 150px; 
            font-weight: 600; 
            color: #666;
            background: #f9f9f9;
        }}
        .info-table .value {{ color: #333; }}
        .info-table code {{ 
            background: #f0f0f0; 
            padding: 2px 6px; 
            border-radius: 3px; 
            font-family: monospace;
        }}
        .info-table-grid {{ width: 100%; border-collapse: collapse; margin-bottom: 15px; }}
        .info-table-grid td {{ padding: 10px 15px; border-bottom: 1px solid #eee; }}
        .info-table-grid .label {{ 
            width: 120px; 
            font-weight: 600; 
            color: #666;
            background: #f9f9f9;
        }}
        .info-table-grid .value {{ color: #333; width: 30%; }}
        .stats-table {{ width: 100%; border-collapse: collapse; margin-bottom: 15px; text-align: center; }}
        .stats-table th {{ 
            padding: 12px 15px; 
            background: #f9f9f9; 
            font-weight: 600; 
            color: #666;
            border-bottom: 2px solid #eee;
        }}
        .stats-table td {{ 
            padding: 15px; 
            font-size: 18px; 
            font-weight: 500; 
            color: #333;
            border-bottom: 1px solid #eee;
        }}
        .file-table {{ width: 100%; border-collapse: collapse; }}
        .file-table th, .file-table td {{ 
            padding: 10px 15px; 
            text-align: left; 
            border-bottom: 1px solid #eee; 
        }}
        .file-table th {{ background: #f9f9f9; font-weight: 600; color: #666; }}
        .file-link {{ 
            color: #667eea; 
            text-decoration: none; 
            padding: 4px 10px;
            border: 1px solid #667eea;
            border-radius: 4px;
            font-size: 12px;
        }}
        .file-link:hover {{ background: #667eea; color: white; }}
        .file-name-link {{
            color: #667eea;
            text-decoration: none;
            font-family: monospace;
            font-weight: 500;
        }}
        .file-name-link:hover {{ text-decoration: underline; }}
        .file-name-disabled {{
            color: #999;
            font-family: monospace;
        }}
        .file-na {{ color: #999; font-size: 12px; }}
        .path-code {{ 
            font-size: 11px; 
            color: #666; 
            word-break: break-all;
            background: #f5f5f5;
            padding: 2px 6px;
            border-radius: 3px;
        }}
        .path-clickable {{
            cursor: pointer;
        }}
        .path-clickable:hover {{
            background: #e8e8e8;
            color: #333;
        }}
        .path-cell {{ max-width: 500px; }}
        .content-preview {{
            background: #f9f9f9;
            border: 1px solid #eee;
            border-radius: 4px;
            padding: 10px;
            font-family: monospace;
            font-size: 12px;
            white-space: pre-wrap;
            word-break: break-all;
            max-height: 150px;
            overflow-y: auto;
            margin: 0;
        }}
        .open-folder-btn {{
            background: #f0f0f0;
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 2px 8px;
            font-size: 12px;
            cursor: pointer;
            margin-left: 8px;
        }}
        .open-folder-btn:hover {{ background: #e0e0e0; }}
        .copy-toast {{
            position: fixed;
            bottom: 20px;
            right: 20px;
            background: #333;
            color: white;
            padding: 10px 20px;
            border-radius: 4px;
            display: none;
            z-index: 1000;
        }}
        .mode-badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 12px;
            font-weight: 600;
            text-transform: uppercase;
        }}
        .mode-geneid {{ background: #e3f2fd; color: #1976d2; }}
        .mode-location {{ background: #f3e5f5; color: #7b1fa2; }}
        .mode-sequence {{ background: #e8f5e9; color: #388e3c; }}
        .visualization-section {{ background: #fafafa; }}
        .viz-container {{ 
            padding: 20px;
            background: white;
            border: 1px solid #eee;
            border-radius: 8px;
            overflow-x: auto;
        }}
        .viz-container svg {{ display: block; margin: 0 auto; width: 100%; height: auto; }}
        .viz-legend-note {{
            margin-top: 15px;
            padding: 10px 15px;
            background: #f5f5f5;
            border-radius: 4px;
            font-size: 13px;
            color: #666;
        }}
        .viz-legend-note p {{ margin: 0; }}
    </style>
</head>
<body>
<div class="container">
"""

    def _get_html_footer(self):
        """生成 HTML 尾部"""
        return r"""
<div id="copyToast" class="copy-toast">已复制到剪贴板</div>
</div>
<script>
// 国际化翻译数据
const translations = {
    // 标题和模块名
    '初始化': 'Initialization',
    '输入文件': 'Input Files',
    '结果文件': 'Output Files',
    '比对可视化': 'Alignment Visualization',
    '宏观比对图': 'Alignment Overview',
    '局部比对图': 'Local Alignment',
    // 表头
    '名称': 'Name',
    '描述': 'Description',
    '路径': 'Path',
    '文件名': 'Filename',
    // 标签
    '参考基因组': 'Reference Genome',
    '查询基因组': 'Query Genome',
    '基因 ID': 'Gene ID',
    'Identity 阈值': 'Identity Threshold',
    '区域': 'Region',
    '序列 ID': 'Sequence ID',
    '序列长度': 'Sequence Length',
    '最小比对长度': 'Min Alignment Length',
    '上游延伸': 'Upstream Extension',
    '下游延伸': 'Downstream Extension',
    '链方向': 'Strand',
    '正链 (+)': 'Forward (+)',
    '正链': 'Forward',
    '负链 (-)': 'Reverse (-)',
    '负链': 'Reverse',
    '注释': 'Annotation',
    '参考基因组注释': 'Reference Annotation',
    '输入序列': 'Input Sequence',
    '选择比对区域：': 'Select region:',
    '比对长度': 'Alignment Length',
    // 描述文本
    'FASTA 序列文件': 'FASTA Sequence File',
    'GFF/GFF3 注释文件': 'GFF/GFF3 Annotation File',
    '用户提供的 FASTA 序列文件': 'User-provided FASTA Sequence File',
    '参考基因组提取的基因序列': 'Gene Sequence Extracted from Reference',
    '基因注释文件（相对坐标）': 'Gene Annotation (Relative Coordinates)',
    'BLAST 比对结果': 'BLAST Alignment Result',
    '比对坐标文件': 'Alignment Coordinates File',
    'SNP/Indel 变异文件': 'SNP/Indel Variant File',
    '提取的区域序列': 'Extracted Region Sequence',
    // 状态文本
    '文件不存在': 'File not found',
    '未生成': 'Not generated',
    // 提示文本
    '点击复制路径': 'Click to copy path',
    '点击查看文件': 'Click to view file',
    '路径已复制到剪贴板': 'Path copied to clipboard',
    '点击比对区域（蓝色块）可跳转到对应的局部比对图': 'Click alignment region (blue block) to jump to the corresponding local alignment'
};

// 反向翻译映射 (英文 -> 中文)
const reverseTranslations = {};
for (const [zh, en] of Object.entries(translations)) {
    reverseTranslations[en] = zh;
}

// 按键长度倒序，优先匹配更长词组
const translationEntries = Object.entries(translations).sort((a, b) => b[0].length - a[0].length);
const reverseEntries = Object.entries(reverseTranslations).sort((a, b) => b[0].length - a[0].length);

let currentLang = 'zh';

function replaceWithSpacing(text, from, to) {
    const idx = text.indexOf(from);
    if (idx === -1) return text;
    const before = text.slice(0, idx);
    const after = text.slice(idx + from.length);
    const left = before.slice(-1);
    const right = after.slice(0, 1);
    const needsLeftSpace = left && /[A-Za-z0-9]/.test(left) && !/^\s/.test(to);
    const needsRightSpace = right && /[A-Za-z0-9]/.test(right) && !/\s$/.test(to);
    const middle = (needsLeftSpace ? ' ' : '') + to + (needsRightSpace ? ' ' : '');
    return before + middle + after;
}

function translateText(text, toLang) {
    if (toLang === 'en') {
        // 精确匹配
        if (translations[text]) return translations[text];
        // 处理带括号的文本，如 "参考基因组 (AGIS-1.0)"
        for (const [zh, en] of translationEntries) {
            if (text.startsWith(zh + ' (')) {
                const suffix = text.substring(zh.length);
                return en + suffix;
            }
            if (text.includes(zh)) {
                return replaceWithSpacing(text, zh, en);
            }
        }
    } else {
        // 精确匹配
        if (reverseTranslations[text]) return reverseTranslations[text];
        // 处理带括号的文本
        for (const [en, zh] of reverseEntries) {
            if (text.startsWith(en + ' (')) {
                const suffix = text.substring(en.length);
                return zh + suffix;
            }
            if (text.includes(en)) {
                return replaceWithSpacing(text, en, zh);
            }
        }
    }
    return text;
}

function switchLang(lang) {
    currentLang = lang;
    
    // 更新按钮状态
    document.querySelectorAll('.lang-btn').forEach(btn => {
        btn.classList.toggle('active', btn.dataset.lang === lang);
    });
    
    // 更新带有 data-zh/data-en 属性的元素
    document.querySelectorAll('[data-zh][data-en]').forEach(el => {
        el.textContent = el.dataset[lang];
    });
    
    // 更新所有需要翻译的文本元素
    const selectors = [
        'h2',           // 模块标题
        'h3',           // 子标题
        'th',           // 表头
        'td.label',     // 标签列
        '.info-table-grid td.value', // 参数值
        '.file-na',     // 状态文本
        '.file-table td'  // 文件表格（名称/描述）
    ];
    
    document.querySelectorAll(selectors.join(', ')).forEach(el => {
        // 跳过包含 code 或 a 标签的元素（路径和链接）
        if (el.querySelector('code, a')) return;
        // 跳过路径单元格
        if (el.classList.contains('path-cell')) return;
        
        const text = el.textContent.trim();
        if (text) {
            const translated = translateText(text, lang);
            if (translated !== text) {
                el.textContent = translated;
            }
        }
    });
    
    // 更新 title 属性
    document.querySelectorAll('[title]').forEach(el => {
        const title = el.getAttribute('title');
        const translated = translateText(title, lang);
        if (translated !== title) {
            el.setAttribute('title', translated);
        }
    });
}

function copyToClipboard(text) {
    navigator.clipboard.writeText(text).then(function() {
        showToast(currentLang === 'zh' ? '路径已复制到剪贴板' : 'Path copied to clipboard');
    }).catch(function(err) {
        // 降级方案
        var temp = document.createElement('input');
        temp.value = text;
        document.body.appendChild(temp);
        temp.select();
        document.execCommand('copy');
        document.body.removeChild(temp);
        showToast(currentLang === 'zh' ? '路径已复制到剪贴板' : 'Path copied to clipboard');
    });
}

function showToast(msg) {
    var toast = document.getElementById('copyToast');
    toast.textContent = msg;
    toast.style.display = 'block';
    setTimeout(function() {
        toast.style.display = 'none';
    }, 2000);
}
</script>
</body>
</html>
"""


class GeneIDVisualizer(BaseVisualizer):
    """Gene ID 模式可视化"""

    def visualize(
        self,
        result,
        ref_genome="",
        qry_genome="",
        identity=90,
        input_source=None,
        genome_files=None,
        legacy_only=False,
        report_suffix=None,
    ):
        """
        生成 Gene ID 模式的报告
        
        result 包含:
        - id: 基因 ID
        - fasta: 基因序列文件
        - gff: 基因注释文件（相对坐标）
        - blast_xml: BLAST 比对结果
        - coords: 比对坐标文件
        - snps: SNP/Indel 文件
        - prefix: 输出前缀
        - extraction_info: 提取信息（可选，包含 gene_rel_start, gene_rel_end, strand, upstream, downstream）
        
        genome_files: 使用的基因组文件路径
        - ref_fasta: 参考基因组 FASTA
        - ref_gff: 参考基因组 GFF
        - qry_fasta: 查询基因组 FASTA
        - qry_gff: 查询基因组 GFF（可选，用于提取比对区域的注释）
        """
        gene_id = result["id"]
        prefix = result.get("prefix", os.path.join(self.output_dir, gene_id))
        extraction_info = result.get("extraction_info", None)
        print(f"[INFO] Gene ID 可视化: {gene_id}")

        # 生成 HTML 报告
        html = self._get_html_header(gene_id)
        
        # Header（标题栏）
        html += self._generate_header("Gene ID 模式", "mode-geneid")
        
        # 模块0：初始化
        params = {
            "基因 ID": gene_id,
            "Identity 阈值": f"{identity}%",
            "最小比对长度": f"{self.min_aln_len} bp"
        }
        # 如果有上下游延伸信息，添加到参数中
        if extraction_info:
            if extraction_info.get("upstream", 0) > 0:
                params["上游延伸"] = f"{extraction_info['upstream']} bp"
            if extraction_info.get("downstream", 0) > 0:
                params["下游延伸"] = f"{extraction_info['downstream']} bp"
            if extraction_info.get("strand"):
                params["链方向"] = "正链 (+)" if extraction_info["strand"] == "+" else "负链 (-)"
        
        html += self._generate_mode_section("Gene ID", ref_genome, qry_genome, params, genome_files)
        
        # 输入文件模块（基因组文件）
        html += self._generate_input_files_section(ref_genome, qry_genome, genome_files)
        
        # 比对可视化模块 - 优先使用 LINKVIEW
        fasta_file = result.get("fasta", "")
        snps_file = result.get("snps", "")
        coords_file = result.get("coords", "")
        gff_file = result.get("gff", "")  # linkview.gff3 (相对坐标，来自参考基因组)
        
        # 尝试从查询基因组 GFF 提取比对区域的注释
        qry_gff_file = None
        if genome_files and genome_files.get('qry_gff') and coords_file:
            qry_gff_file = self._extract_aligned_region_gff(
                genome_files.get('qry_gff'), coords_file, gene_id, self.output_dir, track_name="qry"
            )
        
        # 结果统计模块已移至宏观图下方（_embed_overview_section 中）
        # html += self._generate_stats_section(fasta_file, coords_file, snps_file)
        
        # 尝试使用 LINKVIEW，失败则回退到内置 SVG
        # 传递 extraction_info 用于绘制基因区域大括号
        html += self._generate_linkview_visualization_section(
            gene_id, fasta_file, coords_file, snps_file, gff_file, qry_gff_file, extraction_info
        )
        
        # 结果文件模块
        output_files = [
            (f"{gene_id}.fasta", result.get("fasta", ""), "参考基因组提取的基因序列"),
            (f"{gene_id}.gff3", result.get("gff", ""), "基因注释文件（相对坐标）"),
            (f"{gene_id}.blast.xml", result.get("blast_xml", ""), "BLAST 比对结果"),
            (f"{gene_id}.coords", result.get("coords", ""), "比对坐标文件"),
            (f"{gene_id}.snps", result.get("snps", ""), "SNP/Indel 变异文件"),
        ]
        # 如果生成了查询基因组区域 GFF，也添加到结果文件
        if qry_gff_file:
            output_files.append((f"{gene_id}.qry.linkview.gff3", qry_gff_file, "查询基因组比对区域注释（绝对坐标）"))
        
        html += self._generate_output_section(output_files)
        
        html += self._get_html_footer()
        
        # 保存报告
        target_suffix = report_suffix
        if target_suffix is None and len(result.get("query_results") or []) > 1:
            target_suffix = self._legacy_report_suffix_for_result(result, genome_files)
        report_file = self._legacy_report_path(gene_id, target_suffix)
        with open(report_file, "w", encoding="utf-8") as f:
            f.write(html)
        
        print(f"[INFO] 已生成报告: {report_file}")
        if legacy_only or not self._has_multiple_query_reports(result):
            return report_file
        legacy_report = self._build_legacy_report_map(
            result, "gene_id", report_file, ref_genome, qry_genome, identity, input_source, genome_files
        )
        report_index = self._finalize_single_query_report(
            result,
            "gene_id",
            ref_genome,
            qry_genome,
            identity,
            genome_files,
            legacy_report,
            output_files,
        )
        print(f"[INFO] 已生成新框架报告: {report_index}")
        return report_index


class LocationVisualizer(BaseVisualizer):
    """Location 模式可视化"""

    def visualize(
        self,
        result,
        ref_genome="",
        qry_genome="",
        identity=90,
        input_source=None,
        genome_files=None,
        legacy_only=False,
        report_suffix=None,
    ):
        """
        生成 Location 模式的报告
        
        result 包含:
        - id: 区域名称
        - location: 原始位置字符串 (chr:start-end)
        - fasta: 区域序列文件
        - blast_xml: BLAST 比对结果
        - coords: 比对坐标文件
        - snps: SNP/Indel 文件
        - prefix: 输出前缀
        
        genome_files: 使用的基因组文件路径
        - ref_gff: 参考基因组 GFF（可选，用于提取区域内的基因注释）
        - qry_gff: 查询基因组 GFF（可选，用于提取比对区域的注释）
        """
        loc_id = result["id"]
        location = result.get("location", loc_id.replace("_", ":").replace(":", ":", 1).replace("_", "-"))
        prefix = result.get("prefix", os.path.join(self.output_dir, loc_id))
        print(f"[INFO] Location 可视化: {loc_id}")

        # 生成 HTML 报告
        html = self._get_html_header(loc_id)
        
        # Header（标题栏）
        html += self._generate_header("Location 模式", "mode-location")
        
        # 模块0：初始化
        params = {
            "区域": location,
            "Identity 阈值": f"{identity}%",
            "最小比对长度": f"{self.min_aln_len} bp"
        }
        html += self._generate_mode_section("Location", ref_genome, qry_genome, params, genome_files)
        
        # 输入文件模块（基因组文件）
        html += self._generate_input_files_section(ref_genome, qry_genome, genome_files)
        
        # 比对可视化模块
        fasta_file = result.get("fasta", "")
        snps_file = result.get("snps", "")
        coords_file = result.get("coords", "")
        
        # 尝试从参考基因组 GFF 提取区域内的基因注释（上方轨道）
        ref_gff_file = None
        if genome_files and genome_files.get('ref_gff'):
            ref_gff = genome_files.get('ref_gff')
            chr_name, start, end = self._parse_location(location)
            if chr_name and start and end:
                ref_gff_file = self._extract_region_gff(ref_gff, loc_id, chr_name, start, end, self.output_dir)
        
        # 尝试从查询基因组 GFF 提取比对区域的注释（下方轨道）
        qry_gff_file = None
        if genome_files and genome_files.get('qry_gff') and coords_file:
            qry_gff_file = self._extract_aligned_region_gff(
                genome_files.get('qry_gff'), coords_file, loc_id, self.output_dir, track_name="qry"
            )
        
        # 结果统计模块已移至宏观图下方（_embed_overview_section 中）
        # html += self._generate_stats_section(fasta_file, coords_file, snps_file)
        
        # 比对可视化（如果有 GFF 则显示基因结构，否则只显示 SNP/Indel）
        # 分别传递 ref_gff（上方轨道，相对坐标）和 qry_gff（下方轨道，绝对坐标）
        html += self._generate_linkview_visualization_section(
            loc_id, fasta_file, coords_file, snps_file, ref_gff_file, qry_gff_file
        )
        
        # 结果文件模块
        output_files = [
            (f"{loc_id}.fasta", result.get("fasta", ""), "提取的区域序列"),
            (f"{loc_id}.blast.xml", result.get("blast_xml", ""), "BLAST 比对结果"),
            (f"{loc_id}.coords", result.get("coords", ""), "比对坐标文件"),
            (f"{loc_id}.snps", result.get("snps", ""), "SNP/Indel 变异文件"),
        ]
        # 如果生成了区域 GFF，也添加到结果文件
        if ref_gff_file:
            output_files.append((f"{loc_id}.linkview.gff3", ref_gff_file, "参考基因组区域注释（相对坐标）"))
        if qry_gff_file:
            output_files.append((f"{loc_id}.qry.linkview.gff3", qry_gff_file, "查询基因组比对区域注释（绝对坐标）"))
        
        html += self._generate_output_section(output_files)
        
        html += self._get_html_footer()
        
        # 保存报告
        target_suffix = report_suffix
        if target_suffix is None and len(result.get("query_results") or []) > 1:
            target_suffix = self._legacy_report_suffix_for_result(result, genome_files)
        report_file = self._legacy_report_path(loc_id, target_suffix)
        with open(report_file, "w", encoding="utf-8") as f:
            f.write(html)
        
        print(f"[INFO] 已生成报告: {report_file}")
        if legacy_only or not self._has_multiple_query_reports(result):
            return report_file
        legacy_report = self._build_legacy_report_map(
            result, "location", report_file, ref_genome, qry_genome, identity, input_source, genome_files
        )
        report_index = self._finalize_single_query_report(
            result,
            "location",
            ref_genome,
            qry_genome,
            identity,
            genome_files,
            legacy_report,
            output_files,
        )
        print(f"[INFO] 已生成新框架报告: {report_index}")
        return report_index


class SequenceVisualizer(BaseVisualizer):
    """Sequence 模式可视化"""

    def visualize(
        self,
        result,
        ref_genome="",
        qry_genome=None,
        identity=90,
        input_source=None,
        genome_files=None,
        legacy_only=False,
        report_suffix=None,
    ):
        """
        生成 Sequence 模式的报告
        
        Sequence 模式：用户提供序列，比对到参考基因组
        - 不需要查询基因组
        - genome_files 只需要 ref_fasta
        
        result 包含:
        - id: 序列 ID
        - fasta: 序列文件
        - sequence: 原始序列内容
        - blast_xml: BLAST 比对结果
        - coords: 比对坐标文件
        - snps: SNP/Indel 文件
        - prefix: 输出前缀
        
        genome_files: 使用的基因组文件路径
        - ref_fasta: 参考基因组 FASTA（用于 BLAST 比对）
        - ref_gff: 参考基因组 GFF（可选，用于提取比对区域的注释）
        """
        seq_id = result["id"]
        sequence = result.get("sequence", "")
        input_fasta = result.get("fasta", "")
        prefix = result.get("prefix", os.path.join(self.output_dir, seq_id))
        print(f"[INFO] Sequence 可视化: {seq_id}")

        # 生成 HTML 报告
        html = self._get_html_header(seq_id)
        
        # Header（标题栏）
        html += self._generate_header("Sequence 模式", "mode-sequence")
        
        # 模块0：初始化（Sequence 模式不需要查询基因组）
        seq_len = len(sequence) if sequence else 0
        params = {
            "序列 ID": seq_id,
            "序列长度": f"{seq_len:,} bp" if seq_len else "未知",
            "Identity 阈值": f"{identity}%",
            "最小比对长度": f"{self.min_aln_len} bp"
        }
        html += self._generate_mode_section("Sequence", ref_genome, None, params, genome_files)
        
        # 输入文件模块（包含参考基因组和用户输入的序列文件）
        seq_genome_files = None
        if genome_files:
            seq_genome_files = {
                'ref_fasta': genome_files.get('ref_fasta', ''),
                'ref_gff': genome_files.get('ref_gff', '')
            }
        
        # 添加用户输入的序列文件到输入文件列表
        extra_input_files = []
        if input_fasta:
            extra_input_files.append(("输入序列", "用户提供的 FASTA 序列文件", input_fasta))
        
        html += self._generate_input_files_section_with_extra(ref_genome, None, seq_genome_files, extra_input_files)
        
        # 比对可视化模块
        fasta_file = result.get("fasta", "")
        snps_file = result.get("snps", "")
        coords_file = result.get("coords", "")
        
        # 尝试从参考基因组 GFF 提取比对区域的注释（下方轨道，参考基因组）
        ref_gff_file = None
        if genome_files and genome_files.get('ref_gff') and coords_file:
            ref_gff_file = self._extract_aligned_region_gff(
                genome_files.get('ref_gff'), coords_file, seq_id, self.output_dir, track_name="ref"
            )
        
        # 结果统计模块（已移至宏观图下方的区域统计卡片，见任务5）
        # html += self._generate_stats_section(fasta_file, coords_file, snps_file)
        
        # 比对可视化（如果有 GFF 则显示基因结构，否则只显示 SNP/Indel）
        # Sequence 模式：上方轨道是用户输入序列（无注释），下方轨道是参考基因组（有注释，绝对坐标）
        html += self._generate_linkview_visualization_section(
            seq_id, fasta_file, coords_file, snps_file, None, ref_gff_file
        )
        
        # 结果文件模块（不包含输入序列文件）
        output_files = [
            (f"{seq_id}.blast.xml", result.get("blast_xml", ""), "BLAST 比对结果"),
            (f"{seq_id}.coords", result.get("coords", ""), "比对坐标文件"),
            (f"{seq_id}.snps", result.get("snps", ""), "SNP/Indel 变异文件"),
        ]
        # 如果生成了参考基因组区域 GFF，也添加到结果文件
        if ref_gff_file:
            output_files.append((f"{seq_id}.ref.linkview.gff3", ref_gff_file, "参考基因组比对区域注释（绝对坐标）"))
        
        html += self._generate_output_section(output_files)
        
        html += self._get_html_footer()
        
        # 保存报告
        target_suffix = report_suffix
        if target_suffix is None and len(result.get("query_results") or []) > 1:
            target_suffix = self._legacy_report_suffix_for_result(result, genome_files)
        report_file = self._legacy_report_path(seq_id, target_suffix)
        with open(report_file, "w", encoding="utf-8") as f:
            f.write(html)
        
        print(f"[INFO] 已生成报告: {report_file}")
        if legacy_only or not self._has_multiple_query_reports(result):
            return report_file
        legacy_report = self._build_legacy_report_map(
            result, "sequence", report_file, ref_genome, qry_genome, identity, input_source, genome_files
        )
        report_index = self._finalize_single_query_report(
            result,
            "sequence",
            ref_genome,
            qry_genome,
            identity,
            genome_files,
            legacy_report,
            output_files,
        )
        print(f"[INFO] 已生成新框架报告: {report_index}")
        return report_index
