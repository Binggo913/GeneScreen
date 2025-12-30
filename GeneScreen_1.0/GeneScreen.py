#!/usr/bin/env python3
"""
GeneScreen - 基因组比对与变异分析工具

支持三种输入模式：
1. Gene ID 模式: 从参考基因组提取指定基因，与目标基因组比对
2. Location 模式: 从参考基因组提取指定区域，与目标基因组比对
3. Sequence 模式: 用户提供序列，与参考基因组比对
"""

import os
import subprocess
import argparse

from GeneScreenVisualizer import GeneIDVisualizer, LocationVisualizer, SequenceVisualizer
from GenomeManager import GenomeManager


# ======================= 工具函数 =======================
def run_cmd(cmd, error_msg="命令执行失败"):
    """执行 shell 命令，失败时打印错误信息"""
    print(f"[CMD] {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"[ERROR] {error_msg}")
        return False
    return True


# 全局基因组管理器
genome_manager = GenomeManager()


# ======================= 序列提取模块 =======================
class SequenceExtractor:
    """从参考基因组提取序列"""
    
    def __init__(self, genome, annotation, output_dir):
        self.genome = genome
        self.annotation = annotation
        self.output_dir = output_dir
    
    def extract_by_gene_id(self, gene_id):
        """通过 Gene ID 提取序列"""
        print(f"[INFO] 提取基因 {gene_id} 的序列...")
        
        # 从注释文件提取基因信息
        original_gff = os.path.join(self.output_dir, f"{gene_id}.original.gff3")
        linkview_gff = os.path.join(self.output_dir, f"{gene_id}.linkview.gff3")
        
        awk_cmd = f"awk '$9 ~ /{gene_id}/' {self.annotation} > {original_gff}"
        if not run_cmd(awk_cmd, f"提取 {gene_id} 注释失败"):
            return None, None
        
        # 解析 GFF，生成相对坐标版本
        gene_info = self._parse_gff(original_gff, linkview_gff, gene_id)
        if not gene_info:
            return None, None
        
        # 提取序列
        fasta_file = self._extract_fasta(gene_id, gene_info)
        return fasta_file, linkview_gff
    
    def extract_by_location(self, chrom, start, end, name=None):
        """通过染色体坐标提取序列"""
        loc_name = name or f"{chrom}_{start}_{end}"
        print(f"[INFO] 提取区域 {chrom}:{start}-{end} 的序列...")
        
        # 创建 BED 文件
        bed_file = os.path.join(self.output_dir, f"{loc_name}.bed")
        with open(bed_file, 'w') as f:
            f.write(f"{chrom}\t{start}\t{end}\t{loc_name}\n")
        
        # 提取序列
        fasta_file = os.path.join(self.output_dir, f"{loc_name}.fasta")
        cmd = f"bedtools getfasta -fi {self.genome} -bed {bed_file} -fo {fasta_file} -name"
        if not run_cmd(cmd, f"提取 {loc_name} 序列失败"):
            return None
        
        print(f"[INFO] 已生成序列文件: {fasta_file}")
        return fasta_file
    
    def _parse_gff(self, original_gff, linkview_gff, gene_id):
        """解析 GFF 文件，返回基因坐标信息"""
        start_offset = None
        gene_info = {}
        
        with open(original_gff, 'r') as f_in, open(linkview_gff, 'w') as f_out:
            for line in f_in:
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                if start_offset is None:
                    start_offset = int(fields[3])
                
                # 写入相对坐标版本
                rel_start = int(fields[3]) - start_offset + 1
                rel_end = int(fields[4]) - start_offset + 1
                f_out.write(f"{gene_id}\t{fields[1]}\t{fields[2]}\t{rel_start}\t{rel_end}\t{fields[5]}\t{fields[6]}\t{fields[7]}\t{fields[8]}\n")
                
                if fields[2] == "gene":
                    gene_info = {
                        "chrom": fields[0],
                        "start": fields[3],
                        "end": fields[4]
                    }
        
        if not gene_info:
            print(f"[ERROR] 未找到 {gene_id} 的 gene 记录")
            return None
        return gene_info
    
    def _extract_fasta(self, gene_id, gene_info):
        """使用 bedtools 提取 FASTA 序列"""
        bed_file = os.path.join(self.output_dir, f"{gene_id}.bed")
        fasta_file = os.path.join(self.output_dir, f"{gene_id}.fasta")
        
        with open(bed_file, 'w') as f:
            f.write(f"{gene_info['chrom']}\t{gene_info['start']}\t{gene_info['end']}\t{gene_id}\n")
        
        cmd = f"bedtools getfasta -fi {self.genome} -bed {bed_file} -fo {fasta_file} -name"
        if not run_cmd(cmd, f"提取 {gene_id} 序列失败"):
            return None
        
        print(f"[INFO] 已生成序列文件: {fasta_file}")
        return fasta_file


# ======================= MUMmer 比对模块 =======================
class MummerAligner:
    """MUMmer 基因组比对"""
    
    def __init__(self, output_dir):
        self.output_dir = output_dir
    
    def align(self, reference, query, prefix_name):
        """执行比对流程"""
        prefix = os.path.join(self.output_dir, prefix_name)
        
        # nucmer 比对
        print(f"[INFO] 运行 nucmer 比对...")
        delta_file = f"{prefix}.delta"
        cmd = f"nucmer -p {prefix} {reference} {query} -t 2"
        if not run_cmd(cmd, "nucmer 比对失败"):
            return None
        
        # delta-filter 过滤
        filtered_file = f"{prefix}.filter"
        cmd = f"delta-filter -1 -i 90 -l 100 {delta_file} > {filtered_file}"
        if not run_cmd(cmd, "delta-filter 失败"):
            return None
        
        # 生成 coords 和 snps
        coords_file = f"{prefix}.coords"
        snps_file = f"{prefix}.snps"
        
        cmd = f"show-coords -r -c -l {filtered_file} > {coords_file}"
        if not run_cmd(cmd, "生成 coords 失败"):
            return None
        
        cmd = f"show-snps -C -H -T -r -l {filtered_file} > {snps_file}"
        if not run_cmd(cmd, "生成 snps 失败"):
            return None
        
        print(f"[INFO] 比对完成，结果文件: {coords_file}, {snps_file}")
        return {
            "delta": delta_file,
            "filtered": filtered_file,
            "coords": coords_file,
            "snps": snps_file,
            "prefix": prefix
        }


# ======================= 输入处理器 =======================
class GeneIDProcessor:
    """Gene ID 模式处理器"""
    
    def __init__(self, ref_genome, ref_annotation, query_genome, output_dir):
        self.extractor = SequenceExtractor(ref_genome, ref_annotation, output_dir)
        self.aligner = MummerAligner(output_dir)
        self.query_genome = query_genome
        self.output_dir = output_dir
    
    def process(self, gene_id):
        print(f"\n{'='*50}")
        print(f"[INFO] 处理基因: {gene_id}")
        print(f"{'='*50}")
        
        # 提取基因序列
        fasta_file, gff_file = self.extractor.extract_by_gene_id(gene_id)
        if not fasta_file:
            return None
        
        # 比对
        result = self.aligner.align(self.query_genome, fasta_file, gene_id)
        if not result:
            return None
        
        result["fasta"] = fasta_file
        result["gff"] = gff_file
        result["id"] = gene_id
        
        # 可视化
        visualizer = GeneIDVisualizer(self.output_dir)
        visualizer.visualize(result)
        
        return result


class LocationProcessor:
    """Location 模式处理器"""
    
    def __init__(self, ref_genome, query_genome, output_dir):
        self.extractor = SequenceExtractor(ref_genome, None, output_dir)
        self.aligner = MummerAligner(output_dir)
        self.query_genome = query_genome
        self.output_dir = output_dir
    
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
        
        # 可视化
        visualizer = LocationVisualizer(self.output_dir)
        visualizer.visualize(result)
        
        return result


class SequenceProcessor:
    """Sequence 模式处理器"""
    
    def __init__(self, ref_genome, output_dir):
        self.aligner = MummerAligner(output_dir)
        self.ref_genome = ref_genome
        self.output_dir = output_dir
    
    def process(self, seq_file):
        print(f"\n{'='*50}")
        print(f"[INFO] 处理序列文件: {seq_file}")
        print(f"{'='*50}")
        
        # 确保 FASTA 格式正确
        seq_id, processed_file = self._prepare_fasta(seq_file)
        
        # 比对（序列作为 query，参考基因组作为 reference）
        result = self.aligner.align(self.ref_genome, processed_file, seq_id)
        if not result:
            return None
        
        result["fasta"] = processed_file
        result["id"] = seq_id
        
        # 可视化
        visualizer = SequenceVisualizer(self.output_dir)
        visualizer.visualize(result)
        
        return result
    
    def _prepare_fasta(self, seq_file):
        """确保 FASTA 文件有 ID 行"""
        with open(seq_file, 'r') as f:
            lines = f.readlines()
        
        # 检查是否有 ID 行
        has_id = any(line.startswith('>') for line in lines)
        
        if has_id:
            seq_id = None
            for line in lines:
                if line.startswith('>'):
                    seq_id = line.strip().lstrip('>').split()[0]
                    break
            out_path = os.path.join(self.output_dir, f"{seq_id}.fasta")
            with open(out_path, 'w') as f:
                f.writelines(lines)
            return seq_id, out_path
        else:
            # 添加默认 ID
            seq_id = "query_seq"
            out_path = os.path.join(self.output_dir, f"{seq_id}.fasta")
            with open(out_path, 'w') as f:
                f.write(f">{seq_id}\n")
                for line in lines:
                    if line.strip():
                        f.write(line)
            return seq_id, out_path


# ======================= 主程序 =======================
def main():
    parser = argparse.ArgumentParser(
        description="GeneScreen - 基因组比对与变异分析工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # Gene ID 模式
  python GeneScreen.py -ref Nippon -qry ZS97 -gid LOC_Os06g10990 -o output/

  # Gene ID 列表模式
  python GeneScreen.py -ref Nippon -qry ZS97 -gidl gene_list.txt -o output/

  # Location 模式
  python GeneScreen.py -ref Nippon -qry ZS97 -loc positions.txt -o output/

  # Sequence 模式
  python GeneScreen.py -ref Nippon -seq query.fasta -o output/
        """
    )
    
    # 基因组参数
    parser.add_argument("-ref", required=True, help="参考基因组（内置名称或路径）")
    parser.add_argument("-qry", help="目标基因组（Gene ID/Location 模式必需）")
    parser.add_argument("-ra", help="参考基因组注释文件（使用自定义路径时需要）")
    
    # 输入模式（互斥）
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-gid", help="单个 Gene ID")
    input_group.add_argument("-gidl", help="Gene ID 列表文件")
    input_group.add_argument("-loc", help="位置文件（格式: chrom\\tstart\\tend）")
    input_group.add_argument("-seq", help="序列文件（FASTA 格式）")
    
    # 输出
    parser.add_argument("-o", "--output", required=True, help="输出目录")
    
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    
    # 获取基因组路径
    ref_genome, ref_annotation = genome_manager.get(args.ref)
    if not ref_genome:
        # 提示用户下载或添加
        print(f"\n基因组 '{args.ref}' 不可用。")
        print("可以使用以下命令管理基因组:")
        print(f"  python GenomeManager.py search {args.ref}  # 搜索 IGV 基因组")
        print(f"  python GenomeManager.py download <id>      # 下载 IGV 基因组")
        print(f"  python GenomeManager.py add <name> <fasta> # 添加自定义基因组")
        return
    
    # 如果提供了 -ra 参数，覆盖注释文件
    if args.ra:
        ref_annotation = args.ra
    
    # 根据模式处理
    if args.seq:
        # Sequence 模式
        processor = SequenceProcessor(ref_genome, args.output)
        processor.process(args.seq)
    
    elif args.gid or args.gidl:
        # Gene ID 模式
        if not args.qry:
            raise ValueError("Gene ID 模式需要指定 -qry 目标基因组")
        if not ref_annotation:
            raise ValueError("Gene ID 模式需要参考基因组注释文件")
        
        qry_genome, _ = genome_manager.get(args.qry)
        if not qry_genome:
            print(f"[ERROR] 目标基因组 '{args.qry}' 不可用")
            return
        processor = GeneIDProcessor(ref_genome, ref_annotation, qry_genome, args.output)
        
        gene_ids = []
        if args.gid:
            gene_ids = [args.gid]
        else:
            with open(args.gidl) as f:
                gene_ids = [line.strip() for line in f if line.strip()]
        
        for gene_id in gene_ids:
            processor.process(gene_id)
    
    elif args.loc:
        # Location 模式
        if not args.qry:
            raise ValueError("Location 模式需要指定 -qry 目标基因组")
        
        qry_genome, _ = genome_manager.get(args.qry)
        if not qry_genome:
            print(f"[ERROR] 目标基因组 '{args.qry}' 不可用")
            return
        processor = LocationProcessor(ref_genome, qry_genome, args.output)
        
        with open(args.loc) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                    name = parts[3] if len(parts) > 3 else None
                    processor.process(chrom, start, end, name)
    
    print(f"\n[INFO] 处理完成，结果保存在: {args.output}")


if __name__ == "__main__":
    main()
