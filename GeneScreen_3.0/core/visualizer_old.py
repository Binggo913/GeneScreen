#!/usr/bin/env python3
"""
GeneScreen 3.0 - 可视化与报告模块

整合自 GeneScreen 2.0 的 GeneScreenVisualizer.py:
- LinkviewVisualizer: 使用 LINKVIEW.py 生成比对可视化图
- AlignmentVisualizer: 内置 SVG 生成器（LINKVIEW 不可用时的回退）
- GeneIDVisualizer / LocationVisualizer / SequenceVisualizer: 三种模式的报告生成

改动（相比 2.0）：
- 移除 CLI 相关代码
- 简化 HTML 模板（保留核心功能）
- 添加类型注解
"""

import os
import json
import base64
import subprocess
import datetime
from pathlib import Path
from typing import Optional, Dict, List, Any, Tuple


class LinkviewVisualizer:
    """LINKVIEW 可视化生成器"""
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        self.colors = {
            'snp': 'orange',
            'indel': 'blue',
        }
    
    def generate_k_file(self, gene_id: str, fasta_file: str, coords_file: str) -> Optional[str]:
        """生成 .k 文件（定义显示区域）"""
        k_file = os.path.join(self.output_dir, f"{gene_id}.k")
        min_aln_len = 100
        
        try:
            # 获取基因序列长度
            gene_len = 0
            with open(fasta_file, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        gene_len += len(line.strip())
            
            # 从 coords 文件解析比对区域
            chr_alignments = {}
            
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
                                aln_len = int(parts[4]) if len(parts) > 4 else abs(ref_end - ref_start) + 1
                                
                                if aln_len < min_aln_len:
                                    continue
                                    
                                ref_name = parts[7].split()[0] if parts[7] else "chr"
                                chr_name = self._extract_chr_name(ref_name)
                                
                                if chr_name not in chr_alignments:
                                    chr_alignments[chr_name] = []
                                chr_alignments[chr_name].append((min(ref_start, ref_end), max(ref_start, ref_end), aln_len))
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
                kf.write(f"{gene_id}:1:{gene_len}\n")
                if main_chr and chr_alignments[main_chr]:
                    alns = chr_alignments[main_chr]
                    min_start = min(a[0] for a in alns)
                    max_end = max(a[1] for a in alns)
                    kf.write(f"{main_chr}:{min_start}:{max_end}\n")
            
            print(f"[INFO] 已生成 .k 文件: {k_file}")
            return k_file
        except Exception as e:
            print(f"[ERROR] 生成 .k 文件失败: {e}")
            return None
    
    def generate_hl_file(self, gene_id: str, snps_file: str) -> Optional[str]:
        """生成 .hl 文件（高亮标记）"""
        hl_file = os.path.join(self.output_dir, f"{gene_id}.hl")
        
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
                                    
                                    color = self.colors['snp'] if var_type == 'SNP' else self.colors['indel']
                                    
                                    ref_name = parts[5].split()[0] if len(parts) > 5 else "ref"
                                    qry_name = parts[6] if len(parts) > 6 else gene_id
                                    chr_name = self._extract_chr_name(ref_name)
                                    
                                    if var_type == 'SNP':
                                        hf.write(f"{qry_name}\t{qry_pos-1}\t{qry_pos}\t{color}\n")
                                        hf.write(f"{chr_name}\t{ref_pos-1}\t{ref_pos}\t{color}\n")
                                    elif var_type == 'INS':
                                        hf.write(f"{qry_name}\t{qry_pos-1}\t{qry_pos}\t{color}\n")
                                    elif var_type == 'DEL':
                                        hf.write(f"{chr_name}\t{ref_pos-1}\t{ref_pos}\t{color}\n")
                                except (ValueError, IndexError):
                                    continue
            
            print(f"[INFO] 已生成 .hl 文件: {hl_file}")
            return hl_file
        except Exception as e:
            print(f"[ERROR] 生成 .hl 文件失败: {e}")
            return None
    
    def _extract_chr_name(self, ref_name: str) -> str:
        """从 ref_name 中提取染色体名称"""
        if not ref_name:
            return "chr"
        return ref_name.split()[0] if ref_name else "chr"
    
    def run_linkview(
        self,
        gene_id: str,
        coords_file: str,
        gff_file: Optional[str],
        k_file: str,
        hl_file: str,
        output_prefix: str,
        fasta_file: str
    ) -> Optional[str]:
        """运行 LINKVIEW.py 生成可视化图"""
        output_svg = f"{output_prefix}.svg"
        output_png = f"{output_prefix}.png"
        
        # 生成 LINKVIEW 专用输入文件
        linkview_input = self._generate_linkview_input(gene_id, coords_file, fasta_file)
        if not linkview_input:
            return None
        
        # 获取 LINKVIEW.py 路径（假设在同目录或 PATH 中）
        script_dir = os.path.dirname(os.path.abspath(__file__))
        linkview_path = os.path.join(script_dir, "..", "..", "GeneScreen_2.0", "bin", "LINKVIEW.py")
        
        if not os.path.exists(linkview_path):
            # 尝试从 PATH 中找
            linkview_path = "LINKVIEW.py"
        
        cmd = [
            "python", linkview_path,
            "-t", "2",
            linkview_input,
            "-k", k_file,
            "-hl", hl_file,
            "-o", output_prefix,
            "--min_identity", "90",
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
            print(f"[WARNING] LINKVIEW.py 未找到")
            return None
        except Exception as e:
            print(f"[WARNING] LINKVIEW 运行异常: {e}")
            return None
    
    def _generate_linkview_input(self, gene_id: str, coords_file: str, fasta_file: str) -> Optional[str]:
        """生成 LINKVIEW nucmer coords 格式的输入文件"""
        linkview_input = os.path.join(self.output_dir, f"{gene_id}.linkview.coords")
        min_aln_len = 100
        
        qry_len = 0
        with open(fasta_file, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    qry_len += len(line.strip())
        
        try:
            chr_alignments = {}
            
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
                            identity = float(parts[6]) if len(parts) > 6 else 99.0
                            
                            if ref_aln_len < min_aln_len:
                                continue
                            
                            ref_full_name = parts[7] if len(parts) > 7 else "chr"
                            ref_name = self._extract_chr_name(ref_full_name)
                            qry_name = parts[8]
                            
                            if ref_name not in chr_alignments:
                                chr_alignments[ref_name] = []
                            chr_alignments[ref_name].append((
                                ref_start, ref_end, qry_start, qry_end,
                                ref_aln_len, qry_aln_len, identity, qry_name, ref_full_name
                            ))
                        except (ValueError, IndexError):
                            continue
            
            main_chr = None
            max_total_len = 0
            for chr_name, alns in chr_alignments.items():
                total_len = sum(a[4] for a in alns)
                if total_len > max_total_len:
                    max_total_len = total_len
                    main_chr = chr_name
            
            with open(linkview_input, 'w') as lf:
                lf.write("ref.fasta query.fasta\n")
                lf.write("NUCMER\n")
                lf.write("\n")
                lf.write("    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n")
                lf.write("=" * 120 + "\n")
                
                if main_chr and chr_alignments[main_chr]:
                    for data in chr_alignments[main_chr]:
                        ref_start, ref_end, qry_start, qry_end, ref_aln_len, qry_aln_len, identity, qry_name, _ = data
                        ref_len = max(ref_start, ref_end) + 1000
                        cov_r = ref_aln_len / ref_len * 100
                        cov_q = qry_aln_len / qry_len * 100
                        
                        lf.write(f"{ref_start:>8} {ref_end:>8}  | {qry_start:>8} {qry_end:>8}  | "
                                f"{ref_aln_len:>8} {qry_aln_len:>8}  | {identity:>8.2f}  | "
                                f"{ref_len:>8} {qry_len:>8}  | {cov_r:>8.2f} {cov_q:>8.2f}  | "
                                f"{main_chr}\t{qry_name}\n")
            
            return linkview_input
        except Exception as e:
            print(f"[ERROR] 生成 LINKVIEW 输入文件失败: {e}")
            return None
    
    def generate_visualization(
        self,
        gene_id: str,
        fasta_file: str,
        coords_file: str,
        snps_file: str,
        gff_file: Optional[str]
    ) -> Optional[str]:
        """完整的可视化流程"""
        k_file = self.generate_k_file(gene_id, fasta_file, coords_file)
        if not k_file:
            return None
        
        hl_file = self.generate_hl_file(gene_id, snps_file)
        if not hl_file:
            return None
        
        output_prefix = os.path.join(self.output_dir, gene_id)
        return self.run_linkview(gene_id, coords_file, gff_file, k_file, hl_file, output_prefix, fasta_file)



class ReportGenerator:
    """HTML 报告生成器"""
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        self.linkview_viz = LinkviewVisualizer(output_dir)
    
    def generate_report(
        self,
        result: Dict[str, Any],
        mode: str,
        ref_name: str = "",
        qry_name: str = ""
    ) -> str:
        """
        生成 HTML 报告
        
        Args:
            result: 分析结果字典
            mode: 分析模式 (gene_id, location, sequence)
            ref_name: 参考基因组名称
            qry_name: 查询基因组名称
        
        Returns:
            报告文件路径
        """
        result_id = result.get("id", "unknown")
        
        # 解析结果文件
        fasta_file = result.get("fasta", "")
        coords_file = result.get("coords", "")
        snps_file = result.get("snps", "")
        gff_file = result.get("gff", "")
        
        # 统计信息
        stats = self._calculate_stats(fasta_file, coords_file, snps_file)
        
        # 尝试生成可视化
        image_path = None
        if fasta_file and coords_file:
            image_path = self.linkview_viz.generate_visualization(
                result_id, fasta_file, coords_file, snps_file, gff_file
            )
        
        # 生成 HTML
        html = self._generate_html(
            result_id, mode, ref_name, qry_name, stats, result, image_path
        )
        
        # 保存报告
        report_file = os.path.join(self.output_dir, f"{result_id}.report.html")
        with open(report_file, "w", encoding="utf-8") as f:
            f.write(html)
        
        print(f"[INFO] 已生成报告: {report_file}")
        return report_file
    
    def _calculate_stats(
        self,
        fasta_file: str,
        coords_file: str,
        snps_file: str
    ) -> Dict[str, Any]:
        """计算统计信息"""
        stats = {
            "seq_length": 0,
            "alignment_count": 0,
            "snp_count": 0,
            "ins_count": 0,
            "del_count": 0,
            "identity": 0.0
        }
        
        # 序列长度
        if fasta_file and os.path.exists(fasta_file):
            with open(fasta_file, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        stats["seq_length"] += len(line.strip())
        
        # 比对统计
        if coords_file and os.path.exists(coords_file):
            identities = []
            with open(coords_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('[') or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 7:
                        stats["alignment_count"] += 1
                        try:
                            identities.append(float(parts[6]))
                        except ValueError:
                            pass
            if identities:
                stats["identity"] = sum(identities) / len(identities)
        
        # 变异统计
        if snps_file and os.path.exists(snps_file):
            with open(snps_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('[') or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        var_type = parts[4].upper()
                        if var_type == 'SNP':
                            stats["snp_count"] += 1
                        elif var_type == 'INS':
                            stats["ins_count"] += 1
                        elif var_type == 'DEL':
                            stats["del_count"] += 1
        
        return stats
    
    def _generate_html(
        self,
        result_id: str,
        mode: str,
        ref_name: str,
        qry_name: str,
        stats: Dict[str, Any],
        result: Dict[str, Any],
        image_path: Optional[str]
    ) -> str:
        """生成 HTML 报告"""
        mode_names = {
            "gene_id": "Gene ID 模式",
            "location": "Location 模式",
            "sequence": "Sequence 模式"
        }
        mode_display = mode_names.get(mode, mode)
        
        # 嵌入图片
        image_html = ""
        if image_path and os.path.exists(image_path):
            if image_path.endswith('.svg'):
                with open(image_path, 'r', encoding='utf-8') as f:
                    image_html = f'<div class="visualization">{f.read()}</div>'
            else:
                with open(image_path, 'rb') as f:
                    img_data = base64.b64encode(f.read()).decode()
                image_html = f'<div class="visualization"><img src="data:image/png;base64,{img_data}" style="max-width:100%;"></div>'
        
        html = f'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GeneScreen 3.0 - {result_id}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #f5f5f5; color: #333; line-height: 1.6; }}
        .container {{ max-width: 1200px; margin: 0 auto; padding: 20px; }}
        .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 12px; margin-bottom: 20px; }}
        .header h1 {{ font-size: 24px; margin-bottom: 10px; }}
        .header .mode {{ background: rgba(255,255,255,0.2); padding: 4px 12px; border-radius: 20px; font-size: 14px; display: inline-block; }}
        .section {{ background: white; border-radius: 12px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }}
        .section h2 {{ font-size: 18px; margin-bottom: 15px; color: #333; border-bottom: 2px solid #667eea; padding-bottom: 8px; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 15px; }}
        .stat-card {{ background: #f8f9fa; padding: 15px; border-radius: 8px; text-align: center; }}
        .stat-card .value {{ font-size: 24px; font-weight: bold; color: #667eea; }}
        .stat-card .label {{ font-size: 12px; color: #666; margin-top: 5px; }}
        .info-table {{ width: 100%; border-collapse: collapse; }}
        .info-table td {{ padding: 10px; border-bottom: 1px solid #eee; }}
        .info-table td:first-child {{ font-weight: 500; color: #666; width: 150px; }}
        .visualization {{ margin: 20px 0; overflow-x: auto; }}
        .visualization svg {{ max-width: 100%; height: auto; }}
        .file-list {{ list-style: none; }}
        .file-list li {{ padding: 10px; background: #f8f9fa; margin-bottom: 8px; border-radius: 6px; font-family: monospace; font-size: 13px; }}
        .footer {{ text-align: center; padding: 20px; color: #999; font-size: 12px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>GeneScreen 3.0 分析报告</h1>
            <span class="mode">{mode_display}</span>
        </div>
        
        <div class="section">
            <h2>分析参数</h2>
            <table class="info-table">
                <tr><td>分析 ID</td><td>{result_id}</td></tr>
                <tr><td>参考基因组</td><td>{ref_name or "N/A"}</td></tr>
                <tr><td>查询基因组</td><td>{qry_name or "N/A"}</td></tr>
                <tr><td>分析时间</td><td>{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</td></tr>
            </table>
        </div>
        
        <div class="section">
            <h2>结果统计</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="value">{stats["seq_length"]:,}</div>
                    <div class="label">序列长度 (bp)</div>
                </div>
                <div class="stat-card">
                    <div class="value">{stats["identity"]:.2f}%</div>
                    <div class="label">平均 Identity</div>
                </div>
                <div class="stat-card">
                    <div class="value">{stats["snp_count"]}</div>
                    <div class="label">SNP 数量</div>
                </div>
                <div class="stat-card">
                    <div class="value">{stats["ins_count"]}</div>
                    <div class="label">插入 (INS)</div>
                </div>
                <div class="stat-card">
                    <div class="value">{stats["del_count"]}</div>
                    <div class="label">缺失 (DEL)</div>
                </div>
            </div>
        </div>
        
        <div class="section">
            <h2>比对可视化</h2>
            {image_html if image_html else '<p style="color:#999;">可视化图片生成失败或 LINKVIEW 不可用</p>'}
        </div>
        
        <div class="section">
            <h2>输出文件</h2>
            <ul class="file-list">
                <li>{result.get("fasta", "N/A")}</li>
                <li>{result.get("coords", "N/A")}</li>
                <li>{result.get("snps", "N/A")}</li>
                <li>{result.get("blast_xml", "N/A")}</li>
            </ul>
        </div>
        
        <div class="footer">
            <p>Generated by GeneScreen 3.0 | {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        </div>
    </div>
</body>
</html>'''
        
        return html


# 便捷函数
def generate_report(
    result: Dict[str, Any],
    output_dir: str,
    mode: str = "gene_id",
    ref_name: str = "",
    qry_name: str = ""
) -> str:
    """生成分析报告的便捷函数"""
    generator = ReportGenerator(output_dir)
    return generator.generate_report(result, mode, ref_name, qry_name)
