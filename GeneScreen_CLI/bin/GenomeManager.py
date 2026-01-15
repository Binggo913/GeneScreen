#!/usr/bin/env python3
"""
GenomeManager - 基因组库管理模块 (SQLite 版本)

支持多种来源：
1. 用户自定义基因组（本地路径）
2. IGV 公共基因组
3. Ensembl Plants 植物基因组
4. 已下载基因组的管理（查看/删除）

改动（相比原版）：
- 使用 SQLite 数据库替代 JSON 配置文件
- 通过 database.py 的 Database 类管理数据
"""

import os
import json
import shutil
import ssl
import gzip
import urllib.request
from pathlib import Path
from typing import Optional, List, Dict, Any

from database import get_database, Database


# 数据源 URL
IGV_GENOMES_URL = "https://igv.org/genomes/genomes.json"
ENSEMBL_PLANTS_URL = "https://rest.ensembl.org/info/species?division=EnsemblPlants"
ENSEMBL_FTP_BASE = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current"

# 默认缓存目录
DEFAULT_CACHE_DIR = Path.home() / ".genescreen" / "genomes"

# SSL 上下文（跳过证书验证，解决 Windows 兼容问题）
SSL_CONTEXT = ssl.create_default_context()
SSL_CONTEXT.check_hostname = False
SSL_CONTEXT.verify_mode = ssl.CERT_NONE


class GenomeManager:
    """基因组管理器"""

    def __init__(self, cache_dir: Optional[Path] = None, db: Optional[Database] = None):
        self.cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # 使用 SQLite 数据库
        self.db = db or get_database()

        # 缓存（在线数据源）
        self._igv_genomes = None
        self._ensembl_species = None

    # ==================== 基因组获取 ====================

    def get(self, name_or_path: str, annotation_source: Optional[str] = None):
        """
        获取基因组路径

        参数可以是：
        - 数据库中的基因组名称
        - 本地文件路径

        Args:
            name_or_path: 基因组名称或路径
            annotation_source: 注释版本来源（可选，默认使用 version1 或第一个）

        返回: (genome_path, annotation_path) 或 (None, None)
        """
        # 1. 优先从数据库查找（避免当前目录同名文件/目录干扰）
        genome = self.db.get_genome(name_or_path)
        if genome:
            # 获取注释版本
            if annotation_source:
                ann = self.db.get_annotation_by_source(genome["id"], annotation_source)
                if not ann:
                    print(f"[WARNING] 未找到注释版本 '{annotation_source}'，使用默认注释")
                    ann = self.db.get_default_annotation(genome["id"])
            else:
                ann = self.db.get_default_annotation(genome["id"])
            
            annotation_path = ann["annotation_path"] if ann else genome.get("annotation_path")
            return genome["fasta_path"], annotation_path

        # 2. 数据库没找到，检查是否是本地文件路径
        if os.path.exists(name_or_path):
            return name_or_path, None

        print(f"[ERROR] 未找到基因组: {name_or_path}")
        print("  使用 'python GenomeManager.py search <keyword>' 搜索可用基因组")
        return None, None

    # ==================== 自定义基因组 ====================

    def add_custom(
        self,
        name: str,
        genome_path: str,
        annotation_path: Optional[str] = None
    ) -> bool:
        """添加自定义基因组，支持 .gz 压缩文件"""
        if not os.path.exists(genome_path):
            print(f"[ERROR] 基因组文件不存在: {genome_path}")
            return False

        if annotation_path and not os.path.exists(annotation_path):
            print(f"[ERROR] 注释文件不存在: {annotation_path}")
            return False

        # 检查是否已存在
        if self.db.genome_exists(name):
            print(f"[ERROR] 基因组 '{name}' 已存在")
            return False

        # 处理 .gz 文件
        final_genome_path = genome_path
        final_annotation_path = annotation_path

        if genome_path.endswith(".gz"):
            genome_dir = self.cache_dir / name
            genome_dir.mkdir(parents=True, exist_ok=True)
            final_genome_path = str(genome_dir / f"{name}.fa")
            print(f"[INFO] 解压基因组: {genome_path}")
            self._gunzip_copy(genome_path, final_genome_path)

        if annotation_path and annotation_path.endswith(".gz"):
            genome_dir = self.cache_dir / name
            genome_dir.mkdir(parents=True, exist_ok=True)
            ext = ".gff3" if "gff" in annotation_path.lower() else ".gtf"
            final_annotation_path = str(genome_dir / f"{name}{ext}")
            print(f"[INFO] 解压注释: {annotation_path}")
            self._gunzip_copy(annotation_path, final_annotation_path)

        # 检查并创建索引
        self._ensure_fasta_index(final_genome_path)

        # 保存到数据库
        self.db.add_genome(
            name=name,
            fasta_path=os.path.abspath(final_genome_path),
            annotation_path=os.path.abspath(final_annotation_path) if final_annotation_path else None,
            source="custom"
        )
        print(f"[INFO] 已添加自定义基因组: {name}")
        return True

    def add_annotation(self, genome_name: str, annotation_path: str, source: Optional[str] = None) -> bool:
        """为已有基因组添加注释版本"""
        genome = self.db.get_genome(genome_name)
        if not genome:
            print(f"[ERROR] 未找到基因组: {genome_name}")
            return False

        if not os.path.exists(annotation_path):
            print(f"[ERROR] 注释文件不存在: {annotation_path}")
            return False

        # 解析 source（从文件名或使用默认值）
        if not source:
            source = self._parse_annotation_source(Path(annotation_path).name) or "version1"

        result = self.db.add_annotation(genome["id"], source, os.path.abspath(annotation_path))
        if result == -1:
            print(f"[ERROR] 注释版本 '{source}' 已存在")
            return False

        print(f"[INFO] 已添加注释版本: {genome_name} -> {source}")
        return True

    def list_annotations(self, genome_name: str) -> List[Dict[str, Any]]:
        """列出基因组的所有注释版本"""
        genome = self.db.get_genome(genome_name)
        if not genome:
            return []
        return self.db.get_annotations(genome["id"])

    @staticmethod
    def _parse_annotation_source(filename: str) -> Optional[str]:
        """从文件名解析注释来源，格式: {source}.{species}.{ext}"""
        parts = filename.split(".")
        if len(parts) >= 3:
            # 检查最后一个是否是扩展名
            ext = parts[-1].lower()
            if ext in ("gff", "gff3", "gtf", "gz"):
                if ext == "gz" and len(parts) >= 4:
                    # xxx.gff3.gz
                    return parts[0] if parts[0] not in ("gff", "gff3", "gtf") else None
                return parts[0] if parts[0] not in ("gff", "gff3", "gtf") else None
        return None

    def _gunzip_copy(self, gz_file: str, out_file: str):
        """解压 gzip 文件（保留原文件）"""
        with gzip.open(gz_file, "rb") as f_in:
            with open(out_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    def _gunzip(self, gz_file: str, out_file: str):
        """解压 gzip 文件（删除原文件）"""
        self._gunzip_copy(gz_file, out_file)
        os.remove(gz_file)

    def remove_custom(self, name: str) -> bool:
        """移除自定义基因组（只删除配置，不删除文件）"""
        genome = self.db.get_genome(name)
        if not genome:
            print(f"[ERROR] 未找到基因组: {name}")
            return False
        
        if genome.get("source") != "custom":
            print(f"[ERROR] '{name}' 不是自定义基因组，请使用 remove_downloaded()")
            return False

        self.db.delete_genome(name)
        print(f"[INFO] 已移除自定义基因组: {name}")
        return True

    # ==================== Ensembl Plants ====================

    def _fetch_ensembl_species(self) -> List[Dict]:
        """获取 Ensembl Plants 物种列表"""
        if self._ensembl_species is not None:
            return self._ensembl_species

        print("[INFO] 获取 Ensembl Plants 物种列表...")
        try:
            req = urllib.request.Request(
                ENSEMBL_PLANTS_URL,
                headers={"Accept": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=30, context=SSL_CONTEXT) as response:
                data = json.loads(response.read().decode())
                self._ensembl_species = data.get("species", [])
            return self._ensembl_species
        except Exception as e:
            print(f"[ERROR] 获取 Ensembl Plants 列表失败: {e}")
            return []

    def search_ensembl(self, keyword: str) -> List[Dict]:
        """搜索 Ensembl Plants 物种"""
        species_list = self._fetch_ensembl_species()
        keyword = keyword.lower()
        results = []
        for sp in species_list:
            name = sp.get("name") or ""
            display_name = sp.get("display_name") or ""
            common_name = sp.get("common_name") or ""
            if (
                keyword in name.lower()
                or keyword in display_name.lower()
                or keyword in common_name.lower()
            ):
                results.append({
                    "id": name,
                    "display_name": display_name,
                    "common_name": common_name,
                    "assembly": sp.get("assembly") or "",
                })
        return results

    def download_ensembl(self, species_name: str) -> bool:
        """从 Ensembl Plants 下载基因组"""
        species_list = self._fetch_ensembl_species()

        # 查找物种
        species_info = None
        for sp in species_list:
            if sp.get("name") == species_name:
                species_info = sp
                break

        if not species_info:
            print(f"[ERROR] Ensembl Plants 中未找到: {species_name}")
            return False

        # 检查是否已存在
        if self.db.genome_exists(species_name):
            print(f"[INFO] 基因组 '{species_name}' 已存在")
            return True

        # 创建目录
        genome_dir = self.cache_dir / species_name
        genome_dir.mkdir(parents=True, exist_ok=True)

        assembly = species_info.get("assembly") or ""
        release = species_info.get("release") or "59"
        species_cap = species_name[0].upper() + species_name[1:]

        # 下载 FASTA
        fasta_filename = f"{species_cap}.{assembly}.dna.toplevel.fa.gz"
        fasta_url = f"{ENSEMBL_FTP_BASE}/fasta/{species_name}/dna/{fasta_filename}"
        fasta_gz = genome_dir / fasta_filename
        fasta_file = genome_dir / f"{species_name}.fa"

        print(f"[INFO] 下载 FASTA: {fasta_url}")
        if not self._download_file(fasta_url, str(fasta_gz)):
            # 尝试 dna_sm 格式
            fasta_filename_alt = f"{species_cap}.{assembly}.dna_sm.toplevel.fa.gz"
            fasta_url_alt = f"{ENSEMBL_FTP_BASE}/fasta/{species_name}/dna/{fasta_filename_alt}"
            print(f"[INFO] 尝试备用 URL: {fasta_url_alt}")
            fasta_gz = genome_dir / fasta_filename_alt
            if not self._download_file(fasta_url_alt, str(fasta_gz)):
                return False

        # 解压 FASTA
        print("[INFO] 解压 FASTA...")
        self._gunzip(str(fasta_gz), str(fasta_file))

        # 创建索引
        self._ensure_fasta_index(str(fasta_file))

        # 下载 GFF3（使用 {source}.{species}.{ext} 命名）
        annotation_file = None
        for rel in [release, "62", "61", "60", "59", "58", "57"]:
            gff_filename = f"{species_cap}.{assembly}.{rel}.gff3.gz"
            gff_url = f"{ENSEMBL_FTP_BASE}/gff3/{species_name}/{gff_filename}"
            gff_gz = genome_dir / gff_filename
            gff_file = genome_dir / f"ensembl_plants.{species_name}.gff3"  # 新命名格式

            print(f"[INFO] 下载 GFF3: {gff_url}")
            if self._download_file(gff_url, str(gff_gz)):
                print("[INFO] 解压 GFF3...")
                self._gunzip(str(gff_gz), str(gff_file))
                annotation_file = gff_file
                break

            # 尝试 abinitio 版本
            gff_filename_ab = f"{species_cap}.{assembly}.{rel}.abinitio.gff3.gz"
            gff_url_ab = f"{ENSEMBL_FTP_BASE}/gff3/{species_name}/{gff_filename_ab}"
            gff_gz_ab = genome_dir / gff_filename_ab

            if self._download_file(gff_url_ab, str(gff_gz_ab)):
                print("[INFO] 解压 GFF3 (abinitio)...")
                self._gunzip(str(gff_gz_ab), str(gff_file))
                annotation_file = gff_file
                break

        if not annotation_file:
            print("[WARN] GFF3 注释文件下载失败，基因组仍可使用")

        # 保存到数据库
        genome_id = self.db.add_genome(
            name=species_name,
            fasta_path=str(fasta_file),
            annotation_path=str(annotation_file) if annotation_file else None,
            display_name=species_info.get("display_name") or "",
            source="ensembl_plants",
            assembly=assembly,
            species=species_info.get("display_name") or species_name
        )
        
        # 添加到注释版本表
        if annotation_file:
            self.db.add_annotation(genome_id, "ensembl_plants", str(annotation_file))

        print(f"[INFO] 基因组 {species_name} 下载完成")
        return True

    # ==================== IGV 基因组 ====================

    def _fetch_igv_genomes(self) -> List[Dict]:
        """获取 IGV 基因组列表"""
        if self._igv_genomes is not None:
            return self._igv_genomes

        print("[INFO] 获取 IGV 基因组列表...")
        try:
            with urllib.request.urlopen(IGV_GENOMES_URL, timeout=30, context=SSL_CONTEXT) as response:
                self._igv_genomes = json.loads(response.read().decode())
            return self._igv_genomes
        except Exception as e:
            print(f"[ERROR] 获取 IGV 基因组列表失败: {e}")
            return []

    def _find_igv_genome(self, genome_id: str) -> Optional[Dict]:
        """在 IGV 列表中查找基因组"""
        genomes = self._fetch_igv_genomes()
        for g in genomes:
            if g.get("id") == genome_id:
                return g
        return None

    def search_igv(self, keyword: str) -> List[Dict]:
        """搜索 IGV 基因组"""
        genomes = self._fetch_igv_genomes()
        keyword = keyword.lower()
        results = []
        for g in genomes:
            if keyword in g.get("id", "").lower() or keyword in g.get("name", "").lower():
                results.append({
                    "id": g.get("id"),
                    "name": g.get("name"),
                    "description": g.get("description", ""),
                })
        return results

    def download_igv(self, genome_id: str) -> bool:
        """从 IGV 下载基因组"""
        igv_genome = self._find_igv_genome(genome_id)
        if not igv_genome:
            print(f"[ERROR] IGV 中未找到基因组: {genome_id}")
            return False

        # 检查是否已存在
        if self.db.genome_exists(genome_id):
            print(f"[INFO] 基因组 '{genome_id}' 已存在")
            return True

        # 创建基因组目录
        genome_dir = self.cache_dir / genome_id
        genome_dir.mkdir(parents=True, exist_ok=True)

        # 下载 FASTA
        fasta_url = igv_genome.get("fastaURL")
        if not fasta_url:
            print(f"[ERROR] 基因组 {genome_id} 没有 FASTA URL")
            return False

        fasta_file = genome_dir / f"{genome_id}.fa"
        print(f"[INFO] 下载 FASTA: {fasta_url}")
        if not self._download_file(fasta_url, str(fasta_file)):
            return False

        # 创建索引
        self._ensure_fasta_index(str(fasta_file))

        # 下载注释（如果有，使用 {source}.{species}.{ext} 命名）
        annotation_file = None
        tracks = igv_genome.get("tracks", [])
        for track in tracks:
            if track.get("format") in ["gff3", "gff", "gtf"]:
                ann_url = track.get("url")
                if ann_url:
                    ext = track.get("format", "gff3")
                    annotation_file = genome_dir / f"igv.{genome_id}.{ext}"  # 新命名格式
                    print(f"[INFO] 下载注释: {ann_url}")
                    self._download_file(ann_url, str(annotation_file))
                    break

        # 保存到数据库
        db_genome_id = self.db.add_genome(
            name=genome_id,
            fasta_path=str(fasta_file),
            annotation_path=str(annotation_file) if annotation_file else None,
            display_name=igv_genome.get("name") or genome_id,
            source="igv",
            description=igv_genome.get("description")
        )
        
        # 添加到注释版本表
        if annotation_file:
            self.db.add_annotation(db_genome_id, "igv", str(annotation_file))

        print(f"[INFO] 基因组 {genome_id} 下载完成")
        return True

    # ==================== 通用方法 ====================

    def _download_file(self, url: str, dest: str) -> bool:
        """下载文件"""
        try:
            print(f"  -> {dest}")
            req = urllib.request.Request(url, headers={
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            })
            with urllib.request.urlopen(req, timeout=300, context=SSL_CONTEXT) as response:
                with open(dest, 'wb') as f:
                    shutil.copyfileobj(response, f)
            return True
        except Exception as e:
            print(f"[ERROR] 下载失败: {e}")
            return False

    def _ensure_fasta_index(self, fasta_path: str):
        """确保 FASTA 文件有索引"""
        fai_path = f"{fasta_path}.fai"
        if not os.path.exists(fai_path):
            print(f"[INFO] 创建 FASTA 索引: {fasta_path}")
            from pyfaidx import Fasta
            Fasta(fasta_path)

    # ==================== 基因组管理 ====================

    def list_all(self) -> Dict[str, List[str]]:
        """列出所有可用基因组（按来源分组）"""
        genomes = self.db.list_genomes()
        result = {"custom": [], "downloaded": []}
        
        for g in genomes:
            source = g.get("source", "custom")
            if source == "custom":
                result["custom"].append(g["name"])
            else:
                result["downloaded"].append(g["name"])
        
        return result

    def list_genomes(self) -> List[Dict[str, Any]]:
        """列出所有基因组（完整信息）"""
        return self.db.list_genomes()

    def remove_downloaded(self, genome_id: str) -> bool:
        """删除已下载的基因组（删除文件和配置）"""
        genome = self.db.get_genome(genome_id)
        if not genome:
            print(f"[ERROR] 未找到基因组: {genome_id}")
            return False

        if genome.get("source") == "custom":
            print(f"[ERROR] '{genome_id}' 是自定义基因组，请使用 remove_custom()")
            return False

        # 删除文件
        genome_dir = self.cache_dir / genome_id
        if genome_dir.exists():
            shutil.rmtree(genome_dir)
            print(f"[INFO] 已删除文件: {genome_dir}")

        # 删除数据库记录
        self.db.delete_genome(genome_id)
        print(f"[INFO] 已删除基因组: {genome_id}")
        return True

    def info(self, name: str) -> Optional[Dict[str, Any]]:
        """获取基因组详细信息"""
        genome = self.db.get_genome(name)
        if genome:
            return genome

        # 检查 IGV（在线）
        igv_genome = self._find_igv_genome(name)
        if igv_genome:
            return {"type": "igv_available", **igv_genome}

        return None

    def search(self, keyword: str) -> Dict[str, List[Dict]]:
        """搜索基因组（IGV + Ensembl Plants）"""
        return {
            "igv": self.search_igv(keyword),
            "ensembl": self.search_ensembl(keyword)
        }

    def download(self, genome_id: str) -> bool:
        """下载基因组（自动识别来源）"""
        # 先检查 IGV
        if self._find_igv_genome(genome_id):
            return self.download_igv(genome_id)
        
        # 再检查 Ensembl Plants
        species_list = self._fetch_ensembl_species()
        found = any(sp.get("name") == genome_id for sp in species_list)
        if found:
            return self.download_ensembl(genome_id)
        
        print(f"[ERROR] 未找到基因组: {genome_id}")
        return False


# ==================== CLI 接口 ====================


def main():
    """命令行接口"""
    import argparse

    parser = argparse.ArgumentParser(description="GeneScreen 基因组管理工具 (SQLite 版本)")
    subparsers = parser.add_subparsers(dest="command", help="子命令")

    # list - 列出基因组
    subparsers.add_parser("list", help="列出所有可用基因组")

    # search - 搜索基因组
    search_parser = subparsers.add_parser("search", help="搜索基因组（IGV + Ensembl Plants）")
    search_parser.add_argument("keyword", help="搜索关键词")

    # download - 下载基因组
    download_parser = subparsers.add_parser("download", help="下载基因组")
    download_parser.add_argument("genome_id", help="基因组 ID")

    # add - 添加自定义基因组
    add_parser = subparsers.add_parser("add", help="添加自定义基因组")
    add_parser.add_argument("name", help="基因组名称")
    add_parser.add_argument("genome", help="基因组 FASTA 文件路径")
    add_parser.add_argument("-a", "--annotation", help="注释文件路径")

    # remove - 删除基因组
    remove_parser = subparsers.add_parser("remove", help="删除基因组")
    remove_parser.add_argument("name", help="基因组名称")

    # info - 查看基因组信息
    info_parser = subparsers.add_parser("info", help="查看基因组信息")
    info_parser.add_argument("name", help="基因组名称")

    # add-annotation - 添加注释版本
    add_ann_parser = subparsers.add_parser("add-annotation", help="为基因组添加注释版本")
    add_ann_parser.add_argument("genome", help="基因组名称")
    add_ann_parser.add_argument("annotation", help="注释文件路径")
    add_ann_parser.add_argument("-s", "--source", help="注释来源标识（默认从文件名解析或 version1）")

    args = parser.parse_args()
    manager = GenomeManager()

    if args.command == "list":
        all_genomes = manager.list_genomes()
        
        # 按来源分组显示
        custom_genomes = [g for g in all_genomes if g.get("source") == "custom"]
        downloaded_genomes = [g for g in all_genomes if g.get("source") != "custom"]
        
        def _format_annotations(genome):
            """格式化注释版本信息"""
            anns = manager.list_annotations(genome["name"])
            if anns:
                sources = [a["source"] for a in anns]
                return f" [{', '.join(sources)}]"
            elif genome.get("annotation_path"):
                return " [version1]"
            return ""
        
        print("\n=== 自定义基因组 ===")
        if custom_genomes:
            for g in custom_genomes:
                ann_info = _format_annotations(g)
                print(f"  {g['name']}: {g['fasta_path']}{ann_info}")
        else:
            print("  (无)")
            
        print(f"\n=== 已下载基因组 ({manager.cache_dir}) ===")
        if downloaded_genomes:
            for g in downloaded_genomes:
                source = g.get("source", "unknown")
                ann_info = _format_annotations(g)
                print(f"  {g['name']} [{source}]: {g['fasta_path']}{ann_info}")
        else:
            print("  (无)")

    elif args.command == "search":
        results = manager.search(args.keyword)
        
        if results["igv"]:
            print(f"\n[IGV] 找到 {len(results['igv'])} 个基因组:")
            for r in results["igv"][:10]:
                print(f"  {r['id']}: {r['name']}")
            if len(results["igv"]) > 10:
                print(f"  ... 还有 {len(results['igv']) - 10} 个结果")

        if results["ensembl"]:
            print(f"\n[Ensembl Plants] 找到 {len(results['ensembl'])} 个物种:")
            for r in results["ensembl"][:20]:
                print(f"  {r['id']}: {r['display_name']} ({r['common_name']}) - {r['assembly']}")
            if len(results["ensembl"]) > 20:
                print(f"  ... 还有 {len(results['ensembl']) - 20} 个结果")

        if not results["igv"] and not results["ensembl"]:
            print(f"\n未找到匹配 '{args.keyword}' 的基因组")
        else:
            print(f"\n下载示例: python GenomeManager.py download <genome_id>")

    elif args.command == "download":
        manager.download(args.genome_id)

    elif args.command == "add":
        manager.add_custom(args.name, args.genome, args.annotation)

    elif args.command == "remove":
        genome = manager.db.get_genome(args.name)
        if genome:
            if genome.get("source") == "custom":
                manager.remove_custom(args.name)
            else:
                manager.remove_downloaded(args.name)
        else:
            print(f"[ERROR] 未找到基因组: {args.name}")

    elif args.command == "info":
        info = manager.info(args.name)
        if info:
            # 添加注释版本信息
            anns = manager.list_annotations(args.name)
            if anns:
                info["annotations"] = anns
            print(json.dumps(info, indent=2, ensure_ascii=False, default=str))
        else:
            print(f"未找到基因组: {args.name}")

    elif args.command == "add-annotation":
        manager.add_annotation(args.genome, args.annotation, args.source)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
