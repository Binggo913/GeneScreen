#!/usr/bin/env python3
"""
GenomeManager - 基因组库管理模块

支持三种来源：
1. 用户自定义基因组（本地路径）
2. IGV 公共基因组（自动下载）
3. 已下载基因组的管理（查看/删除）
"""

import os
import json
import shutil
import urllib.request
from pathlib import Path


# IGV 基因组列表 URL
IGV_GENOMES_URL = "https://igv.org/genomes/genomes.json"

# 默认缓存目录
DEFAULT_CACHE_DIR = Path.home() / ".genescreen" / "genomes"


class GenomeManager:
    """基因组管理器"""
    
    def __init__(self, cache_dir=None):
        self.cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # 本地配置文件
        self.config_file = self.cache_dir / "genomes.json"
        self.config = self._load_config()
        
        # IGV 基因组缓存
        self._igv_genomes = None
    
    # ==================== 配置管理 ====================
    
    def _load_config(self):
        """加载本地配置"""
        if self.config_file.exists():
            with open(self.config_file, 'r') as f:
                return json.load(f)
        return {"custom": {}, "downloaded": {}}
    
    def _save_config(self):
        """保存配置"""
        with open(self.config_file, 'w') as f:
            json.dump(self.config, f, indent=2)
    
    # ==================== 基因组获取 ====================
    
    def get(self, name_or_path):
        """
        获取基因组路径
        
        参数可以是：
        - 自定义基因组名称
        - IGV 基因组 ID
        - 本地文件路径
        
        返回: (genome_path, annotation_path) 或 None
        """
        # 1. 检查是否是本地路径
        if os.path.exists(name_or_path):
            return name_or_path, None
        
        # 2. 检查自定义基因组
        if name_or_path in self.config["custom"]:
            entry = self.config["custom"][name_or_path]
            return entry["genome"], entry.get("annotation")
        
        # 3. 检查已下载的 IGV 基因组
        if name_or_path in self.config["downloaded"]:
            entry = self.config["downloaded"][name_or_path]
            return entry["genome"], entry.get("annotation")
        
        # 4. 尝试从 IGV 下载
        igv_genome = self._find_igv_genome(name_or_path)
        if igv_genome:
            print(f"[INFO] 基因组 '{name_or_path}' 未在本地找到，是否从 IGV 下载？")
            return None, None  # 需要用户确认下载
        
        print(f"[ERROR] 未找到基因组: {name_or_path}")
        return None, None
    
    # ==================== 自定义基因组 ====================
    
    def add_custom(self, name, genome_path, annotation_path=None):
        """添加自定义基因组"""
        if not os.path.exists(genome_path):
            print(f"[ERROR] 基因组文件不存在: {genome_path}")
            return False
        
        if annotation_path and not os.path.exists(annotation_path):
            print(f"[ERROR] 注释文件不存在: {annotation_path}")
            return False
        
        self.config["custom"][name] = {
            "genome": os.path.abspath(genome_path),
            "annotation": os.path.abspath(annotation_path) if annotation_path else None,
            "source": "custom"
        }
        self._save_config()
        print(f"[INFO] 已添加自定义基因组: {name}")
        return True
    
    def remove_custom(self, name):
        """移除自定义基因组（只删除配置，不删除文件）"""
        if name in self.config["custom"]:
            del self.config["custom"][name]
            self._save_config()
            print(f"[INFO] 已移除自定义基因组: {name}")
            return True
        print(f"[ERROR] 未找到自定义基因组: {name}")
        return False
    
    # ==================== IGV 基因组 ====================
    
    def _fetch_igv_genomes(self):
        """获取 IGV 基因组列表"""
        if self._igv_genomes is not None:
            return self._igv_genomes
        
        print("[INFO] 获取 IGV 基因组列表...")
        try:
            with urllib.request.urlopen(IGV_GENOMES_URL, timeout=30) as response:
                self._igv_genomes = json.loads(response.read().decode())
            return self._igv_genomes
        except Exception as e:
            print(f"[ERROR] 获取 IGV 基因组列表失败: {e}")
            return []
    
    def _find_igv_genome(self, genome_id):
        """在 IGV 列表中查找基因组"""
        genomes = self._fetch_igv_genomes()
        for g in genomes:
            if g.get("id") == genome_id:
                return g
        return None
    
    def search_igv(self, keyword):
        """搜索 IGV 基因组"""
        genomes = self._fetch_igv_genomes()
        keyword = keyword.lower()
        results = []
        for g in genomes:
            if keyword in g.get("id", "").lower() or keyword in g.get("name", "").lower():
                results.append({
                    "id": g.get("id"),
                    "name": g.get("name"),
                    "description": g.get("description", "")
                })
        return results
    
    def download_igv(self, genome_id):
        """从 IGV 下载基因组"""
        igv_genome = self._find_igv_genome(genome_id)
        if not igv_genome:
            print(f"[ERROR] IGV 中未找到基因组: {genome_id}")
            return False
        
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
        if not self._download_file(fasta_url, fasta_file):
            return False
        
        # 下载索引（如果有）
        index_url = igv_genome.get("indexURL")
        if index_url:
            index_file = genome_dir / f"{genome_id}.fa.fai"
            print(f"[INFO] 下载索引: {index_url}")
            self._download_file(index_url, index_file)
        
        # 下载注释（如果有）
        annotation_file = None
        tracks = igv_genome.get("tracks", [])
        for track in tracks:
            if track.get("format") in ["gff3", "gff", "gtf"]:
                ann_url = track.get("url")
                if ann_url:
                    ext = track.get("format", "gff3")
                    annotation_file = genome_dir / f"{genome_id}.{ext}"
                    print(f"[INFO] 下载注释: {ann_url}")
                    self._download_file(ann_url, annotation_file)
                    break
        
        # 保存配置
        self.config["downloaded"][genome_id] = {
            "genome": str(fasta_file),
            "annotation": str(annotation_file) if annotation_file else None,
            "source": "igv",
            "igv_id": genome_id
        }
        self._save_config()
        
        print(f"[INFO] 基因组 {genome_id} 下载完成")
        return True
    
    def _download_file(self, url, dest):
        """下载文件"""
        try:
            print(f"  -> {dest}")
            urllib.request.urlretrieve(url, dest)
            return True
        except Exception as e:
            print(f"[ERROR] 下载失败: {e}")
            return False
    
    # ==================== 基因组管理 ====================
    
    def list_all(self):
        """列出所有可用基因组"""
        result = {
            "custom": list(self.config["custom"].keys()),
            "downloaded": list(self.config["downloaded"].keys())
        }
        return result
    
    def remove_downloaded(self, genome_id):
        """删除已下载的基因组（删除文件和配置）"""
        if genome_id not in self.config["downloaded"]:
            print(f"[ERROR] 未找到已下载的基因组: {genome_id}")
            return False
        
        # 删除文件
        genome_dir = self.cache_dir / genome_id
        if genome_dir.exists():
            shutil.rmtree(genome_dir)
            print(f"[INFO] 已删除文件: {genome_dir}")
        
        # 删除配置
        del self.config["downloaded"][genome_id]
        self._save_config()
        
        print(f"[INFO] 已删除基因组: {genome_id}")
        return True
    
    def info(self, name):
        """获取基因组详细信息"""
        if name in self.config["custom"]:
            return {"type": "custom", **self.config["custom"][name]}
        if name in self.config["downloaded"]:
            return {"type": "downloaded", **self.config["downloaded"][name]}
        
        # 检查 IGV
        igv_genome = self._find_igv_genome(name)
        if igv_genome:
            return {"type": "igv_available", **igv_genome}
        
        return None


# ==================== CLI 接口 ====================

def main():
    """命令行接口"""
    import argparse
    
    parser = argparse.ArgumentParser(description="GeneScreen 基因组管理工具")
    subparsers = parser.add_subparsers(dest="command", help="子命令")
    
    # list - 列出基因组
    subparsers.add_parser("list", help="列出所有可用基因组")
    
    # search - 搜索 IGV 基因组
    search_parser = subparsers.add_parser("search", help="搜索 IGV 基因组")
    search_parser.add_argument("keyword", help="搜索关键词")
    
    # download - 下载 IGV 基因组
    download_parser = subparsers.add_parser("download", help="从 IGV 下载基因组")
    download_parser.add_argument("genome_id", help="IGV 基因组 ID")
    
    # add - 添加自定义基因组
    add_parser = subparsers.add_parser("add", help="添加自定义基因组")
    add_parser.add_argument("name", help="基因组名称")
    add_parser.add_argument("genome", help="基因组 FASTA 文件路径")
    add_parser.add_argument("-a", "--annotation", help="注释文件路径")
    
    # remove - 删除基因组
    remove_parser = subparsers.add_parser("remove", help="删除基因组")
    remove_parser.add_argument("name", help="基因组名称")
    remove_parser.add_argument("--keep-files", action="store_true", help="保留文件（仅删除配置）")
    
    # info - 查看基因组信息
    info_parser = subparsers.add_parser("info", help="查看基因组信息")
    info_parser.add_argument("name", help="基因组名称")
    
    args = parser.parse_args()
    manager = GenomeManager()
    
    if args.command == "list":
        genomes = manager.list_all()
        print("\n=== 自定义基因组 ===")
        for name in genomes["custom"]:
            info = manager.config["custom"][name]
            print(f"  {name}: {info['genome']}")
        print(f"\n=== 已下载基因组 ({manager.cache_dir}) ===")
        for name in genomes["downloaded"]:
            info = manager.config["downloaded"][name]
            print(f"  {name}: {info['genome']}")
        if not genomes["custom"] and not genomes["downloaded"]:
            print("  (无)")
    
    elif args.command == "search":
        results = manager.search_igv(args.keyword)
        print(f"\n找到 {len(results)} 个基因组:")
        for r in results[:20]:  # 最多显示 20 个
            print(f"  {r['id']}: {r['name']}")
        if len(results) > 20:
            print(f"  ... 还有 {len(results) - 20} 个结果")
    
    elif args.command == "download":
        manager.download_igv(args.genome_id)
    
    elif args.command == "add":
        manager.add_custom(args.name, args.genome, args.annotation)
    
    elif args.command == "remove":
        if args.name in manager.config["custom"]:
            manager.remove_custom(args.name)
        else:
            manager.remove_downloaded(args.name)
    
    elif args.command == "info":
        info = manager.info(args.name)
        if info:
            print(json.dumps(info, indent=2, ensure_ascii=False))
        else:
            print(f"未找到基因组: {args.name}")
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
