#!/usr/bin/env python3
"""
GeneScreen 1.0 - 基因组库管理模块

支持多种来源：
1. 用户自定义基因组（本地路径）
2. IGV 公共基因组
3. Ensembl Plants 植物基因组
4. 已下载基因组的管理（查看/删除）

改动（相比 CLI 版本）：
- 使用 SQLite 数据库替代 JSON 配置文件
- 通过 core/database.py 的 Database 类管理数据
"""

import os
import shutil
import ssl
import json
import gzip
import urllib.request
import threading
import queue
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any

from .database import get_database, Database
from .config import get_db_dir
from .gene_id_utils import (
    GENE_ID_LIST_SUFFIX,
    LEGACY_LIST_SUFFIX,
    extract_gene_ids_from_annotation,
    read_gene_id_list,
    write_gene_id_list,
    count_gene_id_list,
)


# 数据源 URL
IGV_GENOMES_URL = "https://igv.org/genomes/genomes.json"
ENSEMBL_PLANTS_URL = "https://rest.ensembl.org/info/species?division=EnsemblPlants"
ENSEMBL_FTP_BASE = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current"

# 文件扩展名
FASTA_EXTS = ('.fa', '.fasta', '.fa.gz', '.fasta.gz')
ANN_EXTS = ('.gff', '.gff3', '.gtf', '.gff.gz', '.gff3.gz', '.gtf.gz')
IGNORE_MARKER = ".genescreen_ignore"

# SSL 上下文（跳过证书验证，解决 Windows 兼容问题）
SSL_CONTEXT = ssl.create_default_context()
SSL_CONTEXT.check_hostname = False
SSL_CONTEXT.verify_mode = ssl.CERT_NONE


class GenomeManager:
    """基因组管理器"""

    def __init__(self, cache_dir: Optional[Path] = None, db: Optional[Database] = None):
        """
        初始化基因组管理器
        
        Args:
            cache_dir: 基因组文件缓存目录
            db: 数据库实例，默认使用全局单例
        """
        self.cache_dir = Path(cache_dir) if cache_dir else Path(get_db_dir())
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # 使用 SQLite 数据库
        self.db = db or get_database()

        # 缓存（在线数据源）
        self._igv_genomes = None
        self._ensembl_species = None
        self._syncing = False

        # Gene ID 列表异步队列
        self._gene_id_queue: "queue.Queue[Tuple[str, str]]" = queue.Queue()
        self._gene_id_jobs = set()
        self._gene_id_lock = threading.Lock()
        self._gene_id_worker = threading.Thread(target=self._gene_id_worker_loop, daemon=True)
        self._gene_id_worker.start()

        self.sync_cache_dir()

    @staticmethod
    def _strip_known_suffix(filename: str, suffixes: Tuple[str, ...]) -> str:
        lower = filename.lower()
        for suffix in sorted(suffixes, key=len, reverse=True):
            if lower.endswith(suffix):
                return filename[:-len(suffix)]
        return Path(filename).stem

    @staticmethod
    def _parse_annotation_name(filename: str) -> Tuple[str, str]:
        """
        解析注释文件名，提取 source 和 species
        
        格式：{source}.{species}.{ext} 或 {species}.{ext}
        返回：(source, species)
        """
        # 去掉扩展名
        base = filename
        for ext in sorted(ANN_EXTS, key=len, reverse=True):
            if base.lower().endswith(ext):
                base = base[:-len(ext)]
                break
        
        parts = base.split('.', 1)
        if len(parts) == 2:
            # 有 source 前缀：igv.ZS97 -> ('igv', 'ZS97')
            return parts[0], parts[1]
        # 无前缀：ZS97 -> ('version1', 'ZS97')
        return 'version1', parts[0]

    def scan_local_genomes(self, root: Path) -> List[Dict[str, Any]]:
        """
        扫描本地基因组，支持多注释版本
        
        返回格式：[{name, fasta, annotations: [{source, path}]}]
        """
        root = Path(root)
        if not root.exists():
            return []

        # 按目录收集注释文件：{dir: {species_lower: [(source, path)]}}
        ann_by_dir: Dict[Path, Dict[str, List[Tuple[str, Path]]]] = {}
        fasta_files: List[Path] = []

        for f in root.rglob("*"):
            if not f.is_file():
                continue
            # 跳过正在下载的目录
            if (f.parent / ".downloading").exists():
                continue
            if self._is_ignored_cache_path(f, root):
                continue
            name_lower = f.name.lower()
            if any(name_lower.endswith(ext) for ext in FASTA_EXTS):
                fasta_files.append(f)
            elif any(name_lower.endswith(ext) for ext in ANN_EXTS):
                source, species = self._parse_annotation_name(f.name)
                ann_by_dir.setdefault(f.parent, {}).setdefault(species.lower(), []).append((source, f))

        genome_map: Dict[Tuple[Path, str], Dict[str, Any]] = {}
        for fasta in fasta_files:
            base = self._strip_known_suffix(fasta.name, FASTA_EXTS)
            key = (fasta.parent, base.lower())
            fasta_path = str(fasta)
            
            # 获取该目录下匹配的所有注释
            ann_list = ann_by_dir.get(fasta.parent, {}).get(base.lower(), [])
            annotations = [{"source": src, "path": str(p)} for src, p in ann_list]
            
            existing = genome_map.get(key)
            if existing:
                # 优先非压缩 FASTA
                if existing["fasta"].lower().endswith(".gz") and not fasta_path.lower().endswith(".gz"):
                    genome_map[key]["fasta"] = fasta_path
                continue
            
            genome_map[key] = {
                "name": base,
                "fasta": fasta_path,
                "annotations": annotations
            }

        return list(genome_map.values())

    def _is_ignored_cache_path(self, path: Path, root: Path) -> bool:
        """Return True when a cached genome directory was explicitly removed from the library."""
        root = root.resolve()
        current = path.parent.resolve()
        while True:
            if (current / IGNORE_MARKER).exists():
                return True
            if current == root:
                return False
            if root not in current.parents:
                return False
            current = current.parent

    def mark_cache_dir_ignored(self, genome: Dict[str, Any]) -> None:
        """Prevent sync_cache_dir from re-importing a cached genome whose files remain on disk."""
        paths = [
            genome.get("fasta_path") or "",
            genome.get("annotation_path") or "",
            genome.get("gene_ids_path") or "",
        ]
        for raw_path in paths:
            if not raw_path:
                continue
            path = Path(raw_path)
            parent = path if path.is_dir() else path.parent
            try:
                resolved = parent.resolve()
                if resolved != self.cache_dir.resolve() and self.cache_dir.resolve() in resolved.parents:
                    (resolved / IGNORE_MARKER).write_text("Removed from GeneScreen genome library.\n", encoding="utf-8")
                    return
            except OSError:
                continue

    def sync_cache_dir(self) -> Tuple[int, int]:
        """同步缓存目录，支持多注释版本"""
        if self._syncing:
            return 0, 0
        self._syncing = True
        genomes = self.scan_local_genomes(self.cache_dir)
        added = 0
        updated = 0

        try:
            for genome in genomes:
                name = genome["name"]
                fasta_path = os.path.abspath(genome["fasta"])
                annotations = genome.get("annotations", [])
                existing = self.db.get_genome(name)
                
                if existing:
                    genome_id = existing["id"]
                    updates = {}
                    existing_fasta = existing.get("fasta_path")
                    preferred_fasta = existing_fasta or fasta_path
                    
                    if existing_fasta:
                        existing_is_gz = existing_fasta.lower().endswith(".gz")
                        new_is_gz = fasta_path.lower().endswith(".gz")
                        if not os.path.exists(existing_fasta):
                            preferred_fasta = fasta_path
                        elif existing_is_gz and not new_is_gz:
                            preferred_fasta = fasta_path
                        elif not existing_is_gz and new_is_gz:
                            preferred_fasta = existing_fasta
                        else:
                            preferred_fasta = fasta_path

                    if existing_fasta != preferred_fasta:
                        updates["fasta_path"] = preferred_fasta
                    
                    # 同步注释版本到 genome_annotations 表
                    for ann in annotations:
                        source = ann["source"]
                        ann_path = os.path.abspath(ann["path"])
                        existing_ann = self.db.get_annotation_by_source(genome_id, source)
                        if not existing_ann:
                            self.db.add_annotation(genome_id, source, ann_path)
                            updated += 1
                        self.enqueue_gene_id_list(name, ann_path, source)
                    
                    if updates:
                        effective_fasta = updates.get("fasta_path") or existing_fasta or fasta_path
                        if effective_fasta:
                            self._ensure_fasta_index(effective_fasta)
                        self.db.update_genome(name, **updates)
                        updated += 1
                    continue

                # 新增基因组
                # 取第一个注释作为默认（兼容旧逻辑）
                default_ann = annotations[0]["path"] if annotations else None
                if self._add_genome_with_annotations(name, fasta_path, annotations):
                    added += 1
        finally:
            self._syncing = False

        if added or updated:
            print(f"[INFO] 基因组库同步完成: 新增 {added} 个, 更新 {updated} 个")
        return added, updated

    def _add_genome_with_annotations(
        self,
        name: str,
        fasta_path: str,
        annotations: List[Dict[str, str]],
        source: str = "custom",
        species: Optional[str] = None,
        description: Optional[str] = None
    ) -> bool:
        """添加基因组及其注释版本"""
        if not os.path.exists(fasta_path):
            print(f"[ERROR] 基因组文件不存在: {fasta_path}")
            return False
        
        if self.db.genome_exists(name):
            print(f"[ERROR] 基因组 '{name}' 已存在")
            return False
        
        # 创建索引
        self._ensure_fasta_index(fasta_path)
        
        # 取第一个注释作为默认（兼容旧字段）
        default_ann = annotations[0]["path"] if annotations else None
        
        # 添加基因组
        genome_id = self.db.add_genome(
            name=name,
            fasta_path=os.path.abspath(fasta_path),
            annotation_path=os.path.abspath(default_ann) if default_ann else None,
            source=source,
            species=species,
            description=description
        )
        
        # 添加所有注释版本
        for ann in annotations:
            ann_source = ann["source"]
            ann_path = os.path.abspath(ann["path"])
            if os.path.exists(ann_path):
                self.db.add_annotation(genome_id, ann_source, ann_path)
                self.enqueue_gene_id_list(name, ann_path, ann_source)
        
        print(f"[INFO] 已添加基因组: {name} (注释版本: {len(annotations)})")
        return True

    # ==================== 基因组获取 ====================

    def get(self, name_or_path: str) -> Tuple[Optional[str], Optional[str]]:
        """
        获取基因组路径

        参数可以是：
        - 数据库中的基因组名称
        - 本地文件路径

        Returns:
            (genome_path, annotation_path) 或 (None, None)
        """
        # 1. 优先从数据库查找（避免当前目录同名文件/目录干扰）
        genome = self.db.get_genome(name_or_path)
        if genome:
            return genome["fasta_path"], genome.get("annotation_path")

        # 2. 数据库没找到，检查是否是本地文件路径
        if os.path.exists(name_or_path):
            return name_or_path, None

        print(f"[ERROR] 未找到基因组: {name_or_path}")
        return None, None

    # ==================== 自定义基因组 ====================

    def add_custom(
        self,
        name: str,
        genome_path: str,
        annotation_path: Optional[str] = None,
        species: Optional[str] = None,
        description: Optional[str] = None
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
            if not self._gunzip_copy(genome_path, final_genome_path):
                return False

        if annotation_path and annotation_path.endswith(".gz"):
            genome_dir = self.cache_dir / name
            genome_dir.mkdir(parents=True, exist_ok=True)
            ext = ".gff3" if "gff" in annotation_path.lower() else ".gtf"
            final_annotation_path = str(genome_dir / f"{name}{ext}")
            print(f"[INFO] 解压注释: {annotation_path}")
            if not self._gunzip_copy(annotation_path, final_annotation_path):
                final_annotation_path = None

        # 检查并创建索引
        self._ensure_fasta_index(final_genome_path)

        # 保存到数据库
        genome_id = self.db.add_genome(
            name=name,
            fasta_path=os.path.abspath(final_genome_path),
            annotation_path=os.path.abspath(final_annotation_path) if final_annotation_path else None,
            source="custom",
            species=species,
            description=description
        )
        
        # 添加注释版本记录
        if final_annotation_path:
            self.db.add_annotation(genome_id, "version1", os.path.abspath(final_annotation_path))
            self.enqueue_gene_id_list(name, final_annotation_path, "version1")
        
        print(f"[INFO] 已添加自定义基因组: {name}")
        return True

    def _gunzip_copy(self, gz_file: str, out_file: str) -> bool:
        """解压 gzip 文件（保留原文件）"""
        tmp_file = f"{out_file}.tmp"
        try:
            with gzip.open(gz_file, "rb") as f_in:
                with open(tmp_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.replace(tmp_file, out_file)
            return True
        except (OSError, EOFError) as e:
            if os.path.exists(tmp_file):
                try:
                    os.remove(tmp_file)
                except OSError:
                    pass
            print(f"[ERROR] 解压失败: {gz_file} ({e})")
            return False

    def _gunzip(self, gz_file: str, out_file: str) -> bool:
        """解压 gzip 文件（删除原文件）"""
        if not self._gunzip_copy(gz_file, out_file):
            return False
        os.remove(gz_file)
        return True

    def remove_custom(self, name: str) -> bool:
        """移除自定义基因组（只删除配置，不删除文件）"""
        genome = self.db.get_genome(name)
        if not genome:
            print(f"[ERROR] 未找到基因组: {name}")
            return False
        
        if genome.get("source") != "custom":
            print(f"[ERROR] '{name}' 不是自定义基因组，请使用 remove_downloaded()")
            return False

        if not self.db.delete_genome(name):
            print(f"[ERROR] 删除基因组库记录失败: {name}")
            return False
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

    def download_ensembl(self, species_name: str, progress_callback=None, cancel_callback=None) -> bool:
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
            if progress_callback:
                progress_callback(100, "已存在")
            return True

        # 创建目录并放置下载中标记
        genome_dir = self.cache_dir / species_name
        genome_dir.mkdir(parents=True, exist_ok=True)
        downloading_marker = genome_dir / ".downloading"
        downloading_marker.touch()

        try:
            assembly = species_info.get("assembly") or ""
            release = species_info.get("release") or "59"
            species_cap = species_name[0].upper() + species_name[1:]

            if progress_callback:
                progress_callback(5, "下载FASTA")
            if cancel_callback and cancel_callback():
                raise InterruptedError("下载已取消")

            # 下载 FASTA
            fasta_filename = f"{species_cap}.{assembly}.dna.toplevel.fa.gz"
            fasta_url = f"{ENSEMBL_FTP_BASE}/fasta/{species_name}/dna/{fasta_filename}"
            fasta_gz = genome_dir / fasta_filename
            fasta_file = genome_dir / f"{species_name}.fa"

            print(f"[INFO] 下载 FASTA: {fasta_url}")
            if not self._download_file(fasta_url, str(fasta_gz), cancel_callback=cancel_callback):
                # 尝试 dna_sm 格式
                fasta_filename_alt = f"{species_cap}.{assembly}.dna_sm.toplevel.fa.gz"
                fasta_url_alt = f"{ENSEMBL_FTP_BASE}/fasta/{species_name}/dna/{fasta_filename_alt}"
                print(f"[INFO] 尝试备用 URL: {fasta_url_alt}")
                fasta_gz = genome_dir / fasta_filename_alt
                if not self._download_file(fasta_url_alt, str(fasta_gz), cancel_callback=cancel_callback):
                    return False

            if progress_callback:
                progress_callback(30, "解压FASTA")
            if cancel_callback and cancel_callback():
                raise InterruptedError("下载已取消")

            # 解压 FASTA
            print("[INFO] 解压 FASTA...")
            if not self._gunzip(str(fasta_gz), str(fasta_file)):
                return False

            if progress_callback:
                progress_callback(45, "构建索引")
            if cancel_callback and cancel_callback():
                raise InterruptedError("下载已取消")

            # 创建索引
            self._ensure_fasta_index(str(fasta_file))

            if progress_callback:
                progress_callback(55, "下载GFF")
            if cancel_callback and cancel_callback():
                raise InterruptedError("下载已取消")

            # 下载 GFF3 - 命名格式：ensembl_plants.{species_name}.gff3
            annotation_file = None
            ann_source = "ensembl_plants"
            for rel in [release, "62", "61", "60", "59", "58", "57"]:
                gff_filename = f"{species_cap}.{assembly}.{rel}.gff3.gz"
                gff_url = f"{ENSEMBL_FTP_BASE}/gff3/{species_name}/{gff_filename}"
                gff_gz = genome_dir / gff_filename
                gff_file = genome_dir / f"{ann_source}.{species_name}.gff3"

                print(f"[INFO] 下载 GFF3: {gff_url}")
                if self._download_file(gff_url, str(gff_gz), cancel_callback=cancel_callback):
                    if progress_callback:
                        progress_callback(80, "解压GFF")
                    if cancel_callback and cancel_callback():
                        raise InterruptedError("下载已取消")
                    print("[INFO] 解压 GFF3...")
                    if not self._gunzip(str(gff_gz), str(gff_file)):
                        return False
                    annotation_file = gff_file
                    break

                # 尝试 abinitio 版本
                gff_filename_ab = f"{species_cap}.{assembly}.{rel}.abinitio.gff3.gz"
                gff_url_ab = f"{ENSEMBL_FTP_BASE}/gff3/{species_name}/{gff_filename_ab}"
                gff_gz_ab = genome_dir / gff_filename_ab

                if self._download_file(gff_url_ab, str(gff_gz_ab), cancel_callback=cancel_callback):
                    if progress_callback:
                        progress_callback(80, "解压GFF")
                    if cancel_callback and cancel_callback():
                        raise InterruptedError("下载已取消")
                    print("[INFO] 解压 GFF3 (abinitio)...")
                    if not self._gunzip(str(gff_gz_ab), str(gff_file)):
                        return False
                    annotation_file = gff_file
                    break

            if not annotation_file:
                print("[WARN] GFF3 注释文件下载失败，基因组仍可使用")

            if progress_callback:
                progress_callback(90, "保存数据")
            if cancel_callback and cancel_callback():
                raise InterruptedError("下载已取消")

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
            
            # 添加注释版本记录
            if annotation_file:
                self.db.add_annotation(genome_id, ann_source, str(annotation_file))
                self.enqueue_gene_id_list(species_name, str(annotation_file), ann_source)

            if progress_callback:
                progress_callback(100, "完成")

            print(f"[INFO] 基因组 {species_name} 下载完成")
            return True
        
        finally:
            # 确保标记被删除
            if downloading_marker.exists():
                downloading_marker.unlink()

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

    def download_igv(self, genome_id: str, progress_callback=None, cancel_callback=None) -> bool:
        """从 IGV 下载基因组"""
        igv_genome = self._find_igv_genome(genome_id)
        if not igv_genome:
            print(f"[ERROR] IGV 中未找到基因组: {genome_id}")
            return False

        # 检查是否已存在
        if self.db.genome_exists(genome_id):
            print(f"[INFO] 基因组 '{genome_id}' 已存在")
            if progress_callback:
                progress_callback(100, "已存在")
            return True

        # 创建基因组目录并放置下载中标记
        genome_dir = self.cache_dir / genome_id
        genome_dir.mkdir(parents=True, exist_ok=True)
        downloading_marker = genome_dir / ".downloading"
        downloading_marker.touch()

        try:
            # 检查是否有注释文件
            has_annotation = False
            tracks = igv_genome.get("tracks", [])
            for track in tracks:
                if track.get("format") in ["gff3", "gff", "gtf"]:
                    has_annotation = True
                    break

            # 下载 FASTA
            fasta_url = igv_genome.get("fastaURL")
            if not fasta_url:
                print(f"[ERROR] 基因组 {genome_id} 没有 FASTA URL")
                return False

            if progress_callback:
                progress_callback(5, "下载FASTA")
            if cancel_callback and cancel_callback():
                raise InterruptedError("下载已取消")
            
            fasta_file = genome_dir / f"{genome_id}.fa"
            print(f"[INFO] 下载 FASTA: {fasta_url}")
            if not self._download_file(fasta_url, str(fasta_file), cancel_callback=cancel_callback):
                return False

            if progress_callback:
                progress_callback(40, "构建索引")
            if cancel_callback and cancel_callback():
                raise InterruptedError("下载已取消")
            
            # 创建索引
            self._ensure_fasta_index(str(fasta_file))

            if progress_callback:
                progress_callback(60, "索引完成")
            if cancel_callback and cancel_callback():
                raise InterruptedError("下载已取消")

            # 下载注释（如果有）- 命名格式：igv.{genome_id}.{ext}
            annotation_file = None
            ann_source = "igv"
            if has_annotation:
                for track in tracks:
                    if track.get("format") in ["gff3", "gff", "gtf"]:
                        ann_url = track.get("url")
                        if ann_url:
                            if progress_callback:
                                progress_callback(70, "下载注释")
                            ext = track.get("format", "gff3")
                            annotation_file = genome_dir / f"{ann_source}.{genome_id}.{ext}"
                            print(f"[INFO] 下载注释: {ann_url}")
                            self._download_file(ann_url, str(annotation_file), cancel_callback=cancel_callback)
                            break

            if progress_callback:
                progress_callback(90, "保存数据")
            if cancel_callback and cancel_callback():
                raise InterruptedError("下载已取消")

            # 保存到数据库
            genome_db_id = self.db.add_genome(
                name=genome_id,
                fasta_path=str(fasta_file),
                annotation_path=str(annotation_file) if annotation_file else None,
                display_name=igv_genome.get("name") or genome_id,
                source="igv",
                description=igv_genome.get("description")
            )
            
            # 添加注释版本记录
            if annotation_file:
                self.db.add_annotation(genome_db_id, ann_source, str(annotation_file))
                self.enqueue_gene_id_list(genome_id, str(annotation_file), ann_source)

            if progress_callback:
                progress_callback(100, "完成")

            print(f"[INFO] 基因组 {genome_id} 下载完成")
            return True
        
        finally:
            # 确保标记被删除
            if downloading_marker.exists():
                downloading_marker.unlink()

    # ==================== 通用方法 ====================

    def _download_file(self, url: str, dest: str, cancel_callback=None) -> bool:
        """下载文件"""
        try:
            print(f"  -> {dest}")
            # 添加 User-Agent 头，避免 S3 等服务器返回 403
            req = urllib.request.Request(url, headers={
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            })
            with urllib.request.urlopen(req, timeout=30, context=SSL_CONTEXT) as response:
                total = response.headers.get("Content-Length")
                total_size = int(total) if total and total.isdigit() else 0
                downloaded = 0
                with open(dest, 'wb') as f:
                    while True:
                        if cancel_callback and cancel_callback():
                            raise InterruptedError("下载已取消")
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                if total_size and downloaded < total_size:
                    print("[WARNING] 下载可能不完整，文件大小不匹配")
            return True
        except InterruptedError:
            try:
                if os.path.exists(dest):
                    os.remove(dest)
            except OSError:
                pass
            raise
        except Exception as e:
            print(f"[ERROR] 下载失败: {e}")
            try:
                if os.path.exists(dest):
                    os.remove(dest)
            except OSError:
                pass
            return False

    def _ensure_fasta_index(self, fasta_path: str):
        """确保 FASTA 文件有索引（pyfaidx 会自动创建）"""
        fai_path = f"{fasta_path}.fai"
        if not os.path.exists(fai_path):
            if fasta_path.lower().endswith(".gz") and not self._is_bgzf(fasta_path):
                print(f"[WARNING] 跳过索引: 非 BGZF 压缩 FASTA {fasta_path}")
                return
            print(f"[INFO] 创建 FASTA 索引: {fasta_path}")
            try:
                from pyfaidx import Fasta
                Fasta(fasta_path)
            except Exception as e:
                print(f"[WARNING] 创建 FASTA 索引失败: {fasta_path} ({e})")

    @staticmethod
    def _is_bgzf(path: str) -> bool:
        try:
            with open(path, "rb") as handle:
                header = handle.read(4)
            return header == b"\x1f\x8b\x08\x04"
        except OSError:
            return False

    def enqueue_gene_id_list(self, name: str, annotation_path: Optional[str], source: str = "version1") -> None:
        """将 gene id 列表生成任务加入队列"""
        if not annotation_path or not os.path.exists(annotation_path):
            return
        list_path = self._gene_id_list_path(name, source)
        try:
            ann_mtime = os.path.getmtime(annotation_path)
            if list_path.exists() and list_path.stat().st_size > 0:
                if list_path.stat().st_mtime >= ann_mtime:
                    return
        except OSError:
            return

        job_key = f"{name}:{source}"
        with self._gene_id_lock:
            if job_key in self._gene_id_jobs:
                return
            self._gene_id_jobs.add(job_key)
        self._gene_id_queue.put((name, annotation_path, source))

    def get_annotations(self, name: str) -> List[Dict[str, Any]]:
        """获取基因组的所有注释版本"""
        genome = self.db.get_genome(name)
        if not genome:
            return []
        return self.db.get_annotations(genome["id"])

    def get_gene_id_list_path(self, name: str, source: str = "version1") -> str:
        """获取 gene id 列表文件路径"""
        return str(self._gene_id_list_path(name, source))

    def _gene_id_worker_loop(self) -> None:
        while True:
            item = self._gene_id_queue.get()
            # 兼容旧格式 (name, path) 和新格式 (name, path, source)
            if len(item) == 2:
                name, annotation_path = item
                source = "version1"
            else:
                name, annotation_path, source = item
            
            job_key = f"{name}:{source}"
            try:
                gene_ids_path, gene_id_count = self._ensure_gene_id_list(name, annotation_path, source)
                if gene_ids_path:
                    # 更新 genome_annotations 表
                    genome = self.db.get_genome(name)
                    if genome:
                        ann = self.db.get_annotation_by_source(genome["id"], source)
                        if ann:
                            self.db.update_annotation(
                                ann["id"],
                                gene_ids_path=gene_ids_path,
                                gene_id_count=gene_id_count
                            )
            except Exception as e:
                print(f"[WARNING] Gene ID 列表生成失败: {annotation_path} ({e})")
            finally:
                with self._gene_id_lock:
                    self._gene_id_jobs.discard(job_key)
                self._gene_id_queue.task_done()

    def _gene_id_list_path(self, name: str, source: str = "version1") -> Path:
        """生成 gene id 列表文件路径：{source}.{name}.gene_ids.txt"""
        genome_dir = self.cache_dir / name
        genome_dir.mkdir(parents=True, exist_ok=True)
        return genome_dir / f"{source}.{name}.gene_ids.txt"

    def _ensure_gene_id_list(
        self,
        name: str,
        annotation_path: Optional[str],
        source: str = "version1"
    ) -> Tuple[Optional[str], Optional[int]]:
        """确保 gene id 列表存在且最新"""
        if not annotation_path or not os.path.exists(annotation_path):
            return None, None

        list_path = self._gene_id_list_path(name, source)
        # 兼容旧格式
        legacy_path = self.cache_dir / name / f"{name}{GENE_ID_LIST_SUFFIX}"
        legacy_txt = self.cache_dir / name / f"{name}{LEGACY_LIST_SUFFIX}"
        
        try:
            ann_mtime = os.path.getmtime(annotation_path)
            if list_path.exists() and list_path.stat().st_size > 0:
                if list_path.stat().st_mtime >= ann_mtime:
                    return str(list_path), count_gene_id_list(list_path)
            # 兼容旧 json 格式
            if not list_path.exists() and legacy_path.exists() and legacy_path.stat().st_size > 0:
                if legacy_path.stat().st_mtime >= ann_mtime:
                    ids = read_gene_id_list(legacy_path)
                    if ids:
                        write_gene_id_list(ids, list_path)
                        return str(list_path), len(ids)
            # 兼容旧 txt 格式
            if not list_path.exists() and legacy_txt.exists() and legacy_txt.stat().st_size > 0:
                if legacy_txt.stat().st_mtime >= ann_mtime:
                    ids = read_gene_id_list(legacy_txt)
                    if ids:
                        write_gene_id_list(ids, list_path)
                        return str(list_path), len(ids)
        except OSError:
            pass

        try:
            gene_ids = extract_gene_ids_from_annotation(annotation_path)
            if not gene_ids:
                return None, None
            write_gene_id_list(gene_ids, list_path)
            return str(list_path), len(gene_ids)
        except Exception as e:
            print(f"[WARNING] 生成 Gene ID 列表失败: {annotation_path} ({e})")
            return None, None

    def _update_gene_id_list(self, name: str, annotation_path: Optional[str]) -> None:
        gene_ids_path, gene_id_count = self._ensure_gene_id_list(name, annotation_path)
        if gene_ids_path:
            self.db.update_genome(name, notify=False, gene_ids_path=gene_ids_path, gene_id_count=gene_id_count)

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

        # 先删除数据库记录，确保基因组库立即生效；文件清理由后续步骤完成。
        if not self.db.delete_genome(genome_id, notify=False):
            print(f"[ERROR] 删除基因组库记录失败: {genome_id}")
            return False

        # 删除文件
        genome_dir = self.cache_dir / genome_id
        if genome_dir.exists():
            try:
                shutil.rmtree(genome_dir)
                print(f"[INFO] 已删除文件: {genome_dir}")
            except Exception as e:
                self.mark_cache_dir_ignored(genome)
                print(f"[WARNING] 已删除基因组库记录，但文件删除失败: {genome_dir} ({e})")
        else:
            self.mark_cache_dir_ignored(genome)
        self.db._notify_genomes_changed()
        self.db._notify_history_changed()
        print(f"[INFO] 已删除基因组: {genome_id}")
        return True

    def info(self, name: str) -> Optional[Dict[str, Any]]:
        """获取基因组详细信息"""
        # 先从数据库查找
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

    def download(self, genome_id: str, progress_callback=None, cancel_callback=None) -> bool:
        """下载基因组（自动识别来源）"""
        # 先检查 IGV
        if self._find_igv_genome(genome_id):
            return self.download_igv(genome_id, progress_callback, cancel_callback)
        
        # 再检查 Ensembl Plants
        species_list = self._fetch_ensembl_species()
        found = any(sp.get("name") == genome_id for sp in species_list)
        if found:
            return self.download_ensembl(genome_id, progress_callback, cancel_callback)
        
        print(f"[ERROR] 未找到基因组: {genome_id}")
        return False


# 全局实例
_manager_instance: Optional[GenomeManager] = None


def get_genome_manager() -> GenomeManager:
    """获取全局 GenomeManager 实例（单例模式）"""
    global _manager_instance
    expected_cache_dir = Path(get_db_dir())
    if _manager_instance is None or _manager_instance.cache_dir != expected_cache_dir:
        _manager_instance = GenomeManager(cache_dir=expected_cache_dir)
    return _manager_instance
