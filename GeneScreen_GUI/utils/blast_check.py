"""
GeneScreen 1.0 - BLAST+ 检测工具

检测系统是否安装了 BLAST+ 并获取版本信息
"""
import subprocess
import re
import json
import time
from pathlib import Path
from typing import Optional, Tuple
from . import run_subprocess

# 缓存文件路径
_CACHE_FILE = Path.home() / ".genescreen" / "blast_cache.json"
_CACHE_TTL = 86400  # 缓存有效期 24 小时

# 内存缓存
_blast_cache = {"checked": False, "available": False, "version": None}


def _load_file_cache():
    """从文件加载缓存"""
    try:
        if _CACHE_FILE.exists():
            data = json.loads(_CACHE_FILE.read_text())
            if time.time() - data.get("timestamp", 0) < _CACHE_TTL:
                return data
    except:
        pass
    return None


def _save_file_cache():
    """保存缓存到文件"""
    try:
        _CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "available": _blast_cache["available"],
            "version": _blast_cache["version"],
            "timestamp": time.time()
        }
        _CACHE_FILE.write_text(json.dumps(data))
    except:
        pass


def _check_blast_once():
    """检查 BLAST+ 并缓存结果"""
    if _blast_cache["checked"]:
        return
    
    # 先尝试读文件缓存
    file_cache = _load_file_cache()
    if file_cache:
        _blast_cache["available"] = file_cache.get("available", False)
        _blast_cache["version"] = file_cache.get("version")
        _blast_cache["checked"] = True
        return
    
    # 文件缓存无效，调用 subprocess
    try:
        result = run_subprocess(['blastn', '-version'], timeout=10)
        if result.returncode == 0:
            _blast_cache["available"] = True
            match = re.search(r'blastn:\s*(\S+)', result.stdout)
            if match:
                _blast_cache["version"] = match.group(1)
            else:
                lines = result.stdout.strip().split('\n')
                if lines:
                    _blast_cache["version"] = lines[0].replace('blastn: ', '').strip()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    _blast_cache["checked"] = True
    _save_file_cache()


def check_blast() -> bool:
    """检查 BLAST+ 是否已安装并可用"""
    _check_blast_once()
    return _blast_cache["available"]


def get_blast_version() -> Optional[str]:
    """获取 BLAST+ 版本号"""
    _check_blast_once()
    return _blast_cache["version"]


def check_blast_installation() -> Tuple[bool, str]:
    """检查 BLAST+ 安装情况并返回提示信息"""
    if not check_blast():
        return False, "未找到 blastn 命令，请确认 BLAST+ 已安装并加入 PATH。"
    version = get_blast_version()
    if version:
        return True, f"blastn {version}"
    return True, "已检测到 blastn，但无法解析版本号。"


def check_makeblastdb() -> bool:
    """检查 makeblastdb 是否可用"""
    try:
        result = run_subprocess(['makeblastdb', '-version'], timeout=10)
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def get_blast_info() -> Tuple[bool, Optional[str], bool]:
    """获取 BLAST+ 完整信息"""
    blastn_ok = check_blast()
    version = get_blast_version() if blastn_ok else None
    makeblastdb_ok = check_makeblastdb()
    return blastn_ok, version, makeblastdb_ok


def ensure_blast_db(fasta_path: str) -> bool:
    """确保 BLAST 数据库存在，如果不存在则创建"""
    import os

    db_files = [f"{fasta_path}.{ext}" for ext in ['nin', 'nhr', 'nsq']]
    if all(os.path.exists(f) for f in db_files):
        return True

    print(f"[INFO] 创建 BLAST 数据库: {fasta_path}")
    try:
        result = run_subprocess(
            ['makeblastdb', '-in', fasta_path, '-dbtype', 'nucl'],
            timeout=300
        )
        if result.returncode != 0:
            print(f"[ERROR] 创建 BLAST 数据库失败: {result.stderr}")
            return False
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f"[ERROR] 创建 BLAST 数据库失败: {e}")
        return False
