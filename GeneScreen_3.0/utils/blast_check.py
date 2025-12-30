"""
GeneScreen 3.0 - BLAST+ 检测工具

检测系统是否安装了 BLAST+ 并获取版本信息
"""
import subprocess
import re
from typing import Optional, Tuple
from . import run_subprocess


def check_blast() -> bool:
    """检查 BLAST+ 是否已安装并可用"""
    try:
        result = run_subprocess(['blastn', '-version'], timeout=10)
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def get_blast_version() -> Optional[str]:
    """获取 BLAST+ 版本号"""
    try:
        result = run_subprocess(['blastn', '-version'], timeout=10)
        if result.returncode == 0:
            match = re.search(r'blastn:\s*(\S+)', result.stdout)
            if match:
                return match.group(1)
            lines = result.stdout.strip().split('\n')
            if lines:
                return lines[0].replace('blastn: ', '').strip()
        return None
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None


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
