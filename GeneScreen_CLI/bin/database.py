"""
GeneScreen 2.0 - SQLite 数据库管理模块

管理基因组库配置
"""
import sqlite3
from pathlib import Path
from contextlib import contextmanager
from typing import Optional, List, Dict, Any


# 默认数据库路径
DEFAULT_DB_PATH = Path.home() / ".genescreen" / "genescreen.db"


class Database:
    """SQLite 数据库管理器"""

    def __init__(self, db_path: Optional[Path] = None):
        self.db_path = Path(db_path) if db_path else DEFAULT_DB_PATH
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_tables()

    @contextmanager
    def connection(self):
        """获取数据库连接的上下文管理器"""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()

    def _init_tables(self):
        """初始化数据库表"""
        with self.connection() as conn:
            conn.executescript('''
                -- 基因组表
                CREATE TABLE IF NOT EXISTS genomes (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    name TEXT UNIQUE NOT NULL,
                    display_name TEXT,
                    fasta_path TEXT NOT NULL,
                    annotation_path TEXT,
                    source TEXT DEFAULT 'custom',
                    assembly TEXT,
                    species TEXT,
                    description TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                );

                -- 索引
                CREATE INDEX IF NOT EXISTS idx_genomes_name ON genomes(name);
                CREATE INDEX IF NOT EXISTS idx_genomes_source ON genomes(source);
            ''')

    # ==================== 基因组 CRUD ====================

    def add_genome(
        self,
        name: str,
        fasta_path: str,
        annotation_path: Optional[str] = None,
        display_name: Optional[str] = None,
        source: str = "custom",
        assembly: Optional[str] = None,
        species: Optional[str] = None,
        description: Optional[str] = None
    ) -> int:
        """添加基因组，返回新增记录的 ID"""
        with self.connection() as conn:
            cursor = conn.execute('''
                INSERT INTO genomes (
                    name, fasta_path, annotation_path, display_name,
                    source, assembly, species, description
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                name, fasta_path, annotation_path, display_name or name,
                source, assembly, species, description
            ))
            return cursor.lastrowid

    def get_genome(self, name: str) -> Optional[Dict[str, Any]]:
        """根据名称获取基因组"""
        with self.connection() as conn:
            cursor = conn.execute(
                'SELECT * FROM genomes WHERE name = ?', (name,)
            )
            row = cursor.fetchone()
            return dict(row) if row else None

    def list_genomes(self) -> List[Dict[str, Any]]:
        """列出所有基因组"""
        with self.connection() as conn:
            cursor = conn.execute(
                'SELECT * FROM genomes ORDER BY created_at DESC'
            )
            return [dict(row) for row in cursor.fetchall()]

    def list_by_source(self, source: str) -> List[Dict[str, Any]]:
        """按来源列出基因组"""
        with self.connection() as conn:
            cursor = conn.execute(
                'SELECT * FROM genomes WHERE source = ? ORDER BY name', (source,)
            )
            return [dict(row) for row in cursor.fetchall()]

    def update_genome(self, name: str, **kwargs) -> bool:
        """更新基因组信息"""
        if not kwargs:
            return False
        fields = ', '.join(f'{k} = ?' for k in kwargs.keys())
        values = list(kwargs.values()) + [name]
        with self.connection() as conn:
            cursor = conn.execute(
                f'UPDATE genomes SET {fields}, updated_at = CURRENT_TIMESTAMP WHERE name = ?',
                values
            )
            return cursor.rowcount > 0

    def delete_genome(self, name: str) -> bool:
        """删除基因组"""
        with self.connection() as conn:
            cursor = conn.execute(
                'DELETE FROM genomes WHERE name = ?', (name,)
            )
            return cursor.rowcount > 0

    def genome_exists(self, name: str) -> bool:
        """检查基因组是否存在"""
        with self.connection() as conn:
            cursor = conn.execute(
                'SELECT 1 FROM genomes WHERE name = ?', (name,)
            )
            return cursor.fetchone() is not None


# 全局数据库实例
_db_instance: Optional[Database] = None


def get_database() -> Database:
    """获取全局数据库实例（单例模式）"""
    global _db_instance
    if _db_instance is None:
        _db_instance = Database()
    return _db_instance
