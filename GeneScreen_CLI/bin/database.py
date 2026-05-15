"""
GeneScreen 1.0 - SQLite 数据库管理模块

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
        # 启用外键约束（SQLite 默认关闭）
        conn.execute("PRAGMA foreign_keys = ON")
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

                -- 注释版本表（一对多：一个基因组可有多个注释版本）
                CREATE TABLE IF NOT EXISTS genome_annotations (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    genome_id INTEGER NOT NULL,
                    source TEXT NOT NULL,
                    annotation_path TEXT NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (genome_id) REFERENCES genomes(id) ON DELETE CASCADE,
                    UNIQUE(genome_id, source)
                );

                -- 索引
                CREATE INDEX IF NOT EXISTS idx_genomes_name ON genomes(name);
                CREATE INDEX IF NOT EXISTS idx_genomes_source ON genomes(source);
                CREATE INDEX IF NOT EXISTS idx_annotations_genome ON genome_annotations(genome_id);
            ''')
            self._migrate_annotations(conn)

    def _migrate_annotations(self, conn):
        """迁移历史单注释记录到 genome_annotations 表"""
        cursor = conn.execute('''
            SELECT g.id, g.annotation_path
            FROM genomes g
            WHERE g.annotation_path IS NOT NULL AND g.annotation_path != ''
              AND NOT EXISTS (
                  SELECT 1 FROM genome_annotations ga WHERE ga.genome_id = g.id
              )
        ''')
        rows = cursor.fetchall()
        for row in rows:
            conn.execute('''
                INSERT INTO genome_annotations (genome_id, source, annotation_path)
                VALUES (?, 'version1', ?)
            ''', (row[0], row[1]))
        if rows:
            print(f"[INFO] 已迁移 {len(rows)} 条注释记录到 genome_annotations 表")

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

    def get_genome_by_id(self, genome_id: int) -> Optional[Dict[str, Any]]:
        """根据 ID 获取基因组"""
        with self.connection() as conn:
            cursor = conn.execute(
                'SELECT * FROM genomes WHERE id = ?', (genome_id,)
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

    # ==================== 注释版本 CRUD ====================

    def add_annotation(
        self,
        genome_id: int,
        source: str,
        annotation_path: str
    ) -> int:
        """添加注释版本记录，若已存在返回 -1"""
        try:
            with self.connection() as conn:
                cursor = conn.execute('''
                    INSERT INTO genome_annotations (genome_id, source, annotation_path)
                    VALUES (?, ?, ?)
                ''', (genome_id, source, annotation_path))
                return cursor.lastrowid
        except sqlite3.IntegrityError:
            return -1

    def get_annotations(self, genome_id: int) -> List[Dict[str, Any]]:
        """获取基因组的所有注释版本"""
        with self.connection() as conn:
            cursor = conn.execute(
                'SELECT * FROM genome_annotations WHERE genome_id = ? ORDER BY created_at',
                (genome_id,)
            )
            return [dict(row) for row in cursor.fetchall()]

    def get_annotation_by_source(self, genome_id: int, source: str) -> Optional[Dict[str, Any]]:
        """根据 genome_id 和 source 获取注释版本"""
        with self.connection() as conn:
            cursor = conn.execute(
                'SELECT * FROM genome_annotations WHERE genome_id = ? AND source = ?',
                (genome_id, source)
            )
            row = cursor.fetchone()
            return dict(row) if row else None

    def get_default_annotation(self, genome_id: int) -> Optional[Dict[str, Any]]:
        """获取基因组的默认注释（优先 version1，否则第一个）"""
        annotations = self.get_annotations(genome_id)
        if not annotations:
            return None
        for ann in annotations:
            if ann['source'] == 'version1':
                return ann
        return annotations[0]

    def delete_annotation(self, annotation_id: int) -> bool:
        """删除注释版本记录"""
        with self.connection() as conn:
            cursor = conn.execute(
                'DELETE FROM genome_annotations WHERE id = ?', (annotation_id,)
            )
            return cursor.rowcount > 0

    def update_annotation(self, annotation_id: int, **kwargs) -> bool:
        """更新注释版本记录"""
        if not kwargs:
            return False
        fields = ', '.join(f'{k} = ?' for k in kwargs.keys())
        values = list(kwargs.values()) + [annotation_id]
        with self.connection() as conn:
            cursor = conn.execute(
                f'UPDATE genome_annotations SET {fields} WHERE id = ?',
                values
            )
            return cursor.rowcount > 0


# 全局数据库实例
_db_instance: Optional[Database] = None


def get_database() -> Database:
    """获取全局数据库实例（单例模式）"""
    global _db_instance
    if _db_instance is None:
        _db_instance = Database()
    return _db_instance
