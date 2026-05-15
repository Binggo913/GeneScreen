"""
GeneScreen 1.0 - SQLite 数据库管理模块

管理基因组库和分析历史记录
"""
import sqlite3
from pathlib import Path
from contextlib import contextmanager
from datetime import datetime
from typing import Optional, List, Dict, Any, Callable

from .config import get_db_path


# 默认数据库路径
DEFAULT_DB_PATH = get_db_path()


class Database:
    """SQLite 数据库管理器"""

    def __init__(self, db_path: Optional[Path] = None):
        """
        初始化数据库

        Args:
            db_path: 数据库文件路径，默认为 ~/.genescreen/genomes/genescreen.db
        """
        self.db_path = Path(db_path) if db_path else DEFAULT_DB_PATH
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._history_listeners: List[Callable[[], None]] = []
        self._genome_listeners: List[Callable[[], None]] = []
        self._init_tables()

    def add_history_listener(self, callback: Callable[[], None]) -> None:
        if callback not in self._history_listeners:
            self._history_listeners.append(callback)

    def remove_history_listener(self, callback: Callable[[], None]) -> None:
        if callback in self._history_listeners:
            self._history_listeners.remove(callback)

    def _notify_history_changed(self) -> None:
        for callback in list(self._history_listeners):
            try:
                callback()
            except Exception as e:
                print(f"[WARNING] 历史记录监听器异常: {e}")

    def add_genome_listener(self, callback: Callable[[], None]) -> None:
        if callback not in self._genome_listeners:
            self._genome_listeners.append(callback)

    def remove_genome_listener(self, callback: Callable[[], None]) -> None:
        if callback in self._genome_listeners:
            self._genome_listeners.remove(callback)

    def _notify_genomes_changed(self) -> None:
        for callback in list(self._genome_listeners):
            try:
                callback()
            except Exception as e:
                print(f"[WARNING] 基因组监听器异常: {e}")

    @contextmanager
    def connection(self):
        """获取数据库连接的上下文管理器"""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        # 每次连接都启用外键约束（SQLite 默认关闭）
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
                    gene_ids_path TEXT,
                    gene_id_count INTEGER,
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
                    gene_ids_path TEXT,
                    gene_id_count INTEGER,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (genome_id) REFERENCES genomes(id) ON DELETE CASCADE,
                    UNIQUE(genome_id, source)
                );

                -- 分析历史表
                CREATE TABLE IF NOT EXISTS analysis_history (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    mode TEXT NOT NULL,
                    ref_genome_id INTEGER,
                    qry_genome_id INTEGER,
                    input_value TEXT,
                    identity REAL DEFAULT 90,
                    output_dir TEXT,
                    report_path TEXT,
                    status TEXT DEFAULT 'pending',
                    snp_count INTEGER,
                    indel_count INTEGER,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (ref_genome_id) REFERENCES genomes(id),
                    FOREIGN KEY (qry_genome_id) REFERENCES genomes(id)
                );

                -- 索引
                CREATE INDEX IF NOT EXISTS idx_genomes_name ON genomes(name);
                CREATE INDEX IF NOT EXISTS idx_history_created ON analysis_history(created_at);
                CREATE INDEX IF NOT EXISTS idx_annotations_genome ON genome_annotations(genome_id);
            ''')
            self._ensure_genome_columns(conn)
            self._ensure_history_columns(conn)
            self._migrate_annotations(conn)

    def _ensure_genome_columns(self, conn):
        """确保基因组表包含新增字段"""
        cursor = conn.execute("PRAGMA table_info(genomes)")
        columns = {row[1] for row in cursor.fetchall()}
        if "gene_ids_path" not in columns:
            conn.execute("ALTER TABLE genomes ADD COLUMN gene_ids_path TEXT")
        if "gene_id_count" not in columns:
            conn.execute("ALTER TABLE genomes ADD COLUMN gene_id_count INTEGER")

    def _ensure_history_columns(self, conn):
        """确保历史表包含新增字段"""
        cursor = conn.execute("PRAGMA table_info(analysis_history)")
        columns = {row[1] for row in cursor.fetchall()}
        if "report_path" not in columns:
            conn.execute("ALTER TABLE analysis_history ADD COLUMN report_path TEXT")

    def _migrate_annotations(self, conn):
        """迁移历史单注释记录到 genome_annotations 表"""
        # 查找 genomes 表中有 annotation_path 但未迁移到 genome_annotations 的记录
        cursor = conn.execute('''
            SELECT g.id, g.annotation_path, g.gene_ids_path, g.gene_id_count
            FROM genomes g
            WHERE g.annotation_path IS NOT NULL AND g.annotation_path != ''
              AND NOT EXISTS (
                  SELECT 1 FROM genome_annotations ga WHERE ga.genome_id = g.id
              )
        ''')
        rows = cursor.fetchall()
        for row in rows:
            conn.execute('''
                INSERT INTO genome_annotations (genome_id, source, annotation_path, gene_ids_path, gene_id_count)
                VALUES (?, 'version1', ?, ?, ?)
            ''', (row[0], row[1], row[2], row[3]))
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
        """
        添加基因组

        Returns:
            新增记录的 ID
        """
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
            genome_id = cursor.lastrowid
        self._notify_genomes_changed()
        return genome_id

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

    def update_genome(self, name: str, notify: bool = True, **kwargs) -> bool:
        """
        更新基因组信息

        Args:
            name: 基因组名称
            **kwargs: 要更新的字段

        Returns:
            是否更新成功
        """
        if not kwargs:
            return False

        # 构建 UPDATE 语句
        fields = ', '.join(f'{k} = ?' for k in kwargs.keys())
        values = list(kwargs.values()) + [name]

        with self.connection() as conn:
            cursor = conn.execute(
                f'UPDATE genomes SET {fields}, updated_at = CURRENT_TIMESTAMP WHERE name = ?',
                values
            )
            updated = cursor.rowcount > 0
        if updated and notify:
            self._notify_genomes_changed()
        return updated

    def delete_genome(self, name: str) -> bool:
        """删除基因组"""
        with self.connection() as conn:
            cursor = conn.execute(
                'DELETE FROM genomes WHERE name = ?', (name,)
            )
            deleted = cursor.rowcount > 0
        if deleted:
            self._notify_genomes_changed()
        return deleted

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
        annotation_path: str,
        gene_ids_path: Optional[str] = None,
        gene_id_count: Optional[int] = None
    ) -> int:
        """添加注释版本记录，若已存在返回 -1"""
        try:
            with self.connection() as conn:
                cursor = conn.execute('''
                    INSERT INTO genome_annotations (genome_id, source, annotation_path, gene_ids_path, gene_id_count)
                    VALUES (?, ?, ?, ?, ?)
                ''', (genome_id, source, annotation_path, gene_ids_path, gene_id_count))
                ann_id = cursor.lastrowid
            self._notify_genomes_changed()
            return ann_id
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

    def get_annotation_by_id(self, annotation_id: int) -> Optional[Dict[str, Any]]:
        """根据 ID 获取注释版本"""
        with self.connection() as conn:
            cursor = conn.execute(
                'SELECT * FROM genome_annotations WHERE id = ?', (annotation_id,)
            )
            row = cursor.fetchone()
            return dict(row) if row else None

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
            updated = cursor.rowcount > 0
        if updated:
            self._notify_genomes_changed()
        return updated

    def delete_annotation(self, annotation_id: int) -> bool:
        """删除注释版本记录"""
        with self.connection() as conn:
            cursor = conn.execute(
                'DELETE FROM genome_annotations WHERE id = ?', (annotation_id,)
            )
            deleted = cursor.rowcount > 0
        if deleted:
            self._notify_genomes_changed()
        return deleted

    # ==================== 分析历史 CRUD ====================

    def add_history(
        self,
        mode: str,
        ref_genome_id: Optional[int] = None,
        qry_genome_id: Optional[int] = None,
        input_value: Optional[str] = None,
        identity: float = 90,
        output_dir: Optional[str] = None,
        report_path: Optional[str] = None,
        status: str = "pending"
    ) -> int:
        """
        添加分析历史记录

        Returns:
            新增记录的 ID
        """
        with self.connection() as conn:
            cursor = conn.execute('''
                INSERT INTO analysis_history (
                    mode, ref_genome_id, qry_genome_id, input_value,
                    identity, output_dir, report_path, status
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                mode, ref_genome_id, qry_genome_id, input_value,
                identity, output_dir, report_path, status
            ))
            history_id = cursor.lastrowid
        self._notify_history_changed()
        return history_id

    def update_history(self, history_id: int, **kwargs) -> bool:
        """更新分析历史"""
        if not kwargs:
            return False

        fields = ', '.join(f'{k} = ?' for k in kwargs.keys())
        values = list(kwargs.values()) + [history_id]

        with self.connection() as conn:
            cursor = conn.execute(
                f'UPDATE analysis_history SET {fields} WHERE id = ?',
                values
            )
            updated = cursor.rowcount > 0
        if updated:
            self._notify_history_changed()
        return updated

    def get_history(self, limit: int = 50) -> List[Dict[str, Any]]:
        """
        获取分析历史

        Args:
            limit: 返回记录数量限制

        Returns:
            历史记录列表
        """
        with self.connection() as conn:
            cursor = conn.execute('''
                SELECT
                    h.*,
                    r.name as ref_name,
                    q.name as qry_name
                FROM analysis_history h
                LEFT JOIN genomes r ON h.ref_genome_id = r.id
                LEFT JOIN genomes q ON h.qry_genome_id = q.id
                ORDER BY h.created_at DESC
                LIMIT ?
            ''', (limit,))
            return [dict(row) for row in cursor.fetchall()]

    def get_history_report_paths(self) -> List[str]:
        """获取所有历史记录的报告路径"""
        with self.connection() as conn:
            cursor = conn.execute('''
                SELECT report_path
                FROM analysis_history
                WHERE report_path IS NOT NULL AND report_path != ''
            ''')
            return [row["report_path"] for row in cursor.fetchall() if row["report_path"]]

    def get_history_by_id(self, history_id: int) -> Optional[Dict[str, Any]]:
        """根据 ID 获取历史记录"""
        with self.connection() as conn:
            cursor = conn.execute('''
                SELECT
                    h.*,
                    r.name as ref_name,
                    q.name as qry_name
                FROM analysis_history h
                LEFT JOIN genomes r ON h.ref_genome_id = r.id
                LEFT JOIN genomes q ON h.qry_genome_id = q.id
                WHERE h.id = ?
            ''', (history_id,))
            row = cursor.fetchone()
            return dict(row) if row else None

    def delete_history(self, history_id: int) -> bool:
        """删除历史记录"""
        with self.connection() as conn:
            cursor = conn.execute(
                'DELETE FROM analysis_history WHERE id = ?', (history_id,)
            )
            deleted = cursor.rowcount > 0
        if deleted:
            self._notify_history_changed()
        return deleted

    def clear_history(self) -> int:
        """清空所有历史记录"""
        with self.connection() as conn:
            cursor = conn.execute('DELETE FROM analysis_history')
            count = cursor.rowcount
        if count:
            self._notify_history_changed()
        return count


# 全局数据库实例
_db_instance: Optional[Database] = None


def get_database() -> Database:
    """获取全局数据库实例（单例模式）"""
    global _db_instance
    expected_path = Path(get_db_path())
    if _db_instance is None or _db_instance.db_path != expected_path:
        _db_instance = Database(expected_path)
    return _db_instance
