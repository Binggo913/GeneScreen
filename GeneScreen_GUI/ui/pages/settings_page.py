"""
GeneScreen 1.0 - 设置页面
"""
from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QGroupBox, QFormLayout, QFileDialog, QMessageBox
)
from PySide6.QtCore import Qt

from core.config import get_db_dir, set_db_dir, get_output_dir, set_output_dir
from core.database import get_database
from core.genome_manager import get_genome_manager


class SettingsPage(QWidget):
    """设置页面"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._init_ui()

    def _init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 30, 30, 30)
        layout.setSpacing(20)
        layout.setAlignment(Qt.AlignTop)

        title = QLabel("设置")
        title.setProperty("role", "pageTitle")
        layout.addWidget(title)

        desc = QLabel("配置应用路径与存储位置")
        desc.setProperty("role", "pageDesc")
        layout.addWidget(desc)

        db_group = QGroupBox("数据库")
        db_layout = QFormLayout(db_group)

        path_layout = QHBoxLayout()
        self.db_path_input = QLineEdit()
        self.db_path_input.setPlaceholderText("选择数据库存储目录")
        current_path = str(get_db_dir())
        self.db_path_input.setText(current_path)
        self._last_saved_path = current_path
        path_layout.addWidget(self.db_path_input, 1)
        self.db_path_input.editingFinished.connect(self._save_settings)

        browse_btn = QPushButton("浏览...")
        browse_btn.clicked.connect(self._browse_db_dir)
        path_layout.addWidget(browse_btn)

        db_layout.addRow("数据库路径:", path_layout)

        hint = QLabel("默认路径: 用户目录/.genescreen/genomes/")
        hint.setProperty("role", "muted")
        db_layout.addRow("", hint)

        layout.addWidget(db_group)

        output_group = QGroupBox("结果输出")
        output_layout = QFormLayout(output_group)

        output_path_layout = QHBoxLayout()
        self.output_path_input = QLineEdit()
        self.output_path_input.setPlaceholderText("选择结果输出目录")
        current_output_path = str(get_output_dir())
        self.output_path_input.setText(current_output_path)
        self._last_saved_output_path = current_output_path
        output_path_layout.addWidget(self.output_path_input, 1)
        self.output_path_input.editingFinished.connect(self._save_output_settings)

        output_browse_btn = QPushButton("浏览...")
        output_browse_btn.clicked.connect(self._browse_output_dir)
        output_path_layout.addWidget(output_browse_btn)

        output_layout.addRow("结果输出路径:", output_path_layout)

        output_hint = QLabel("默认路径: 数据库同级/GeneScreenOutput/")
        output_hint.setProperty("role", "muted")
        output_layout.addRow("", output_hint)

        layout.addWidget(output_group)

    def _browse_db_dir(self):
        current = self.db_path_input.text().strip()
        start_dir = current if current else str(get_db_dir())
        directory = QFileDialog.getExistingDirectory(self, "选择数据库存储目录", start_dir)
        if directory:
            self.db_path_input.setText(directory)
            self._save_settings()

    def _save_settings(self):
        path_text = self.db_path_input.text().strip()
        if not path_text or path_text == self._last_saved_path:
            return

        db_dir = Path(path_text).expanduser()
        try:
            db_dir.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            QMessageBox.warning(self, "失败", f"无法创建目录: {e}")
            return

        set_db_dir(db_dir)
        get_database()
        get_genome_manager()
        self._last_saved_path = str(db_dir)
        QMessageBox.information(self, "已保存", "数据库路径已更新，建议重启应用以确保所有页面生效。")

    def _browse_output_dir(self):
        current = self.output_path_input.text().strip()
        start_dir = current if current else str(get_output_dir())
        directory = QFileDialog.getExistingDirectory(self, "选择结果输出目录", start_dir)
        if directory:
            self.output_path_input.setText(directory)
            self._save_output_settings()

    def _save_output_settings(self):
        path_text = self.output_path_input.text().strip()
        if not path_text or path_text == self._last_saved_output_path:
            return

        output_dir = Path(path_text).expanduser()
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            QMessageBox.warning(self, "失败", f"无法创建目录: {e}")
            return

        set_output_dir(output_dir)
        self._last_saved_output_path = str(output_dir)
        QMessageBox.information(self, "已保存", "结果输出路径已更新。")
