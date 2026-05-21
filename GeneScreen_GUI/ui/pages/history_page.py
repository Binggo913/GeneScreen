"""
GeneScreen 1.0 - 历史记录页面

查看和管理分析历史
"""
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QHeaderView, QMessageBox,
    QCheckBox, QDialog
)
from PySide6.QtCore import Qt, QUrl, QSize
from PySide6.QtGui import QDesktopServices, QColor, QIcon, QPalette
from datetime import datetime, timezone
import os
import re
import shutil
from pathlib import Path

from core import get_database
from core.config import get_output_dir
from ui.utils.icon_utils import draw_sidebar_icon
from ui.widgets.table_checkbox import (
    configure_table_checkbox,
    configure_table_checkbox_column,
    make_table_checkbox_cell,
    position_header_checkbox,
)


class HistoryPage(QWidget):
    """历史记录页面"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._db = get_database()
        self._syncing_history = False
        self._init_ui()
        self._load_history()
        self._db.add_history_listener(self._on_history_changed)
    
    def _init_ui(self):
        """初始化 UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 30, 30, 30)
        layout.setSpacing(20)
        
        # 标题
        title = QLabel("历史记录")
        title.setProperty("role", "pageTitle")
        layout.addWidget(title)
        
        desc = QLabel("查看和管理分析历史记录")
        desc.setProperty("role", "pageDesc")
        layout.addWidget(desc)
        
        panel = QWidget()
        panel_layout = QVBoxLayout(panel)
        panel_layout.setContentsMargins(0, 0, 0, 0)
        panel_layout.setSpacing(10)

        # 工具栏
        toolbar = QHBoxLayout()
        toolbar.addStretch()

        refresh_btn = QPushButton("刷新")
        refresh_btn.setFixedSize(50, 24)
        refresh_btn.setStyleSheet("font-size: 11px; padding: 2px 6px;")
        refresh_btn.clicked.connect(self._load_history)
        toolbar.addWidget(refresh_btn)
        
        delete_btn = QPushButton("删除")
        delete_btn.setFixedSize(50, 24)
        delete_btn.setStyleSheet("font-size: 11px; padding: 2px 6px;")
        delete_btn.clicked.connect(self._delete_selected)
        toolbar.addWidget(delete_btn)
        
        panel_layout.addLayout(toolbar)
        
        # 历史表格
        self.history_table = QTableWidget()
        self.history_table.setColumnCount(7)
        self.history_table.setHorizontalHeaderLabels([
            "", "ID", "模式", "时间", "状态", "报告", "结果路径"
        ])
        self.history_table.setIconSize(QSize(18, 18))
        configure_table_checkbox_column(self.history_table)
        self.history_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.history_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.history_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self.history_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeToContents)
        self.history_table.horizontalHeader().setSectionResizeMode(5, QHeaderView.Stretch)
        self.history_table.horizontalHeader().setSectionResizeMode(6, QHeaderView.Stretch)
        self.history_table.verticalHeader().setVisible(False)
        self.history_table.setSelectionMode(QTableWidget.NoSelection)
        self.history_table.setAlternatingRowColors(True)
        self.history_table.cellClicked.connect(self._on_cell_clicked)
        panel_layout.addWidget(self.history_table)
        layout.addWidget(panel)

        self.history_select_all = QCheckBox(self.history_table.horizontalHeader())
        configure_table_checkbox(self.history_select_all)
        self.history_select_all.setTristate(False)
        self.history_select_all.stateChanged.connect(self._toggle_history_all)
        header = self.history_table.horizontalHeader()
        self._position_history_header_checkbox()
        header.sectionResized.connect(self._position_history_header_checkbox)
        header.sectionMoved.connect(self._position_history_header_checkbox)
        
    def _format_path_display(self, path: str) -> str:
        if not path:
            return "-"
        name = os.path.basename(os.path.normpath(path))
        return name or path

    def _format_created_at(self, created_at: str) -> str:
        if not created_at:
            return ""
        try:
            dt = datetime.strptime(created_at[:19], "%Y-%m-%d %H:%M:%S")
            dt = dt.replace(tzinfo=timezone.utc).astimezone()
            return dt.strftime("%Y-%m-%d %H:%M:%S")
        except ValueError:
            return created_at[:19]

    def _make_checkbox_cell(self, checked: bool = False):
        checkbox = configure_table_checkbox(QCheckBox())
        checkbox.setChecked(checked)
        checkbox.stateChanged.connect(self._update_history_header_checkbox)
        checkbox_widget = make_table_checkbox_cell(checkbox)
        return checkbox_widget, checkbox

    def _get_checked_record_ids(self):
        checked_ids = []
        for row in range(self.history_table.rowCount()):
            checkbox_widget = self.history_table.cellWidget(row, 0)
            if not checkbox_widget:
                continue
            checkbox = checkbox_widget.findChild(QCheckBox)
            if checkbox and checkbox.isChecked():
                record = self.history_table.item(row, 1).data(Qt.UserRole)
                if record:
                    checked_ids.append(record.get("id"))
        return checked_ids

    def _position_history_header_checkbox(self):
        position_header_checkbox(self.history_table, self.history_select_all)

    def _toggle_history_all(self, state):
        checked = self.history_select_all.isChecked()
        for row in range(self.history_table.rowCount()):
            checkbox_widget = self.history_table.cellWidget(row, 0)
            if checkbox_widget:
                checkbox = checkbox_widget.findChild(QCheckBox)
                if checkbox:
                    checkbox.blockSignals(True)
                    checkbox.setChecked(checked)
                    checkbox.blockSignals(False)
        self._update_history_header_checkbox()

    def _update_history_header_checkbox(self):
        total = self.history_table.rowCount()
        if total == 0:
            self.history_select_all.setChecked(False)
            return
        checked = sum(
            1
            for row in range(total)
            if (w := self.history_table.cellWidget(row, 0))
            and (cb := w.findChild(QCheckBox))
            and cb.isChecked()
        )
        self.history_select_all.blockSignals(True)
        self.history_select_all.setChecked(checked == total)
        self.history_select_all.blockSignals(False)

    def _load_history(self):
        """加载历史记录"""
        if self._syncing_history:
            return
        self._syncing_history = True
        self._sync_history_from_output_dir()
        self._syncing_history = False
        checked_ids = set(self._get_checked_record_ids())
        history = self._db.get_history(limit=100)
        
        self.history_table.setRowCount(len(history))
        
        for i, record in enumerate(history):
            checkbox_widget, checkbox = self._make_checkbox_cell(
                record.get("id") in checked_ids
            )
            self.history_table.setCellWidget(i, 0, checkbox_widget)

            mode = record.get("mode", "")
            display_id = record.get("input_value") or str(record.get("id", ""))
            if mode == "location":
                match = re.match(r"^(\w+)_([0-9]+)_([0-9]+)$", str(display_id))
                if match:
                    display_id = f"{match.group(1)}:{match.group(2)}-{match.group(3)}"
            self.history_table.setItem(i, 1, QTableWidgetItem(str(display_id)))
            self.history_table.item(i, 1).setFlags(
                self.history_table.item(i, 1).flags() & ~Qt.ItemIsEditable
            )
            
            # 模式
            mode_display = {
                "gene_id": "Gene ID",
                "location": "Location",
                "sequence": "Sequence"
            }.get(mode, mode)
            mode_item = QTableWidgetItem(mode_display)
            icon_name = {
                "gene_id": "id",
                "location": "pin",
                "sequence": "code"
            }.get(mode)
            if icon_name:
                icon_color = QColor("#606266")
                if self.palette().color(QPalette.Base).lightness() < 128:
                    icon_color = QColor("#c7c9cc")
                pixmap = draw_sidebar_icon(
                    icon_name,
                    icon_color,
                    size=18,
                    widget=self
                )
                mode_item.setIcon(QIcon(pixmap))
            self.history_table.setItem(i, 2, mode_item)
            self.history_table.item(i, 2).setFlags(
                self.history_table.item(i, 2).flags() & ~Qt.ItemIsEditable
            )
            
            # 时间
            created_at = self._format_created_at(record.get("created_at", ""))
            self.history_table.setItem(i, 3, QTableWidgetItem(created_at))
            self.history_table.item(i, 3).setFlags(
                self.history_table.item(i, 3).flags() & ~Qt.ItemIsEditable
            )

            status_item = self._make_status_item(record.get("status", ""))
            self.history_table.setItem(i, 4, status_item)

            report_path = record.get("report_path", "") or "-"
            report_display = self._format_path_display(report_path)
            report_item = QTableWidgetItem(report_display)
            report_item.setToolTip(report_path)
            if report_path and report_path != "-":
                report_item.setForeground(QColor("#1d4ed8"))
                report_font = report_item.font()
                report_font.setUnderline(True)
                report_item.setFont(report_font)
            else:
                report_item.setForeground(QColor("#999999"))
            report_item.setFlags(report_item.flags() & ~Qt.ItemIsEditable)
            self.history_table.setItem(i, 5, report_item)

            output_dir = record.get("output_dir", "") or "-"
            output_display = self._format_path_display(output_dir)
            output_item = QTableWidgetItem(output_display)
            output_item.setToolTip(output_dir)
            if output_dir and output_dir != "-":
                output_item.setForeground(QColor("#1d4ed8"))
                output_font = output_item.font()
                output_font.setUnderline(True)
                output_item.setFont(output_font)
            else:
                output_item.setForeground(QColor("#999999"))
            output_item.setFlags(output_item.flags() & ~Qt.ItemIsEditable)
            self.history_table.setItem(i, 6, output_item)
            
            # 存储完整记录
            self.history_table.item(i, 1).setData(Qt.UserRole, record)
        self._update_history_header_checkbox()

    def _on_cell_clicked(self, row: int, column: int):
        if column == 0:
            checkbox_widget = self.history_table.cellWidget(row, 0)
            if checkbox_widget:
                checkbox = checkbox_widget.findChild(QCheckBox)
                if checkbox:
                    checkbox.setChecked(not checkbox.isChecked())
            return
        if column in (5, 6):
            self._open_result(row, column)

    def _make_status_item(self, status: str) -> QTableWidgetItem:
        status = (status or "pending").lower()
        label_map = {
            "completed": "✓ 已完成",
            "running": "Running",
            "reporting": "Reporting",
            "pending": "Pending",
            "failed": "Failed",
        }
        item = QTableWidgetItem(label_map.get(status, status))
        color_map = {
            "completed": QColor("#16a34a"),
            "running": QColor("#2563eb"),
            "reporting": QColor("#7c3aed"),
            "pending": QColor("#8a6d3b"),
            "failed": QColor("#dc2626"),
        }
        item.setForeground(color_map.get(status, QColor("#606266")))
        item.setFlags(item.flags() & ~Qt.ItemIsEditable)
        return item

    def _on_history_changed(self):
        if not self._syncing_history:
            self._load_history()
    
    def _delete_selected(self):
        """删除选中的记录"""
        checked_records = self._get_checked_records()
        if not checked_records:
            QMessageBox.warning(self, "提示", "请先勾选要删除的记录")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("确认删除")
        dialog.setMinimumWidth(360)
        if self.windowIcon():
            dialog.setWindowIcon(self.windowIcon())
        layout = QVBoxLayout(dialog)
        msg_label = QLabel(f"确定要删除选中的 {len(checked_records)} 条记录吗？")
        layout.addWidget(msg_label)

        delete_files_checkbox = QCheckBox("同时删除对应文件")
        delete_files_checkbox.setChecked(True)
        layout.addWidget(delete_files_checkbox)

        note_label = QLabel("注意：删除后无法恢复")
        note_label.setStyleSheet("color: #e74c3c; font-size: 12px;")
        layout.addWidget(note_label)

        layout.addSpacing(10)
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        cancel_btn = QPushButton("取消")
        cancel_btn.clicked.connect(dialog.reject)
        btn_layout.addWidget(cancel_btn)
        confirm_btn = QPushButton("删除")
        confirm_btn.setStyleSheet("background-color: #e74c3c; color: white;")
        confirm_btn.clicked.connect(dialog.accept)
        btn_layout.addWidget(confirm_btn)
        layout.addLayout(btn_layout)

        if dialog.exec() == QDialog.Accepted:
            delete_files = delete_files_checkbox.isChecked()
            deleting_ids = {
                int(record.get("id"))
                for record in checked_records
                if record.get("id") is not None
            }
            deleted_paths = set()
            for record in checked_records:
                history_id = record.get("id")
                files_deleted = False
                if delete_files:
                    files_deleted = self._delete_output_artifacts(record, deleting_ids, deleted_paths)
                if history_id:
                    self._db.delete_history(
                        history_id,
                        suppress_report=not delete_files or not files_deleted
                    )
            self._load_history()
    
    def _open_result(self, row: int, column: int):
        """打开结果目录"""
        record = self.history_table.item(row, 1).data(Qt.UserRole)

        report_path = record.get("report_path", "")
        output_dir = record.get("output_dir", "")

        if column == 5:
            if report_path and os.path.exists(report_path):
                QDesktopServices.openUrl(QUrl.fromLocalFile(report_path))
                return
            QMessageBox.warning(self, "提示", "报告文件不存在或未记录")
            return

        if column == 6:
            if output_dir and os.path.exists(output_dir):
                QDesktopServices.openUrl(QUrl.fromLocalFile(output_dir))
                return
            QMessageBox.warning(self, "提示", "输出目录不存在或未记录")
            return


    def closeEvent(self, event):
        if self._db:
            self._db.remove_history_listener(self._on_history_changed)
        super().closeEvent(event)

    def _get_checked_records(self):
        records = []
        for row in range(self.history_table.rowCount()):
            checkbox_widget = self.history_table.cellWidget(row, 0)
            if not checkbox_widget:
                continue
            checkbox = checkbox_widget.findChild(QCheckBox)
            if checkbox and checkbox.isChecked():
                record = self.history_table.item(row, 1).data(Qt.UserRole)
                if record:
                    records.append(record)
        return records

    @staticmethod
    def _canonical_path(path: str) -> str:
        return os.path.normcase(os.path.normpath(path))

    def _delete_output_artifacts(self, record: dict, deleting_ids: set, deleted_paths: set) -> bool:
        """删除当前历史项对应文件；共享路径仍被其他历史项引用时保留。"""
        output_dir = record.get("output_dir", "") or ""
        report_path = record.get("report_path", "") or ""
        files_deleted = False
        history_id = record.get("id")
        exclude_ids = sorted(deleting_ids)

        try:
            if output_dir:
                normalized_output = self._canonical_path(output_dir)
                path = Path(output_dir)
                shared_output = self._db.count_history_with_output_dir(
                    output_dir,
                    exclude_ids=exclude_ids
                ) > 0
                if normalized_output not in deleted_paths and path.exists() and not shared_output:
                    if path.is_dir():
                        shutil.rmtree(path)
                    else:
                        path.unlink()
                    deleted_paths.add(normalized_output)
                    return True

            if report_path:
                normalized_report = self._canonical_path(report_path)
                path = Path(report_path)
                shared_report = self._db.count_history_with_report_path(
                    report_path,
                    exclude_ids=[int(history_id)] if history_id is not None else []
                ) > 0
                if normalized_report not in deleted_paths and path.exists() and not shared_report:
                    path.unlink()
                    deleted_paths.add(normalized_report)
                    files_deleted = True
        except Exception:
            return False
        return files_deleted

    def _sync_history_from_output_dir(self) -> None:
        output_root = Path(get_output_dir())
        if not output_root.exists():
            return
        existing_reports = set(
            os.path.normcase(os.path.normpath(p))
            for p in self._db.get_history_report_paths()
        )
        deleted_reports = set(
            os.path.normcase(os.path.normpath(p))
            for p in self._db.get_deleted_history_report_paths()
        )
        mode_map = {
            "Gene_ID": "gene_id",
            "Location": "location",
            "Sequence": "sequence"
        }
        for folder, mode in mode_map.items():
            mode_dir = output_root / folder
            if not mode_dir.exists():
                continue
            for run_dir in mode_dir.iterdir():
                if not run_dir.is_dir():
                    continue
                new_report = run_dir / "report" / "index.html"
                if new_report.exists():
                    report_files = [new_report]
                else:
                    report_files = list((run_dir / "report").glob("*.report.html"))
                    report_files.extend(run_dir.glob("*.report.html"))
                for report_file in report_files:
                    report_path = os.path.normcase(os.path.normpath(str(report_file)))
                    if report_path in existing_reports or report_path in deleted_reports:
                        continue
                    if report_file.name == "index.html" and report_file.parent.name == "report":
                        input_value = run_dir.name
                    else:
                        filename = report_file.name
                        input_value = filename[:-len(".report.html")] if filename.endswith(".report.html") else report_file.stem
                    if mode == "location":
                        match = re.match(r"^(\w+)_([0-9]+)_([0-9]+)$", input_value)
                        if match:
                            input_value = f"{match.group(1)}:{match.group(2)}-{match.group(3)}"
                    self._db.add_history(
                        mode=mode,
                        input_value=input_value,
                        output_dir=str(run_dir),
                        report_path=str(report_file),
                        status="completed"
                    )
