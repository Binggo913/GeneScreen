"""
GeneScreen 1.0 - 基因组管理页面

管理本地基因组库：添加、删除、搜索、下载
"""
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTabWidget,
    QLabel, QPushButton, QLineEdit, QTableWidget,
    QTableWidgetItem, QHeaderView, QFileDialog, QMessageBox,
    QGroupBox, QFormLayout, QTextEdit, QCheckBox, QStyledItemDelegate,
    QStyle, QStyleOptionViewItem, QApplication, QDialog, QScrollArea,
    QProgressBar, QFrame
)
from PySide6.QtCore import Qt, QThread, Signal, QTimer, QRectF, QSize
from PySide6.QtGui import QTextDocument, QTextOption

from pathlib import Path
import shutil

from core import get_genome_manager, get_database
from ui.widgets.table_checkbox import (
    TABLE_CHECKBOX_ROW_MIN_HEIGHT,
    configure_table_checkbox_column,
    create_table_checkbox,
    make_table_checkbox_cell,
    position_header_checkbox,
)


class NameEditDelegate(QStyledItemDelegate):
    """名称列编辑委托：编辑框完全填充单元格"""
    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        editor.setFrame(False)
        editor.setStyleSheet("""
            QLineEdit {
                border: 1px solid #667eea;
                border-radius: 3px;
                padding: 2px 4px;
                margin: 0px;
                background: transparent;
                color: palette(text);
            }
        """)
        return editor
    
    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)


class DownloadThread(QThread):
    """下载线程，支持进度报告和取消"""
    progress = Signal(int, str)  # (进度百分比, 阶段描述)
    finished = Signal(bool, str)
    
    def __init__(self, genome_id: str, parent=None):
        super().__init__(parent)
        self.genome_id = genome_id
        self._cancelled = False
        self._genome_dir = None
    
    def cancel(self):
        self._cancelled = True
    
    def run(self):
        try:
            manager = get_genome_manager()
            
            # 设置进度回调
            def progress_callback(percent, stage):
                if self._cancelled:
                    raise InterruptedError("下载已取消")
                self.progress.emit(percent, stage)
            
            def cancel_callback() -> bool:
                return self._cancelled
            
            # 记录基因组目录用于取消时清理
            self._genome_dir = manager.cache_dir / self.genome_id
            
            success = manager.download(
                self.genome_id,
                progress_callback=progress_callback,
                cancel_callback=cancel_callback
            )
            
            if self._cancelled:
                self._cleanup()
                self.finished.emit(False, "已取消")
            elif success:
                self.finished.emit(True, f"下载完成: {self.genome_id}")
            else:
                self.finished.emit(False, f"下载失败: {self.genome_id}")
        except InterruptedError:
            self._cleanup()
            self.finished.emit(False, "已取消")
        except Exception as e:
            self.finished.emit(False, f"下载出错: {str(e)}")
    
    def _cleanup(self):
        """清理下载的文件"""
        if self._genome_dir and self._genome_dir.exists():
            try:
                shutil.rmtree(self._genome_dir)
            except:
                pass


class SearchThread(QThread):
    """搜索线程"""
    finished = Signal(dict)
    
    def __init__(self, keyword: str, parent=None):
        super().__init__(parent)
        self.keyword = keyword
    
    def run(self):
        try:
            manager = get_genome_manager()
            results = manager.search(self.keyword)
            self.finished.emit(results)
        except Exception as e:
            self.finished.emit({"error": str(e)})


class GenomeManagerPage(QWidget):
    """基因组管理页面"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.search_thread = None
        self.download_thread = None
        self._gene_list_pending = []
        self._gene_list_timer = QTimer(self)
        self._gene_list_timer.setInterval(1000)
        self._gene_list_timer.timeout.connect(self._check_gene_list_ready)
        self._init_ui()
        self._load_genomes()

    def closeEvent(self, event):
        """关闭时等待线程完成"""
        if self.search_thread and self.search_thread.isRunning():
            self.search_thread.wait(1000)
        if self.download_thread and self.download_thread.isRunning():
            self.download_thread.wait(1000)
        super().closeEvent(event)
    
    def _init_ui(self):
        """初始化 UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 30, 30, 30)
        layout.setSpacing(20)
        
        title = QLabel("基因组管理")
        title.setProperty("role", "pageTitle")
        layout.addWidget(title)
        
        desc = QLabel("管理本地基因组库：添加本地基因组、搜索和下载在线基因组")
        desc.setProperty("role", "pageDesc")
        layout.addWidget(desc)
        
        tabs = QTabWidget()
        tabs.addTab(self._create_local_tab(), "📁 基因组库")
        tabs.addTab(self._create_search_tab(), "🔍 添加在线基因组")
        tabs.addTab(self._create_add_tab(), "➕ 添加本地基因组")
        layout.addWidget(tabs)
    
    def _create_local_tab(self) -> QWidget:
        """创建基因组库标签页"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        toolbar = QHBoxLayout()
        toolbar.addStretch()
        
        refresh_btn = QPushButton("刷新")
        refresh_btn.setFixedSize(50, 24)
        refresh_btn.setStyleSheet("font-size: 11px; padding: 2px 6px;")
        refresh_btn.clicked.connect(self._load_genomes)
        toolbar.addWidget(refresh_btn)

        edit_btn = QPushButton("编辑")
        edit_btn.setFixedSize(50, 24)
        edit_btn.setStyleSheet("font-size: 11px; padding: 2px 6px;")
        edit_btn.clicked.connect(self._edit_selected)
        toolbar.addWidget(edit_btn)

        delete_btn = QPushButton("删除")
        delete_btn.setFixedSize(50, 24)
        delete_btn.setStyleSheet("font-size: 11px; padding: 2px 6px;")
        delete_btn.clicked.connect(self._delete_selected)
        toolbar.addWidget(delete_btn)
        
        layout.addLayout(toolbar)
        
        self.genome_table = QTableWidget()
        self.genome_table.setColumnCount(5)
        self.genome_table.setHorizontalHeaderLabels(["", "名称", "来源", "基因组文件", "注释文件"])
        self.genome_table.verticalHeader().setVisible(False)
        
        header = self.genome_table.horizontalHeader()
        configure_table_checkbox_column(self.genome_table)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.Stretch)
        header.setSectionResizeMode(4, QHeaderView.Stretch)
        self.genome_table.setSelectionMode(QTableWidget.NoSelection)
        self.genome_table.setAlternatingRowColors(True)
        self.genome_table.setWordWrap(True)
        self.genome_table.setTextElideMode(Qt.ElideNone)

        self.genome_select_all = create_table_checkbox(self.genome_table.horizontalHeader())
        self.genome_select_all.setTristate(False)
        self.genome_select_all.stateChanged.connect(self._toggle_genome_all)
        self._position_genome_header_checkbox()
        header.sectionResized.connect(self._position_genome_header_checkbox)
        header.sectionMoved.connect(self._position_genome_header_checkbox)
        
        name_delegate = NameEditDelegate(self.genome_table)
        self.genome_table.setItemDelegateForColumn(1, name_delegate)
        self.genome_table.cellClicked.connect(self._on_cell_clicked)
        self.genome_table.itemChanged.connect(self._on_item_changed)
        
        layout.addWidget(self.genome_table)
        return widget

    def _on_search_table_resize(self, event):
        QTableWidget.resizeEvent(self.search_table, event)
        self._update_search_status_position()

    def _update_search_status_position(self):
        header_height = self.search_table.horizontalHeader().height()
        self.search_status.setGeometry(
            0, header_height,
            self.search_table.width(),
            self.search_table.height() - header_height
        )
        self.search_status.raise_()

    def _create_search_tab(self) -> QWidget:
        """创建添加在线基因组标签页"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        search_layout = QHBoxLayout()
        search_layout.setContentsMargins(0, 0, 0, 0)
        search_layout.setSpacing(10)

        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("输入关键词搜索 (如: rice, arabidopsis, human)")
        self.search_input.setMinimumWidth(350)
        self.search_input.returnPressed.connect(self._search_online)
        search_layout.addWidget(self.search_input)

        search_btn = QPushButton("搜索")
        search_btn.setFixedSize(50, 24)
        search_btn.setStyleSheet("font-size: 11px; padding: 2px 6px;")
        search_btn.clicked.connect(self._search_online)
        search_layout.addWidget(search_btn)

        download_selected_btn = QPushButton("下载选中")
        download_selected_btn.setFixedSize(70, 24)
        download_selected_btn.setStyleSheet("font-size: 11px; padding: 2px 6px;")
        download_selected_btn.clicked.connect(self._download_selected_genomes)
        search_layout.addWidget(download_selected_btn)

        search_layout.addStretch()
        layout.addLayout(search_layout)

        self.search_table = QTableWidget()
        self.search_table.setColumnCount(5)
        self.search_table.setHorizontalHeaderLabels(["", "来源", "ID", "名称", "描述"])
        self.search_table.verticalHeader().setVisible(False)
        header = self.search_table.horizontalHeader()
        configure_table_checkbox_column(self.search_table)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.Stretch)
        header.setSectionResizeMode(4, QHeaderView.Stretch)
        self.search_table.setSelectionMode(QTableWidget.NoSelection)
        self.search_table.setAlternatingRowColors(True)
        self.search_table.setWordWrap(True)
        self.search_table.setTextElideMode(Qt.ElideNone)
        layout.addWidget(self.search_table)

        self.search_status = QLabel("请输入关键词进行搜索", self.search_table)
        self.search_status.setAlignment(Qt.AlignCenter)
        self.search_status.setProperty("role", "searchStatus")
        self.search_status.setAttribute(Qt.WA_TransparentForMouseEvents)
        self.search_table.resizeEvent = self._on_search_table_resize

        self.search_select_all = create_table_checkbox(self.search_table.horizontalHeader())
        self.search_select_all.setTristate(False)
        self.search_select_all.stateChanged.connect(self._toggle_search_all)
        self._position_search_header_checkbox()
        header.sectionResized.connect(self._position_search_header_checkbox)
        header.sectionMoved.connect(self._position_search_header_checkbox)
        
        # 下载状态面板
        download_group = QGroupBox("下载状态")
        download_group.setMaximumHeight(200)
        download_group_layout = QVBoxLayout(download_group)
        download_group_layout.setContentsMargins(4, 2, 4, 2)
        download_group_layout.setSpacing(0)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QScrollArea.NoFrame)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        
        self.download_list_widget = QWidget()
        self.download_list_layout = QVBoxLayout(self.download_list_widget)
        self.download_list_layout.setContentsMargins(0, 0, 0, 0)
        self.download_list_layout.setSpacing(2)
        self.download_list_layout.addStretch()
        
        scroll_area.setWidget(self.download_list_widget)
        download_group_layout.addWidget(scroll_area)
        layout.addWidget(download_group)
        
        self.download_queue = []
        self.current_download = None
        self.download_items = {}
        self._gene_list_pending = []

        QTimer.singleShot(0, self._update_search_status_position)
        return widget

    def _create_add_tab(self) -> QWidget:
        """创建添加本地基因组标签页"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # 单个添加面板
        form_group = QGroupBox("添加本地基因组")
        form_layout = QFormLayout(form_group)
        
        # 名称
        self.add_name = QLineEdit()
        self.add_name.setPlaceholderText("默认从基因组文件名获取")
        form_layout.addRow("名称:", self.add_name)
        
        # 基因组文件（必须）
        fasta_layout = QHBoxLayout()
        self.add_fasta = QLineEdit()
        self.add_fasta.setPlaceholderText("选择基因组文件 (.fa, .fasta, .fa.gz)")
        self.add_fasta.textChanged.connect(self._on_fasta_changed)
        fasta_layout.addWidget(self.add_fasta, 1)
        fasta_btn = QPushButton("浏览...")
        fasta_btn.setProperty("secondary", True)
        fasta_btn.setProperty("compactAction", True)
        fasta_btn.setMinimumWidth(96)
        fasta_btn.setFixedHeight(38)
        fasta_btn.clicked.connect(self._browse_fasta)
        fasta_layout.addWidget(fasta_btn)
        form_layout.addRow("基因组文件:", fasta_layout)
        
        # 注释文件（可选）
        gff_layout = QHBoxLayout()
        self.add_gff = QLineEdit()
        self.add_gff.setPlaceholderText("选择注释文件 (.gff, .gff3, .gtf) [可选]")
        gff_layout.addWidget(self.add_gff, 1)
        gff_btn = QPushButton("浏览...")
        gff_btn.setProperty("secondary", True)
        gff_btn.setProperty("compactAction", True)
        gff_btn.setMinimumWidth(96)
        gff_btn.setFixedHeight(38)
        gff_btn.clicked.connect(self._browse_gff)
        gff_layout.addWidget(gff_btn)
        form_layout.addRow("注释文件:", gff_layout)
        
        layout.addWidget(form_group)
        
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        add_btn = QPushButton("➕ 添加基因组")
        add_btn.setMinimumWidth(150)
        add_btn.clicked.connect(self._add_custom_genome)
        btn_layout.addWidget(add_btn)
        layout.addLayout(btn_layout)
        
        # 批量添加面板
        batch_group = QGroupBox("批量添加本地基因组")
        batch_layout = QVBoxLayout(batch_group)
        
        batch_desc = QLabel("选择文件夹（含子目录），自动扫描基因组文件 (.fa, .fasta)，同名注释文件 (.gff, .gff3, .gtf) 会自动关联")
        batch_desc.setProperty("role", "muted")
        batch_desc.setWordWrap(True)
        batch_layout.addWidget(batch_desc)
        
        dir_layout = QHBoxLayout()
        self.batch_dir = QLineEdit()
        self.batch_dir.setPlaceholderText("选择包含基因组文件的文件夹")
        dir_layout.addWidget(self.batch_dir, 1)
        dir_btn = QPushButton("浏览...")
        dir_btn.setProperty("secondary", True)
        dir_btn.setProperty("compactAction", True)
        dir_btn.setMinimumWidth(96)
        dir_btn.setFixedHeight(38)
        dir_btn.clicked.connect(self._browse_batch_dir)
        dir_layout.addWidget(dir_btn)
        batch_layout.addLayout(dir_layout)
        
        # 扫描结果预览
        self.batch_preview = QLabel("")
        self.batch_preview.setProperty("role", "muted")
        batch_layout.addWidget(self.batch_preview)
        
        batch_btn_layout = QHBoxLayout()
        batch_btn_layout.addStretch()
        batch_add_btn = QPushButton("➕ 批量添加")
        batch_add_btn.setMinimumWidth(100)
        batch_add_btn.clicked.connect(self._add_batch_genomes)
        batch_btn_layout.addWidget(batch_add_btn)
        batch_layout.addLayout(batch_btn_layout)
        
        layout.addWidget(batch_group)
        
        layout.addStretch()
        return widget

    def _browse_batch_dir(self):
        """浏览批量添加的文件夹"""
        dir_path = QFileDialog.getExistingDirectory(self, "选择文件夹")
        if dir_path:
            self.batch_dir.setText(dir_path)
            self._scan_batch_dir()

    def _scan_batch_dir(self):
        """扫描文件夹中的基因组文件"""
        dir_path = self.batch_dir.text().strip()
        if not dir_path:
            self.batch_preview.setText("")
            self._batch_genomes = []
            return
        
        folder = Path(dir_path)
        if not folder.exists() or not folder.is_dir():
            self.batch_preview.setText("文件夹不存在")
            self._batch_genomes = []
            return
        
        manager = get_genome_manager()
        self._batch_genomes = manager.scan_local_genomes(folder)
        
        # 更新预览（统计注释版本数）
        if self._batch_genomes:
            total_anns = sum(len(g.get("annotations", [])) for g in self._batch_genomes)
            with_ann = sum(1 for g in self._batch_genomes if g.get("annotations"))
            self.batch_preview.setText(
                f"找到 {len(self._batch_genomes)} 个基因组，"
                f"其中 {with_ann} 个有注释文件（共 {total_anns} 个注释版本）"
            )
        else:
            self.batch_preview.setText("未找到基因组文件")

    def _add_batch_genomes(self):
        """批量添加基因组（支持多注释版本）"""
        if not hasattr(self, '_batch_genomes') or not self._batch_genomes:
            QMessageBox.warning(self, "提示", "请先扫描文件夹")
            return
        
        manager = get_genome_manager()
        success_count = 0
        fail_count = 0
        
        for genome in self._batch_genomes:
            annotations = genome.get("annotations", [])
            success = manager._add_genome_with_annotations(
                genome["name"], 
                genome["fasta"], 
                annotations
            )
            if success:
                success_count += 1
            else:
                fail_count += 1
        
        self._load_genomes()
        self._batch_genomes = []
        self.batch_dir.clear()
        self.batch_preview.setText("")
        
        if fail_count == 0:
            QMessageBox.information(self, "成功", f"已添加 {success_count} 个基因组")
        else:
            QMessageBox.warning(self, "完成", f"成功添加 {success_count} 个，失败 {fail_count} 个（可能名称已存在）")

    def _on_fasta_changed(self, text: str):
        """基因组文件路径变化时，自动填充名称"""
        if text and not self.add_name.text().strip():
            # 从文件路径提取文件名（不含扩展名）
            filename = Path(text).stem
            # 去除 .fa 等后缀
            if filename.endswith('.fa'):
                filename = filename[:-3]
            elif filename.endswith('.fasta'):
                filename = filename[:-6]
            self.add_name.setText(filename)
    
    def _on_cell_clicked(self, row: int, column: int):
        if column == 0:
            checkbox_widget = self.genome_table.cellWidget(row, 0)
            if checkbox_widget:
                checkbox = checkbox_widget.findChild(QCheckBox)
                if checkbox:
                    checkbox.setChecked(not checkbox.isChecked())
    
    def _on_item_changed(self, item: QTableWidgetItem):
        if item.column() == 1:
            old_name = item.data(Qt.UserRole)
            new_display_name = item.text().strip()
            if old_name and new_display_name:
                db = get_database()
                db.update_genome(old_name, display_name=new_display_name)
    
    def _load_genomes(self):
        self.genome_table.blockSignals(True)
        db = get_database()
        genomes = db.list_genomes()
        self.genome_table.setRowCount(len(genomes))
        
        for i, genome in enumerate(genomes):
            checkbox = create_table_checkbox()
            checkbox_widget = make_table_checkbox_cell(checkbox)
            self.genome_table.setCellWidget(i, 0, checkbox_widget)
            checkbox.stateChanged.connect(self._update_genome_header_checkbox)
            
            name_item = QTableWidgetItem(genome.get("display_name") or genome["name"])
            name_item.setData(Qt.UserRole, genome["name"])
            name_item.setFlags(name_item.flags() | Qt.ItemIsEditable)
            self.genome_table.setItem(i, 1, name_item)
            
            source_item = QTableWidgetItem(genome.get("source", "custom"))
            source_item.setFlags(source_item.flags() & ~Qt.ItemIsEditable)
            self.genome_table.setItem(i, 2, source_item)
            
            fasta_path = genome.get("fasta_path", "")
            fasta_display = Path(fasta_path).name if fasta_path else "-"
            fasta_item = QTableWidgetItem(fasta_display)
            fasta_item.setFlags(fasta_item.flags() & ~Qt.ItemIsEditable)
            fasta_item.setToolTip(fasta_path)
            self.genome_table.setItem(i, 3, fasta_item)

            # 获取所有注释版本
            annotations = db.get_annotations(genome["id"])
            if annotations:
                # 显示所有 source，用换行分隔
                sources = [ann["source"] for ann in annotations]
                ann_display = "\n".join(sources)
                # tooltip 显示完整路径
                tooltip_lines = [f"{ann['source']}: {ann['annotation_path']}" for ann in annotations]
                ann_tooltip = "\n".join(tooltip_lines)
            else:
                ann_display = "-"
                ann_tooltip = "无注释文件"
            
            ann_item = QTableWidgetItem(ann_display)
            ann_item.setFlags(ann_item.flags() & ~Qt.ItemIsEditable)
            ann_item.setToolTip(ann_tooltip)
            self.genome_table.setItem(i, 4, ann_item)
            
            # 根据注释版本数调整行高
            if len(annotations) > 1:
                self.genome_table.setRowHeight(
                    i, max(TABLE_CHECKBOX_ROW_MIN_HEIGHT, 20 * len(annotations))
                )
        
        self.genome_table.blockSignals(False)
        self._update_genome_header_checkbox()

    def _position_genome_header_checkbox(self):
        position_header_checkbox(self.genome_table, self.genome_select_all)

    def _toggle_genome_all(self, state):
        checked = self.genome_select_all.isChecked()
        for row in range(self.genome_table.rowCount()):
            checkbox_widget = self.genome_table.cellWidget(row, 0)
            if checkbox_widget:
                checkbox = checkbox_widget.findChild(QCheckBox)
                if checkbox:
                    checkbox.blockSignals(True)
                    checkbox.setChecked(checked)
                    checkbox.blockSignals(False)
        self._update_genome_header_checkbox()

    def _update_genome_header_checkbox(self):
        total = self.genome_table.rowCount()
        if total == 0:
            self.genome_select_all.setChecked(False)
            return
        checked = sum(
            1
            for row in range(total)
            if (w := self.genome_table.cellWidget(row, 0))
            and (cb := w.findChild(QCheckBox))
            and cb.isChecked()
        )
        self.genome_select_all.blockSignals(True)
        self.genome_select_all.setChecked(checked == total)
        self.genome_select_all.blockSignals(False)

    def _delete_selected(self):
        selected_names = self._get_selected_names()
        if not selected_names:
            QMessageBox.warning(self, "提示", "请先勾选要删除的基因组")
            return
        
        # 统计注释版本数量
        db = get_database()
        total_annotations = 0
        for name in selected_names:
            genome = db.get_genome(name)
            if genome:
                anns = db.get_annotations(genome["id"])
                total_annotations += len(anns)
        
        dialog = QDialog(self)
        dialog.setWindowTitle("确认删除")
        dialog.setMinimumWidth(350)
        
        layout = QVBoxLayout(dialog)
        msg_label = QLabel(f"确定要删除 {len(selected_names)} 个基因组吗？")
        layout.addWidget(msg_label)
        
        # 显示注释版本影响
        if total_annotations > 0:
            impact_label = QLabel(f"将同时删除 {total_annotations} 个注释版本及其 Gene ID 列表")
            impact_label.setStyleSheet("color: #f39c12; font-size: 12px;")
            layout.addWidget(impact_label)
        
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
            db = get_database()
            delete_files = delete_files_checkbox.isChecked()
            manager = get_genome_manager()
            cache_dir = Path(manager.cache_dir).resolve()
            deleted_any = False
            for name in selected_names:
                genome = db.get_genome(name)
                deleted = db.delete_genome(name, notify=False)
                if not deleted:
                    continue
                deleted_any = True
                if delete_files:
                    if genome and not self._delete_genome_files(genome, cache_dir):
                        manager.mark_cache_dir_ignored(genome)
                elif genome:
                    manager.mark_cache_dir_ignored(genome)
            if deleted_any:
                db._notify_genomes_changed()
                db._notify_history_changed()
            self._load_genomes()

    def _get_selected_names(self) -> list:
        selected_names = []
        for row in range(self.genome_table.rowCount()):
            checkbox_widget = self.genome_table.cellWidget(row, 0)
            if checkbox_widget:
                checkbox = checkbox_widget.findChild(QCheckBox)
                if checkbox and checkbox.isChecked():
                    name_item = self.genome_table.item(row, 1)
                    if name_item:
                        selected_names.append(name_item.data(Qt.UserRole))
        return selected_names

    def _delete_genome_files(self, genome: dict, cache_dir: Path) -> bool:
        paths = []
        fasta_path = genome.get("fasta_path", "") or ""
        ann_path = genome.get("annotation_path", "") or ""
        gene_ids_path = genome.get("gene_ids_path", "") or ""
        if fasta_path:
            paths.append(Path(fasta_path))
        if ann_path:
            paths.append(Path(ann_path))
        if gene_ids_path:
            paths.append(Path(gene_ids_path))

        existing = [p for p in paths if p.exists()]
        if not existing:
            return True

        delete_target = None
        for path in existing:
            if path.is_dir():
                delete_target = path
                break

        if delete_target is None:
            parents = [p.resolve().parent for p in existing]
            if parents and all(p == parents[0] for p in parents):
                common_parent = parents[0]
                if common_parent != cache_dir and cache_dir in common_parent.parents:
                    delete_target = common_parent

        if delete_target is None:
            for path in existing:
                resolved = path.resolve()
                parent = resolved.parent
                if cache_dir in resolved.parents and parent != cache_dir:
                    delete_target = parent
                    break

        if delete_target:
            try:
                shutil.rmtree(delete_target)
                return True
            except Exception:
                return False

        ok = True
        for path in existing:
            try:
                path.unlink()
            except Exception:
                ok = False
        return ok

    def _edit_selected(self):
        selected_names = self._get_selected_names()
        if len(selected_names) != 1:
            QMessageBox.warning(self, "提示", "请先勾选一个基因组进行编辑")
            return

        db = get_database()
        genome = db.get_genome(selected_names[0])
        if not genome:
            QMessageBox.warning(self, "提示", "未找到该基因组记录")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("编辑基因组")
        dialog.setMinimumWidth(520)

        layout = QVBoxLayout(dialog)
        form = QFormLayout()
        layout.addLayout(form)

        name_input = QLineEdit(genome.get("display_name") or genome["name"])
        form.addRow("名称:", name_input)

        # 基因组文件
        fasta_layout = QHBoxLayout()
        fasta_input = QLineEdit(genome.get("fasta_path", ""))
        fasta_input.setReadOnly(True)
        fasta_layout.addWidget(fasta_input, 1)
        fasta_btn = QPushButton("浏览")
        fasta_btn.setProperty("secondary", True)
        fasta_btn.setProperty("compactAction", True)
        fasta_btn.setProperty("smallAction", True)
        fasta_btn.setFixedWidth(64)
        fasta_btn.setFixedHeight(28)
        def _pick_fasta():
            start_dir = str(Path(fasta_input.text()).parent) if fasta_input.text().strip() else ""
            path, _ = QFileDialog.getOpenFileName(dialog, "选择 FASTA 文件", start_dir, "FASTA (*.fa *.fasta *.fa.gz)")
            if path:
                fasta_input.setText(path)
        fasta_btn.clicked.connect(_pick_fasta)
        fasta_layout.addWidget(fasta_btn)
        form.addRow("基因组文件:", fasta_layout)

        # 注释文件列表
        ann_group = QGroupBox("注释文件")
        ann_group_layout = QVBoxLayout(ann_group)
        ann_group_layout.setSpacing(8)
        ann_group_layout.setContentsMargins(10, 10, 10, 10)
        
        # 存储注释行的容器
        ann_rows_widget = QWidget()
        ann_rows_layout = QVBoxLayout(ann_rows_widget)
        ann_rows_layout.setContentsMargins(0, 0, 0, 0)
        ann_rows_layout.setSpacing(6)
        
        # 按钮样式
        small_btn_style = "font-size: 11px; padding: 2px 6px; min-width: 36px;"
        del_btn_style = "font-size: 11px; padding: 2px 6px; min-width: 36px; background-color: #e74c3c;"
        
        def _create_ann_row(ann_id, source, ann_path):
            """创建一个注释版本行（来源 + 路径输入框 + 浏览 + 删除）"""
            row_widget = QWidget()
            row_layout = QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setSpacing(6)
            
            # 来源标签
            source_label = QLabel(source)
            source_label.setFixedWidth(100)
            source_label.setStyleSheet("font-size: 12px;")
            row_layout.addWidget(source_label)
            
            # 路径输入框（只读，显示完整路径）
            path_input = QLineEdit(ann_path)
            path_input.setReadOnly(True)
            path_input.setProperty("ann_id", ann_id)
            row_layout.addWidget(path_input, 1)
            
            # 浏览按钮
            browse_btn = QPushButton("浏览")
            browse_btn.setProperty("secondary", True)
            browse_btn.setProperty("compactAction", True)
            browse_btn.setProperty("smallAction", True)
            browse_btn.setFixedWidth(64)
            browse_btn.setFixedHeight(28)
            def _browse():
                start_dir = str(Path(path_input.text()).parent) if path_input.text().strip() else ""
                path, _ = QFileDialog.getOpenFileName(dialog, "选择注释文件", start_dir, "注释 (*.gff *.gff3 *.gtf *.gff.gz)")
                if path:
                    path_input.setText(path)
                    db.update_annotation(ann_id, annotation_path=path)
            browse_btn.clicked.connect(_browse)
            row_layout.addWidget(browse_btn)
            
            # 删除按钮
            del_btn = QPushButton("删除")
            del_btn.setStyleSheet(del_btn_style)
            del_btn.setFixedWidth(45)
            def _delete():
                db.delete_annotation(ann_id)
                _refresh_ann_list()
            del_btn.clicked.connect(_delete)
            row_layout.addWidget(del_btn)
            
            return row_widget
        
        def _create_new_ann_row():
            """创建新增注释行（只有新增按钮）"""
            row_widget = QWidget()
            row_layout = QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setSpacing(6)
            
            row_layout.addStretch()
            
            # 新增按钮（点击后选择文件并直接添加）
            add_btn = QPushButton("+ 新增注释")
            add_btn.setStyleSheet(small_btn_style)
            add_btn.setFixedWidth(80)
            def _add_new():
                path, _ = QFileDialog.getOpenFileName(dialog, "选择注释文件", "", "注释 (*.gff *.gff3 *.gtf *.gff.gz)")
                if not path:
                    return
                # 自动解析来源
                from core.genome_manager import GenomeManager
                parsed_source, _ = GenomeManager._parse_annotation_name(Path(path).name)
                source = parsed_source or "version1"
                
                result = db.add_annotation(genome["id"], source, path)
                if result == -1:
                    QMessageBox.warning(dialog, "提示", f"来源 '{source}' 已存在")
                    return
                _refresh_ann_list()
                get_genome_manager().enqueue_gene_id_list(genome["name"], path, source)
            add_btn.clicked.connect(_add_new)
            row_layout.addWidget(add_btn)
            
            return row_widget
        
        def _refresh_ann_list():
            """刷新注释列表"""
            # 清空现有行
            while ann_rows_layout.count():
                item = ann_rows_layout.takeAt(0)
                if item.widget():
                    item.widget().deleteLater()
            
            # 添加现有注释行
            anns = db.get_annotations(genome["id"])
            for ann in anns:
                row = _create_ann_row(ann["id"], ann["source"], ann["annotation_path"])
                ann_rows_layout.addWidget(row)
            
            # 添加新增行
            new_row = _create_new_ann_row()
            ann_rows_layout.addWidget(new_row)
        
        # 初始化
        _refresh_ann_list()
        
        ann_group_layout.addWidget(ann_rows_widget)
        layout.addWidget(ann_group)

        # 按钮行
        btn_row = QHBoxLayout()
        btn_row.addStretch()
        cancel_btn = QPushButton("取消")
        cancel_btn.clicked.connect(dialog.reject)
        btn_row.addWidget(cancel_btn)
        save_btn = QPushButton("保存")
        btn_row.addWidget(save_btn)
        layout.addLayout(btn_row)

        def _save():
            display_name = name_input.text().strip()
            fasta_path = fasta_input.text().strip()
            if not display_name:
                QMessageBox.warning(dialog, "提示", "名称不能为空")
                return
            if not fasta_path:
                QMessageBox.warning(dialog, "提示", "基因组文件不能为空")
                return
            db.update_genome(genome["name"], display_name=display_name, fasta_path=fasta_path)
            dialog.accept()

        save_btn.clicked.connect(_save)
        if dialog.exec():
            self._load_genomes()

    def _search_online(self):
        keyword = self.search_input.text().strip()
        if not keyword:
            self.search_status.setText("请输入关键词进行搜索")
            self.search_status.show()
            self.search_table.setRowCount(0)
            return

        if self.search_thread and self.search_thread.isRunning():
            return

        self.search_status.setText("搜索中...")
        self.search_status.show()
        self.search_table.setRowCount(0)

        self.search_thread = SearchThread(keyword)
        self.search_thread.finished.connect(self._on_search_finished)
        self.search_thread.start()

    def _on_search_finished(self, results: dict):
        if "error" in results:
            self.search_status.setText(f"搜索出错: {results['error']}")
            self.search_status.show()
            return

        rows = []
        for r in results.get("igv", []):
            rows.append({"source": "igv", "id": r.get("id", ""), "name": r.get("name", ""), "desc": r.get("description", "") or "-"})
        for r in results.get("ensembl", []):
            desc_parts = [p for p in [r.get("common_name", ""), r.get("assembly", "")] if p]
            rows.append({"source": "ensembl_plants", "id": r.get("id", ""), "name": r.get("display_name", "") or r.get("id", ""), "desc": " / ".join(desc_parts) if desc_parts else "-"})

        if not rows:
            self.search_status.setText("未找到匹配结果")
            self.search_status.show()
            self.search_table.setRowCount(0)
            return

        self.search_status.hide()
        self.search_table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            checkbox = create_table_checkbox()
            checkbox_widget = make_table_checkbox_cell(checkbox)
            self.search_table.setCellWidget(i, 0, checkbox_widget)
            checkbox.stateChanged.connect(self._update_search_header_checkbox)

            for col, key in enumerate(["source", "id", "name", "desc"], 1):
                item = QTableWidgetItem(row[key])
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.search_table.setItem(i, col, item)


    def _position_search_header_checkbox(self):
        position_header_checkbox(self.search_table, self.search_select_all)

    def _toggle_search_all(self, state):
        checked = self.search_select_all.isChecked()
        for row in range(self.search_table.rowCount()):
            checkbox_widget = self.search_table.cellWidget(row, 0)
            if checkbox_widget:
                checkbox = checkbox_widget.findChild(QCheckBox)
                if checkbox:
                    checkbox.blockSignals(True)
                    checkbox.setChecked(checked)
                    checkbox.blockSignals(False)

    def _update_search_header_checkbox(self):
        total = self.search_table.rowCount()
        if total == 0:
            self.search_select_all.setChecked(False)
            return
        checked = sum(1 for row in range(total) if (w := self.search_table.cellWidget(row, 0)) and (cb := w.findChild(QCheckBox)) and cb.isChecked())
        self.search_select_all.blockSignals(True)
        self.search_select_all.setChecked(checked == total)
        self.search_select_all.blockSignals(False)

    def _download_selected_genomes(self):
        selected = []
        for row in range(self.search_table.rowCount()):
            checkbox_widget = self.search_table.cellWidget(row, 0)
            if checkbox_widget:
                checkbox = checkbox_widget.findChild(QCheckBox)
                if checkbox and checkbox.isChecked():
                    genome_id = self.search_table.item(row, 2).text()
                    genome_name = self.search_table.item(row, 3).text()
                    selected.append({"id": genome_id, "name": genome_name})
        
        if not selected:
            QMessageBox.warning(self, "提示", "请先勾选要下载的基因组")
            return
        
        for item in selected:
            exists = any(q["id"] == item["id"] for q in self.download_queue) or (self.current_download and self.current_download["id"] == item["id"])
            if not exists:
                self.download_queue.append({"id": item["id"], "name": item["name"], "status": "待下载"})
        
        self._update_download_list()
        self._process_download_queue()

    def _process_download_queue(self):
        if self.current_download or not self.download_queue:
            return
        
        self.current_download = self.download_queue.pop(0)
        self.current_download["status"] = "下载中"
        self.current_download["progress"] = 0
        self.current_download["stage"] = "准备中"
        self._update_download_list()
        
        self.download_thread = DownloadThread(self.current_download["id"])
        self.download_thread.progress.connect(self._on_download_progress)
        self.download_thread.finished.connect(self._on_download_finished)
        self.download_thread.start()

    def _on_download_progress(self, percent: int, stage: str):
        """更新下载进度"""
        if self.current_download:
            if self.current_download.get("_cancel_requested"):
                return
            self.current_download["progress"] = percent
            self.current_download["stage"] = stage
            # 直接更新进度条，避免重建整个列表
            if "_progress_bar" in self.current_download:
                self.current_download["_progress_bar"].setValue(percent)
            if "_stage_label" in self.current_download:
                self.current_download["_stage_label"].setText(stage)

    def _update_download_list(self):
        while self.download_list_layout.count() > 1:
            item = self.download_list_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        # 只显示待下载和下载中的条目
        all_items = []
        if self.current_download and self.current_download.get("status") in {"下载中", "待下载"}:
            all_items.append(self.current_download)
        all_items.extend([q for q in self.download_queue if q.get("status") in {"下载中", "待下载"}])
        all_items.extend(self._gene_list_pending)
        
        self.download_items = {}
        
        for idx, item in enumerate(all_items):
            # 添加分隔线（除了第一个）
            if idx > 0:
                line = QFrame()
                line.setFrameShape(QFrame.HLine)
                line.setFixedHeight(1)
                line.setProperty("role", "divider")
                self.download_list_layout.insertWidget(self.download_list_layout.count() - 1, line)
            
            row_widget = QWidget()
            row_layout = QVBoxLayout(row_widget)
            row_layout.setContentsMargins(4, 4, 4, 4)
            row_layout.setSpacing(4)
            
            # 第一行：名称、状态、按钮
            top_layout = QHBoxLayout()
            top_layout.setSpacing(8)
            
            name_label = QLabel(f"{item['name']} ({item['id']})")
            top_layout.addWidget(name_label, 1)
            
            status = item["status"]
            status_label = QLabel(status)
            status_label.setFixedWidth(60)
            status_label.setProperty("role", "downloadStatus")
            status_label.setProperty("state", "active" if status in {"下载中", "生成中"} else "pending")
            top_layout.addWidget(status_label)
            
            # 取消按钮（下载中和待下载都可以取消）
            if status == "下载中":
                cancel_btn = QPushButton("取消")
                cancel_btn.setFixedSize(40, 20)
                cancel_btn.setStyleSheet("font-size: 10px; padding: 0; color: red;")
                cancel_btn.clicked.connect(lambda _=False: self._cancel_current_download())
                top_layout.addWidget(cancel_btn)
            elif status == "待下载":
                cancel_btn = QPushButton("取消")
                cancel_btn.setFixedSize(40, 20)
                cancel_btn.setStyleSheet("font-size: 10px; padding: 0;")
                gid = item["id"]
                cancel_btn.clicked.connect(lambda _=False, g=gid: self._cancel_pending_download(g))
                top_layout.addWidget(cancel_btn)
            
            row_layout.addLayout(top_layout)
            
            # 第二行：进度条（下载中/生成中显示）
            if status in {"下载中", "生成中"}:
                progress_layout = QHBoxLayout()
                progress_layout.setSpacing(8)
                
                progress_bar = QProgressBar()
                progress_bar.setFixedHeight(16)
                if status == "生成中":
                    progress_bar.setRange(0, 0)
                else:
                    progress_bar.setRange(0, 100)
                    progress_bar.setValue(item.get("progress", 0))
                progress_bar.setProperty("role", "downloadProgress")
                progress_layout.addWidget(progress_bar)
                
                stage_label = QLabel(item.get("stage", ""))
                stage_label.setFixedWidth(80)
                stage_label.setProperty("role", "muted")
                progress_layout.addWidget(stage_label)
                
                row_layout.addLayout(progress_layout)
                
                # 保存进度条引用以便更新
                item["_progress_bar"] = progress_bar
                item["_stage_label"] = stage_label
            
            self.download_list_layout.insertWidget(self.download_list_layout.count() - 1, row_widget)
            self.download_items[item["id"]] = item

    def _cancel_current_download(self):
        """取消当前下载"""
        if self.download_thread and self.download_thread.isRunning():
            self.download_thread.cancel()
        if self.current_download:
            # 取消请求中仍显示为“下载中”，避免阻塞队列显示
            self.current_download["_cancel_requested"] = True
            self.current_download["stage"] = "取消中"
            self._update_download_list()

    def _cancel_pending_download(self, genome_id: str):
        self.download_queue = [q for q in self.download_queue if q["id"] != genome_id]
        if genome_id in self.download_items:
            del self.download_items[genome_id]
        self._update_download_list()
    
    def _on_download_finished(self, success: bool, message: str):
        finished_item = self.current_download
        if self.current_download:
            # 清理临时引用
            self.current_download.pop("_progress_bar", None)
            self.current_download.pop("_stage_label", None)
            self.current_download.pop("_cancel_requested", None)
            self.current_download = None
        
        if success:
            if finished_item:
                self._enqueue_gene_list_status(finished_item)
            self._load_genomes()
        self._update_download_list()
        self._process_download_queue()

    def _enqueue_gene_list_status(self, item: dict) -> None:
        """将 gene list 生成任务加入状态显示"""
        genome_id = item.get("id")
        if not genome_id:
            return
        if any(p.get("id") == genome_id for p in self._gene_list_pending):
            return
        
        db = get_database()
        genome = db.get_genome(genome_id)
        if not genome:
            return
        
        # 获取所有注释版本
        annotations = db.get_annotations(genome["id"])
        if not annotations:
            return
        
        display_name = genome.get("display_name") or item.get("name") or genome_id
        gm = get_genome_manager()
        
        # 为每个注释版本创建 gene list 生成任务
        for ann in annotations:
            source = ann["source"]
            list_path = Path(gm.get_gene_id_list_path(genome_id, source))
            try:
                if list_path.exists() and list_path.stat().st_size > 0:
                    continue
            except OSError:
                continue
            
            task_id = f"{genome_id}:{source}"
            if any(p.get("id") == task_id for p in self._gene_list_pending):
                continue
            
            self._gene_list_pending.append({
                "id": task_id,
                "genome_id": genome_id,
                "source": source,
                "name": f"{display_name} [{source}]",
                "status": "生成中",
                "stage": "生成gene list",
                "gene_list_path": str(list_path)
            })
        
        if self._gene_list_pending and not self._gene_list_timer.isActive():
            self._gene_list_timer.start()

    def _check_gene_list_ready(self) -> None:
        if not self._gene_list_pending:
            self._gene_list_timer.stop()
            return
        remaining = []
        for item in self._gene_list_pending:
            path = item.get("gene_list_path")
            if not path:
                continue
            list_path = Path(path)
            try:
                if list_path.exists() and list_path.stat().st_size > 0:
                    continue
            except OSError:
                pass
            remaining.append(item)
        if len(remaining) != len(self._gene_list_pending):
            self._gene_list_pending = remaining
            self._update_download_list()
        if not self._gene_list_pending:
            self._gene_list_timer.stop()

    def _browse_fasta(self):
        path, _ = QFileDialog.getOpenFileName(self, "选择 FASTA 文件", "", "FASTA 文件 (*.fa *.fasta *.fa.gz *.fasta.gz);;所有文件 (*)")
        if path:
            self.add_fasta.setText(path)
    
    def _browse_gff(self):
        path, _ = QFileDialog.getOpenFileName(self, "选择注释文件", "", "注释文件 (*.gff *.gff3 *.gtf *.gff.gz *.gff3.gz);;所有文件 (*)")
        if path:
            self.add_gff.setText(path)
    
    def _add_custom_genome(self):
        name = self.add_name.text().strip()
        fasta = self.add_fasta.text().strip()
        gff = self.add_gff.text().strip() or None
        
        if not fasta:
            QMessageBox.warning(self, "提示", "请选择基因组文件")
            return
        
        # 如果名称为空，从文件名获取
        if not name:
            filename = Path(fasta).stem
            if filename.endswith('.fa'):
                filename = filename[:-3]
            elif filename.endswith('.fasta'):
                filename = filename[:-6]
            name = filename
        
        if not name:
            QMessageBox.warning(self, "提示", "请输入基因组名称")
            return
        
        manager = get_genome_manager()
        success = manager.add_custom(name, fasta, gff)
        
        if success:
            QMessageBox.information(self, "成功", f"已添加基因组: {name}")
            self._load_genomes()
            self.add_name.clear()
            self.add_fasta.clear()
            self.add_gff.clear()
        else:
            QMessageBox.warning(self, "失败", "添加基因组失败，请检查文件路径或名称是否已存在")
