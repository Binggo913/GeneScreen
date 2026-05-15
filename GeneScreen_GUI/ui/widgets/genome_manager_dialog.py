"""
GeneScreen 1.0 - 基因组管理对话框

管理本地基因组库：添加、删除、搜索、下载
"""
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTabWidget,
    QWidget, QLabel, QPushButton, QLineEdit, QTableWidget,
    QTableWidgetItem, QHeaderView, QFileDialog, QMessageBox,
    QGroupBox, QFormLayout, QProgressBar, QTextEdit,
    QStyledItemDelegate, QStyle, QStyleOptionViewItem, QApplication
)
from PySide6.QtCore import Qt, QThread, Signal, QRectF, QSize
from PySide6.QtGui import QTextDocument, QTextOption

from core import get_genome_manager, get_database


class WrapAnywhereDelegate(QStyledItemDelegate):
    """表格单元格文本自动换行（路径等无空格字符串）"""

    def paint(self, painter, option, index):
        options = QStyleOptionViewItem(option)
        self.initStyleOption(options, index)

        text = options.text
        options.text = ""

        style = options.widget.style() if options.widget else QApplication.style()
        style.drawControl(QStyle.CE_ItemViewItem, options, painter, options.widget)

        doc = QTextDocument()
        doc.setDefaultFont(options.font)
        text_option = doc.defaultTextOption()
        text_option.setWrapMode(QTextOption.WrapAnywhere)
        doc.setDefaultTextOption(text_option)
        doc.setPlainText(text)

        text_rect = style.subElementRect(QStyle.SE_ItemViewItemText, options, options.widget)
        doc.setTextWidth(text_rect.width())

        painter.save()
        painter.translate(text_rect.topLeft())
        doc.drawContents(painter, QRectF(0, 0, text_rect.width(), text_rect.height()))
        painter.restore()

    def sizeHint(self, option, index):
        options = QStyleOptionViewItem(option)
        self.initStyleOption(options, index)

        doc = QTextDocument()
        doc.setDefaultFont(options.font)
        text_option = doc.defaultTextOption()
        text_option.setWrapMode(QTextOption.WrapAnywhere)
        doc.setDefaultTextOption(text_option)
        doc.setPlainText(options.text)

        width = options.rect.width()
        if width <= 0 and options.widget:
            width = options.widget.columnWidth(index.column())
        if width <= 0:
            width = 200

        doc.setTextWidth(width)
        size = doc.size()
        return QSize(int(size.width()), int(size.height()))


class DownloadThread(QThread):
    """下载线程"""
    progress = Signal(str)
    finished = Signal(bool, str)
    
    def __init__(self, genome_id: str, parent=None):
        super().__init__(parent)
        self.genome_id = genome_id
    
    def run(self):
        try:
            manager = get_genome_manager()
            self.progress.emit(f"正在下载 {self.genome_id}...")
            success = manager.download(self.genome_id)
            if success:
                self.finished.emit(True, f"下载完成: {self.genome_id}")
            else:
                self.finished.emit(False, f"下载失败: {self.genome_id}")
        except Exception as e:
            self.finished.emit(False, f"下载出错: {str(e)}")


class GenomeManagerDialog(QDialog):
    """基因组管理对话框"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("基因组管理")
        self.setMinimumSize(800, 600)
        self.resize(900, 650)
        
        # 继承父窗口图标
        if parent and parent.windowIcon():
            self.setWindowIcon(parent.windowIcon())
        
        self._init_ui()
        self._load_genomes()
    
    def _init_ui(self):
        """初始化 UI"""
        layout = QVBoxLayout(self)
        layout.setSpacing(15)
        
        # 标签页
        tabs = QTabWidget()
        
        # 本地基因组
        tabs.addTab(self._create_local_tab(), "📁 本地基因组")
        
        # 在线搜索
        tabs.addTab(self._create_search_tab(), "🔍 在线搜索")
        
        # 添加自定义
        tabs.addTab(self._create_add_tab(), "➕ 添加自定义")
        
        layout.addWidget(tabs)
        
        # 底部按钮
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        
        close_btn = QPushButton("关闭")
        close_btn.clicked.connect(self.accept)
        btn_layout.addWidget(close_btn)
        
        layout.addLayout(btn_layout)
    
    def _create_local_tab(self) -> QWidget:
        """创建本地基因组标签页"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # 工具栏
        toolbar = QHBoxLayout()
        
        refresh_btn = QPushButton("🔄 刷新")
        refresh_btn.clicked.connect(self._load_genomes)
        toolbar.addWidget(refresh_btn)
        
        toolbar.addStretch()
        
        delete_btn = QPushButton("🗑️ 删除选中")
        delete_btn.clicked.connect(self._delete_selected)
        toolbar.addWidget(delete_btn)
        
        layout.addLayout(toolbar)
        
        # 表格
        self.genome_table = QTableWidget()
        self.genome_table.setColumnCount(5)
        self.genome_table.setHorizontalHeaderLabels(["名称", "来源", "FASTA 路径", "注释文件", "添加时间"])
        self.genome_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.genome_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.genome_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.genome_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.Stretch)
        self.genome_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeToContents)
        self.genome_table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.genome_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.genome_table.setAlternatingRowColors(True)
        self.genome_table.setWordWrap(True)
        self.genome_table.setTextElideMode(Qt.ElideNone)
        self.genome_table.setItemDelegateForColumn(2, WrapAnywhereDelegate(self.genome_table))
        self.genome_table.setItemDelegateForColumn(3, WrapAnywhereDelegate(self.genome_table))
        layout.addWidget(self.genome_table)
        
        return widget
    
    def _create_search_tab(self) -> QWidget:
        """创建在线搜索标签页"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # 搜索框
        search_layout = QHBoxLayout()
        
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("输入关键词搜索 (如: rice, arabidopsis, human)")
        self.search_input.returnPressed.connect(self._search_online)
        search_layout.addWidget(self.search_input, 1)
        
        search_btn = QPushButton("搜索")
        search_btn.clicked.connect(self._search_online)
        search_layout.addWidget(search_btn)
        
        layout.addLayout(search_layout)
        
        # 搜索结果
        self.search_results = QTextEdit()
        self.search_results.setReadOnly(True)
        self.search_results.setPlaceholderText("搜索结果将显示在这里...")
        layout.addWidget(self.search_results)
        
        # 下载区域
        download_layout = QHBoxLayout()
        
        self.download_input = QLineEdit()
        self.download_input.setPlaceholderText("输入基因组 ID 进行下载")
        download_layout.addWidget(self.download_input, 1)
        
        download_btn = QPushButton("下载")
        download_btn.clicked.connect(self._download_genome)
        download_layout.addWidget(download_btn)
        
        layout.addLayout(download_layout)
        
        # 进度
        self.download_progress = QLabel("")
        layout.addWidget(self.download_progress)
        
        return widget
    
    def _create_add_tab(self) -> QWidget:
        """创建添加自定义标签页"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # 表单
        form_group = QGroupBox("添加自定义基因组")
        form_layout = QFormLayout(form_group)
        
        # 名称
        self.add_name = QLineEdit()
        self.add_name.setPlaceholderText("基因组名称 (如: my_genome)")
        form_layout.addRow("名称:", self.add_name)
        
        # FASTA 文件
        fasta_layout = QHBoxLayout()
        self.add_fasta = QLineEdit()
        self.add_fasta.setPlaceholderText("FASTA 文件路径 (.fa, .fasta, .fa.gz)")
        fasta_layout.addWidget(self.add_fasta, 1)
        
        fasta_btn = QPushButton("浏览...")
        fasta_btn.clicked.connect(self._browse_fasta)
        fasta_layout.addWidget(fasta_btn)
        form_layout.addRow("FASTA:", fasta_layout)
        
        # 注释文件
        gff_layout = QHBoxLayout()
        self.add_gff = QLineEdit()
        self.add_gff.setPlaceholderText("注释文件路径 (.gff, .gff3, .gtf) [可选]")
        gff_layout.addWidget(self.add_gff, 1)
        
        gff_btn = QPushButton("浏览...")
        gff_btn.clicked.connect(self._browse_gff)
        gff_layout.addWidget(gff_btn)
        form_layout.addRow("注释:", gff_layout)
        
        # 物种
        self.add_species = QLineEdit()
        self.add_species.setPlaceholderText("物种名称 [可选]")
        form_layout.addRow("物种:", self.add_species)
        
        # 描述
        self.add_desc = QLineEdit()
        self.add_desc.setPlaceholderText("描述信息 [可选]")
        form_layout.addRow("描述:", self.add_desc)
        
        layout.addWidget(form_group)
        
        # 添加按钮
        add_btn = QPushButton("添加基因组")
        add_btn.clicked.connect(self._add_custom_genome)
        layout.addWidget(add_btn)
        
        layout.addStretch()
        
        return widget
    
    def _load_genomes(self):
        """加载基因组列表"""
        db = get_database()
        genomes = db.list_genomes()
        
        self.genome_table.setRowCount(len(genomes))
        
        for i, genome in enumerate(genomes):
            self.genome_table.setItem(i, 0, QTableWidgetItem(genome.get("display_name") or genome["name"]))
            self.genome_table.setItem(i, 1, QTableWidgetItem(genome.get("source", "custom")))
            self.genome_table.setItem(i, 2, QTableWidgetItem(genome.get("fasta_path", "")))
            self.genome_table.setItem(i, 3, QTableWidgetItem(genome.get("annotation_path", "") or ""))
            self.genome_table.setItem(i, 4, QTableWidgetItem(genome.get("created_at", "")[:10] if genome.get("created_at") else ""))
            
            # 存储名称用于删除
            self.genome_table.item(i, 0).setData(Qt.UserRole, genome["name"])
        self.genome_table.resizeRowsToContents()
    
    def _delete_selected(self):
        """删除选中的基因组"""
        selected = self.genome_table.selectedItems()
        if not selected:
            QMessageBox.warning(self, "提示", "请先选择要删除的基因组")
            return
        
        row = selected[0].row()
        name = self.genome_table.item(row, 0).data(Qt.UserRole)
        
        reply = QMessageBox.question(
            self, "确认删除",
            f"确定要删除基因组 '{name}' 吗？\n\n注意：这只会删除数据库记录，不会删除文件。",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            db = get_database()
            db.delete_genome(name)
            self._load_genomes()
    
    def _search_online(self):
        """在线搜索基因组"""
        keyword = self.search_input.text().strip()
        if not keyword:
            return
        
        self.search_results.setText("搜索中...")
        
        manager = get_genome_manager()
        results = manager.search(keyword)
        
        output = []
        
        # IGV 结果
        igv_results = results.get("igv", [])
        if igv_results:
            output.append(f"=== IGV 基因组 ({len(igv_results)} 个) ===\n")
            for r in igv_results[:10]:
                output.append(f"  {r['id']}: {r['name']}")
            if len(igv_results) > 10:
                output.append(f"  ... 还有 {len(igv_results) - 10} 个结果")
            output.append("")
        
        # Ensembl Plants 结果
        ensembl_results = results.get("ensembl", [])
        if ensembl_results:
            output.append(f"=== Ensembl Plants ({len(ensembl_results)} 个) ===\n")
            for r in ensembl_results[:15]:
                output.append(f"  {r['id']}: {r['display_name']} ({r.get('common_name', '')}) - {r.get('assembly', '')}")
            if len(ensembl_results) > 15:
                output.append(f"  ... 还有 {len(ensembl_results) - 15} 个结果")
            output.append("")
        
        if not igv_results and not ensembl_results:
            output.append(f"未找到匹配 '{keyword}' 的基因组")
        else:
            output.append("\n提示: 复制上面的 ID 到下方输入框进行下载")
        
        self.search_results.setText("\n".join(output))
    
    def _download_genome(self):
        """下载基因组"""
        genome_id = self.download_input.text().strip()
        if not genome_id:
            QMessageBox.warning(self, "提示", "请输入基因组 ID")
            return
        
        self.download_progress.setText(f"正在下载 {genome_id}...")
        
        self.download_thread = DownloadThread(genome_id)
        self.download_thread.progress.connect(lambda msg: self.download_progress.setText(msg))
        self.download_thread.finished.connect(self._on_download_finished)
        self.download_thread.start()
    
    def _on_download_finished(self, success: bool, message: str):
        """下载完成回调"""
        self.download_progress.setText(message)
        if success:
            self._load_genomes()
            QMessageBox.information(self, "下载完成", message)
        else:
            QMessageBox.warning(self, "下载失败", message)
    
    def _browse_fasta(self):
        """浏览 FASTA 文件"""
        path, _ = QFileDialog.getOpenFileName(
            self, "选择 FASTA 文件", "",
            "FASTA 文件 (*.fa *.fasta *.fa.gz *.fasta.gz);;所有文件 (*)"
        )
        if path:
            self.add_fasta.setText(path)
    
    def _browse_gff(self):
        """浏览 GFF 文件"""
        path, _ = QFileDialog.getOpenFileName(
            self, "选择注释文件", "",
            "注释文件 (*.gff *.gff3 *.gtf *.gff.gz *.gff3.gz);;所有文件 (*)"
        )
        if path:
            self.add_gff.setText(path)
    
    def _add_custom_genome(self):
        """添加自定义基因组"""
        name = self.add_name.text().strip()
        fasta = self.add_fasta.text().strip()
        gff = self.add_gff.text().strip() or None
        species = self.add_species.text().strip() or None
        desc = self.add_desc.text().strip() or None
        
        if not name:
            QMessageBox.warning(self, "提示", "请输入基因组名称")
            return
        
        if not fasta:
            QMessageBox.warning(self, "提示", "请选择 FASTA 文件")
            return
        
        manager = get_genome_manager()
        success = manager.add_custom(name, fasta, gff, species, desc)
        
        if success:
            QMessageBox.information(self, "成功", f"已添加基因组: {name}")
            self._load_genomes()
            # 清空表单
            self.add_name.clear()
            self.add_fasta.clear()
            self.add_gff.clear()
            self.add_species.clear()
            self.add_desc.clear()
        else:
            QMessageBox.warning(self, "失败", "添加基因组失败，请检查文件路径")
