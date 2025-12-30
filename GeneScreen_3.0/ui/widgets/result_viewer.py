"""
GeneScreen 3.0 - 结果查看器组件

显示分析结果：统计信息、变异列表、可视化图
"""
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QHeaderView, QGroupBox,
    QTabWidget, QTextBrowser, QFileDialog, QMessageBox
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QDesktopServices
from PySide6.QtCore import QUrl
import os


class ResultViewer(QWidget):
    """结果查看器"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._init_ui()
        self.result = None
    
    def _init_ui(self):
        """初始化 UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # 标签页
        self.tabs = QTabWidget()
        
        # 统计信息
        self.tabs.addTab(self._create_stats_tab(), "📊 统计")
        
        # 变异列表
        self.tabs.addTab(self._create_variants_tab(), "🧬 变异")
        
        # 文件列表
        self.tabs.addTab(self._create_files_tab(), "📁 文件")
        
        layout.addWidget(self.tabs)
    
    def _create_stats_tab(self) -> QWidget:
        """创建统计标签页"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # 统计卡片
        cards_layout = QHBoxLayout()
        
        self.stat_cards = {}
        for key, label in [
            ("seq_length", "序列长度"),
            ("identity", "平均 Identity"),
            ("snp_count", "SNP"),
            ("ins_count", "插入"),
            ("del_count", "缺失")
        ]:
            card = self._create_stat_card(label, "0")
            self.stat_cards[key] = card
            cards_layout.addWidget(card)
        
        layout.addLayout(cards_layout)
        layout.addStretch()
        
        return widget
    
    def _create_stat_card(self, label: str, value: str) -> QGroupBox:
        """创建统计卡片"""
        card = QGroupBox()
        card.setStyleSheet("""
            QGroupBox {
                background: white;
                border: 1px solid #e0e0e0;
                border-radius: 8px;
                padding: 15px;
            }
        """)
        
        layout = QVBoxLayout(card)
        layout.setAlignment(Qt.AlignCenter)
        
        value_label = QLabel(value)
        value_label.setObjectName("statValue")
        value_label.setStyleSheet("font-size: 28px; font-weight: bold; color: #667eea;")
        value_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(value_label)
        
        name_label = QLabel(label)
        name_label.setStyleSheet("font-size: 12px; color: #666;")
        name_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(name_label)
        
        card.value_label = value_label
        return card
    
    def _create_variants_tab(self) -> QWidget:
        """创建变异标签页"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # 变异表格
        self.variants_table = QTableWidget()
        self.variants_table.setColumnCount(6)
        self.variants_table.setHorizontalHeaderLabels([
            "类型", "参考位置", "查询位置", "参考碱基", "变异碱基", "长度"
        ])
        self.variants_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.variants_table.setAlternatingRowColors(True)
        self.variants_table.setSelectionBehavior(QTableWidget.SelectRows)
        layout.addWidget(self.variants_table)
        
        return widget
    
    def _create_files_tab(self) -> QWidget:
        """创建文件标签页"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # 文件表格
        self.files_table = QTableWidget()
        self.files_table.setColumnCount(3)
        self.files_table.setHorizontalHeaderLabels(["文件名", "路径", "操作"])
        self.files_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.files_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.files_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.files_table.setAlternatingRowColors(True)
        layout.addWidget(self.files_table)
        
        # 打开目录按钮
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        
        self.open_dir_btn = QPushButton("📂 打开输出目录")
        self.open_dir_btn.clicked.connect(self._open_output_dir)
        btn_layout.addWidget(self.open_dir_btn)
        
        layout.addLayout(btn_layout)
        
        return widget
    
    def set_result(self, result: dict):
        """设置结果数据"""
        self.result = result
        self._update_stats()
        self._update_variants()
        self._update_files()
    
    def _update_stats(self):
        """更新统计信息"""
        if not self.result:
            return
        
        # 读取 coords 文件计算统计
        coords_file = self.result.get("coords", "")
        snps_file = self.result.get("snps", "")
        fasta_file = self.result.get("fasta", "")
        
        # 序列长度
        seq_length = 0
        if fasta_file and os.path.exists(fasta_file):
            with open(fasta_file, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        seq_length += len(line.strip())
        self.stat_cards["seq_length"].value_label.setText(f"{seq_length:,}")
        
        # Identity
        identity = 0.0
        if coords_file and os.path.exists(coords_file):
            identities = []
            with open(coords_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('[') or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 7:
                        try:
                            identities.append(float(parts[6]))
                        except ValueError:
                            pass
            if identities:
                identity = sum(identities) / len(identities)
        self.stat_cards["identity"].value_label.setText(f"{identity:.1f}%")
        
        # 变异统计
        snp_count = 0
        ins_count = 0
        del_count = 0
        if snps_file and os.path.exists(snps_file):
            with open(snps_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('[') or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        var_type = parts[4].upper()
                        if var_type == 'SNP':
                            snp_count += 1
                        elif var_type == 'INS':
                            ins_count += 1
                        elif var_type == 'DEL':
                            del_count += 1
        
        self.stat_cards["snp_count"].value_label.setText(str(snp_count))
        self.stat_cards["ins_count"].value_label.setText(str(ins_count))
        self.stat_cards["del_count"].value_label.setText(str(del_count))
    
    def _update_variants(self):
        """更新变异列表"""
        self.variants_table.setRowCount(0)
        
        if not self.result:
            return
        
        snps_file = self.result.get("snps", "")
        if not snps_file or not os.path.exists(snps_file):
            return
        
        variants = []
        with open(snps_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('[') or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 5:
                    variants.append({
                        'ref_pos': parts[0],
                        'ref_base': parts[1],
                        'alt_base': parts[2],
                        'qry_pos': parts[3],
                        'type': parts[4].upper()
                    })
        
        self.variants_table.setRowCount(len(variants))
        for i, v in enumerate(variants):
            self.variants_table.setItem(i, 0, QTableWidgetItem(v['type']))
            self.variants_table.setItem(i, 1, QTableWidgetItem(v['ref_pos']))
            self.variants_table.setItem(i, 2, QTableWidgetItem(v['qry_pos']))
            self.variants_table.setItem(i, 3, QTableWidgetItem(v['ref_base']))
            self.variants_table.setItem(i, 4, QTableWidgetItem(v['alt_base']))
            
            # 计算长度
            if v['type'] == 'SNP':
                length = "1"
            elif v['type'] == 'INS':
                length = str(len(v['alt_base']))
            elif v['type'] == 'DEL':
                length = str(len(v['ref_base']))
            else:
                length = "-"
            self.variants_table.setItem(i, 5, QTableWidgetItem(length))
    
    def _update_files(self):
        """更新文件列表"""
        self.files_table.setRowCount(0)
        
        if not self.result:
            return
        
        files = [
            ("FASTA", self.result.get("fasta", "")),
            ("Coords", self.result.get("coords", "")),
            ("SNPs", self.result.get("snps", "")),
            ("BLAST XML", self.result.get("blast_xml", "")),
            ("GFF", self.result.get("gff", "")),
        ]
        
        row = 0
        for name, path in files:
            if path and os.path.exists(path):
                self.files_table.insertRow(row)
                self.files_table.setItem(row, 0, QTableWidgetItem(name))
                self.files_table.setItem(row, 1, QTableWidgetItem(path))
                
                open_btn = QPushButton("打开")
                open_btn.clicked.connect(lambda checked, p=path: self._open_file(p))
                self.files_table.setCellWidget(row, 2, open_btn)
                
                row += 1
    
    def _open_file(self, path: str):
        """打开文件"""
        QDesktopServices.openUrl(QUrl.fromLocalFile(path))
    
    def _open_output_dir(self):
        """打开输出目录"""
        if self.result:
            # 从任意文件路径获取目录
            for key in ["fasta", "coords", "snps"]:
                path = self.result.get(key, "")
                if path and os.path.exists(path):
                    dir_path = os.path.dirname(path)
                    QDesktopServices.openUrl(QUrl.fromLocalFile(dir_path))
                    return
        
        QMessageBox.warning(self, "提示", "无法找到输出目录")
