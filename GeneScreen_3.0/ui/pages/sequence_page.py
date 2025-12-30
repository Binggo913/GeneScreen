"""
GeneScreen 3.0 - Sequence 模式页面

输入序列，与参考基因组进行比对
"""
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QTextEdit, QGroupBox, QFormLayout,
    QSpinBox, QProgressBar, QMessageBox, QFileDialog,
    QAbstractSpinBox
)
from PySide6.QtCore import QThread, Signal, Qt
from datetime import datetime
import os
from pathlib import Path

from ui.widgets.genome_selector import GenomeSelector
from core import SequenceProcessor, get_database
from ui.widgets.report_worker import ReportWorker
from core.config import get_output_dir


class SequenceAnalysisThread(QThread):
    """分析线程"""
    progress = Signal(str)
    finished = Signal(bool, dict, str)
    
    def __init__(self, processor, sequence: str, seq_id: str, parent=None):
        super().__init__(parent)
        self.processor = processor
        self.sequence = sequence
        self.seq_id = seq_id
    
    def run(self):
        self.progress.emit(f"处理序列 {self.seq_id}...")
        try:
            result = self.processor.process_sequence_text(self.sequence, self.seq_id)
            if result:
                self.finished.emit(True, result, "分析完成")
            else:
                self.finished.emit(False, {}, "分析失败")
        except Exception as e:
            self.finished.emit(False, {}, f"分析出错: {str(e)}")


class SequencePage(QWidget):
    """Sequence 模式页面"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.history_id = None
        self._default_seq_id = "query_seq"
        self._seq_id_manual = False
        self._last_auto_seq_id = ""
        self._auto_output_dir = ""
        self._report_worker = None
        self._pending_finish_message = ""
        self._last_report_path = ""
        self._init_ui()
    
    def _init_ui(self):
        """初始化 UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 30, 30, 30)
        layout.setSpacing(20)
        
        # 标题
        title = QLabel("Sequence 模式")
        title.setProperty("role", "pageTitle")
        layout.addWidget(title)
        
        desc = QLabel("输入序列，与参考基因组进行比对，检测 SNP 和 Indel 变异")
        desc.setProperty("role", "pageDesc")
        layout.addWidget(desc)
        
        # 基因组选择（只需要参考基因组）
        genome_group = QGroupBox("参考基因组")
        genome_layout = QVBoxLayout(genome_group)
        self.genome_selector = GenomeSelector("参考基因组", show_manage_btn=False)
        genome_layout.addWidget(self.genome_selector)
        layout.addWidget(genome_group)
        
        # 输入区域
        input_group = QGroupBox("输入序列")
        input_layout = QVBoxLayout(input_group)
        
        # 序列 ID
        id_layout = QHBoxLayout()
        id_label = QLabel("序列 ID:")
        id_label.setMinimumWidth(80)
        id_layout.addWidget(id_label)
        
        self.seq_id_input = QLineEdit()
        self.seq_id_input.setPlaceholderText("输入序列 ID (如: my_sequence)")
        self.seq_id_input.setText(self._default_seq_id)
        self.seq_id_input.textEdited.connect(self._on_seq_id_edited)
        id_layout.addWidget(self.seq_id_input, 1)
        input_layout.addLayout(id_layout)
        
        # 序列输入
        input_layout.addSpacing(10)
        seq_label = QLabel("序列 (FASTA 格式或纯序列):")
        seq_label.setProperty("role", "muted")
        input_layout.addWidget(seq_label)
        
        self.sequence_input = QTextEdit()
        self.sequence_input.setPlaceholderText(">my_sequence\nATGCATGCATGC...\n\n或直接输入序列:\nATGCATGCATGC...")
        self.sequence_input.setMinimumHeight(150)
        self.sequence_input.textChanged.connect(self._on_sequence_text_changed)
        input_layout.addWidget(self.sequence_input)
        
        # 从文件加载
        file_layout = QHBoxLayout()
        file_layout.addStretch()
        
        load_btn = QPushButton("📂 从文件加载")
        load_btn.clicked.connect(self._load_from_file)
        file_layout.addWidget(load_btn)
        input_layout.addLayout(file_layout)
        
        layout.addWidget(input_group)
        
        # 参数设置
        param_group = QGroupBox("参数设置")
        param_group.setObjectName("paramGroup")
        param_layout = QFormLayout(param_group)
        param_layout.setLabelAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        
        spin_height = 30

        # Identity 阈值 + 最小比对长度
        self.identity_input = QSpinBox()
        self.identity_input.setRange(70, 100)
        self.identity_input.setValue(90)
        self.identity_input.setFixedHeight(spin_height)
        self.identity_input.setProperty("paramInput", True)
        self.identity_input.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.min_aln_len_input = QSpinBox()
        self.min_aln_len_input.setRange(1, 1000000)
        self.min_aln_len_input.setValue(100)
        self.min_aln_len_input.setFixedHeight(spin_height)
        self.min_aln_len_input.setProperty("paramInput", True)
        self.min_aln_len_input.setButtonSymbols(QAbstractSpinBox.NoButtons)

        def make_param_label(text: str) -> QLabel:
            label = QLabel(text)
            label.setFixedHeight(spin_height)
            label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
            return label

        row1 = QHBoxLayout()
        row1.addWidget(make_param_label("Identity 阈值:"))
        row1.addWidget(self.identity_input)
        row1.addWidget(make_param_label("%"))
        row1.addSpacing(20)
        row1.addWidget(make_param_label("最小比对长度:"))
        row1.addWidget(self.min_aln_len_input)
        row1.addWidget(make_param_label("bp"))
        row1.addStretch()
        param_layout.addRow(row1)
        
        # 输出目录
        output_layout = QHBoxLayout()
        self.output_dir = QLineEdit()
        self.output_dir.setPlaceholderText("选择输出目录")
        self._set_default_output_dir()
        self.output_dir.setProperty("paramInput", True)
        self.output_dir.setFixedHeight(spin_height)
        output_layout.addWidget(self.output_dir, 1)
        
        browse_btn = QPushButton("浏览...")
        browse_btn.clicked.connect(self._browse_output)
        output_layout.addWidget(browse_btn)
        param_layout.addRow("输出目录:", output_layout)
        
        layout.addWidget(param_group)
        
        # 运行按钮
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        
        self.run_btn = QPushButton("🚀 开始分析")
        self.run_btn.setMinimumWidth(150)
        self.run_btn.setMinimumHeight(45)
        self.run_btn.clicked.connect(self._run_analysis)
        btn_layout.addWidget(self.run_btn)
        
        layout.addLayout(btn_layout)
        
        # 进度
        self.progress_label = QLabel("")
        self.progress_label.setProperty("role", "muted")
        layout.addWidget(self.progress_label)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        
        layout.addStretch()
    
    def _browse_output(self):
        """浏览输出目录"""
        path = QFileDialog.getExistingDirectory(self, "选择输出目录")
        if path:
            self.output_dir.setText(path)

    def _default_output_dir(self) -> str:
        base = get_output_dir()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M")
        return os.path.join(str(base), "Sequence", timestamp)

    def _set_default_output_dir(self):
        self._auto_output_dir = self._default_output_dir()
        self.output_dir.setText(self._auto_output_dir)

    def _resolve_output_dir(self) -> str:
        output_dir = self.output_dir.text().strip()
        if not output_dir or output_dir == self._auto_output_dir:
            self._set_default_output_dir()
            output_dir = self._auto_output_dir
        return output_dir

    def showEvent(self, event):
        super().showEvent(event)
        current = self.output_dir.text().strip()
        if not current or current == self._auto_output_dir:
            self._set_default_output_dir()
    
    def _load_from_file(self):
        """从文件加载序列"""
        path, _ = QFileDialog.getOpenFileName(
            self, "选择序列文件", "",
            "FASTA 文件 (*.fa *.fasta *.fna);;文本文件 (*.txt);;所有文件 (*)"
        )
        if path:
            try:
                with open(path, 'r') as f:
                    content = f.read()
                self.sequence_input.setText(content)
                self._update_seq_id_from_text(content, force=True)
            except Exception as e:
                QMessageBox.warning(self, "错误", f"读取文件失败: {str(e)}")

    def _on_seq_id_edited(self, text: str) -> None:
        if text.strip():
            if text.strip() != self._last_auto_seq_id:
                self._seq_id_manual = True
        else:
            self._seq_id_manual = False

    def _on_sequence_text_changed(self) -> None:
        text = self.sequence_input.toPlainText()
        self._update_seq_id_from_text(text)

    def _extract_seq_id_from_text(self, text: str) -> str:
        for line in text.splitlines():
            line = line.strip()
            if not line.startswith('>'):
                continue
            header = line[1:].strip()
            if not header:
                continue
            return header.split()[0]
        return ""

    def _update_seq_id_from_text(self, text: str, force: bool = False) -> None:
        seq_id = self._extract_seq_id_from_text(text)
        if not seq_id:
            if not self.seq_id_input.text().strip():
                self.seq_id_input.setText(self._default_seq_id)
            return
        current = self.seq_id_input.text().strip()
        if force or not self._seq_id_manual or current in ("", self._default_seq_id, self._last_auto_seq_id):
            self.seq_id_input.setText(seq_id)
            self._last_auto_seq_id = seq_id
            self._seq_id_manual = False
    
    def _parse_sequence(self, text: str) -> str:
        """解析序列文本，返回纯序列"""
        lines = text.strip().split('\n')
        sequence = ""
        for line in lines:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line
        return sequence
    
    def _run_analysis(self):
        """运行分析"""
        # 获取基因组
        ref_genome = self.genome_selector.get_selected_genome()
        
        if not ref_genome.get("fasta_path"):
            QMessageBox.warning(self, "提示", "请选择参考基因组")
            return
        
        # 获取序列
        seq_text = self.sequence_input.toPlainText().strip()
        if not seq_text:
            QMessageBox.warning(self, "提示", "请输入序列")
            return
        
        sequence = self._parse_sequence(seq_text)
        if len(sequence) < 50:
            QMessageBox.warning(self, "提示", "序列太短 (至少 50 bp)")
            return
        
        self._update_seq_id_from_text(seq_text)
        seq_id = self.seq_id_input.text().strip() or self._default_seq_id
        
        # 创建输出目录
        output_dir = self._resolve_output_dir()
        if not output_dir:
            QMessageBox.warning(self, "提示", "请选择输出目录")
            return
        
        os.makedirs(output_dir, exist_ok=True)
        
        # 创建处理器
        identity = self.identity_input.value()
        min_aln_len = self.min_aln_len_input.value()
        db = get_database()
        self.history_id = db.add_history(
            mode="sequence",
            ref_genome_id=ref_genome.get("id"),
            qry_genome_id=None,
            input_value=seq_id,
            identity=identity,
            output_dir=output_dir,
            status="running"
        )

        processor = SequenceProcessor(
            ref_genome=ref_genome["fasta_path"],
            output_dir=output_dir,
            ref_name=ref_genome.get("name", ""),
            identity=identity,
            ref_gff=ref_genome.get("annotation_path"),
            min_aln_len=min_aln_len
        )
        
        # 启动分析线程
        self.run_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        self.analysis_thread = SequenceAnalysisThread(processor, sequence, seq_id)
        self.analysis_thread.progress.connect(lambda msg: self.progress_label.setText(msg))
        self.analysis_thread.finished.connect(self._on_analysis_finished)
        self.analysis_thread.start()
    
    def _on_analysis_finished(self, success: bool, result: dict, message: str):
        """分析完成回调"""
        self.run_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.progress_label.setText(message)
        
        if success and result:
            output_dir = self.output_dir.text().strip()
            ref_name = self.genome_selector.get_selected_genome().get("name", "")
            self._pending_finish_message = message
            self._report_worker = ReportWorker(result, output_dir, "sequence", ref_name, "", self)
            self._report_worker.report_ready.connect(self._on_report_ready)
            self._report_worker.start()
        else:
            if self.history_id:
                db = get_database()
                db.update_history(self.history_id, status="failed")
                self.history_id = None
            QMessageBox.warning(self, "分析失败", message)

    def _on_report_ready(self, success: bool, report_path: str, error: str):
        if self._report_worker:
            self._report_worker.deleteLater()
            self._report_worker = None
        if success:
            if self.history_id:
                db = get_database()
                db.update_history(
                    self.history_id,
                    status="completed",
                    report_path=report_path,
                    input_value=Path(report_path).stem
                )
                self.history_id = None
            self._last_report_path = report_path
            if self._pending_finish_message:
                QMessageBox.information(
                    self, "分析完成",
                    f"{self._pending_finish_message}\n\n报告已保存到:\n{report_path}"
                )
                self._pending_finish_message = ""
        else:
            if self.history_id:
                db = get_database()
                db.update_history(self.history_id, status="failed")
                self.history_id = None
            QMessageBox.warning(self, "分析失败", error or "报告生成失败")
