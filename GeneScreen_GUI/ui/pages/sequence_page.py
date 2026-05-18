"""
GeneScreen 1.0 - Sequence 模式页面

输入序列，与参考基因组进行比对
"""
from PySide6.QtWidgets import (
    QWidget, QHBoxLayout, QLabel,
    QPushButton, QTextEdit, QFormLayout,
    QSpinBox, QProgressBar, QMessageBox, QFileDialog,
    QAbstractSpinBox, QLineEdit, QCheckBox
)
from PySide6.QtCore import QThread, Signal, Qt
from datetime import datetime
import os
import re

from ui.widgets.genome_selector import GenomePairSelector
from ui.widgets.analysis_layout import (
    configure_pairwise_limit_controls,
    create_card,
    create_scroll_content,
)
from core import SequenceProcessor, get_database
from ui.widgets.report_worker import ReportWorker
from core.config import get_output_dir


class SequenceAnalysisThread(QThread):
    """分析线程 - 支持多序列批量处理"""
    progress = Signal(str)
    finished = Signal(bool, dict, str)
    item_finished = Signal(str, object, str)  # seq_id, result, error
    
    def __init__(self, processor, sequences: list, output_dir_builder=None, parent=None):
        super().__init__(parent)
        self.processor = processor
        self.sequences = sequences  # [{"seq_id": "xxx", "sequence": "ATGC..."}, ...]
        self.output_dir_builder = output_dir_builder
    
    def run(self):
        results = []
        for i, seq in enumerate(self.sequences):
            seq_id = seq["seq_id"]
            self.progress.emit(f"处理序列 {seq_id} ({i+1}/{len(self.sequences)})...")
            try:
                # 设置输出目录
                if self.output_dir_builder:
                    output_dir = self.output_dir_builder(seq)
                    if output_dir:
                        os.makedirs(output_dir, exist_ok=True)
                        self.processor.set_output_dir(output_dir)
                
                result = self.processor.process_sequence_text(seq["sequence"], seq_id)
                if result and "output_dir" not in result:
                    result["output_dir"] = getattr(self.processor, "output_dir", "")
                
                if result:
                    self.item_finished.emit(seq_id, result, "")
                    results.append(result)
                else:
                    self.item_finished.emit(seq_id, None, "未产生结果")
            except Exception as e:
                self.progress.emit(f"处理 {seq_id} 出错: {str(e)}")
                self.item_finished.emit(seq_id, None, str(e))
        
        if results:
            self.finished.emit(True, results[-1], f"完成 {len(results)}/{len(self.sequences)} 个序列")
        else:
            self.finished.emit(False, {}, "所有序列处理失败")


class SequencePage(QWidget):
    """Sequence 模式页面"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.history_id = None
        self._auto_output_dir = ""
        self._batch_mode = False
        self._history_map = {}
        self._output_dir_map = {}
        self._report_queue = []
        self._report_worker = None
        self._active_report_job = None
        self._pending_finish_message = ""
        self._last_report_path = ""
        self._init_ui()
    
    def _init_ui(self):
        """初始化 UI"""
        layout = create_scroll_content(self)
        
        # 标题
        title = QLabel("Sequence 模式")
        title.setProperty("role", "pageTitle")
        layout.addWidget(title)
        
        desc = QLabel("输入序列，与查询基因组进行比对，检测 SNP 和 Indel 变异")
        desc.setProperty("role", "pageDesc")
        layout.addWidget(desc)
        
        self.genome_selector = GenomePairSelector(
            show_manage_btn=False,
            multi_query=True,
            show_reference=False,
            use_cards=True,
            query_title="查询基因组"
        )
        layout.addWidget(self.genome_selector)
        
        # 输入区域
        input_group, input_layout = create_card("Seq 序列")
        
        # 序列输入（支持多序列 FASTA）
        seq_label = QLabel("序列 (支持多序列 FASTA 格式):")
        seq_label.setProperty("role", "muted")
        input_layout.addWidget(seq_label)
        
        # 从文件加载
        file_layout = QHBoxLayout()
        file_layout.addStretch()
        
        load_btn = QPushButton("📂 从文件加载")
        load_btn.setMinimumWidth(150)
        load_btn.setFixedHeight(36)
        load_btn.clicked.connect(self._load_from_file)
        file_layout.addWidget(load_btn)
        input_layout.addLayout(file_layout)

        self.sequence_input = QTextEdit()
        self.sequence_input.setPlaceholderText(
            ">seq1\nATGCATGCATGC...\n>seq2\nGCTAGCTAGCTA...\n\n"
            "或直接输入纯序列 (自动命名为 query_seq_1, query_seq_2...)"
        )
        self.sequence_input.setMinimumHeight(170)
        input_layout.addWidget(self.sequence_input)
        
        layout.addWidget(input_group)
        
        # 参数设置
        param_group, param_layout = create_card("参数设置", QFormLayout)
        param_group.setObjectName("paramGroup")
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

        self.upstream_input = QSpinBox()
        self.upstream_input.setRange(0, 1000000)
        self.upstream_input.setValue(0)
        self.upstream_input.setFixedHeight(spin_height)
        self.upstream_input.setProperty("paramInput", True)
        self.upstream_input.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.downstream_input = QSpinBox()
        self.downstream_input.setRange(0, 1000000)
        self.downstream_input.setValue(0)
        self.downstream_input.setFixedHeight(spin_height)
        self.downstream_input.setProperty("paramInput", True)
        self.downstream_input.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.candidate_limit_input = QSpinBox()
        self.candidate_limit_input.setRange(1, 1000)
        self.candidate_limit_input.setValue(3)
        self.candidate_limit_input.setFixedHeight(spin_height)
        self.candidate_limit_input.setProperty("paramInput", True)
        self.candidate_limit_input.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.pairwise_all_input = QCheckBox("全量候选")
        configure_pairwise_limit_controls(self.candidate_limit_input, self.pairwise_all_input)

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
        row1.addSpacing(20)
        row1.addWidget(make_param_label("查询上游延伸:"))
        row1.addWidget(self.upstream_input)
        row1.addWidget(make_param_label("bp"))
        row1.addSpacing(20)
        row1.addWidget(make_param_label("查询下游延伸:"))
        row1.addWidget(self.downstream_input)
        row1.addWidget(make_param_label("bp"))
        row1.addStretch()
        param_layout.addRow(row1)

        row2 = QHBoxLayout()
        row2.addWidget(make_param_label("Pairwise Top-N:"))
        row2.addWidget(self.candidate_limit_input)
        row2.addSpacing(20)
        row2.addWidget(self.pairwise_all_input)
        row2.addStretch()
        param_layout.addRow(row2)
        
        # 输出目录
        output_layout = QHBoxLayout()
        self.output_dir = QLineEdit()
        self.output_dir.setPlaceholderText("选择输出目录")
        self._set_default_output_dir()
        self.output_dir.setProperty("paramInput", True)
        self.output_dir.setFixedHeight(spin_height)
        output_layout.addWidget(self.output_dir, 1)
        
        browse_btn = QPushButton("浏览...")
        browse_btn.setProperty("secondary", True)
        browse_btn.setProperty("compactAction", True)
        browse_btn.setMinimumWidth(96)
        browse_btn.setFixedHeight(38)
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
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
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
        """从文件加载序列（支持多序列 FASTA）"""
        path, _ = QFileDialog.getOpenFileName(
            self, "选择序列文件", "",
            "FASTA 文件 (*.fa *.fasta *.fna);;文本文件 (*.txt);;所有文件 (*)"
        )
        if path:
            try:
                with open(path, 'r') as f:
                    content = f.read()
                self.sequence_input.setText(content)
            except Exception as e:
                QMessageBox.warning(self, "错误", f"读取文件失败: {str(e)}")

    def _parse_multi_fasta(self, text: str) -> list:
        """
        解析多序列 FASTA 文本
        返回: [{"seq_id": "xxx", "sequence": "ATGC..."}, ...]
        纯序列自动命名为 query_seq_1, query_seq_2...
        """
        sequences = []
        lines = text.strip().split('\n')
        current_id = None
        current_seq = []
        auto_idx = 0
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # 保存上一个序列
                if current_seq:
                    seq_str = ''.join(current_seq)
                    if seq_str:
                        if current_id is None:
                            auto_idx += 1
                            current_id = f"query_seq_{auto_idx}"
                        sequences.append({"seq_id": current_id, "sequence": seq_str})
                # 开始新序列
                header = line[1:].strip()
                current_id = header.split()[0] if header else None
                current_seq = []
            else:
                # 序列行
                current_seq.append(line)
        
        # 保存最后一个序列
        if current_seq:
            seq_str = ''.join(current_seq)
            if seq_str:
                if current_id is None:
                    auto_idx += 1
                    current_id = f"query_seq_{auto_idx}"
                sequences.append({"seq_id": current_id, "sequence": seq_str})
        
        return sequences
    
    def _run_analysis(self):
        """运行分析"""
        # 获取基因组
        qry_genome = self.genome_selector.get_qry_genome()
        qry_genomes = self.genome_selector.get_qry_genomes()
        
        if not qry_genome.get("fasta_path"):
            QMessageBox.warning(self, "提示", "请选择查询基因组")
            return
        
        # 解析序列
        seq_text = self.sequence_input.toPlainText().strip()
        if not seq_text:
            QMessageBox.warning(self, "提示", "请输入序列")
            return
        
        sequences = self._parse_multi_fasta(seq_text)
        if not sequences:
            QMessageBox.warning(self, "提示", "未检测到有效序列")
            return
        
        # 检查序列长度
        for seq in sequences:
            if len(seq["sequence"]) < 50:
                QMessageBox.warning(self, "提示", f"序列 {seq['seq_id']} 太短 (至少 50 bp)")
                return
        
        # 创建输出目录
        output_dir = self._resolve_output_dir()
        if not output_dir:
            QMessageBox.warning(self, "提示", "请选择输出目录")
            return
        
        os.makedirs(output_dir, exist_ok=True)
        multi_mode = len(sequences) > 1
        self._batch_mode = multi_mode
        
        # 创建处理器
        identity = self.identity_input.value()
        min_aln_len = self.min_aln_len_input.value()
        upstream = self.upstream_input.value()
        downstream = self.downstream_input.value()
        pairwise_all = self.pairwise_all_input.isChecked()
        candidate_limit = None if pairwise_all else self.candidate_limit_input.value()
        db = get_database()
        self._history_map = {}
        self._output_dir_map = {}
        self._last_report_path = ""
        
        # 为每个序列创建历史记录和输出目录
        for seq in sequences:
            seq_id = seq["seq_id"]
            safe_name = self._sanitize_path_segment(seq_id) or "sequence"
            item_output_dir = output_dir if not multi_mode else os.path.join(output_dir, safe_name)
            os.makedirs(item_output_dir, exist_ok=True)
            history_id = db.add_history(
                mode="sequence",
                ref_genome_id=None,
                qry_genome_id=qry_genome.get("id"),
                input_value=seq_id,
                identity=identity,
                output_dir=item_output_dir,
                status="running"
            )
            self._history_map[seq_id] = history_id
            self._output_dir_map[seq_id] = item_output_dir

        processor = SequenceProcessor(
            ref_genome=qry_genome["fasta_path"],
            output_dir=output_dir,
            ref_name=qry_genome.get("name", ""),
            identity=identity,
            ref_gff=qry_genome.get("annotation_path"),
            min_aln_len=min_aln_len,
            query_genomes=qry_genomes,
            upstream=upstream,
            downstream=downstream,
            candidate_limit=candidate_limit,
            pairwise_all=pairwise_all
        )
        
        # 启动分析线程
        self.run_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        output_dir_builder = None
        if multi_mode:
            def output_dir_builder(seq: dict) -> str:
                safe_name = self._sanitize_path_segment(seq["seq_id"]) or "sequence"
                return os.path.join(output_dir, safe_name)
        
        self.analysis_thread = SequenceAnalysisThread(processor, sequences, output_dir_builder=output_dir_builder)
        self.analysis_thread.progress.connect(lambda msg: self.progress_label.setText(msg))
        self.analysis_thread.item_finished.connect(self._on_item_finished)
        self.analysis_thread.finished.connect(self._on_analysis_finished)
        self.analysis_thread.start()
    
    def _on_item_finished(self, seq_id: str, result: object, error: str):
        """单个序列分析完成回调"""
        history_id = self._history_map.get(seq_id)
        if not history_id:
            return
        db = get_database()
        if result:
            output_dir = self._output_dir_map.get(seq_id, "") or result.get("output_dir", "")
            ref_name = self.genome_selector.get_qry_genome().get("name", "")
            job = {
                "history_id": history_id,
                "result": result,
                "output_dir": output_dir,
                "mode": "sequence",
                "ref_name": ref_name,
                "qry_name": ""
            }
            self._enqueue_report_job(job)
        else:
            db.update_history(history_id, status="failed")

    def _enqueue_report_job(self, job: dict):
        """添加报告生成任务到队列"""
        self._report_queue.append(job)
        self._process_report_queue()

    def _process_report_queue(self):
        """处理报告生成队列"""
        if self._report_worker and self._report_worker.isRunning():
            return
        if not self._report_queue:
            if not self._batch_mode and self._pending_finish_message and self._last_report_path:
                QMessageBox.information(
                    self, "分析完成",
                    f"{self._pending_finish_message}\n\n报告已保存到:\n{self._last_report_path}"
                )
                self._pending_finish_message = ""
            return
        self._active_report_job = self._report_queue.pop(0)
        job = self._active_report_job
        self._report_worker = ReportWorker(
            job["result"],
            job["output_dir"],
            job["mode"],
            job["ref_name"],
            job["qry_name"],
            self
        )
        self._report_worker.report_ready.connect(self._on_report_ready)
        self._report_worker.start()
    
    def _on_analysis_finished(self, success: bool, result: dict, message: str):
        """分析完成回调"""
        self.run_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.progress_label.setText(message)
        
        if success:
            if self._batch_mode:
                QMessageBox.information(self, "分析完成", f"{message}\n\n报告生成中...")
            else:
                if self._report_worker or self._report_queue:
                    self._pending_finish_message = message
                elif self._last_report_path:
                    QMessageBox.information(
                        self, "分析完成",
                        f"{message}\n\n报告已保存到:\n{self._last_report_path}"
                    )
                else:
                    QMessageBox.information(self, "分析完成", message)
        else:
            QMessageBox.warning(self, "分析失败", message)

    def _on_report_ready(self, success: bool, report_path: str, error: str):
        """报告生成完成回调"""
        job = self._active_report_job
        self._active_report_job = None
        if not job:
            return
        history_id = job.get("history_id")
        output_dir = job.get("output_dir", "")
        db = get_database()
        if success:
            db.update_history(
                history_id,
                status="completed",
                report_path=report_path,
                output_dir=output_dir
            )
            if not self._batch_mode:
                self._last_report_path = report_path
        else:
            db.update_history(history_id, status="failed")
        self._report_worker.deleteLater()
        self._report_worker = None
        self._process_report_queue()

    @staticmethod
    def _sanitize_path_segment(text: str) -> str:
        """清理路径中的非法字符"""
        return re.sub(r'[<>:"/\\\\|?*]', "_", text).strip()
