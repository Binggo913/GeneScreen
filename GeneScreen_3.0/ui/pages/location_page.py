"""
GeneScreen 3.0 - Location 模式页面

输入染色体坐标，提取指定区域并与目标基因组比对
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
import re

from ui.widgets.genome_selector import GenomePairSelector
from core import LocationProcessor, get_genome_manager, get_database
from ui.widgets.report_worker import ReportWorker
from core.config import get_output_dir


class LocationAnalysisThread(QThread):
    """分析线程"""
    progress = Signal(str)
    finished = Signal(bool, dict, str)
    item_finished = Signal(str, object, str)  # loc_key, result, error
    
    def __init__(self, processor, locations: list, output_dir_builder=None, parent=None):
        super().__init__(parent)
        self.processor = processor
        self.locations = locations
        self.output_dir_builder = output_dir_builder
    
    def run(self):
        results = []
        for i, loc in enumerate(self.locations):
            loc_label = loc.get("label") or f"{loc['chrom']}:{loc['start']}-{loc['end']}"
            self.progress.emit(f"处理 {loc_label} ({i+1}/{len(self.locations)})...")
            try:
                loc_key = loc_label
                loc_name = loc.get("file_tag") or f"{loc['chrom']}_{loc['start']}_{loc['end']}"
                if self.output_dir_builder:
                    output_dir = self.output_dir_builder(loc)
                    if output_dir:
                        os.makedirs(output_dir, exist_ok=True)
                        if hasattr(self.processor, "set_output_dir"):
                            self.processor.set_output_dir(output_dir)
                        else:
                            self.processor.output_dir = output_dir
                result = self.processor.process(loc['chrom'], loc['start'], loc['end'], loc_name)
                if result and "output_dir" not in result:
                    result["output_dir"] = getattr(self.processor, "output_dir", "")
                if result:
                    self.item_finished.emit(loc_key, result, "")
                else:
                    self.item_finished.emit(loc_key, None, "未产生结果")
                if result:
                    results.append(result)
            except Exception as e:
                self.progress.emit(f"处理出错: {str(e)}")
                self.item_finished.emit(loc_label, None, str(e))
        
        if results:
            self.finished.emit(True, results[-1], f"完成 {len(results)}/{len(self.locations)} 个区域")
        else:
            self.finished.emit(False, {}, "所有区域处理失败")


class LocationPage(QWidget):
    """Location 模式页面"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.history_id = None
        self._auto_output_dir = ""
        self._syncing_inputs = False
        self._batch_mode = False
        self._history_map = {}
        self._output_dir_map = {}
        self._last_report_path = ""
        self._report_queue = []
        self._report_worker = None
        self._active_report_job = None
        self._pending_finish_message = ""
        self._init_ui()
    
    def _init_ui(self):
        """初始化 UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 30, 30, 30)
        layout.setSpacing(20)
        
        # 标题
        title = QLabel("Location 模式")
        title.setProperty("role", "pageTitle")
        layout.addWidget(title)
        
        desc = QLabel("输入染色体坐标，提取指定区域并与目标基因组比对，检测 SNP 和 Indel 变异")
        desc.setProperty("role", "pageDesc")
        layout.addWidget(desc)
        
        # 基因组选择
        genome_group = QGroupBox("基因组选择")
        genome_layout = QVBoxLayout(genome_group)
        self.genome_selector = GenomePairSelector(show_manage_btn=False)
        genome_layout.addWidget(self.genome_selector)
        layout.addWidget(genome_group)
        
        # 输入区域
        input_group = QGroupBox("输入")
        input_layout = QVBoxLayout(input_group)
        
        # 单个位置输入
        loc_layout = QHBoxLayout()
        loc_label = QLabel("位置:")
        loc_label.setMinimumWidth(80)
        loc_layout.addWidget(loc_label)
        
        self.location_input = QLineEdit()
        self.location_input.setPlaceholderText("输入位置 (如: Chr1:1000000-1050000)")
        self.location_input.textChanged.connect(self._on_single_location_changed)
        loc_layout.addWidget(self.location_input, 1)
        input_layout.addLayout(loc_layout)
        
        # 批量输入
        input_layout.addSpacing(10)
        batch_label = QLabel("或批量输入 (每行一个位置，格式: Chr1:100-200):")
        batch_label.setProperty("role", "muted")
        input_layout.addWidget(batch_label)
        
        self.batch_input = QTextEdit()
        self.batch_input.setPlaceholderText("Chr1:1000000-1050000\nChr2:2000000-2100000\n...")
        self.batch_input.setMaximumHeight(100)
        self.batch_input.textChanged.connect(self._on_batch_locations_changed)
        input_layout.addWidget(self.batch_input)
        
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
        return os.path.join(str(base), "Location", timestamp)

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

    def _on_single_location_changed(self, text: str):
        if self._syncing_inputs:
            return
        if text.strip():
            self._clear_batch_input()

    def _on_batch_locations_changed(self):
        if self._syncing_inputs:
            return
        if self.batch_input.toPlainText().strip():
            self._clear_single_input()

    def _clear_batch_input(self):
        if not self.batch_input.toPlainText().strip():
            return
        self._syncing_inputs = True
        self.batch_input.blockSignals(True)
        self.batch_input.clear()
        self.batch_input.blockSignals(False)
        self._syncing_inputs = False

    def _clear_single_input(self):
        if not self.location_input.text().strip():
            return
        self._syncing_inputs = True
        self.location_input.blockSignals(True)
        self.location_input.clear()
        self.location_input.blockSignals(False)
        self._syncing_inputs = False
    
    def _parse_location(self, text: str) -> dict:
        """解析位置字符串"""
        # 格式1: Chr1:1000000-1050000
        match = re.match(r'(\w+):(\d+)-(\d+)', text.strip())
        if match:
            chrom = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))
            return {
                'chrom': chrom,
                'start': start,
                'end': end,
                'label': f"{chrom}:{start}-{end}",
                'file_tag': f"{chrom}_{start}_{end}",
                'name': None
            }
        return None
    
    def _parse_batch_line(self, line: str) -> dict:
        """解析批量输入行"""
        return self._parse_location(line)
    
    def _run_analysis(self):
        """运行分析"""
        # 获取基因组
        ref_genome = self.genome_selector.get_ref_genome()
        qry_genome = self.genome_selector.get_qry_genome()
        
        if not ref_genome.get("fasta_path"):
            QMessageBox.warning(self, "提示", "请选择参考基因组")
            return
        
        if not qry_genome.get("fasta_path"):
            QMessageBox.warning(self, "提示", "请选择查询基因组")
            return
        
        # 获取位置
        locations = []
        single_loc = self.location_input.text().strip()
        if single_loc:
            loc = self._parse_location(single_loc)
            if loc:
                locations.append(loc)
        
        batch_text = self.batch_input.toPlainText().strip()
        if batch_text:
            for line in batch_text.split('\n'):
                if line.strip():
                    loc = self._parse_batch_line(line)
                    if loc:
                        locations.append(loc)
        
        if not locations:
            QMessageBox.warning(self, "提示", "请输入至少一个有效位置")
            return
        
        # 创建输出目录
        output_dir = self._resolve_output_dir()
        if not output_dir:
            QMessageBox.warning(self, "提示", "请选择输出目录")
            return
        
        os.makedirs(output_dir, exist_ok=True)
        multi_mode = len(locations) > 1
        self._batch_mode = multi_mode
        
        # 创建处理器
        identity = self.identity_input.value()
        min_aln_len = self.min_aln_len_input.value()
        db = get_database()
        self._history_map = {}
        self._output_dir_map = {}
        self._last_report_path = ""
        for loc in locations:
            loc_label = loc.get("label") or f"{loc['chrom']}:{loc['start']}-{loc['end']}"
            loc_tag = loc.get("file_tag") or f"{loc['chrom']}_{loc['start']}_{loc['end']}"
            safe_name = self._sanitize_path_segment(loc_tag) or "location"
            item_output_dir = output_dir if not multi_mode else f"{output_dir}_{safe_name}"
            os.makedirs(item_output_dir, exist_ok=True)
            history_id = db.add_history(
                mode="location",
                ref_genome_id=ref_genome.get("id"),
                qry_genome_id=qry_genome.get("id"),
                input_value=loc_label,
                identity=identity,
                output_dir=item_output_dir,
                status="running"
            )
            self._history_map[loc_label] = history_id
            self._output_dir_map[loc_label] = item_output_dir

        processor = LocationProcessor(
            ref_genome=ref_genome["fasta_path"],
            query_genome=qry_genome["fasta_path"],
            output_dir=output_dir,
            ref_name=ref_genome.get("name", ""),
            qry_name=qry_genome.get("name", ""),
            identity=identity,
            ref_gff=ref_genome.get("annotation_path"),
            qry_gff=qry_genome.get("annotation_path"),
            min_aln_len=min_aln_len
        )
        
        # 启动分析线程
        self.run_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        output_dir_builder = None
        if multi_mode:
            def output_dir_builder(loc: dict) -> str:
                name = loc.get("file_tag") or f"{loc['chrom']}_{loc['start']}_{loc['end']}"
                safe_name = self._sanitize_path_segment(name) or "location"
                return f"{output_dir}_{safe_name}"
        self.analysis_thread = LocationAnalysisThread(processor, locations, output_dir_builder=output_dir_builder)
        self.analysis_thread.progress.connect(lambda msg: self.progress_label.setText(msg))
        self.analysis_thread.item_finished.connect(self._on_item_finished)
        self.analysis_thread.finished.connect(self._on_analysis_finished)
        self.analysis_thread.start()
    
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

    def _on_item_finished(self, loc_key: str, result: object, error: str):
        history_id = self._history_map.get(loc_key)
        if not history_id:
            return
        db = get_database()
        if result:
            output_dir = self._output_dir_map.get(loc_key, "") or result.get("output_dir", "")
            ref_name = self.genome_selector.get_ref_genome().get("name", "")
            qry_name = self.genome_selector.get_qry_genome().get("name", "")
            job = {
                "history_id": history_id,
                "result": result,
                "output_dir": output_dir,
                "mode": "location",
                "ref_name": ref_name,
                "qry_name": qry_name
            }
            self._enqueue_report_job(job)
        else:
            db.update_history(history_id, status="failed")

    def _enqueue_report_job(self, job: dict):
        self._report_queue.append(job)
        self._process_report_queue()

    def _process_report_queue(self):
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

    def _on_report_ready(self, success: bool, report_path: str, error: str):
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
        return re.sub(r'[<>:"/\\\\|?*]', "_", text).strip()
