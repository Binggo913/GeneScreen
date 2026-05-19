"""
GeneScreen 1.0 - Location 模式页面

输入染色体坐标，提取指定区域并与目标基因组比对
"""
from PySide6.QtWidgets import (
    QWidget, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QTextEdit, QFormLayout,
    QSpinBox, QMessageBox, QFileDialog,
    QAbstractSpinBox, QCheckBox
)
from PySide6.QtCore import QThread, Signal, Qt
from datetime import datetime
import os
import re

from ui.widgets.genome_selector import GenomePairSelector
from ui.widgets.analysis_layout import (
    configure_auto_growing_text_edit,
    configure_pairwise_limit_controls,
    create_card,
    create_scroll_content,
)
from core import LocationProcessor, get_genome_manager, get_database
from core.task_manager import AnalysisTask, get_analysis_task_manager
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
        layout = create_scroll_content(self)
        
        # 标题
        title = QLabel("Location 模式")
        title.setProperty("role", "pageTitle")
        layout.addWidget(title)
        
        desc = QLabel("输入染色体坐标，提取指定区域并与目标基因组比对，检测 SNP 和 Indel 变异")
        desc.setProperty("role", "pageDesc")
        layout.addWidget(desc)
        
        self.genome_selector = GenomePairSelector(
            show_manage_btn=False,
            multi_query=True,
            use_cards=True,
            ref_title="参考",
            query_title="查询"
        )
        layout.addWidget(self.genome_selector)
        
        # 输入区域
        input_group, input_layout = create_card("参考上的位置")
        
        # 位置输入说明
        loc_label = QLabel("输入位置 (每行一个，格式: Chr1:1000000-1050000):")
        loc_label.setProperty("role", "muted")
        input_layout.addWidget(loc_label)
        
        # 批量输入
        self.batch_input = QTextEdit()
        self.batch_input.setPlaceholderText("Chr1:1000000-1050000\nChr2:2000000-2100000\n...")
        configure_auto_growing_text_edit(self.batch_input, min_lines=2)
        input_layout.addWidget(self.batch_input)
        
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
        
        self.run_btn = QPushButton("🚀 提交分析")
        self.run_btn.setProperty("primaryAction", True)
        self.run_btn.setMinimumWidth(150)
        self.run_btn.setMinimumHeight(45)
        self.run_btn.clicked.connect(self._run_analysis)
        btn_layout.addWidget(self.run_btn)
        
        layout.addLayout(btn_layout)
        
        # 进度
        self.progress_label = QLabel("")
        self.progress_label.setProperty("role", "muted")
        layout.addWidget(self.progress_label)
        
        layout.addStretch()
    
    def _browse_output(self):
        """浏览输出目录"""
        path = QFileDialog.getExistingDirectory(self, "选择输出目录")
        if path:
            self.output_dir.setText(path)

    def _default_output_dir(self) -> str:
        base = get_output_dir()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
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
        qry_genomes = self.genome_selector.get_qry_genomes()
        
        if not ref_genome.get("fasta_path"):
            QMessageBox.warning(self, "提示", "请选择参考基因组")
            return
        
        if not qry_genome.get("fasta_path"):
            QMessageBox.warning(self, "提示", "请选择查询基因组")
            return
        
        # 获取位置
        locations = []
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
        upstream = self.upstream_input.value()
        downstream = self.downstream_input.value()
        pairwise_all = self.pairwise_all_input.isChecked()
        candidate_limit = None if pairwise_all else self.candidate_limit_input.value()
        db = get_database()
        self._history_map = {}
        self._output_dir_map = {}
        self._last_report_path = ""
        for loc in locations:
            loc_label = loc.get("label") or f"{loc['chrom']}:{loc['start']}-{loc['end']}"
            loc_tag = loc.get("file_tag") or f"{loc['chrom']}_{loc['start']}_{loc['end']}"
            safe_name = self._sanitize_path_segment(loc_tag) or "location"
            item_output_dir = output_dir if not multi_mode else os.path.join(output_dir, safe_name)
            os.makedirs(item_output_dir, exist_ok=True)
            history_id = db.add_history(
                mode="location",
                ref_genome_id=ref_genome.get("id"),
                qry_genome_id=qry_genome.get("id"),
                input_value=loc_label,
                identity=identity,
                output_dir=item_output_dir,
                status="pending"
            )
            self._history_map[loc_label] = history_id
            self._output_dir_map[loc_label] = item_output_dir
        history_map = dict(self._history_map)
        output_dir_map = dict(self._output_dir_map)

        processor = LocationProcessor(
            ref_genome=ref_genome["fasta_path"],
            query_genome=qry_genome["fasta_path"],
            output_dir=output_dir,
            ref_name=ref_genome.get("name", ""),
            qry_name=qry_genome.get("name", ""),
            identity=identity,
            ref_gff=ref_genome.get("annotation_path"),
            qry_gff=qry_genome.get("annotation_path"),
            min_aln_len=min_aln_len,
            upstream=upstream,
            downstream=downstream,
            query_genomes=qry_genomes,
            candidate_limit=candidate_limit,
            pairwise_all=pairwise_all
        )
        
        # 提交到后台任务队列
        self.progress_label.setText("已提交后台队列")
        
        output_dir_builder = None
        if multi_mode:
            def output_dir_builder(loc: dict) -> str:
                name = loc.get("file_tag") or f"{loc['chrom']}_{loc['start']}_{loc['end']}"
                safe_name = self._sanitize_path_segment(name) or "location"
                loc_label = loc.get("label") or f"{loc['chrom']}:{loc['start']}-{loc['end']}"
                return output_dir_map.get(loc_label, os.path.join(output_dir, safe_name))
        self.analysis_thread = LocationAnalysisThread(processor, locations, output_dir_builder=output_dir_builder)
        self.analysis_thread.history_map = history_map
        self.analysis_thread.output_dir_map = output_dir_map
        self.analysis_thread.batch_mode = multi_mode
        self.analysis_thread.ref_name = ref_genome.get("name", "")
        self.analysis_thread.qry_name = qry_genome.get("name", "")
        self.analysis_thread.progress.connect(lambda msg: self.progress_label.setText(msg))
        self.analysis_thread.item_finished.connect(self._on_item_finished)
        self.analysis_thread.finished.connect(self._on_analysis_finished)
        history_ids = list(history_map.values())

        def on_started():
            self.progress_label.setText("后台任务执行中...")

        get_analysis_task_manager().submit(
            AnalysisTask(
                thread=self.analysis_thread,
                history_ids=history_ids,
                label=f"Location: {len(locations)} item(s)",
                on_started=on_started,
            )
        )
        QMessageBox.information(self, "提交成功", "分析任务已提交到后台队列，可在历史记录中查看状态。")
    
    def _on_analysis_finished(self, success: bool, result: dict, message: str):
        """分析完成回调"""
        self.progress_label.setText(message)
        
        if not success:
            QMessageBox.warning(self, "分析失败", message)

    def _on_item_finished(self, loc_key: str, result: object, error: str):
        sender = self.sender()
        history_map = getattr(sender, "history_map", self._history_map)
        output_dir_map = getattr(sender, "output_dir_map", self._output_dir_map)
        history_id = history_map.get(loc_key)
        if not history_id:
            return
        db = get_database()
        if result:
            output_dir = output_dir_map.get(loc_key, "") or result.get("output_dir", "")
            ref_name = getattr(sender, "ref_name", "") or self.genome_selector.get_ref_genome().get("name", "")
            qry_name = getattr(sender, "qry_name", "") or self.genome_selector.get_qry_genome().get("name", "")
            job = {
                "history_id": history_id,
                "result": result,
                "output_dir": output_dir,
                "mode": "location",
                "ref_name": ref_name,
                "qry_name": qry_name,
                "batch_mode": getattr(sender, "batch_mode", self._batch_mode),
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
            if not job.get("batch_mode", self._batch_mode):
                self._last_report_path = report_path
        else:
            db.update_history(history_id, status="failed")
        self._report_worker.deleteLater()
        self._report_worker = None
        self._process_report_queue()

    @staticmethod
    def _sanitize_path_segment(text: str) -> str:
        return re.sub(r'[<>:"/\\\\|?*]', "_", text).strip()
