"""
GeneScreen 1.0 - Gene ID 模式页面

输入基因 ID，从参考基因组提取序列并与目标基因组比对
"""
from PySide6.QtWidgets import (
    QWidget, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QTextEdit, QFormLayout,
    QSpinBox, QProgressBar, QMessageBox, QFileDialog,
    QComboBox, QCompleter, QAbstractSpinBox, QListWidget, QListWidgetItem,
    QCheckBox, QSizePolicy
)
from PySide6.QtCore import (
    QThread, Signal, Qt, QTimer, QAbstractListModel,
    QModelIndex, QStringListModel, QEvent
)
from datetime import datetime
import os
import re

from ui.widgets.genome_selector import GenomePairSelector
from ui.widgets.analysis_layout import create_card, create_scroll_content
from core import GeneIDProcessor, get_database, get_genome_manager
from ui.widgets.report_worker import ReportWorker
from core.gene_id_utils import load_gene_ids
from core.config import get_output_dir


class AnalysisThread(QThread):
    """分析线程"""
    progress = Signal(str)
    finished = Signal(bool, dict, str)  # success, result, message
    item_finished = Signal(str, object, str)  # gene_id, result, error
    
    def __init__(self, processor, gene_ids: list, output_dir_builder=None, parent=None):
        super().__init__(parent)
        self.processor = processor
        self.gene_ids = gene_ids
        self.output_dir_builder = output_dir_builder
    
    def run(self):
        results = []
        for i, gene_id in enumerate(self.gene_ids):
            self.progress.emit(f"处理 {gene_id} ({i+1}/{len(self.gene_ids)})...")
            try:
                if self.output_dir_builder:
                    output_dir = self.output_dir_builder(gene_id)
                    if output_dir:
                        os.makedirs(output_dir, exist_ok=True)
                        if hasattr(self.processor, "set_output_dir"):
                            self.processor.set_output_dir(output_dir)
                        else:
                            self.processor.output_dir = output_dir
                result = self.processor.process(gene_id)
                if result and "output_dir" not in result:
                    result["output_dir"] = getattr(self.processor, "output_dir", "")
                if result:
                    self.item_finished.emit(gene_id, result, "")
                else:
                    self.item_finished.emit(gene_id, None, "未产生结果")
                if result:
                    results.append(result)
            except Exception as e:
                self.progress.emit(f"处理 {gene_id} 出错: {str(e)}")
                self.item_finished.emit(gene_id, None, str(e))
        
        if results:
            self.finished.emit(True, results[-1], f"完成 {len(results)}/{len(self.gene_ids)} 个基因")
        else:
            self.finished.emit(False, {}, "所有基因处理失败")


class GeneIdLoadThread(QThread):
    """从注释文件加载 Gene ID 列表"""
    finished = Signal(str, list, str)  # annotation_path, ids, error

    def __init__(self, annotation_path: str, parent=None):
        super().__init__(parent)
        self.annotation_path = annotation_path

    def run(self):
        try:
            ids = load_gene_ids(self.annotation_path)
            self.finished.emit(self.annotation_path, ids, "")
        except Exception as e:
            self.finished.emit(self.annotation_path, [], str(e))


class LazyGeneIdModel(QAbstractListModel):
    """延迟加载 Gene ID，避免一次性填充造成卡顿。"""

    def __init__(self, items=None, chunk_size: int = 800, initial_chunk: int = 200, parent=None):
        super().__init__(parent)
        self._items = list(items or [])
        self._chunk_size = max(50, chunk_size)
        self._initial_chunk = max(1, initial_chunk)
        self._loaded = min(len(self._items), self._initial_chunk)

    def rowCount(self, parent=QModelIndex()):
        if parent.isValid():
            return 0
        return self._loaded

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid() or not (0 <= index.row() < self._loaded):
            return None
        if role in (Qt.DisplayRole, Qt.EditRole):
            return self._items[index.row()]
        return None

    def canFetchMore(self, parent=QModelIndex()):
        if parent.isValid():
            return False
        return self._loaded < len(self._items)

    def fetchMore(self, parent=QModelIndex()):
        if parent.isValid() or self._loaded >= len(self._items):
            return
        remaining = len(self._items) - self._loaded
        to_load = min(self._chunk_size, remaining)
        start = self._loaded
        end = self._loaded + to_load - 1
        self.beginInsertRows(QModelIndex(), start, end)
        self._loaded += to_load
        self.endInsertRows()

    def set_items(self, items):
        self.beginResetModel()
        self._items = list(items or [])
        self._loaded = min(len(self._items), self._initial_chunk)
        self.endResetModel()


class GeneIDPage(QWidget):
    """Gene ID 模式页面"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.history_id = None
        self._gene_id_cache = {}
        self._gene_id_threads = []
        self._pending_gene_id_path = ""
        self._gene_id_poll_timer = None
        self._gene_id_poll_path = ""
        self._suppress_popup = False
        self._auto_output_dir = ""
        self._gene_id_popup_model = QStringListModel([], self)
        self._gene_id_completer = None
        self._gene_id_all_ids = []
        self._gene_id_all_ids_lower = []
        self._gene_id_all_ids_set = set()
        self._last_filter_text = ""
        self._syncing_inputs = False
        self._batch_mode = False
        self._history_map = {}
        self._output_dir_map = {}
        self._last_report_path = ""
        self._report_queue = []
        self._report_worker = None
        self._active_report_job = None
        self._pending_finish_message = ""
        # 防抖定时器
        self._debounce_timer = None
        self._init_ui()
    
    def _init_ui(self):
        """初始化 UI"""
        layout = create_scroll_content(self)
        
        # 标题
        title = QLabel("Gene ID 模式")
        title.setProperty("role", "pageTitle")
        layout.addWidget(title)
        
        desc = QLabel("输入基因 ID，从参考基因组提取序列并与目标基因组比对，检测 SNP 和 Indel 变异")
        desc.setProperty("role", "pageDesc")
        layout.addWidget(desc)
        
        self.genome_selector = GenomePairSelector(show_manage_btn=False, show_ref_annotation=True,
                                                   multi_query=True,
                                                   ref_annotation_label="参考注释版本",
                                                   use_cards=True,
                                                   ref_title="参考",
                                                   query_title="查询")
        self.genome_selector.ref_changed.connect(self._on_ref_genome_changed)
        self.genome_selector.ref_annotation_changed.connect(self._on_ref_annotation_changed)
        layout.addWidget(self.genome_selector)
        
        # 输入区域
        input_group, input_layout = create_card("Gene ID")
        
        # Gene ID 输入（用 QLineEdit + 独立弹窗列表）
        id_layout = QHBoxLayout()
        id_label = QLabel("Gene ID:")
        id_label.setMinimumWidth(80)
        id_layout.addWidget(id_label)
        
        self.gene_id_input = QLineEdit()
        self.gene_id_input.setMinimumHeight(36)
        self.gene_id_input.setProperty("paramInput", True)
        self.gene_id_input.setPlaceholderText("请选择参考基因组的gene id")
        self.gene_id_input.textEdited.connect(self._on_gene_id_text_edited)
        self.gene_id_input.textChanged.connect(self._on_single_gene_id_changed)
        self.gene_id_input.installEventFilter(self)
        id_layout.addWidget(self.gene_id_input, 1)
        input_layout.addLayout(id_layout)
        
        # 独立的选择列表（作为子控件，不是独立窗口）
        self._popup_list = QListWidget(self)
        self._popup_list.setMaximumHeight(300)
        self._popup_list.setFocusPolicy(Qt.NoFocus)
        self._popup_list.itemClicked.connect(self._on_popup_item_clicked)
        self._popup_list.hide()
        
        self._set_gene_id_loading_state("请选择参考基因组的gene id")
        
        # 批量输入
        input_layout.addSpacing(10)
        
        batch_header = QHBoxLayout()
        batch_label = QLabel("或批量输入 (每行一个 ID):")
        batch_label.setProperty("role", "muted")
        batch_header.addWidget(batch_label)
        batch_header.addStretch()

        import_gene_ids_btn = QPushButton("从文件导入")
        import_gene_ids_btn.setProperty("secondary", True)
        import_gene_ids_btn.setProperty("compactAction", True)
        import_gene_ids_btn.setMinimumWidth(132)
        import_gene_ids_btn.setFixedHeight(38)
        import_gene_ids_btn.clicked.connect(self._import_gene_ids_from_file)
        batch_header.addWidget(import_gene_ids_btn)
        
        # Gene list 路径显示（可点击打开）
        self.gene_list_path_label = QLabel("")
        self.gene_list_path_label.setProperty("role", "muted")
        self.gene_list_path_label.setStyleSheet("color: #667eea;")
        self.gene_list_path_label.setCursor(Qt.PointingHandCursor)
        self.gene_list_path_label.mousePressEvent = self._open_gene_list_file
        self.gene_list_path_label.hide()
        batch_header.addWidget(self.gene_list_path_label)
        
        input_layout.addLayout(batch_header)
        
        self.batch_input = QTextEdit()
        self.batch_input.setPlaceholderText("每行一个 ID")
        self.batch_input.setFixedHeight(self._batch_gene_id_input_height(2))
        self.batch_input.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.batch_input.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.batch_input.textChanged.connect(self._on_batch_gene_ids_changed)
        input_layout.addWidget(self.batch_input)
        
        layout.addWidget(input_group)
        
        # 参数设置
        param_group, param_layout = create_card("参数设置", QFormLayout)
        param_group.setObjectName("paramGroup")
        param_layout.setLabelAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        
        spin_height = 30
        spin_width = 100

        # Identity 阈值 + 最小比对长度
        self.identity_input = QSpinBox()
        self.identity_input.setRange(70, 100)
        self.identity_input.setValue(90)
        self.identity_input.setFixedHeight(spin_height)
        self.identity_input.setProperty("paramInput", True)
        self.identity_input.setFixedWidth(spin_width)
        self.identity_input.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.min_aln_len_input = QSpinBox()
        self.min_aln_len_input.setRange(1, 1000000)
        self.min_aln_len_input.setValue(100)
        self.min_aln_len_input.setFixedHeight(spin_height)
        self.min_aln_len_input.setProperty("paramInput", True)
        self.min_aln_len_input.setFixedWidth(spin_width)
        self.min_aln_len_input.setButtonSymbols(QAbstractSpinBox.NoButtons)

        # 上下游延伸
        self.upstream_input = QSpinBox()
        self.upstream_input.setRange(0, 1000000)
        self.upstream_input.setValue(0)
        self.upstream_input.setFixedHeight(spin_height)
        self.upstream_input.setProperty("paramInput", True)
        self.upstream_input.setFixedWidth(spin_width)
        self.upstream_input.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.downstream_input = QSpinBox()
        self.downstream_input.setRange(0, 1000000)
        self.downstream_input.setValue(0)
        self.downstream_input.setFixedHeight(spin_height)
        self.downstream_input.setProperty("paramInput", True)
        self.downstream_input.setFixedWidth(spin_width)
        self.downstream_input.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.candidate_limit_input = QSpinBox()
        self.candidate_limit_input.setRange(1, 1000)
        self.candidate_limit_input.setValue(3)
        self.candidate_limit_input.setFixedHeight(spin_height)
        self.candidate_limit_input.setProperty("paramInput", True)
        self.candidate_limit_input.setFixedWidth(spin_width)
        self.candidate_limit_input.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.pairwise_all_input = QCheckBox("全量候选")

        def make_param_label(text: str) -> QLabel:
            label = QLabel(text)
            label.setFixedHeight(spin_height)
            label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
            return label

        row1 = QHBoxLayout()
        row1.addWidget(make_param_label("Identity 阈值:"))
        row1.addWidget(self.identity_input)
        row1.addWidget(make_param_label("%"))
        row1.addSpacing(16)
        row1.addWidget(make_param_label("最小比对长度:"))
        row1.addWidget(self.min_aln_len_input)
        row1.addWidget(make_param_label("bp"))
        row1.addSpacing(16)
        row1.addWidget(make_param_label("上游延伸:"))
        row1.addWidget(self.upstream_input)
        row1.addWidget(make_param_label("bp"))
        row1.addSpacing(16)
        row1.addWidget(make_param_label("下游延伸:"))
        row1.addWidget(self.downstream_input)
        row1.addWidget(make_param_label("bp"))
        row1.addStretch()
        param_layout.addRow(row1)

        row2 = QHBoxLayout()
        row2.addWidget(make_param_label("Pairwise Top-N:"))
        row2.addWidget(self.candidate_limit_input)
        row2.addSpacing(16)
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
        browse_btn.clicked.connect(self._browse_output)
        output_layout.addWidget(browse_btn)
        param_layout.addRow("输出目录:", output_layout)
        
        layout.addWidget(param_group)
        
        # 运行按钮
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        
        self.run_btn = QPushButton("🚀 开始分析")
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
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        
        layout.addStretch()

    def _batch_gene_id_input_height(self, line_count: int) -> int:
        return self.batch_input.fontMetrics().lineSpacing() * max(2, line_count) + 30

    def _resize_batch_gene_id_input(self):
        if not hasattr(self, "batch_input"):
            return
        line_count = self.batch_input.toPlainText().count("\n") + 1
        self.batch_input.setFixedHeight(self._batch_gene_id_input_height(line_count))
        self.batch_input.updateGeometry()

    def _on_ref_genome_changed(self, name: str, genome: dict):
        """参考基因组改变时，等待注释版本选择"""
        self._stop_gene_id_poll()
        self._current_ref_name = name
        if not name:
            self._set_gene_id_loading_state("请选择参考基因组的gene id", disable=False)
            self._update_gene_list_path_label("")
            return
        # 注释版本会通过 _on_ref_annotation_changed 触发加载
    
    def _on_ref_annotation_changed(self, source: str, ann: dict):
        """参考基因组注释版本改变时，加载对应的 gene id 列表"""
        self._stop_gene_id_poll()
        
        if not source or not ann:
            self._set_gene_id_loading_state("该基因组无注释文件，无法使用 Gene ID 模式", disable=True)
            self._update_gene_list_path_label("")
            return
        
        gene_ids_path = ann.get("gene_ids_path", "")
        ann_path = ann.get("annotation_path", "")
        
        # 更新 gene list 路径显示
        if gene_ids_path:
            self._update_gene_list_path_label(gene_ids_path)
        else:
            # 预测路径
            ref_name = getattr(self, '_current_ref_name', '')
            if ref_name:
                manager = get_genome_manager()
                predicted_path = manager.get_gene_id_list_path(ref_name, source)
                self._update_gene_list_path_label(predicted_path)
        
        if gene_ids_path and os.path.exists(gene_ids_path):
            if gene_ids_path in self._gene_id_cache:
                self._apply_gene_id_list(self._gene_id_cache[gene_ids_path])
                return
            self._load_gene_ids_async(gene_ids_path)
            return
        
        if not ann_path:
            self._set_gene_id_loading_state("注释文件不存在", disable=False)
            return
        
        # 触发 gene list 生成
        ref_name = getattr(self, '_current_ref_name', '')
        if ref_name:
            manager = get_genome_manager()
            manager.enqueue_gene_id_list(ref_name, ann_path, source)
            self._gene_id_poll_path = manager.get_gene_id_list_path(ref_name, source)
            self._set_gene_id_loading_state("Gene ID 列表生成中，可直接输入", disable=False)
            self._start_gene_id_poll()
    
    def _update_gene_list_path_label(self, path: str):
        """更新 gene list 路径显示"""
        if hasattr(self, 'gene_list_path_label'):
            if path:
                self.gene_list_path_label.setText(f"📄 {path}")
                self.gene_list_path_label.setToolTip(path)
                self.gene_list_path_label.show()
            else:
                self.gene_list_path_label.hide()

    def _set_gene_id_loading_state(self, placeholder: str, disable: bool = True):
        self.gene_id_input.setEnabled(not disable)
        self.gene_id_input.setReadOnly(False)
        self.gene_id_input.blockSignals(True)
        self.gene_id_input.clear()
        self.gene_id_input.blockSignals(False)
        self._gene_id_all_ids_set = set()
        self.gene_id_input.setPlaceholderText(placeholder)

    def _load_gene_ids_async(self, annotation_path: str):
        self._pending_gene_id_path = annotation_path
        self._set_gene_id_loading_state("正在加载基因 ID...", disable=True)
        loader = GeneIdLoadThread(annotation_path, self)
        loader.finished.connect(self._on_gene_ids_loaded)
        loader.finished.connect(lambda *_ , thread=loader: self._cleanup_gene_id_thread(thread))
        self._gene_id_threads.append(loader)
        loader.start()

    def _on_gene_ids_loaded(self, annotation_path: str, ids: list, error: str):
        if annotation_path != self._pending_gene_id_path:
            return
        if error:
            self._set_gene_id_loading_state("基因 ID 加载失败", disable=False)
            return
        self._gene_id_cache[annotation_path] = ids
        self._apply_gene_id_list(ids)

    def _apply_gene_id_list(self, ids: list):
        self.gene_id_input.setEnabled(True)
        self.gene_id_input.setReadOnly(False)
        self.gene_id_input.blockSignals(True)
        self.gene_id_input.clear()
        self.gene_id_input.blockSignals(False)
        self.gene_id_input.setPlaceholderText("请选择参考基因组的gene id")
        self._gene_id_all_ids = list(ids)
        self._gene_id_all_ids_lower = [gid.lower() for gid in self._gene_id_all_ids]
        self._gene_id_all_ids_set = set(self._gene_id_all_ids)
        self._last_filter_text = ""
        self._update_batch_placeholder(ids)

    def _on_gene_id_text_edited(self, text: str):
        # 输入时实时更新并显示匹配项
        self._show_gene_id_popup()

    def _update_popup_debounced(self):
        pass

    def _on_gene_id_selected(self, text: str):
        pass

    def _on_popup_item_clicked(self, item):
        """点击弹窗列表项"""
        self.gene_id_input.setText(item.text())
        self._popup_list.hide()
        self.gene_id_input.setFocus()

    def _show_gene_id_popup(self):
        """显示选择框"""
        if not self._gene_id_all_ids:
            self._popup_list.hide()
            return
        
        current_text = self.gene_id_input.text().strip()
        
        # 计算匹配项
        if not current_text:
            matches = self._gene_id_all_ids[:20]
        else:
            needle = current_text.lower()
            matches = []
            for gid, gid_lower in zip(self._gene_id_all_ids, self._gene_id_all_ids_lower):
                if needle in gid_lower:
                    matches.append(gid)
                    if len(matches) >= 20:
                        break
        
        if not matches:
            self._popup_list.hide()
            return
        
        # 更新列表
        self._popup_list.clear()
        self._popup_list.addItems(matches)
        
        # 动态调整高度（每项约25px，最大300px）
        item_height = 25
        content_height = len(matches) * item_height + 4
        actual_height = min(content_height, 300)
        self._popup_list.setFixedHeight(actual_height)
        
        # 定位到输入框下方（相对于父控件）
        pos = self.gene_id_input.mapTo(self, self.gene_id_input.rect().bottomLeft())
        self._popup_list.setFixedWidth(self.gene_id_input.width())
        self._popup_list.move(pos)
        self._popup_list.raise_()  # 置顶
        self._popup_list.show()

    def _update_popup_content_only(self):
        pass

    def _update_popup_items(self):
        pass

    def eventFilter(self, obj, event):
        if obj == self.gene_id_input:
            if event.type() == QEvent.MouseButtonPress:
                QTimer.singleShot(50, self._show_gene_id_popup)
            elif event.type() == QEvent.FocusOut:
                # 延迟隐藏，让点击事件先处理
                QTimer.singleShot(200, self._hide_popup_if_no_focus)
        return super().eventFilter(obj, event)

    def _hide_popup_if_no_focus(self):
        if not self.gene_id_input.hasFocus() and not self._popup_list.hasFocus():
            self._popup_list.hide()

    def _clear_popup_suppress(self) -> None:
        self._suppress_popup = False

    def _on_single_gene_id_changed(self, text: str):
        if self._syncing_inputs:
            return
        if text.strip():
            self._clear_batch_input()
            # 校验 gene id 是否在当前注释版本的列表中
            self._validate_single_gene_id(text.strip())
    
    def _validate_single_gene_id(self, gene_id: str):
        """校验单个 gene id 是否有效"""
        if not gene_id:
            self.gene_id_input.setStyleSheet("")
            return
        if not self._gene_id_all_ids_set:
            # 列表未加载，不校验
            self.gene_id_input.setStyleSheet("")
            return
        if gene_id in self._gene_id_all_ids_set:
            # 有效
            self.gene_id_input.setStyleSheet("")
        else:
            # 无效，标红
            self.gene_id_input.setStyleSheet("border: 1px solid #e74c3c;")

    def _on_batch_gene_ids_changed(self):
        self._resize_batch_gene_id_input()
        if self._syncing_inputs:
            return
        if self.batch_input.toPlainText().strip():
            self._clear_single_input()
            # 校验批量 gene id
            self._validate_batch_gene_ids()
        else:
            self.batch_input.setStyleSheet("")
            self.batch_input.setToolTip("")

    def _clear_batch_input(self):
        if not self.batch_input.toPlainText().strip():
            return
        self._syncing_inputs = True
        self.batch_input.blockSignals(True)
        self.batch_input.clear()
        self.batch_input.blockSignals(False)
        self._syncing_inputs = False
        self._resize_batch_gene_id_input()
        self.batch_input.setStyleSheet("")  # 清除校验样式
    
    def _validate_batch_gene_ids(self):
        """校验批量 gene id"""
        text = self.batch_input.toPlainText().strip()
        if not text:
            self.batch_input.setStyleSheet("")
            return
        if not self._gene_id_all_ids_set:
            # 列表未加载，不校验
            self.batch_input.setStyleSheet("")
            return
        
        ids = [line.strip() for line in text.split('\n') if line.strip()]
        invalid_ids = [gid for gid in ids if gid not in self._gene_id_all_ids_set]
        
        if invalid_ids:
            # 有无效 ID，标黄
            self.batch_input.setStyleSheet("border: 1px solid #f39c12;")
            self.batch_input.setToolTip(f"以下 {len(invalid_ids)} 个 ID 未找到:\n" + "\n".join(invalid_ids[:10]))
        else:
            self.batch_input.setStyleSheet("")
            self.batch_input.setToolTip("")

    def _clear_single_input(self):
        if not self.gene_id_input.text().strip():
            return
        self._syncing_inputs = True
        self.gene_id_input.blockSignals(True)
        self.gene_id_input.clear()
        self.gene_id_input.blockSignals(False)
        self._syncing_inputs = False

    def _update_batch_placeholder(self, ids: list):
        if not ids:
            self.batch_input.setPlaceholderText("每行一个 ID")
            return
        first = ids[0]
        second = ids[1] if len(ids) > 1 else ""
        lines = [first]
        if second:
            lines.append(second)
        lines.append("...")
        self.batch_input.setPlaceholderText("\n".join(lines))

    def _import_gene_ids_from_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "选择 Gene ID 文件",
            "",
            "文本文件 (*.txt *.list *.tsv *.csv);;所有文件 (*)"
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8-sig") as handle:
                gene_ids = [line.strip() for line in handle if line.strip()]
        except Exception as exc:
            QMessageBox.warning(self, "错误", f"读取 Gene ID 文件失败: {exc}")
            return
        if not gene_ids:
            QMessageBox.warning(self, "提示", "文件中未读取到 Gene ID")
            return
        self.batch_input.setPlainText("\n".join(gene_ids))
        self._resize_batch_gene_id_input()

    def _update_completer_matches(self, text: str):
        if text == self._last_filter_text:
            return
        self._last_filter_text = text
        if not self._gene_id_all_ids:
            self._gene_id_popup_model.setStringList([])
            return
        needle = text.strip().lower()
        if not needle:
            matches = self._gene_id_all_ids[:20]
        else:
            matches = []
            for gid, gid_lower in zip(self._gene_id_all_ids, self._gene_id_all_ids_lower):
                if needle in gid_lower:
                    matches.append(gid)
                    if len(matches) >= 20:
                        break

    def _cleanup_gene_id_thread(self, thread: QThread):
        if thread in self._gene_id_threads:
            self._gene_id_threads.remove(thread)
        thread.deleteLater()

    def _start_gene_id_poll(self):
        if not self._gene_id_poll_path:
            return
        if self._gene_id_poll_timer is None:
            self._gene_id_poll_timer = QTimer(self)
            self._gene_id_poll_timer.timeout.connect(self._check_gene_id_list_ready)
        if not self._gene_id_poll_timer.isActive():
            self._gene_id_poll_timer.start(800)

    def _stop_gene_id_poll(self):
        if self._gene_id_poll_timer and self._gene_id_poll_timer.isActive():
            self._gene_id_poll_timer.stop()

    def _check_gene_id_list_ready(self):
        if not self._gene_id_poll_path:
            return
        if os.path.exists(self._gene_id_poll_path) and os.path.getsize(self._gene_id_poll_path) > 0:
            path = self._gene_id_poll_path
            self._stop_gene_id_poll()
            if path in self._gene_id_cache:
                self._apply_gene_id_list(self._gene_id_cache[path])
            else:
                self._load_gene_ids_async(path)

    def closeEvent(self, event):
        self._stop_gene_id_poll()
        for thread in list(self._gene_id_threads):
            thread.wait(2000)
        super().closeEvent(event)
    
    def _browse_output(self):
        """浏览输出目录"""
        path = QFileDialog.getExistingDirectory(self, "选择输出目录")
        if path:
            self.output_dir.setText(path)

    def _default_output_dir(self) -> str:
        base = get_output_dir()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return os.path.join(str(base), "Gene_ID", timestamp)

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
        
        if not ref_genome.get("annotation_path"):
            QMessageBox.warning(self, "提示", "参考基因组缺少注释文件 (GFF)")
            return
        
        # 获取 Gene ID
        gene_ids = []
        single_id = self.gene_id_input.text().strip()
        if single_id:
            gene_ids.append(single_id)
        
        batch_text = self.batch_input.toPlainText().strip()
        if batch_text:
            gene_ids.extend([line.strip() for line in batch_text.split('\n') if line.strip()])
        
        if not gene_ids:
            QMessageBox.warning(self, "提示", "请输入至少一个 Gene ID")
            return
        
        # 批量校验：检查无效 gene id
        if self._gene_id_all_ids_set and len(gene_ids) > 1:
            invalid_ids = [gid for gid in gene_ids if gid not in self._gene_id_all_ids_set]
            if invalid_ids:
                msg = f"以下 {len(invalid_ids)} 个 Gene ID 未在当前注释版本中找到：\n"
                msg += "\n".join(invalid_ids[:10])
                if len(invalid_ids) > 10:
                    msg += f"\n... 还有 {len(invalid_ids) - 10} 个"
                msg += "\n\n是否继续？"
                reply = QMessageBox.question(self, "确认", msg, 
                                             QMessageBox.Yes | QMessageBox.No)
                if reply != QMessageBox.Yes:
                    return
        
        # 创建输出目录
        output_dir = self._resolve_output_dir()
        if not output_dir:
            QMessageBox.warning(self, "提示", "请选择输出目录")
            return
        
        os.makedirs(output_dir, exist_ok=True)
        multi_mode = len(gene_ids) > 1
        self._batch_mode = multi_mode
        
        # 创建处理器
        identity = self.identity_input.value()
        upstream = self.upstream_input.value()
        downstream = self.downstream_input.value()
        min_aln_len = self.min_aln_len_input.value()
        candidate_limit = self.candidate_limit_input.value()
        pairwise_all = self.pairwise_all_input.isChecked()
        db = get_database()
        self._history_map = {}
        self._output_dir_map = {}
        self._last_report_path = ""
        for gene_id in gene_ids:
            safe_id = self._sanitize_path_segment(gene_id) or "gene"
            item_output_dir = output_dir if not multi_mode else os.path.join(output_dir, safe_id)
            os.makedirs(item_output_dir, exist_ok=True)
            history_id = db.add_history(
                mode="gene_id",
                ref_genome_id=ref_genome.get("id"),
                qry_genome_id=qry_genome.get("id"),
                input_value=gene_id,
                identity=identity,
                output_dir=item_output_dir,
                status="running"
            )
            self._history_map[gene_id] = history_id
            self._output_dir_map[gene_id] = item_output_dir

        processor = GeneIDProcessor(
            ref_genome=ref_genome["fasta_path"],
            ref_annotation=ref_genome["annotation_path"],
            query_genome=qry_genome["fasta_path"],
            output_dir=output_dir,
            ref_name=ref_genome.get("name", ""),
            qry_name=qry_genome.get("name", ""),
            identity=identity,
            qry_gff=qry_genome.get("annotation_path"),
            upstream=upstream,
            downstream=downstream,
            min_aln_len=min_aln_len,
            query_genomes=qry_genomes,
            candidate_limit=candidate_limit,
            pairwise_all=pairwise_all
        )
        
        # 启动分析线程
        self.run_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # 不确定进度
        
        output_dir_builder = None
        if multi_mode:
            def output_dir_builder(gene_id: str) -> str:
                return self._output_dir_map.get(gene_id, output_dir)
        self.analysis_thread = AnalysisThread(processor, gene_ids, output_dir_builder=output_dir_builder)
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

    def _on_item_finished(self, gene_id: str, result: object, error: str):
        history_id = self._history_map.get(gene_id)
        if not history_id:
            return
        db = get_database()
        if result:
            output_dir = self._output_dir_map.get(gene_id, "") or result.get("output_dir", "")
            ref_name = self.genome_selector.get_ref_genome().get("name", "")
            qry_name = self.genome_selector.get_qry_genome().get("name", "")
            job = {
                "history_id": history_id,
                "result": result,
                "output_dir": output_dir,
                "mode": "gene_id",
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
    
    def _open_gene_list_file(self, event):
        """打开 gene list 文件"""
        path = self.gene_list_path_label.toolTip()
        if path and os.path.exists(path):
            import subprocess
            import sys
            if sys.platform == 'win32':
                os.startfile(path)
            elif sys.platform == 'darwin':
                subprocess.run(['open', path])
            else:
                subprocess.run(['xdg-open', path])
