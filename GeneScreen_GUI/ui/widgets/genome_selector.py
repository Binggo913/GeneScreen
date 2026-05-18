"""
GeneScreen 1.0 - 基因组选择器组件

用于选择参考基因组和查询基因组的下拉组件
"""
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QComboBox,
    QLabel, QPushButton, QMessageBox, QListWidget, QListWidgetItem,
    QSizePolicy
)
from PySide6.QtCore import Signal, Qt

from core import get_genome_manager, get_database

from .analysis_layout import create_card


class GenomeSelector(QWidget):
    """基因组选择器（含注释版本选择）"""
    
    # 信号：基因组选择改变
    genome_changed = Signal(str, dict)  # (genome_name, genome_info)
    # 信号：注释版本选择改变
    annotation_changed = Signal(str, dict)  # (source, annotation_info)
    
    def __init__(self, label: str = "基因组", show_manage_btn: bool = True, 
                 show_annotation: bool = False, annotation_label: str = "注释版本", parent=None):
        super().__init__(parent)
        self.label_text = label
        self.show_manage_btn = show_manage_btn
        self.show_annotation = show_annotation
        self.annotation_label_text = annotation_label
        self._db = get_database()
        self._init_ui()
        self._load_genomes()
        self._db.add_genome_listener(self._on_genomes_changed)
    
    def _init_ui(self):
        """初始化 UI"""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(8)
        
        # 基因组选择行
        genome_layout = QHBoxLayout()
        genome_layout.setContentsMargins(0, 0, 0, 0)
        genome_layout.setSpacing(10)
        
        # 标签
        label = QLabel(self.label_text)
        label.setMinimumWidth(80)
        label.setProperty("role", "fieldLabel")
        genome_layout.addWidget(label)
        
        # 下拉框
        self.combo = QComboBox()
        self.combo.setMinimumWidth(250)
        self.combo.setMinimumHeight(36)
        self.combo.currentIndexChanged.connect(self._on_selection_changed)
        genome_layout.addWidget(self.combo, 1)
        
        # 管理按钮
        if self.show_manage_btn:
            self.manage_btn = QPushButton("管理")
            self.manage_btn.setProperty("secondary", True)
            self.manage_btn.setFixedWidth(60)
            self.manage_btn.clicked.connect(self._open_manager)
            genome_layout.addWidget(self.manage_btn)
        
        main_layout.addLayout(genome_layout)
        
        # 注释版本选择行（可选）
        if self.show_annotation:
            ann_layout = QHBoxLayout()
            ann_layout.setContentsMargins(0, 0, 0, 0)
            ann_layout.setSpacing(10)
            
            ann_label = QLabel(self.annotation_label_text)
            ann_label.setMinimumWidth(80)
            ann_label.setProperty("role", "fieldLabel")
            ann_layout.addWidget(ann_label)
            
            self.ann_combo = QComboBox()
            self.ann_combo.setMinimumWidth(250)
            self.ann_combo.setMinimumHeight(32)
            self.ann_combo.currentIndexChanged.connect(self._on_annotation_changed)
            ann_layout.addWidget(self.ann_combo, 1)
            
            # 占位，与上面对齐
            if self.show_manage_btn:
                spacer = QWidget()
                spacer.setFixedWidth(60)
                ann_layout.addWidget(spacer)
            
            main_layout.addLayout(ann_layout)
    
    def _load_genomes(self):
        """加载基因组列表"""
        self.combo.clear()
        self.combo.addItem("-- 请选择 --", None)
        get_genome_manager().sync_cache_dir()
        genomes = self._db.list_genomes()
        
        for genome in genomes:
            display_name = genome.get("display_name") or genome["name"]
            source = genome.get("source", "custom")
            
            # 添加来源标识
            if source == "custom":
                item_text = f"📁 {display_name}"
            elif source == "ensembl_plants":
                item_text = f"🌱 {display_name}"
            elif source == "igv":
                item_text = f"🌐 {display_name}"
            else:
                item_text = display_name
            
            self.combo.addItem(item_text, genome)
    
    def _on_selection_changed(self, index: int):
        """选择改变时触发"""
        genome = self.combo.itemData(index)
        if genome:
            self.genome_changed.emit(genome["name"], genome)
            # 加载注释版本
            if self.show_annotation:
                self._load_annotations(genome["id"])
        else:
            self.genome_changed.emit("", {})
            if self.show_annotation:
                self._clear_annotations()
    
    def _load_annotations(self, genome_id: int):
        """加载注释版本列表"""
        if not hasattr(self, 'ann_combo'):
            return
        self.ann_combo.clear()
        annotations = self._db.get_annotations(genome_id)
        if not annotations:
            self.ann_combo.addItem("-- 无注释 --", None)
            self.ann_combo.setEnabled(False)
            self.annotation_changed.emit("", {})
        else:
            self.ann_combo.setEnabled(True)
            for ann in annotations:
                self.ann_combo.addItem(ann["source"], ann)
            # 默认选中第一个
            if annotations:
                self.annotation_changed.emit(annotations[0]["source"], annotations[0])
    
    def _clear_annotations(self):
        """清空注释版本列表"""
        if not hasattr(self, 'ann_combo'):
            return
        self.ann_combo.clear()
        self.ann_combo.addItem("-- 请先选择基因组 --", None)
        self.ann_combo.setEnabled(False)
        self.annotation_changed.emit("", {})
    
    def _on_annotation_changed(self, index: int):
        """注释版本选择改变"""
        if not hasattr(self, 'ann_combo'):
            return
        ann = self.ann_combo.itemData(index)
        if ann:
            self.annotation_changed.emit(ann["source"], ann)
        else:
            self.annotation_changed.emit("", {})
    
    def get_selected_annotation(self) -> dict:
        """获取当前选中的注释版本"""
        if not hasattr(self, 'ann_combo'):
            return {}
        return self.ann_combo.currentData() or {}
    
    def get_selected_annotation_source(self) -> str:
        """获取当前选中的注释版本 source"""
        ann = self.get_selected_annotation()
        return ann.get("source", "")
    
    def has_annotation(self) -> bool:
        """当前基因组是否有注释"""
        if not hasattr(self, 'ann_combo'):
            return False
        return self.ann_combo.isEnabled() and self.ann_combo.currentData() is not None
    
    def _open_manager(self):
        """打开基因组管理对话框"""
        from .genome_manager_dialog import GenomeManagerDialog
        dialog = GenomeManagerDialog(self)
        if dialog.exec():
            self._load_genomes()
    
    def get_selected_genome(self) -> dict:
        """获取当前选中的基因组"""
        return self.combo.currentData() or {}
    
    def get_selected_name(self) -> str:
        """获取当前选中的基因组名称"""
        genome = self.get_selected_genome()
        return genome.get("name", "")
    
    def set_selected_genome(self, name: str):
        """设置选中的基因组"""
        self.set_selected_genome_with_emit(name, emit=True)

    def set_selected_genome_with_emit(self, name: str, emit: bool = True):
        """设置选中的基因组（可选是否触发信号）"""
        if not emit:
            self.combo.blockSignals(True)
        for i in range(self.combo.count()):
            genome = self.combo.itemData(i)
            if genome and genome.get("name") == name:
                self.combo.setCurrentIndex(i)
                if not emit:
                    self.combo.blockSignals(False)
                return
        if not emit:
            self.combo.blockSignals(False)

    def clear_selection(self, emit: bool = True):
        """清空选择"""
        if not emit:
            self.combo.blockSignals(True)
        self.combo.setCurrentIndex(0)
        if not emit:
            self.combo.blockSignals(False)
        if emit:
            self.genome_changed.emit("", {})
    
    def refresh(self):
        """刷新基因组列表"""
        current = self.get_selected_name()
        self._load_genomes()
        if current:
            self.set_selected_genome(current)

    def _on_genomes_changed(self):
        self.refresh()

    def closeEvent(self, event):
        if self._db:
            self._db.remove_genome_listener(self._on_genomes_changed)
        super().closeEvent(event)


class GenomePairSelector(QWidget):
    """基因组对选择器（参考 + 查询）"""
    QUERY_LIST_ROW_HEIGHT = 36
    QUERY_LIST_VERTICAL_PADDING = 18
    
    # 信号
    ref_changed = Signal(str, dict)
    qry_changed = Signal(str, dict)
    ref_annotation_changed = Signal(str, dict)  # 参考基因组注释版本改变
    
    def __init__(self, parent=None, show_manage_btn: bool = True, show_ref_annotation: bool = False,
                 multi_query: bool = False,
                 ref_annotation_label: str = "注释版本",
                 show_reference: bool = True,
                 use_cards: bool = False,
                 ref_title: str = "参考",
                 query_title: str = "查询"):
        super().__init__(parent)
        self.show_manage_btn = show_manage_btn
        self.show_ref_annotation = show_ref_annotation
        self.multi_query = multi_query
        self.ref_annotation_label = ref_annotation_label
        self.show_reference = show_reference
        self.use_cards = use_cards
        self.ref_title = ref_title
        self.query_title = query_title
        self._syncing_selection = False
        self._last_ref_name = ""
        self._last_qry_name = ""
        self.ref_selector = None
        self.qry_selector = None
        self._init_ui()
        if self.multi_query:
            self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
    
    def _init_ui(self):
        """初始化 UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(15)
        
        if self.show_reference:
            ref_parent_layout = layout
            if self.use_cards:
                ref_card, ref_parent_layout = create_card(self.ref_title)
                layout.addWidget(ref_card)

            # 参考基因组（可选显示注释版本）
            self.ref_selector = GenomeSelector(
                "参考基因组",
                show_manage_btn=self.show_manage_btn,
                show_annotation=self.show_ref_annotation,
                annotation_label=self.ref_annotation_label
            )
            self.ref_selector.genome_changed.connect(self._on_ref_changed)
            if self.show_ref_annotation:
                self.ref_selector.annotation_changed.connect(self._on_ref_annotation_changed)
            ref_parent_layout.addWidget(self.ref_selector)
        
        # 查询基因组
        query_parent_layout = layout
        if self.use_cards:
            query_card, query_parent_layout = create_card(self.query_title)
            layout.addWidget(query_card)

        self.qry_selector = GenomeSelector("查询基因组", show_manage_btn=False)
        self.qry_selector.genome_changed.connect(self._on_qry_changed)
        if not self.multi_query:
            query_parent_layout.addWidget(self.qry_selector)
        else:
            query_box = QVBoxLayout()
            query_box.setContentsMargins(0, 0, 0, 0)
            query_box.setSpacing(10)
            query_box.addWidget(self.qry_selector)

            self.selected_qry_label = QLabel("已选查询基因组")
            self.selected_qry_label.setProperty("role", "fieldLabel")

            action_layout = QHBoxLayout()
            action_layout.setContentsMargins(80, 0, 0, 0)
            action_layout.setSpacing(8)
            action_layout.addWidget(self.selected_qry_label)
            action_layout.addStretch()

            self.add_qry_btn = QPushButton("添加")
            self.add_qry_btn.setProperty("secondary", True)
            self.add_qry_btn.setMinimumWidth(88)
            self.add_qry_btn.setFixedHeight(32)
            self.add_qry_btn.clicked.connect(self._add_current_query)
            action_layout.addWidget(self.add_qry_btn)

            self.remove_qry_btn = QPushButton("移除")
            self.remove_qry_btn.setProperty("secondary", True)
            self.remove_qry_btn.setMinimumWidth(88)
            self.remove_qry_btn.setFixedHeight(32)
            self.remove_qry_btn.clicked.connect(self._remove_selected_query)
            action_layout.addWidget(self.remove_qry_btn)
            query_box.addLayout(action_layout)

            self.qry_list = QListWidget()
            self.qry_list.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
            self.qry_list.setAlternatingRowColors(True)
            query_box.addWidget(self.qry_list)
            self._sync_selected_query_section()

            query_parent_layout.addLayout(query_box)

    def _sync_selected_query_section(self):
        if not self.multi_query or not hasattr(self, "qry_list"):
            return
        rows = self.qry_list.count()
        has_selection = rows > 0
        self.selected_qry_label.setVisible(has_selection)
        self.remove_qry_btn.setVisible(has_selection)
        self.qry_list.setVisible(has_selection)
        self.qry_list.setMinimumHeight(
            rows * self.QUERY_LIST_ROW_HEIGHT + self.QUERY_LIST_VERTICAL_PADDING
            if has_selection else 0
        )
        self.qry_list.updateGeometry()
        self.updateGeometry()
    
    def get_ref_genome(self) -> dict:
        """获取参考基因组"""
        if not self.ref_selector:
            return {}
        return self.ref_selector.get_selected_genome()
    
    def get_qry_genome(self) -> dict:
        """获取查询基因组"""
        genomes = self.get_qry_genomes()
        return genomes[0] if genomes else {}

    def get_qry_genomes(self) -> list:
        """获取查询基因组列表。未添加列表时兼容返回当前下拉选择。"""
        if not self.multi_query:
            genome = self.qry_selector.get_selected_genome()
            return [genome] if genome else []
        genomes = []
        if hasattr(self, "qry_list"):
            for index in range(self.qry_list.count()):
                item = self.qry_list.item(index)
                genome = item.data(Qt.UserRole) or {}
                if genome:
                    genomes.append(genome)
        if not genomes:
            genome = self.qry_selector.get_selected_genome()
            if genome:
                genomes.append(genome)
        return genomes
    
    def get_ref_annotation(self) -> dict:
        """获取参考基因组的注释版本"""
        if not self.ref_selector:
            return {}
        return self.ref_selector.get_selected_annotation()
    
    def ref_has_annotation(self) -> bool:
        """参考基因组是否有注释"""
        if not self.ref_selector:
            return False
        return self.ref_selector.has_annotation()
    
    def refresh(self):
        """刷新列表"""
        if self.ref_selector:
            self.ref_selector.refresh()
        if self.qry_selector:
            self.qry_selector.refresh()

    def _on_ref_changed(self, name: str, genome: dict):
        if self._syncing_selection:
            return
        qry_name = self.qry_selector.get_selected_name() if self.qry_selector else ""
        if name and name == qry_name:
            QMessageBox.warning(self, "提示", "参考基因组和查询基因组不能相同")
            self._syncing_selection = True
            if self._last_ref_name:
                self.ref_selector.set_selected_genome_with_emit(self._last_ref_name, emit=False)
            else:
                self.ref_selector.clear_selection(emit=False)
            self._syncing_selection = False
            self.ref_changed.emit("", {})
            return
        self._last_ref_name = name
        self.ref_changed.emit(name, genome)
    
    def _on_ref_annotation_changed(self, source: str, ann: dict):
        """参考基因组注释版本改变"""
        self.ref_annotation_changed.emit(source, ann)

    def _on_qry_changed(self, name: str, genome: dict):
        if self._syncing_selection:
            return
        ref_name = self.ref_selector.get_selected_name() if self.ref_selector else ""
        if name and name == ref_name:
            QMessageBox.warning(self, "提示", "参考基因组和查询基因组不能相同")
            self._syncing_selection = True
            if self._last_qry_name:
                self.qry_selector.set_selected_genome_with_emit(self._last_qry_name, emit=False)
            else:
                self.qry_selector.clear_selection(emit=False)
            self._syncing_selection = False
            self.qry_changed.emit("", {})
            return
        self._last_qry_name = name
        self.qry_changed.emit(name, genome)

    def _add_current_query(self):
        if not self.multi_query:
            return
        genome = self.qry_selector.get_selected_genome()
        name = genome.get("name", "")
        if not name:
            QMessageBox.warning(self, "提示", "请选择要添加的查询基因组")
            return
        ref_name = self.ref_selector.get_selected_name() if self.ref_selector else ""
        if name == ref_name:
            QMessageBox.warning(self, "提示", "参考基因组和查询基因组不能相同")
            return
        for index in range(self.qry_list.count()):
            item = self.qry_list.item(index)
            existing = item.data(Qt.UserRole) or {}
            if existing.get("name") == name:
                QMessageBox.information(self, "提示", "该查询基因组已添加")
                return
        display_name = genome.get("display_name") or name
        annotation_state = "有注释" if genome.get("annotation_path") else "无注释"
        item = QListWidgetItem(f"{display_name} ({annotation_state})")
        item.setData(Qt.UserRole, genome)
        self.qry_list.addItem(item)
        self._sync_selected_query_section()
        self._last_qry_name = name
        self.qry_changed.emit(name, genome)

    def _remove_selected_query(self):
        if not self.multi_query or not hasattr(self, "qry_list"):
            return
        row = self.qry_list.currentRow()
        if row < 0:
            return
        self.qry_list.takeItem(row)
        self._sync_selected_query_section()
        genomes = self.get_qry_genomes()
        if genomes:
            self.qry_changed.emit(genomes[0].get("name", ""), genomes[0])
        else:
            self.qry_changed.emit("", {})
