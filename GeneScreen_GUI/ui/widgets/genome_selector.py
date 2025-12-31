"""
GeneScreen 3.0 - 基因组选择器组件

用于选择参考基因组和查询基因组的下拉组件
"""
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QComboBox,
    QLabel, QPushButton, QFrame, QMessageBox
)
from PySide6.QtCore import Signal, Qt

from core import get_genome_manager, get_database


class GenomeSelector(QWidget):
    """基因组选择器"""
    
    # 信号：基因组选择改变
    genome_changed = Signal(str, dict)  # (genome_name, genome_info)
    
    def __init__(self, label: str = "基因组", show_manage_btn: bool = True, parent=None):
        super().__init__(parent)
        self.label_text = label
        self.show_manage_btn = show_manage_btn
        self._db = get_database()
        self._init_ui()
        self._load_genomes()
        self._db.add_genome_listener(self._on_genomes_changed)
    
    def _init_ui(self):
        """初始化 UI"""
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)
        
        # 标签
        label = QLabel(self.label_text)
        label.setMinimumWidth(80)
        label.setProperty("role", "fieldLabel")
        layout.addWidget(label)
        
        # 下拉框
        self.combo = QComboBox()
        self.combo.setMinimumWidth(250)
        self.combo.setMinimumHeight(36)
        self.combo.currentIndexChanged.connect(self._on_selection_changed)
        layout.addWidget(self.combo, 1)
        
        # 管理按钮
        if self.show_manage_btn:
            self.manage_btn = QPushButton("管理")
            self.manage_btn.setProperty("secondary", True)
            self.manage_btn.setFixedWidth(60)
            self.manage_btn.clicked.connect(self._open_manager)
            layout.addWidget(self.manage_btn)
    
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
        else:
            self.genome_changed.emit("", {})
    
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
    
    # 信号
    ref_changed = Signal(str, dict)
    qry_changed = Signal(str, dict)
    
    def __init__(self, parent=None, show_manage_btn: bool = True):
        super().__init__(parent)
        self.show_manage_btn = show_manage_btn
        self._syncing_selection = False
        self._last_ref_name = ""
        self._last_qry_name = ""
        self._init_ui()
    
    def _init_ui(self):
        """初始化 UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(15)
        
        # 参考基因组
        self.ref_selector = GenomeSelector("参考基因组", show_manage_btn=self.show_manage_btn)
        self.ref_selector.genome_changed.connect(self._on_ref_changed)
        layout.addWidget(self.ref_selector)
        
        # 查询基因组
        self.qry_selector = GenomeSelector("查询基因组", show_manage_btn=False)
        self.qry_selector.genome_changed.connect(self._on_qry_changed)
        layout.addWidget(self.qry_selector)
    
    def get_ref_genome(self) -> dict:
        """获取参考基因组"""
        return self.ref_selector.get_selected_genome()
    
    def get_qry_genome(self) -> dict:
        """获取查询基因组"""
        return self.qry_selector.get_selected_genome()
    
    def refresh(self):
        """刷新列表"""
        self.ref_selector.refresh()
        self.qry_selector.refresh()

    def _on_ref_changed(self, name: str, genome: dict):
        if self._syncing_selection:
            return
        qry_name = self.qry_selector.get_selected_name()
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

    def _on_qry_changed(self, name: str, genome: dict):
        if self._syncing_selection:
            return
        ref_name = self.ref_selector.get_selected_name()
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
