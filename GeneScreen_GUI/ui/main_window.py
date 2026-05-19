"""
GeneScreen 1.0 - 主窗口

PySide6 实现的桌面 GUI 主窗口
"""
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QStackedWidget, QPushButton, QLabel, QFrame,
    QMessageBox, QComboBox, QAbstractScrollArea
)
from PySide6.QtCore import Qt, QSize, QEvent, QPoint, QPointF, QRect
from PySide6.QtGui import QIcon, QColor, QPainter, QPen, QPixmap, QBrush, QPainterPath, QRegion, QPalette

import math
import sys
import ctypes
from pathlib import Path

from utils.blast_check import check_blast_installation
from core.task_manager import get_analysis_task_manager


class SidebarButton(QPushButton):
    """侧边栏按钮"""

    def __init__(self, text: str, parent=None):
        super().__init__(text, parent)
        self.setCheckable(True)
        self.setMinimumHeight(44)
        self.setCursor(Qt.PointingHandCursor)


class TitleBar(QFrame):
    """自定义标题栏"""

    def __init__(self, parent=None, on_double_click=None):
        super().__init__(parent)
        self._drag_offset = None
        self._on_double_click = on_double_click
        self.setMouseTracking(True)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self._drag_offset = event.globalPosition().toPoint() - self.window().frameGeometry().topLeft()
            event.accept()

    def mouseMoveEvent(self, event):
        if event.buttons() & Qt.LeftButton and self._drag_offset is not None:
            window = self.window()
            if window.isMaximized():
                return
            window.move(event.globalPosition().toPoint() - self._drag_offset)
            event.accept()

    def mouseReleaseEvent(self, event):
        self._drag_offset = None
        super().mouseReleaseEvent(event)

    def mouseDoubleClickEvent(self, event):
        if event.button() == Qt.LeftButton and self._on_double_click:
            self._on_double_click()
            event.accept()
        super().mouseDoubleClickEvent(event)


class WindowRoot(QWidget):
    """窗口圆角容器，使用抗锯齿绘制背景与边框"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._bg = QColor("#f8f9fb")
        self._border = QColor("#e6e8eb")
        self._radius = 12
        self.setAutoFillBackground(False)
        self.setAttribute(Qt.WA_TranslucentBackground, True)

    def set_theme(self, bg_color: str, border_color: str, radius: int = 12):
        self._bg = QColor(bg_color)
        self._border = QColor(border_color)
        self._radius = radius
        self._update_mask()
        self.update()

    def _update_mask(self):
        if self._radius <= 0:
            self.clearMask()
            return
        window = self.window()
        if window and (window.isMaximized() or window.isFullScreen()):
            self.clearMask()
            return
        rect = self.rect().adjusted(0, 0, -1, -1)
        path = QPainterPath()
        path.addRoundedRect(rect, self._radius, self._radius)
        self.setMask(QRegion(path.toFillPolygon().toPolygon()))

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._update_mask()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        rect = self.rect().adjusted(0, 0, -1, -1)
        path = QPainterPath()
        path.addRoundedRect(rect, self._radius, self._radius)
        painter.setPen(QPen(self._border, 1))
        painter.setBrush(QBrush(self._bg))
        painter.drawPath(path)


class MainWindow(QMainWindow):
    """GeneScreen 1.0 主窗口"""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("GeneScreen")
        self.setMinimumSize(1200, 800)
        self.resize(1400, 900)
        self._sidebar_expanded = True
        self._sidebar_full_width = 220
        self._sidebar_collapsed_width = 64
        self._title_bar_height = 52
        self._dark_mode = False
        self._resize_margin = 6
        self._resize_cursor_active = False
        app = QApplication.instance()
        self._base_app_stylesheet = app.styleSheet() if app else ""
        self._theme_app_stylesheet_cache = {}
        self.setWindowFlags(self.windowFlags() | Qt.FramelessWindowHint)
        self.setAttribute(Qt.WA_TranslucentBackground, True)
        
        # 检查 BLAST+ 安装
        self._check_dependencies()
        
        # 初始化 UI
        self._init_ui()
        self._apply_styles()
        if app:
            app.installEventFilter(self)
    
    def _check_dependencies(self):
        """检查依赖"""
        blast_ok, blast_msg = check_blast_installation()
        if not blast_ok:
            QMessageBox.warning(
                self,
                "依赖检查",
                f"BLAST+ 未正确安装:\n{blast_msg}\n\n请安装 BLAST+ 后重启程序。"
            )
    
    def _init_ui(self):
        """初始化 UI"""
        # 主容器
        central = WindowRoot()
        self.window_root = central
        central.setObjectName("windowRoot")
        self.setCentralWidget(central)
        outer_layout = QVBoxLayout(central)
        outer_layout.setContentsMargins(1, 1, 1, 1)
        outer_layout.setSpacing(0)

        # 标题栏
        title_bar = self._create_title_bar()
        outer_layout.addWidget(title_bar)

        body_layout = QHBoxLayout()
        body_layout.setContentsMargins(0, 0, 0, 0)
        body_layout.setSpacing(0)

        # 侧边栏
        sidebar = self._create_sidebar()
        body_layout.addWidget(sidebar)
        
        # 内容区域
        content = self._create_content_area()
        body_layout.addWidget(content, 1)

        outer_layout.addLayout(body_layout, 1)

        # 状态栏
        self.statusBar().hide()
    
    def _create_sidebar(self) -> QFrame:
        """创建侧边栏"""
        self.sidebar = QFrame()
        self.sidebar.setObjectName("sidebar")
        self.sidebar.setFixedWidth(self._sidebar_full_width)
        
        layout = QVBoxLayout(self.sidebar)
        layout.setContentsMargins(16, 18, 16, 18)
        layout.setSpacing(8)
        
        # Logo / 标题
        header = QWidget()
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(0, 0, 0, 0)
        header_layout.setSpacing(10)

        self.sidebar_toggle = QPushButton()
        self.sidebar_toggle.setObjectName("sidebarToggle")
        self.sidebar_toggle.setFixedSize(30, 30)
        self.sidebar_toggle.setIcon(self._build_toggle_icon("#8b8f97"))
        self.sidebar_toggle.setIconSize(QSize(22, 16))
        self.sidebar_toggle.setToolTip("展开/收缩侧边栏")
        self.sidebar_toggle.clicked.connect(self._toggle_sidebar)
        header_layout.addWidget(self.sidebar_toggle, 0, Qt.AlignLeft | Qt.AlignVCenter)

        title_layout = QVBoxLayout()
        title_layout.setContentsMargins(0, 0, 0, 0)
        title_layout.setSpacing(0)

        self.sidebar_title = QLabel("导航菜单")
        self.sidebar_title.setObjectName("sidebarTitle")
        self.sidebar_title.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        title_layout.addWidget(self.sidebar_title)
        
        self.sidebar_version = QLabel("")
        self.sidebar_version.setObjectName("sidebarVersion")
        self.sidebar_version.setAlignment(Qt.AlignLeft)
        self.sidebar_version.setVisible(False)
        title_layout.addWidget(self.sidebar_version)

        header_layout.addLayout(title_layout)
        header_layout.addStretch()
        layout.addWidget(header)

        header_line = QFrame()
        header_line.setFrameShape(QFrame.HLine)
        header_line.setObjectName("sidebarHeaderLine")
        layout.addWidget(header_line)
        
        layout.addSpacing(6)
        
        # 导航按钮
        self.nav_buttons = []
        self.sidebar_text_buttons = []
        
        # Gene ID 模式
        btn_gene_id = SidebarButton("Gene ID")
        btn_gene_id.setObjectName("navButton")
        btn_gene_id.setChecked(True)
        btn_gene_id.clicked.connect(lambda: self._switch_page(0))
        layout.addWidget(btn_gene_id)
        self.nav_buttons.append(btn_gene_id)
        self._register_sidebar_button(btn_gene_id, "Gene ID", "id")
        
        # Location 模式
        btn_location = SidebarButton("Location")
        btn_location.setObjectName("navButton")
        btn_location.clicked.connect(lambda: self._switch_page(1))
        layout.addWidget(btn_location)
        self.nav_buttons.append(btn_location)
        self._register_sidebar_button(btn_location, "Location", "pin")
        
        # Sequence 模式
        btn_sequence = SidebarButton("Sequence")
        btn_sequence.setObjectName("navButton")
        btn_sequence.clicked.connect(lambda: self._switch_page(2))
        layout.addWidget(btn_sequence)
        self.nav_buttons.append(btn_sequence)
        self._register_sidebar_button(btn_sequence, "Sequence", "code")
        
        layout.addSpacing(20)
        
        # 分隔线
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setObjectName("sidebarLine")
        layout.addWidget(line)
        
        layout.addSpacing(20)
        
        # 历史记录
        btn_history = SidebarButton("历史记录")
        btn_history.setObjectName("navButton")
        btn_history.clicked.connect(lambda: self._switch_page(3))
        layout.addWidget(btn_history)
        self.nav_buttons.append(btn_history)
        self._register_sidebar_button(btn_history, "历史记录", "clock")
        
        # 基因组管理
        btn_genomes = SidebarButton("基因组管理")
        btn_genomes.setObjectName("navButton")
        btn_genomes.clicked.connect(lambda: self._switch_page(4))
        layout.addWidget(btn_genomes)
        self.nav_buttons.append(btn_genomes)
        self._register_sidebar_button(btn_genomes, "基因组管理", "database")
        
        layout.addStretch()
        
        # 设置按钮
        btn_settings = SidebarButton("设置")
        btn_settings.setObjectName("navButtonSecondary")
        btn_settings.clicked.connect(lambda: self._switch_page(5))
        layout.addWidget(btn_settings)
        self._register_sidebar_button(btn_settings, "设置", "settings")
        self.nav_buttons.append(btn_settings)
        
        return self.sidebar

    def _create_title_bar(self) -> QFrame:
        bar = TitleBar(self, on_double_click=self._toggle_maximize)
        bar.setObjectName("titleBar")
        bar.setFixedHeight(self._title_bar_height)

        layout = QHBoxLayout(bar)
        layout.setContentsMargins(16, 0, 16, 0)
        layout.setSpacing(10)

        icon_label = QLabel()
        icon_label.setObjectName("titleBarIcon")
        icon_label.setPixmap(self._build_app_icon())
        icon_label.setFixedSize(22, 22)
        icon_label.setScaledContents(True)
        icon_label.setAttribute(Qt.WA_TransparentForMouseEvents, True)
        layout.addWidget(icon_label)

        title = QLabel("GeneScreen")
        title.setObjectName("titleBarTitle")
        title.setAttribute(Qt.WA_TransparentForMouseEvents, True)
        layout.addWidget(title)

        layout.addStretch()

        self.theme_btn = QPushButton()
        self.theme_btn.setObjectName("titleButton")
        self.theme_btn.setProperty("titleRole", "theme")
        self.theme_btn.setFixedSize(26, 26)
        self.theme_btn.clicked.connect(self._toggle_theme)
        layout.addWidget(self.theme_btn)

        self.min_btn = QPushButton()
        self.min_btn.setObjectName("titleButton")
        self.min_btn.setProperty("titleRole", "min")
        self.min_btn.setFixedSize(26, 26)
        self.min_btn.clicked.connect(self.showMinimized)
        layout.addWidget(self.min_btn)

        self.max_btn = QPushButton()
        self.max_btn.setObjectName("titleButton")
        self.max_btn.setProperty("titleRole", "max")
        self.max_btn.setFixedSize(26, 26)
        self.max_btn.clicked.connect(self._toggle_maximize)
        layout.addWidget(self.max_btn)

        self.close_btn = QPushButton()
        self.close_btn.setObjectName("titleButton")
        self.close_btn.setProperty("titleRole", "close")
        self.close_btn.setFixedSize(26, 26)
        self.close_btn.clicked.connect(self.close)
        layout.addWidget(self.close_btn)

        self._update_theme_button_icon()
        self._update_maximize_button_icon()
        self._update_title_button_icons()

        return bar

    def _register_sidebar_button(self, button: SidebarButton, label: str, icon_name: str):
        button._label = label
        button._icon_name = icon_name
        button.setText(label)
        button.setToolTip(label)
        button.setIcon(self._build_sidebar_icon(icon_name))
        button.setIconSize(QSize(22, 22))
        button.setProperty("collapsed", False)
        self.sidebar_text_buttons.append(button)

    def _device_pixel_ratio(self) -> float:
        screen = self.windowHandle().screen() if self.windowHandle() else QApplication.primaryScreen()
        return screen.devicePixelRatio() if screen else 1.0

    def _create_pixmap(self, size: int) -> QPixmap:
        dpr = self._device_pixel_ratio()
        pixmap = QPixmap(int(size * dpr), int(size * dpr))
        pixmap.setDevicePixelRatio(dpr)
        pixmap.fill(Qt.transparent)
        return pixmap

    def _build_toggle_icon(self, color: str) -> QIcon:
        size = 22
        pixmap = self._create_pixmap(size)
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.Antialiasing)
        pen = QPen(QColor(color))
        pen.setWidthF(1.8)
        pen.setCapStyle(Qt.RoundCap)
        painter.setPen(pen)
        for y in (7, 11, 15):
            painter.drawLine(4, y, size - 4, y)
        painter.end()
        icon = QIcon()
        icon.addPixmap(pixmap, QIcon.Normal, QIcon.Off)
        return icon

    def _build_sidebar_icon(self, name: str) -> QIcon:
        normal_color = QColor("#606266") if not self._dark_mode else QColor("#c7c9cc")
        normal = self._draw_sidebar_icon(name, normal_color)
        active = self._draw_sidebar_icon(name, QColor("#ffffff"))
        icon = QIcon()
        icon.addPixmap(normal, QIcon.Normal, QIcon.Off)
        icon.addPixmap(active, QIcon.Normal, QIcon.On)
        return icon

    def _draw_sidebar_icon(self, name: str, color: QColor, size: int = 20) -> QPixmap:
        pixmap = self._create_pixmap(size)
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.Antialiasing)
        pen = QPen(color)
        pen.setWidthF(1.8)
        pen.setCapStyle(Qt.RoundCap)
        pen.setJoinStyle(Qt.RoundJoin)
        painter.setPen(pen)

        if name == "grid":
            margin = 3
            gap = 3
            cell = (size - margin * 2 - gap) / 2
            for row in range(2):
                for col in range(2):
                    x = margin + col * (cell + gap)
                    y = margin + row * (cell + gap)
                    painter.drawRect(x, y, cell, cell)
        elif name == "id":
            margin = 3
            width = size - margin * 2
            height = size - margin * 2
            painter.drawRoundedRect(margin, margin + 1, width, height - 2, 3, 3)
            painter.drawEllipse(margin + 3, margin + 4, 4, 4)
            painter.drawLine(margin + 9, margin + 6, size - margin - 3, margin + 6)
            painter.drawLine(margin + 9, margin + 10, size - margin - 4, margin + 10)
        elif name == "pin":
            center_x = size / 2
            center_y = size / 2 - 2
            radius = 4.5
            painter.drawEllipse(center_x - radius, center_y - radius, radius * 2, radius * 2)
            painter.drawLine(center_x, center_y + radius - 0.5, center_x, size - 3)
        elif name == "code":
            left = 5
            right = size - 5
            top = 5
            bottom = size - 5
            mid = size / 2
            painter.drawLine(left + 3, top, left, mid)
            painter.drawLine(left, mid, left + 3, bottom)
            painter.drawLine(right - 3, top, right, mid)
            painter.drawLine(right, mid, right - 3, bottom)
        elif name == "clock":
            radius = 7
            center = size / 2
            painter.drawEllipse(center - radius, center - radius, radius * 2, radius * 2)
            painter.drawLine(center, center, center, center - 4)
            painter.drawLine(center, center, center + 3.5, center + 2)
        elif name == "database":
            margin = 3
            width = size - margin * 2
            top_h = 4
            painter.drawEllipse(margin, margin, width, top_h * 2)
            painter.drawLine(margin, margin + top_h, margin, size - margin - top_h)
            painter.drawLine(margin + width, margin + top_h, margin + width, size - margin - top_h)
            painter.drawArc(margin, size - margin - top_h * 2, width, top_h * 2, 0, -180 * 16)
        else:
            y_positions = (5, 10, 15)
            for i, y in enumerate(y_positions):
                painter.drawLine(4, y, size - 4, y)
                knob_x = 7 if i % 2 == 0 else size - 7
                painter.drawEllipse(knob_x - 2.2, y - 2.2, 4.4, 4.4)

        painter.end()
        return pixmap

    def _build_app_icon(self) -> QPixmap:
        size = 20
        icon_path = Path(__file__).resolve().parent / "resources" / "icons" / "app.ico"
        if icon_path.exists():
            pixmap = QPixmap(str(icon_path))
            if not pixmap.isNull():
                return pixmap.scaled(size, size, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        pixmap = self._create_pixmap(size)
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setPen(Qt.NoPen)
        painter.setBrush(QBrush(QColor("#3b82f6")))
        painter.drawRoundedRect(0, 0, size, size, 5, 5)
        painter.setBrush(QBrush(QColor("#1d4ed8")))
        painter.drawRoundedRect(4, 4, size - 8, size - 8, 3, 3)
        painter.end()
        return pixmap

    def _build_title_icon(self, name: str, color: QColor) -> QIcon:
        size = 14
        pixmap = self._create_pixmap(size)
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.Antialiasing)
        pen = QPen(color)
        pen.setWidthF(1.6)
        pen.setCapStyle(Qt.RoundCap)
        pen.setJoinStyle(Qt.RoundJoin)
        painter.setPen(pen)

        if name == "sun":
            center = size / 2
            radius = 3.5
            painter.drawEllipse(center - radius, center - radius, radius * 2, radius * 2)
            for angle in range(0, 360, 45):
                rad = math.radians(angle)
                x1 = center + math.cos(rad) * (radius + 1.5)
                y1 = center + math.sin(rad) * (radius + 1.5)
                x2 = center + math.cos(rad) * (radius + 4)
                y2 = center + math.sin(rad) * (radius + 4)
                painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
        elif name == "moon":
            painter.setPen(Qt.NoPen)
            painter.setBrush(QBrush(color))
            outer = QPainterPath()
            inner = QPainterPath()
            outer.addEllipse(3, 3, 8, 8)
            inner.addEllipse(6, 2, 8, 8)
            crescent = outer.subtracted(inner)
            painter.drawPath(crescent)
            painter.setPen(pen)
        elif name == "min":
            painter.drawLine(3, size / 2, size - 3, size / 2)
        elif name == "max":
            painter.drawRect(3, 3, size - 6, size - 6)
        elif name == "restore":
            painter.drawRect(4, 5, size - 7, size - 7)
            painter.drawRect(2, 3, size - 7, size - 7)
        else:
            painter.drawLine(4, 4, size - 4, size - 4)
            painter.drawLine(size - 4, 4, 4, size - 4)

        painter.end()
        icon = QIcon()
        icon.addPixmap(pixmap, QIcon.Normal, QIcon.Off)
        return icon

    def _toggle_theme(self):
        self._dark_mode = not self._dark_mode
        self._apply_styles()
        self._refresh_sidebar_icons()
        toggle_color = "#8b8f97" if not self._dark_mode else "#d1d5db"
        self.sidebar_toggle.setIcon(self._build_toggle_icon(toggle_color))
        self._update_title_button_icons()

    def _toggle_maximize(self):
        if self.isMaximized():
            self.showNormal()
        else:
            self.showMaximized()
        self._update_maximize_button_icon()
        if hasattr(self, "window_root"):
            self.window_root._update_mask()

    def _update_theme_button_icon(self):
        icon_color = QColor("#5f6368") if not self._dark_mode else QColor("#e5e7eb")
        icon_name = "sun" if self._dark_mode else "moon"
        self.theme_btn.setIcon(self._build_title_icon(icon_name, icon_color))
        icon_size = 18 if icon_name == "moon" else 14
        self.theme_btn.setIconSize(QSize(icon_size, icon_size))

    def _update_maximize_button_icon(self):
        icon_name = "restore" if self.isMaximized() else "max"
        icon_color = QColor("#1f2937") if not self._dark_mode else QColor("#111827")
        self.max_btn.setIcon(self._build_title_icon(icon_name, icon_color))
        self.max_btn.setIconSize(QSize(14, 14))

    def _update_title_button_icons(self):
        icon_color = QColor("#1f2937") if not self._dark_mode else QColor("#111827")
        self.min_btn.setIcon(self._build_title_icon("min", icon_color))
        self.min_btn.setIconSize(QSize(14, 14))
        self.close_btn.setIcon(self._build_title_icon("close", icon_color))
        self.close_btn.setIconSize(QSize(14, 14))
        self._update_theme_button_icon()
        self._update_maximize_button_icon()

    def _refresh_sidebar_icons(self):
        for btn in self.sidebar_text_buttons:
            icon_name = getattr(btn, "_icon_name", None)
            if icon_name:
                btn.setIcon(self._build_sidebar_icon(icon_name))

    def changeEvent(self, event):
        if event.type() == QEvent.WindowStateChange:
            self._update_maximize_button_icon()
        super().changeEvent(event)

    def closeEvent(self, event):
        task_manager = get_analysis_task_manager()
        if task_manager.has_tasks():
            reply = QMessageBox.warning(
                self,
                "仍有分析任务",
                (
                    f"当前仍有 {task_manager.running_count()} 个 Running 任务、"
                    f"{task_manager.pending_count()} 个 Pending 任务。\n\n"
                    "关闭软件会中断正在执行或排队的分析任务，是否继续关闭？"
                ),
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                event.ignore()
                return
            task_manager.mark_unfinished_failed()
        super().closeEvent(event)

    def showEvent(self, event):
        super().showEvent(event)
        self._apply_windows_rounding()

    def _apply_windows_rounding(self):
        if not sys.platform.startswith("win"):
            return
        try:
            hwnd = int(self.winId())
            DWMWA_WINDOW_CORNER_PREFERENCE = 33
            DWMWCP_ROUND = 2
            preference = ctypes.c_int(DWMWCP_ROUND)
            ctypes.windll.dwmapi.DwmSetWindowAttribute(
                hwnd,
                DWMWA_WINDOW_CORNER_PREFERENCE,
                ctypes.byref(preference),
                ctypes.sizeof(preference),
            )
        except Exception:
            pass

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Wheel:
            combo = self._combo_for_wheel_target(obj)
            if combo and self._suppress_combo_wheel(combo, event):
                return True

        if self.isMaximized() or self.isFullScreen():
            return super().eventFilter(obj, event)

        if event.type() == QEvent.MouseMove:
            if isinstance(obj, QWidget) and (obj is self or self.isAncestorOf(obj)):
                if hasattr(event, "globalPosition"):
                    edges = self._hit_test_edges(event.globalPosition().toPoint())
                    self._update_resize_cursor(edges)
            return super().eventFilter(obj, event)

        if event.type() in (QEvent.Leave, QEvent.HoverLeave):
            if self._resize_cursor_active:
                self.unsetCursor()
                self._resize_cursor_active = False
            return super().eventFilter(obj, event)

        if event.type() == QEvent.MouseButtonPress and event.button() == Qt.LeftButton:
            if isinstance(obj, QWidget) and (obj is self or self.isAncestorOf(obj)):
                edges = self._hit_test_edges(event.globalPosition().toPoint())
                if edges != Qt.Edges():
                    if self.windowHandle():
                        self.windowHandle().startSystemResize(edges)
                        return True
        return super().eventFilter(obj, event)

    def _combo_for_wheel_target(self, obj) -> QComboBox:
        widget = obj if isinstance(obj, QWidget) else None
        while widget:
            if isinstance(widget, QComboBox):
                return widget
            widget = widget.parentWidget()
        return None

    def _suppress_combo_wheel(self, combo: QComboBox, event) -> bool:
        """Disable mouse-wheel selection changes on closed combo boxes."""
        self._scroll_nearest_parent(combo, event)
        event.accept()
        return True

    def _scroll_nearest_parent(self, widget: QWidget, event) -> None:
        parent = widget.parentWidget()
        while parent:
            if isinstance(parent, QAbstractScrollArea):
                scroll_bar = parent.verticalScrollBar()
                pixel_delta = event.pixelDelta().y()
                angle_delta = event.angleDelta().y()
                if pixel_delta:
                    delta = pixel_delta
                elif angle_delta:
                    delta = int(angle_delta / 120 * scroll_bar.singleStep() * 3)
                else:
                    delta = 0
                if delta:
                    scroll_bar.setValue(scroll_bar.value() - delta)
                return
            parent = parent.parentWidget()

    def _hit_test_edges(self, global_pos: QPoint) -> Qt.Edges:
        rect = self.frameGeometry()
        margin = self._resize_margin
        left = abs(global_pos.x() - rect.left()) <= margin
        right = abs(global_pos.x() - rect.right()) <= margin
        top = abs(global_pos.y() - rect.top()) <= margin
        bottom = abs(global_pos.y() - rect.bottom()) <= margin

        edges = Qt.Edges()
        if left:
            edges |= Qt.LeftEdge
        if right:
            edges |= Qt.RightEdge
        if top:
            edges |= Qt.TopEdge
        if bottom:
            edges |= Qt.BottomEdge
        return edges

    def _update_resize_cursor(self, edges: Qt.Edges) -> None:
        if edges == Qt.Edges():
            if self._resize_cursor_active:
                self.unsetCursor()
                self._resize_cursor_active = False
            return
        if (edges & Qt.LeftEdge and edges & Qt.TopEdge) or (edges & Qt.RightEdge and edges & Qt.BottomEdge):
            cursor = Qt.SizeFDiagCursor
        elif (edges & Qt.RightEdge and edges & Qt.TopEdge) or (edges & Qt.LeftEdge and edges & Qt.BottomEdge):
            cursor = Qt.SizeBDiagCursor
        elif edges & (Qt.LeftEdge | Qt.RightEdge):
            cursor = Qt.SizeHorCursor
        else:
            cursor = Qt.SizeVerCursor
        self.setCursor(cursor)
        self._resize_cursor_active = True

    
    def _create_content_area(self) -> QWidget:
        """创建内容区域"""
        content = QWidget()
        content.setObjectName("contentArea")
        layout = QVBoxLayout(content)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # 页面堆栈
        self.page_stack = QStackedWidget()
        self.page_stack.setObjectName("pageStack")
        
        # 懒加载页面（先只创建第一个，其他用占位符）
        from ui.pages import GeneIDPage
        
        self._page_classes = None  # 延迟导入
        self._pages_loaded = [True, False, False, False, False, False]
        
        self.page_stack.addWidget(GeneIDPage())  # 首页立即加载
        for _ in range(5):  # 其他页面用占位符
            self.page_stack.addWidget(QWidget())
        
        layout.addWidget(self.page_stack)
        
        return content
    
    def _ensure_page_loaded(self, index: int):
        """确保页面已加载（懒加载）"""
        if self._pages_loaded[index]:
            return
        
        # 延迟导入页面类
        if self._page_classes is None:
            from ui.pages import GeneIDPage, LocationPage, SequencePage, HistoryPage, GenomeManagerPage, SettingsPage
            self._page_classes = [GeneIDPage, LocationPage, SequencePage, HistoryPage, GenomeManagerPage, SettingsPage]
        
        # 替换占位符为实际页面
        old_widget = self.page_stack.widget(index)
        new_widget = self._page_classes[index]()
        self.page_stack.removeWidget(old_widget)
        old_widget.deleteLater()
        self.page_stack.insertWidget(index, new_widget)
        self._pages_loaded[index] = True
    
    def _switch_page(self, index: int):
        """切换页面"""
        self._ensure_page_loaded(index)
        self.page_stack.setCurrentIndex(index)
        
        # 更新按钮状态
        for i, btn in enumerate(self.nav_buttons):
            btn.setChecked(i == index)

    def _toggle_sidebar(self):
        self._set_sidebar_collapsed(self._sidebar_expanded)

    def _set_sidebar_collapsed(self, collapsed: bool):
        self._sidebar_expanded = not collapsed
        width = self._sidebar_collapsed_width if collapsed else self._sidebar_full_width
        self.sidebar.setFixedWidth(width)
        self.sidebar_title.setVisible(not collapsed)
        self.sidebar_version.setVisible(False)
        for btn in self.sidebar_text_buttons:
            btn.setText("" if collapsed else btn._label)
            btn.setToolTip(btn._label if collapsed else "")
            btn.setProperty("collapsed", collapsed)
            btn.style().unpolish(btn)
            btn.style().polish(btn)

    def _apply_app_palette(
        self,
        theme_bg: str,
        panel_bg: str,
        panel_text: str,
        input_bg: str,
        input_text: str,
        input_placeholder: str,
        input_border: str,
        table_alt: str,
        secondary_button_bg: str,
        secondary_button_text: str,
        readonly_text: str,
    ):
        app = QApplication.instance()
        if not app:
            return

        palette = QPalette()
        palette.setColor(QPalette.Window, QColor(theme_bg))
        palette.setColor(QPalette.WindowText, QColor(panel_text))
        palette.setColor(QPalette.Base, QColor(input_bg))
        palette.setColor(QPalette.AlternateBase, QColor(table_alt))
        palette.setColor(QPalette.ToolTipBase, QColor(panel_bg))
        palette.setColor(QPalette.ToolTipText, QColor(panel_text))
        palette.setColor(QPalette.Text, QColor(input_text))
        palette.setColor(QPalette.Button, QColor(secondary_button_bg))
        palette.setColor(QPalette.ButtonText, QColor(secondary_button_text))
        palette.setColor(QPalette.BrightText, QColor("#ffffff"))
        palette.setColor(QPalette.Link, QColor(secondary_button_text))
        palette.setColor(QPalette.Highlight, QColor("#2f80ff"))
        palette.setColor(QPalette.HighlightedText, QColor("#ffffff"))
        palette.setColor(QPalette.PlaceholderText, QColor(input_placeholder))

        palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor(readonly_text))
        palette.setColor(QPalette.Disabled, QPalette.Text, QColor(readonly_text))
        palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(readonly_text))
        palette.setColor(QPalette.Disabled, QPalette.Base, QColor(input_bg))
        palette.setColor(QPalette.Disabled, QPalette.Button, QColor(input_border))

        app.setPalette(palette)
    
    def _apply_styles(self):
        """应用样式"""
        theme_bg = "#f8f9fb"
        theme_border = "#e6e8eb"
        title_text = "#1f2328"
        sidebar_bg = "#f8f9fb"
        sidebar_border = "#e6e8eb"
        nav_text = "#606266"
        nav_hover_bg = "#f2f4f7"
        nav_hover_text = "#303133"
        content_bg = "#f5f6f8"
        theme_toggle_bg = "#d1d5db"
        theme_toggle_hover = "#c7cbd3"
        panel_bg = "#ffffff"
        panel_border = "#e6e8eb"
        panel_text = "#303133"
        muted_text = "#666666"
        input_bg = "#ffffff"
        input_border = "#dcdfe6"
        input_text = "#303133"
        input_placeholder = "#a0a4ad"
        table_alt = "#f7f8fa"
        tab_bg = "#ffffff"
        tab_text = "#606266"
        tab_selected_bg = "#ffffff"
        tab_selected_text = "#303133"
        tab_hover_bg = "#f2f4f7"
        readonly_text = "#8c9197"
        scrollbar_bg = "#f0f2f5"
        scrollbar_handle = "#c4c8cf"
        scrollbar_handle_hover = "#aeb4bd"
        secondary_button_bg = "#f3f4f6"
        secondary_button_hover = "#e9edf5"
        secondary_button_pressed = "#dde3ee"
        secondary_button_border = "#cfd6e4"
        secondary_button_text = "#2f80ff"
        primary_button_bg = "#2f80ff"
        primary_button_hover = "#1f72e8"
        primary_button_pressed = "#195fc4"
        primary_button_border = "#1f66d1"
        primary_button_text = "#ffffff"
        primary_button_disabled_bg = "#d6dde8"
        primary_button_disabled_text = "#8c98aa"
        primary_button_disabled_border = "#c5cfdd"

        if self._dark_mode:
            theme_bg = "#202124"
            theme_border = "#2b2b2b"
            title_text = "#e5e7eb"
            sidebar_bg = "#1b1c1e"
            sidebar_border = "#2b2b2b"
            nav_text = "#c7c9cc"
            nav_hover_bg = "#2a2b2f"
            nav_hover_text = "#ffffff"
            theme_toggle_bg = "#3a3b3f"
            theme_toggle_hover = "#4a4b50"
            content_bg = "#14161a"
            panel_bg = "#1f2329"
            panel_border = "#2b2f36"
            panel_text = "#e5e7eb"
            muted_text = "#9aa0a6"
            input_bg = "#1b1f25"
            input_border = "#2f343c"
            input_text = "#e5e7eb"
            input_placeholder = "#8b9099"
            table_alt = "#1a1d22"
            tab_bg = "#1f2329"
            tab_text = "#b9bec7"
            tab_selected_bg = "#2a2f36"
            tab_selected_text = "#ffffff"
            tab_hover_bg = "#2b2f36"
            readonly_text = "#9aa0a6"
            scrollbar_bg = "#1b1f25"
            scrollbar_handle = "#3b4048"
            scrollbar_handle_hover = "#4b515b"
            secondary_button_bg = "#242a32"
            secondary_button_hover = "#2d3642"
            secondary_button_pressed = "#354052"
            secondary_button_border = "#3a4657"
            secondary_button_text = "#8ab4ff"
            primary_button_bg = "#2f80ff"
            primary_button_hover = "#4c94ff"
            primary_button_pressed = "#1f66d1"
            primary_button_border = "#5aa0ff"
            primary_button_text = "#ffffff"
            primary_button_disabled_bg = "#303844"
            primary_button_disabled_text = "#7b8594"
            primary_button_disabled_border = "#3d4654"

        if hasattr(self, "window_root"):
            self.window_root.set_theme(theme_bg, theme_border, radius=12)

        self._apply_app_palette(
            theme_bg,
            panel_bg,
            panel_text,
            input_bg,
            input_text,
            input_placeholder,
            input_border,
            table_alt,
            secondary_button_bg,
            secondary_button_text,
            readonly_text,
        )

        theme_stylesheet = f"""
            /* 标题栏 */
            #titleBar {{
                background: {theme_bg};
                border-bottom: 1px solid {theme_border};
            }}

            #titleBarTitle {{
                font-size: 17px;
                font-weight: 700;
                color: {title_text};
            }}

            #titleBarIcon {{
                background: transparent;
            }}

            #titleButton {{
                border: none;
                border-radius: 13px;
            }}

            #titleButton[titleRole="theme"] {{
                background: {theme_toggle_bg};
            }}

            #titleButton[titleRole="theme"]:hover {{
                background: {theme_toggle_hover};
            }}

            #titleButton[titleRole="min"] {{
                background: #f6c344;
            }}

            #titleButton[titleRole="max"] {{
                background: #5ac85a;
            }}

            #titleButton[titleRole="close"] {{
                background: #ff5f56;
            }}

            #titleButton[titleRole="min"]:hover {{
                background: #eab736;
            }}

            #titleButton[titleRole="max"]:hover {{
                background: #4bb94d;
            }}

            #titleButton[titleRole="close"]:hover {{
                background: #ef4b41;
            }}

            #titleButton[titleRole="min"]:pressed {{
                background: #d5a62f;
            }}

            #titleButton[titleRole="max"]:pressed {{
                background: #3ea041;
            }}

            #titleButton[titleRole="close"]:pressed {{
                background: #d9423a;
            }}

            /* 窗口圆角 */
            #windowRoot {{
                background: transparent;
            }}

            /* 侧边栏 */
            #sidebar {{
                background: {sidebar_bg};
                border-right: 1px solid {sidebar_border};
            }}
            
            #sidebarTitle {{
                font-size: 15px;
                font-weight: 600;
                color: {nav_hover_text};
            }}
            
            #sidebarVersion {{
                font-size: 12px;
                color: {nav_text};
            }}

            #sidebarToggle {{
                background: transparent;
                border: none;
                border-radius: 6px;
            }}

            #sidebarToggle:hover {{
                background: {nav_hover_bg};
            }}

            #sidebarToggle:pressed {{
                background: {theme_border};
            }}

            #sidebarHeaderLine {{
                background: {sidebar_border};
                max-height: 1px;
                margin: 10px 0 6px 0;
            }}
            
            #sidebarLine {{
                background: {sidebar_border};
                max-height: 1px;
                margin: 8px 0;
            }}
            
            /* 导航按钮 */
            #navButton, #navButtonSecondary {{
                background: transparent;
                border: none;
                border-radius: 10px;
                color: {nav_text};
                font-size: 14px;
                text-align: left;
                padding: 10px 12px;
            }}

            #navButton[collapsed="true"], #navButtonSecondary[collapsed="true"] {{
                text-align: center;
                padding: 10px 0px;
            }}
            
            #navButton:hover, #navButtonSecondary:hover {{
                background: {nav_hover_bg};
                color: {nav_hover_text};
            }}
            
            #navButton:checked, #navButtonSecondary:checked {{
                background: #2f80ff;
                color: white;
                font-weight: 600;
            }}
            
            /* 内容区域 */
            #contentArea {{
                background: {content_bg};
            }}
            
            #pageStack {{
                background: {content_bg};
            }}

            QLabel {{
                color: {panel_text};
            }}

            QLabel[role="pageTitle"] {{
                font-size: 24px;
                font-weight: 700;
                color: {panel_text};
            }}

            QLabel[role="pageDesc"] {{
                font-size: 14px;
                color: {muted_text};
            }}

            QLabel[role="muted"] {{
                color: {muted_text};
            }}

            QLabel[role="fieldLabel"] {{
                font-weight: 500;
                color: {panel_text};
            }}

            QLabel[role="searchStatus"] {{
                font-size: 18px;
                color: {muted_text};
                background: {panel_bg};
            }}

            QListWidget[role="selectedQueryList"] {{
                background: {input_bg};
                color: {panel_text};
                border: 1px solid {input_border};
                border-radius: 6px;
                padding: 4px;
                outline: 0;
                selection-background-color: {tab_hover_bg};
                selection-color: {panel_text};
            }}

            QListWidget[role="selectedQueryList"]::item {{
                min-height: 26px;
                padding: 4px 8px;
                border-radius: 4px;
            }}

            QListWidget[role="selectedQueryList"]::item:selected,
            QListWidget[role="selectedQueryList"]::item:selected:active,
            QListWidget[role="selectedQueryList"]::item:selected:!active {{
                background: {tab_hover_bg};
                color: {panel_text};
                border: 1px solid #2f80ff;
            }}

            QCheckBox {{
                color: {panel_text};
                min-height: 28px;
                spacing: 6px;
            }}

            QCheckBox::indicator {{
                width: 16px;
                height: 16px;
                border: 1px solid {input_border};
                border-radius: 4px;
                background: {input_bg};
            }}

            QCheckBox::indicator:hover {{
                border-color: {secondary_button_text};
            }}

            QCheckBox::indicator:checked {{
                background: #2f80ff;
                border: 1px solid #2f80ff;
            }}

            QCheckBox::indicator:disabled {{
                background: {secondary_button_bg};
                border-color: {secondary_button_border};
            }}

            QCheckBox#tableCheckBox {{
                min-width: 22px;
                max-width: 22px;
                min-height: 22px;
                max-height: 22px;
                spacing: 0px;
                padding: 0px;
                margin: 0px;
            }}

            QCheckBox#tableCheckBox::indicator {{
                width: 18px;
                height: 18px;
                margin: 2px;
            }}

            QCheckBox#tableCheckBox::indicator:checked {{
                background: #2f80ff;
                border: 1px solid #2f80ff;
            }}

            QCheckBox#tableCheckBox::indicator:disabled {{
                background: {secondary_button_bg};
                border-color: {secondary_button_border};
            }}

            QTabWidget::pane {{
                background: {panel_bg};
                border: 1px solid {panel_border};
                border-radius: 10px;
                top: -1px;
                padding: 8px;
            }}

            QTabBar::tab {{
                background: {tab_bg};
                color: {tab_text};
                border: 1px solid {panel_border};
                border-bottom: none;
                padding: 6px 14px;
                margin-right: 6px;
                border-top-left-radius: 8px;
                border-top-right-radius: 8px;
            }}

            QTabBar::tab:selected {{
                background: {tab_selected_bg};
                color: {tab_selected_text};
                font-weight: 600;
            }}

            QTabBar::tab:hover {{
                background: {tab_hover_bg};
            }}

            /* 子面板与输入控件 */
            QGroupBox {{
                background: {panel_bg};
                border: 1px solid {panel_border};
                border-radius: 10px;
                margin-top: 10px;
            }}

            QGroupBox::title {{
                subcontrol-origin: margin;
                subcontrol-position: top left;
                padding: 0 6px;
                color: {panel_text};
            }}

            QGroupBox QLabel {{
                color: {panel_text};
            }}

            QLineEdit, QTextEdit, QPlainTextEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
                background: {input_bg};
                color: {input_text};
                border: 1px solid {input_border};
                border-radius: 6px;
                padding: 6px 10px;
            }}

            QSpinBox[paramInput="true"], QDoubleSpinBox[paramInput="true"], QLineEdit[paramInput="true"], QComboBox[paramInput="true"] {{
                border-radius: 5px;
                padding: 4px 8px;
                min-height: 28px;
                min-width: 88px;
            }}

            QLineEdit:focus, QTextEdit:focus, QPlainTextEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
                border-color: #2f80ff;
            }}

            QLineEdit[role="readonly"] {{
                background: {input_bg};
                color: {readonly_text};
            }}

            QLineEdit::placeholder {{
                color: {input_placeholder};
            }}

            QPushButton[secondary="true"] {{
                background: {secondary_button_bg};
                color: {secondary_button_text};
                border: 1px solid {secondary_button_border};
            }}

            QPushButton[secondary="true"]:hover {{
                background: {secondary_button_hover};
                border-color: {secondary_button_text};
            }}

            QPushButton[secondary="true"]:pressed {{
                background: {secondary_button_pressed};
            }}

            QPushButton[compactAction="true"] {{
                padding: 6px 14px;
            }}

            QPushButton[smallAction="true"] {{
                font-size: 11px;
                padding: 2px 6px;
            }}

            QPushButton[primaryAction="true"] {{
                background: {primary_button_bg};
                color: {primary_button_text};
                border: 1px solid {primary_button_border};
                border-radius: 8px;
                font-weight: 600;
            }}

            QPushButton[primaryAction="true"]:hover {{
                background: {primary_button_hover};
                border-color: {primary_button_border};
            }}

            QPushButton[primaryAction="true"]:pressed {{
                background: {primary_button_pressed};
                border-color: {primary_button_border};
            }}

            QPushButton[primaryAction="true"]:disabled {{
                background: {primary_button_disabled_bg};
                color: {primary_button_disabled_text};
                border-color: {primary_button_disabled_border};
            }}

            QComboBox QAbstractItemView {{
                background: {panel_bg};
                color: {panel_text};
                border: 1px solid {panel_border};
                outline: 0;
                selection-background-color: #2f80ff;
                selection-color: white;
            }}

            QComboBox QAbstractItemView::item {{
                background: {panel_bg};
                color: {panel_text};
                min-height: 28px;
                padding: 4px 10px;
            }}

            QComboBox QAbstractItemView::item:hover {{
                background: {tab_hover_bg};
                color: {panel_text};
            }}

            QComboBox QAbstractItemView::item:selected {{
                background: #2f80ff;
                color: white;
            }}

            QTableWidget, QTableView {{
                background: {panel_bg};
                color: {panel_text};
                border: 1px solid {panel_border};
                gridline-color: {panel_border};
                alternate-background-color: {table_alt};
            }}

            QTableView::viewport, QTableWidget::viewport {{
                background: {panel_bg};
            }}

            QHeaderView::section {{
                background: {panel_bg};
                color: {panel_text};
                border: 1px solid {panel_border};
                padding: 6px 8px;
            }}

            QTableCornerButton::section {{
                background: {panel_bg};
                border: 1px solid {panel_border};
            }}

            QTableWidget::item:selected, QTableView::item:selected {{
                background: #2f80ff;
                color: white;
            }}

            QScrollArea {{
                background: {panel_bg};
                border: 1px solid {panel_border};
                border-radius: 10px;
            }}

            QScrollArea QWidget {{
                background: transparent;
            }}

            QProgressBar {{
                background: {input_bg};
                color: {muted_text};
                border: 1px solid {panel_border};
                border-radius: 4px;
                text-align: center;
            }}

            QProgressBar::chunk {{
                background: #2f80ff;
                border-radius: 3px;
            }}

            QProgressBar[role="downloadProgress"] {{
                font-size: 10px;
            }}

            QFrame[role="divider"] {{
                background: {panel_border};
            }}

            QLabel[role="downloadStatus"] {{
                font-size: 12px;
                color: {muted_text};
            }}

            QLabel[role="downloadStatus"][state="active"] {{
                color: #2f80ff;
                font-weight: 600;
            }}

            QScrollBar:vertical {{
                background: {scrollbar_bg};
                width: 10px;
                border-radius: 5px;
            }}

            QScrollBar::handle:vertical {{
                background: {scrollbar_handle};
                border-radius: 5px;
                min-height: 30px;
            }}

            QScrollBar::handle:vertical:hover {{
                background: {scrollbar_handle_hover};
            }}

            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
                height: 0;
            }}

            QScrollBar:horizontal {{
                background: {scrollbar_bg};
                height: 10px;
                border-radius: 5px;
            }}

            QScrollBar::handle:horizontal {{
                background: {scrollbar_handle};
                border-radius: 5px;
                min-width: 30px;
            }}

            QScrollBar::handle:horizontal:hover {{
                background: {scrollbar_handle_hover};
            }}

            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
                width: 0;
            }}

            QDialog, QMessageBox {{
                background-color: {panel_bg};
                color: {panel_text};
            }}

            QDialog QWidget, QMessageBox QWidget {{
                background-color: {panel_bg};
                color: {panel_text};
            }}

            QDialog QLabel, QMessageBox QLabel {{
                background-color: {panel_bg};
                color: {panel_text};
            }}

            QMessageBox * {{
                background-color: {panel_bg};
                color: {panel_text};
            }}

            QMessageBox QPushButton, QDialog QPushButton {{
                background-color: {secondary_button_bg};
                color: {secondary_button_text};
                border: 1px solid {secondary_button_border};
                border-radius: 6px;
                padding: 6px 18px;
                min-width: 80px;
            }}

            QMessageBox QPushButton:hover, QDialog QPushButton:hover {{
                background-color: {secondary_button_hover};
                border-color: {secondary_button_text};
            }}

            QMessageBox QPushButton:pressed, QDialog QPushButton:pressed {{
                background-color: {secondary_button_pressed};
            }}
            
            /* 状态栏 */
            QStatusBar {{
                background: #f3f4f6;
                color: #666666;
            }}
        """

        app = QApplication.instance()
        if app:
            theme_key = "dark" if self._dark_mode else "light"
            app_stylesheet = self._theme_app_stylesheet_cache.get(theme_key)
            if app_stylesheet is None:
                app_stylesheet = f"{self._base_app_stylesheet}\n{theme_stylesheet}"
                self._theme_app_stylesheet_cache[theme_key] = app_stylesheet
            if app.styleSheet() != app_stylesheet:
                self.setUpdatesEnabled(False)
                try:
                    app.setStyleSheet(app_stylesheet)
                finally:
                    self.setUpdatesEnabled(True)
                    self.update()
