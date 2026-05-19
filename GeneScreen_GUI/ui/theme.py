"""
Central theme tokens for the GeneScreen desktop UI.

The application stylesheet is installed once and uses Qt palette roles. Theme
switching should only update the palette and a few custom-painted widgets.
"""
from dataclasses import dataclass

from PySide6.QtGui import QColor, QPalette


@dataclass(frozen=True)
class ThemeTokens:
    name: str
    theme_bg: str
    theme_border: str
    panel_bg: str
    panel_border: str
    panel_text: str
    muted_text: str
    sidebar_bg: str
    nav_text: str
    nav_hover_bg: str
    nav_hover_text: str
    input_bg: str
    input_border: str
    input_text: str
    input_placeholder: str
    table_alt: str
    secondary_button_bg: str
    secondary_button_hover: str
    secondary_button_pressed: str
    secondary_button_border: str
    secondary_button_text: str
    readonly_text: str
    scrollbar_bg: str
    scrollbar_handle: str
    scrollbar_handle_hover: str
    theme_toggle_bg: str
    theme_toggle_hover: str
    primary_button_bg: str
    primary_button_hover: str
    primary_button_pressed: str
    primary_button_border: str
    primary_button_text: str
    primary_button_disabled_bg: str
    primary_button_disabled_text: str
    primary_button_disabled_border: str


LIGHT_THEME = ThemeTokens(
    name="light",
    theme_bg="#f8f9fb",
    theme_border="#e6e8eb",
    panel_bg="#ffffff",
    panel_border="#e6e8eb",
    panel_text="#303133",
    muted_text="#666666",
    sidebar_bg="#f8f9fb",
    nav_text="#606266",
    nav_hover_bg="#f2f4f7",
    nav_hover_text="#303133",
    input_bg="#ffffff",
    input_border="#dcdfe6",
    input_text="#303133",
    input_placeholder="#a0a4ad",
    table_alt="#f7f8fa",
    secondary_button_bg="#f3f4f6",
    secondary_button_hover="#e9edf5",
    secondary_button_pressed="#dde3ee",
    secondary_button_border="#cfd6e4",
    secondary_button_text="#2f80ff",
    readonly_text="#8c9197",
    scrollbar_bg="#f0f2f5",
    scrollbar_handle="#c4c8cf",
    scrollbar_handle_hover="#aeb4bd",
    theme_toggle_bg="#d1d5db",
    theme_toggle_hover="#c7cbd3",
    primary_button_bg="#2f80ff",
    primary_button_hover="#1f72e8",
    primary_button_pressed="#195fc4",
    primary_button_border="#1f66d1",
    primary_button_text="#ffffff",
    primary_button_disabled_bg="#d6dde8",
    primary_button_disabled_text="#8c98aa",
    primary_button_disabled_border="#c5cfdd",
)


DARK_THEME = ThemeTokens(
    name="dark",
    theme_bg="#202124",
    theme_border="#2b2b2b",
    panel_bg="#1f2329",
    panel_border="#2b2f36",
    panel_text="#e5e7eb",
    muted_text="#9aa0a6",
    sidebar_bg="#1b1c1e",
    nav_text="#c7c9cc",
    nav_hover_bg="#2a2b2f",
    nav_hover_text="#ffffff",
    input_bg="#1b1f25",
    input_border="#2f343c",
    input_text="#e5e7eb",
    input_placeholder="#8b9099",
    table_alt="#1a1d22",
    secondary_button_bg="#242a32",
    secondary_button_hover="#2d3642",
    secondary_button_pressed="#354052",
    secondary_button_border="#3a4657",
    secondary_button_text="#8ab4ff",
    readonly_text="#9aa0a6",
    scrollbar_bg="#1b1f25",
    scrollbar_handle="#3b4048",
    scrollbar_handle_hover="#4b515b",
    theme_toggle_bg="#3a3b3f",
    theme_toggle_hover="#4a4b50",
    primary_button_bg="#2f80ff",
    primary_button_hover="#4c94ff",
    primary_button_pressed="#1f66d1",
    primary_button_border="#5aa0ff",
    primary_button_text="#ffffff",
    primary_button_disabled_bg="#303844",
    primary_button_disabled_text="#7b8594",
    primary_button_disabled_border="#3d4654",
)


def get_theme_tokens(dark_mode: bool) -> ThemeTokens:
    return DARK_THEME if dark_mode else LIGHT_THEME


def build_palette(tokens: ThemeTokens) -> QPalette:
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(tokens.theme_bg))
    palette.setColor(QPalette.WindowText, QColor(tokens.panel_text))
    palette.setColor(QPalette.Base, QColor(tokens.input_bg))
    palette.setColor(QPalette.AlternateBase, QColor(tokens.table_alt))
    palette.setColor(QPalette.ToolTipBase, QColor(tokens.panel_bg))
    palette.setColor(QPalette.ToolTipText, QColor(tokens.panel_text))
    palette.setColor(QPalette.Text, QColor(tokens.input_text))
    palette.setColor(QPalette.Button, QColor(tokens.secondary_button_bg))
    palette.setColor(QPalette.ButtonText, QColor(tokens.secondary_button_text))
    palette.setColor(QPalette.BrightText, QColor("#ffffff"))
    palette.setColor(QPalette.Link, QColor(tokens.primary_button_hover))
    palette.setColor(QPalette.LinkVisited, QColor(tokens.primary_button_pressed))
    palette.setColor(QPalette.Highlight, QColor(tokens.primary_button_bg))
    palette.setColor(QPalette.HighlightedText, QColor(tokens.primary_button_text))
    palette.setColor(QPalette.PlaceholderText, QColor(tokens.input_placeholder))

    palette.setColor(QPalette.Light, QColor(tokens.nav_hover_bg))
    palette.setColor(QPalette.Midlight, QColor(tokens.secondary_button_hover))
    palette.setColor(QPalette.Mid, QColor(tokens.input_border))
    palette.setColor(QPalette.Dark, QColor(tokens.panel_border))
    palette.setColor(QPalette.Shadow, QColor(tokens.sidebar_bg))

    palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor(tokens.readonly_text))
    palette.setColor(QPalette.Disabled, QPalette.Text, QColor(tokens.readonly_text))
    palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(tokens.readonly_text))
    palette.setColor(QPalette.Disabled, QPalette.Base, QColor(tokens.input_bg))
    palette.setColor(QPalette.Disabled, QPalette.Button, QColor(tokens.input_border))
    return palette


STATIC_THEME_STYLESHEET = """
/* Fixed theme stylesheet. Runtime colors come from QApplication.palette(). */
#titleBar {
    background: palette(window);
    border-bottom: 1px solid palette(dark);
}

#titleBarTitle {
    font-size: 17px;
    font-weight: 700;
    color: palette(window-text);
}

#titleBarIcon {
    background: transparent;
}

#titleButton {
    border: none;
    border-radius: 13px;
}

#titleButton[titleRole="theme"] {
    background: palette(mid);
}

#titleButton[titleRole="theme"]:hover {
    background: palette(midlight);
}

#titleButton[titleRole="min"] {
    background: #f6c344;
}

#titleButton[titleRole="max"] {
    background: #5ac85a;
}

#titleButton[titleRole="close"] {
    background: #ff5f56;
}

#titleButton[titleRole="min"]:hover {
    background: #eab736;
}

#titleButton[titleRole="max"]:hover {
    background: #4bb94d;
}

#titleButton[titleRole="close"]:hover {
    background: #ef4b41;
}

#titleButton[titleRole="min"]:pressed {
    background: #d5a62f;
}

#titleButton[titleRole="max"]:pressed {
    background: #3ea041;
}

#titleButton[titleRole="close"]:pressed {
    background: #d9423a;
}

#windowRoot {
    background: transparent;
}

#sidebar {
    background: palette(shadow);
    border-right: 1px solid palette(dark);
}

#sidebarTitle {
    font-size: 15px;
    font-weight: 600;
    color: palette(window-text);
}

#sidebarVersion {
    font-size: 12px;
    color: palette(text);
}

#sidebarToggle {
    background: transparent;
    border: none;
    border-radius: 6px;
}

#sidebarToggle:hover {
    background: palette(light);
}

#sidebarToggle:pressed {
    background: palette(dark);
}

#sidebarHeaderLine {
    background: palette(dark);
    max-height: 1px;
    margin: 10px 0 6px 0;
}

#sidebarLine {
    background: palette(dark);
    max-height: 1px;
    margin: 8px 0;
}

#navButton, #navButtonSecondary {
    background: transparent;
    border: none;
    border-radius: 10px;
    color: palette(text);
    font-size: 14px;
    text-align: left;
    padding: 10px 12px;
}

#navButton[collapsed="true"], #navButtonSecondary[collapsed="true"] {
    text-align: center;
    padding: 10px 0px;
}

#navButton:hover, #navButtonSecondary:hover {
    background: palette(light);
    color: palette(window-text);
}

#navButton:checked, #navButtonSecondary:checked {
    background: palette(highlight);
    color: palette(highlighted-text);
    font-weight: 600;
}

#contentArea {
    background: palette(window);
}

#pageStack {
    background: palette(window);
}

QLabel {
    color: palette(window-text);
}

QLabel[role="pageTitle"] {
    font-size: 24px;
    font-weight: 700;
    color: palette(window-text);
}

QLabel[role="pageDesc"],
QLabel[role="muted"],
QLabel[role="searchStatus"],
QLabel[role="downloadStatus"] {
    color: palette(placeholder-text);
}

QLabel[role="pageDesc"] {
    font-size: 14px;
}

QLabel[role="fieldLabel"] {
    font-weight: 500;
    color: palette(window-text);
}

QLabel[role="searchStatus"] {
    font-size: 18px;
    background: palette(base);
}

QLabel[role="downloadStatus"] {
    font-size: 12px;
}

QLabel[role="downloadStatus"][state="active"] {
    color: palette(highlight);
    font-weight: 600;
}

QListWidget[role="selectedQueryList"] {
    background: palette(base);
    color: palette(window-text);
    border: 1px solid palette(mid);
    border-radius: 6px;
    padding: 4px;
    outline: 0;
    selection-background-color: palette(light);
    selection-color: palette(window-text);
}

QListWidget[role="selectedQueryList"]::item {
    min-height: 26px;
    padding: 4px 8px;
    border-radius: 4px;
}

QListWidget[role="selectedQueryList"]::item:selected,
QListWidget[role="selectedQueryList"]::item:selected:active,
QListWidget[role="selectedQueryList"]::item:selected:!active {
    background: palette(light);
    color: palette(window-text);
    border: 1px solid palette(highlight);
}

QCheckBox {
    color: palette(window-text);
    min-height: 28px;
    spacing: 6px;
}

QCheckBox::indicator {
    width: 16px;
    height: 16px;
    border: 1px solid palette(mid);
    border-radius: 4px;
    background: palette(base);
}

QCheckBox::indicator:hover {
    border-color: palette(button-text);
}

QCheckBox::indicator:checked {
    background: palette(highlight);
    border: 1px solid palette(highlight);
}

QCheckBox::indicator:disabled {
    background: palette(button);
    border-color: palette(dark);
}

QCheckBox#tableCheckBox {
    min-width: 22px;
    max-width: 22px;
    min-height: 22px;
    max-height: 22px;
    spacing: 0px;
    padding: 0px;
    margin: 0px;
}

QCheckBox#tableCheckBox::indicator {
    width: 18px;
    height: 18px;
    margin: 2px;
}

QTabWidget::pane {
    background: palette(base);
    border: 1px solid palette(dark);
    border-radius: 10px;
    top: -1px;
    padding: 8px;
}

QTabBar::tab {
    background: palette(base);
    color: palette(text);
    border: 1px solid palette(dark);
    border-bottom: none;
    padding: 6px 14px;
    margin-right: 6px;
    border-top-left-radius: 8px;
    border-top-right-radius: 8px;
}

QTabBar::tab:selected {
    background: palette(base);
    color: palette(window-text);
    font-weight: 600;
}

QTabBar::tab:hover {
    background: palette(light);
}

QGroupBox {
    background: palette(base);
    border: 1px solid palette(dark);
    border-radius: 10px;
    margin-top: 10px;
}

QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 0 6px;
    color: palette(window-text);
}

QGroupBox QLabel {
    color: palette(window-text);
}

QLineEdit, QTextEdit, QPlainTextEdit, QSpinBox, QDoubleSpinBox, QComboBox {
    background: palette(base);
    color: palette(text);
    border: 1px solid palette(mid);
    border-radius: 6px;
    padding: 6px 10px;
}

QLineEdit:disabled, QTextEdit:disabled, QPlainTextEdit:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled, QComboBox:disabled {
    background: palette(base);
    color: palette(placeholder-text);
    border: 1px solid palette(mid);
}

QSpinBox[paramInput="true"], QDoubleSpinBox[paramInput="true"], QLineEdit[paramInput="true"], QComboBox[paramInput="true"] {
    border-radius: 5px;
    padding: 4px 8px;
    min-height: 28px;
    min-width: 88px;
}

QLineEdit:focus, QTextEdit:focus, QPlainTextEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {
    border-color: palette(highlight);
}

QLineEdit[role="readonly"] {
    background: palette(base);
    color: palette(placeholder-text);
}

QPushButton[secondary="true"] {
    background: palette(button);
    color: palette(button-text);
    border: 1px solid palette(dark);
}

QPushButton[secondary="true"]:hover {
    background: palette(midlight);
    border-color: palette(button-text);
}

QPushButton[secondary="true"]:pressed {
    background: palette(light);
}

QPushButton[compactAction="true"] {
    padding: 6px 14px;
}

QPushButton[smallAction="true"] {
    font-size: 11px;
    padding: 2px 6px;
}

QPushButton[primaryAction="true"] {
    background: palette(highlight);
    color: palette(highlighted-text);
    border: 1px solid palette(highlight);
    border-radius: 8px;
    font-weight: 600;
}

QPushButton[primaryAction="true"]:hover {
    background: palette(link);
    border-color: palette(link);
}

QPushButton[primaryAction="true"]:pressed {
    background: palette(link-visited);
    border-color: palette(link-visited);
}

QPushButton[primaryAction="true"]:disabled {
    background: palette(button);
    color: palette(placeholder-text);
    border-color: palette(dark);
}

QComboBox QAbstractItemView {
    background: palette(base);
    color: palette(window-text);
    border: 1px solid palette(dark);
    outline: 0;
    selection-background-color: palette(highlight);
    selection-color: palette(highlighted-text);
}

QComboBox QAbstractItemView::item {
    background: palette(base);
    color: palette(window-text);
    min-height: 28px;
    padding: 4px 10px;
}

QComboBox QAbstractItemView::item:hover {
    background: palette(light);
    color: palette(window-text);
}

QComboBox QAbstractItemView::item:selected {
    background: palette(highlight);
    color: palette(highlighted-text);
}

QListWidget[role="geneIdPopup"] {
    background: palette(base);
    color: palette(window-text);
    border: 1px solid palette(dark);
    border-radius: 6px;
    padding: 4px;
    outline: 0;
}

QListWidget[role="geneIdPopup"]::item {
    background: palette(base);
    color: palette(window-text);
    min-height: 28px;
    padding: 4px 10px;
    border-radius: 4px;
}

QListWidget[role="geneIdPopup"]::item:hover,
QListWidget[role="geneIdPopup"]::item:selected,
QListWidget[role="geneIdPopup"]::item:selected:active,
QListWidget[role="geneIdPopup"]::item:selected:!active {
    background: palette(highlight);
    color: palette(highlighted-text);
}

QTableWidget, QTableView {
    background: palette(base);
    color: palette(window-text);
    border: 1px solid palette(dark);
    gridline-color: palette(dark);
    alternate-background-color: palette(alternate-base);
}

QTableView::viewport, QTableWidget::viewport {
    background: palette(base);
}

QHeaderView::section {
    background: palette(base);
    color: palette(window-text);
    border: 1px solid palette(dark);
    padding: 6px 8px;
}

QTableCornerButton::section {
    background: palette(base);
    border: 1px solid palette(dark);
}

QTableWidget::item:selected, QTableView::item:selected {
    background: palette(highlight);
    color: palette(highlighted-text);
}

QScrollArea {
    background: palette(base);
    border: 1px solid palette(dark);
    border-radius: 10px;
}

QScrollArea QWidget {
    background: transparent;
}

QProgressBar {
    background: palette(base);
    color: palette(placeholder-text);
    border: 1px solid palette(dark);
    border-radius: 4px;
    text-align: center;
}

QProgressBar::chunk {
    background: palette(highlight);
    border-radius: 3px;
}

QProgressBar[role="downloadProgress"] {
    font-size: 10px;
}

QFrame[role="divider"] {
    background: palette(dark);
}

QScrollBar:vertical {
    background: palette(base);
    width: 10px;
    border-radius: 5px;
}

QScrollBar::handle:vertical {
    background: palette(mid);
    border-radius: 5px;
    min-height: 30px;
}

QScrollBar::handle:vertical:hover {
    background: palette(midlight);
}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0;
}

QScrollBar:horizontal {
    background: palette(base);
    height: 10px;
    border-radius: 5px;
}

QScrollBar::handle:horizontal {
    background: palette(mid);
    border-radius: 5px;
    min-width: 30px;
}

QScrollBar::handle:horizontal:hover {
    background: palette(midlight);
}

QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
    width: 0;
}

QDialog, QMessageBox {
    background-color: palette(base);
    color: palette(window-text);
}

QDialog QWidget, QMessageBox QWidget {
    background-color: palette(base);
    color: palette(window-text);
}

QDialog QLabel, QMessageBox QLabel {
    background-color: palette(base);
    color: palette(window-text);
}

QMessageBox * {
    background-color: palette(base);
    color: palette(window-text);
}

QMessageBox QPushButton, QDialog QPushButton {
    background-color: palette(button);
    color: palette(button-text);
    border: 1px solid palette(dark);
    border-radius: 6px;
    padding: 6px 18px;
    min-width: 80px;
}

QMessageBox QPushButton:hover, QDialog QPushButton:hover {
    background-color: palette(midlight);
    border-color: palette(button-text);
}

QMessageBox QPushButton:pressed, QDialog QPushButton:pressed {
    background-color: palette(light);
}

QStatusBar {
    background: palette(button);
    color: palette(text);
}

QToolTip {
    background-color: palette(base);
    color: palette(window-text);
    border: 1px solid palette(dark);
    border-radius: 4px;
    padding: 4px 8px;
}
"""
