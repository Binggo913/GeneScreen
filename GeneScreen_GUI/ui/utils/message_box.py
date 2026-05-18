"""Application-themed QMessageBox helpers."""

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication, QMessageBox, QWidget


def _standard_button(name: str):
    if hasattr(QMessageBox, name):
        return getattr(QMessageBox, name)
    return getattr(QMessageBox.StandardButton, name)


def _message_option(name: str):
    if hasattr(QMessageBox, name):
        return getattr(QMessageBox, name)
    return getattr(getattr(QMessageBox, "Option", object), name, None)


def _message_icon(name: str):
    if hasattr(QMessageBox, name):
        return getattr(QMessageBox, name)
    return getattr(QMessageBox.Icon, name)


def _set_non_native(message_box: QMessageBox) -> None:
    option = _message_option("DontUseNativeDialog")
    if option is not None:
        message_box.setOption(option, True)


def _dialog_colors(message_box: QMessageBox) -> dict:
    app = QApplication.instance()
    palette = app.palette() if app else message_box.palette()
    panel = palette.color(QPalette.ToolTipBase)
    text = palette.color(QPalette.ToolTipText)
    button = palette.color(QPalette.Button)
    button_text = palette.color(QPalette.Link)
    border = QColor("#3a4657") if panel.lightness() < 128 else QColor("#cfd6e4")
    hover = QColor("#2d3642") if panel.lightness() < 128 else QColor("#e9edf5")
    pressed = QColor("#354052") if panel.lightness() < 128 else QColor("#dde3ee")
    return {
        "panel": panel.name(),
        "text": text.name(),
        "button": button.name(),
        "button_text": button_text.name(),
        "border": border.name(),
        "hover": hover.name(),
        "pressed": pressed.name(),
    }


def apply_message_box_theme(message_box: QMessageBox) -> QMessageBox:
    _set_non_native(message_box)
    colors = _dialog_colors(message_box)
    message_box.setAttribute(Qt.WA_StyledBackground, True)
    message_box.setAutoFillBackground(True)
    message_box.setStyleSheet(f"""
        QMessageBox {{
            background-color: {colors["panel"]};
            color: {colors["text"]};
        }}
        QMessageBox QWidget {{
            background-color: {colors["panel"]};
            color: {colors["text"]};
        }}
        QMessageBox QLabel {{
            background-color: {colors["panel"]};
            color: {colors["text"]};
        }}
        QMessageBox QPushButton {{
            background-color: {colors["button"]};
            color: {colors["button_text"]};
            border: 1px solid {colors["border"]};
            border-radius: 6px;
            padding: 6px 18px;
            min-width: 80px;
        }}
        QMessageBox QPushButton:hover {{
            background-color: {colors["hover"]};
            border-color: {colors["button_text"]};
        }}
        QMessageBox QPushButton:pressed {{
            background-color: {colors["pressed"]};
        }}
    """)

    palette = message_box.palette()
    palette.setColor(QPalette.Window, QColor(colors["panel"]))
    palette.setColor(QPalette.Base, QColor(colors["panel"]))
    palette.setColor(QPalette.WindowText, QColor(colors["text"]))
    palette.setColor(QPalette.Text, QColor(colors["text"]))
    message_box.setPalette(palette)
    for child in message_box.findChildren(QWidget):
        child.setAttribute(Qt.WA_StyledBackground, True)
        child.setAutoFillBackground(True)
        child.setPalette(palette)
    return message_box


def show_message_box(parent, icon, title, text, buttons=None, default_button=None):
    if buttons is None:
        buttons = _standard_button("Ok")
    message_box = QMessageBox(parent)
    message_box.setIcon(icon)
    message_box.setWindowTitle(title)
    message_box.setText(text)
    message_box.setStandardButtons(buttons)
    if default_button is not None and default_button != _standard_button("NoButton"):
        message_box.setDefaultButton(default_button)
    apply_message_box_theme(message_box)
    return message_box.exec()


def install_themed_message_boxes() -> None:
    if getattr(QMessageBox, "_genescreen_themed_static_methods", False):
        return

    ok = _standard_button("Ok")
    no_button = _standard_button("NoButton")
    yes = _standard_button("Yes")
    no = _standard_button("No")

    def warning(parent, title, text, buttons=ok, defaultButton=no_button):
        return show_message_box(parent, _message_icon("Warning"), title, text, buttons, defaultButton)

    def information(parent, title, text, buttons=ok, defaultButton=no_button):
        return show_message_box(parent, _message_icon("Information"), title, text, buttons, defaultButton)

    def critical(parent, title, text, buttons=ok, defaultButton=no_button):
        return show_message_box(parent, _message_icon("Critical"), title, text, buttons, defaultButton)

    def question(parent, title, text, buttons=(yes | no), defaultButton=no_button):
        return show_message_box(parent, _message_icon("Question"), title, text, buttons, defaultButton)

    QMessageBox.warning = staticmethod(warning)
    QMessageBox.information = staticmethod(information)
    QMessageBox.critical = staticmethod(critical)
    QMessageBox.question = staticmethod(question)
    QMessageBox._genescreen_themed_static_methods = True
