"""
GeneScreen 1.0 - 基因组比对与变异分析工具

主程序入口
"""
import sys
import ctypes
from pathlib import Path

from PySide6.QtWidgets import QApplication, QMessageBox, QLabel
from PySide6.QtGui import QIcon
from PySide6.QtCore import Qt

# 添加项目根目录到 Python 路径
ROOT_DIR = Path(__file__).parent
sys.path.insert(0, str(ROOT_DIR))

from utils.blast_check import check_blast, check_blast_installation
from ui.main_window import MainWindow


def _apply_qt_app_attributes():
    attr = getattr(Qt, "AA_DontUseNativeDialogs", None)
    if attr is None:
        attr = getattr(getattr(Qt, "ApplicationAttribute", object), "AA_DontUseNativeDialogs", None)
    if attr is not None:
        QApplication.setAttribute(attr, True)


def main():
    """主程序入口"""
    if sys.platform.startswith("win"):
        try:
            app_id = "GeneScreen.App"
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(app_id)
        except Exception:
            pass

    # 创建应用
    _apply_qt_app_attributes()
    app = QApplication(sys.argv)
    app.setApplicationName("GeneScreen")
    app.setApplicationVersion("1.0")
    app.setOrganizationName("GeneScreen")

    # 设置应用图标
    icon_path = ROOT_DIR / "ui" / "resources" / "icons" / "app.ico"
    if icon_path.exists():
        icon = QIcon(str(icon_path))
        app.setWindowIcon(icon)

    # 检查 BLAST+ 是否安装
    if not check_blast():
        link = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
        _, blast_msg = check_blast_installation()
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setWindowTitle("BLAST+ 未安装")
        msg.setTextFormat(Qt.RichText)
        msg.setText(
            "GeneScreen 需要 BLAST+ 才能运行。<br><br>"
            f"{blast_msg}<br><br>"
            "请先安装 BLAST+：<br>"
            f'<a href="{link}">{link}</a><br><br>'
            "安装完成后，请确保 blastn 命令可在命令行中使用，"
            "然后重新启动 GeneScreen。"
        )
        copy_btn = msg.addButton("复制链接", QMessageBox.ActionRole)
        msg.addButton(QMessageBox.Ok)
        label = msg.findChild(QLabel, "qt_msgbox_label")
        if label:
            label.setTextInteractionFlags(
                Qt.TextSelectableByMouse
                | Qt.TextSelectableByKeyboard
                | Qt.LinksAccessibleByMouse
                | Qt.LinksAccessibleByKeyboard
            )
            label.setOpenExternalLinks(True)
        msg.exec()
        if msg.clickedButton() == copy_btn:
            QApplication.clipboard().setText(link)
        sys.exit(1)

    # 加载样式表（如果存在）
    style_path = ROOT_DIR / "ui" / "resources" / "styles.qss"
    if style_path.exists():
        with open(style_path, "r", encoding="utf-8") as f:
            app.setStyleSheet(f.read())

    # 创建并显示主窗口
    window = MainWindow()

    if icon_path.exists():
        window.setWindowIcon(icon)
    window.show()

    # 运行应用
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
