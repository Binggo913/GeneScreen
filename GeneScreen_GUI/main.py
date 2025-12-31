"""
GeneScreen 3.0 - 基因组比对与变异分析工具

主程序入口
"""
import sys
import time
_start = time.time()

import ctypes
from pathlib import Path

from PySide6.QtWidgets import QApplication, QMessageBox, QLabel
from PySide6.QtGui import QIcon
from PySide6.QtCore import Qt
print(f"[PERF] PySide6 import: {time.time()-_start:.2f}s")

# 添加项目根目录到 Python 路径
ROOT_DIR = Path(__file__).parent
sys.path.insert(0, str(ROOT_DIR))

_t = time.time()
from utils.blast_check import check_blast, get_blast_version, check_blast_installation
print(f"[PERF] blast_check import: {time.time()-_t:.2f}s")

_t = time.time()
from ui.main_window import MainWindow
print(f"[PERF] MainWindow import: {time.time()-_t:.2f}s")
print(f"[PERF] 总 import 时间: {time.time()-_start:.2f}s")


def main():
    """主程序入口"""
    if sys.platform.startswith("win"):
        try:
            app_id = "GeneScreen.App"
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(app_id)
        except Exception:
            pass

    # 创建应用
    app = QApplication(sys.argv)
    app.setApplicationName("GeneScreen")
    app.setApplicationVersion("1.0.0")
    app.setOrganizationName("GeneScreen")

    # 设置应用图标
    icon_path = ROOT_DIR / "ui" / "resources" / "icons" / "app.ico"
    if icon_path.exists():
        icon = QIcon(str(icon_path))
        app.setWindowIcon(icon)

    # 检查 BLAST+ 是否安装
    _t = time.time()
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
    print(f"[PERF] check_blast: {time.time()-_t:.2f}s")

    # 显示 BLAST 版本信息
    _t = time.time()
    blast_version = get_blast_version()
    print(f"[INFO] 检测到 BLAST+ 版本: {blast_version}")
    print(f"[PERF] get_blast_version: {time.time()-_t:.2f}s")

    # 加载样式表（如果存在）
    _t = time.time()
    style_path = ROOT_DIR / "ui" / "resources" / "styles.qss"
    if style_path.exists():
        with open(style_path, "r", encoding="utf-8") as f:
            app.setStyleSheet(f.read())
    print(f"[PERF] load_stylesheet: {time.time()-_t:.2f}s")

    # 创建并显示主窗口
    _t = time.time()
    window = MainWindow()
    print(f"[PERF] MainWindow(): {time.time()-_t:.2f}s")
    
    _t = time.time()
    if icon_path.exists():
        window.setWindowIcon(icon)
    window.show()
    print(f"[PERF] window.show(): {time.time()-_t:.2f}s")
    print(f"[PERF] 总启动时间: {time.time()-_start:.2f}s")

    # 运行应用
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
