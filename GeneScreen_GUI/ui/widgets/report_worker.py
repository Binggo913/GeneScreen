"""
GeneScreen 3.0 - 报告生成线程
"""
from PySide6.QtCore import QThread, Signal

from core.visualizer import generate_report


class ReportWorker(QThread):
    """后台生成报告，避免阻塞 UI"""
    report_ready = Signal(bool, str, str)  # success, report_path, error

    def __init__(self, result: dict, output_dir: str, mode: str, ref_name: str = "", qry_name: str = "", parent=None):
        super().__init__(parent)
        self.result = result
        self.output_dir = output_dir
        self.mode = mode
        self.ref_name = ref_name
        self.qry_name = qry_name

    def run(self):
        try:
            report_path = generate_report(self.result, self.output_dir, self.mode, self.ref_name, self.qry_name)
            self.report_ready.emit(True, report_path, "")
        except Exception as e:
            self.report_ready.emit(False, "", str(e))
