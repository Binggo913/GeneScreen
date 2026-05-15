# GeneScreen 1.0 Utils
import subprocess
import sys

# Windows 下隐藏 subprocess 弹出的黑窗口
SUBPROCESS_FLAGS = subprocess.CREATE_NO_WINDOW if sys.platform == 'win32' else 0


def run_subprocess(cmd, capture_output=True, text=True, timeout=None, shell=False, **kwargs):
    """
    封装 subprocess.run，Windows 下自动隐藏黑窗口
    """
    return subprocess.run(
        cmd,
        capture_output=capture_output,
        text=text,
        timeout=timeout,
        shell=shell,
        creationflags=SUBPROCESS_FLAGS,
        **kwargs
    )
