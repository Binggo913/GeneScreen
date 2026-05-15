#!/usr/bin/env python
"""
Nuitka 打包脚本 - GeneScreen 1.0
自动识别平台，生成 GitHub Actions 使用的 GUI 构建目录

用法:
  python build.py
"""
import platform
import subprocess
import sys
from pathlib import Path


APP_NAME = "GeneScreen"
APP_VERSION = "1.0"


def build_nuitka():
    """Nuitka 编译"""
    system = platform.system()

    cmd = [
        sys.executable, "-m", "nuitka",
        "--standalone",
        "--enable-plugin=pyside6",
        "--output-dir=dist",
        "--include-data-dir=ui/resources=ui/resources",
        # LINKVIEW.py 作为外部脚本需要单独打包（subprocess 调用）
        "--include-data-files=core/LINKVIEW.py=core/LINKVIEW.py",
        "--include-data-files=core/interval.py=core/interval.py",
        "--assume-yes-for-downloads",
        "--remove-output",
    ]

    if system == "Windows":
        cmd.extend([
            "--windows-console-mode=disable",
            f"--output-filename={APP_NAME}.exe",
        ])
        icon_path = Path("ui/resources/icons/app.ico")
        if icon_path.exists():
            cmd.append(f"--windows-icon-from-ico={icon_path}")
    elif system == "Darwin":
        cmd.extend([
            "--macos-create-app-bundle",
            f"--macos-app-name={APP_NAME}",
            f"--macos-app-version={APP_VERSION}",
        ])
    else:
        cmd.extend([
            "--static-libpython=no",
            f"--output-filename={APP_NAME}",
        ])

    cmd.append("main.py")

    print(f"Platform: {system}")
    print("Starting Nuitka build...")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    print("\nNuitka build completed.")


def main():
    build_nuitka()


if __name__ == "__main__":
    main()
