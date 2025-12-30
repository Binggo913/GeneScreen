#!/usr/bin/env python
"""
Nuitka 打包脚本 - GeneScreen 3.0
自动识别平台，生成对应的可执行文件

用法: python build.py
"""
import subprocess
import sys
import platform

def build():
    is_windows = platform.system() == "Windows"
    
    cmd = [
        sys.executable, "-m", "nuitka",
        "--standalone",
        "--onefile",
        "--enable-plugin=pyside6",
        "--output-dir=dist",
        "--include-data-dir=ui/resources=ui/resources",
        "--assume-yes-for-downloads",
        "--remove-output",
    ]
    
    if is_windows:
        cmd.extend([
            "--windows-console-mode=disable",
            "--windows-icon-from-ico=ui/resources/icons/app.ico",
            "--output-filename=GeneScreen.exe",
        ])
    else:
        cmd.extend([
            "--static-libpython=no",  # micromamba 兼容
            "--output-filename=GeneScreen",
        ])
    
    cmd.append("main.py")
    
    print(f"平台: {platform.system()}")
    print("开始 Nuitka 编译...")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"\n编译完成！输出: dist/{'GeneScreen.exe' if is_windows else 'GeneScreen'}")

if __name__ == "__main__":
    build()
