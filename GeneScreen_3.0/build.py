#!/usr/bin/env python
"""
Nuitka 打包脚本 - GeneScreen 3.0
自动识别平台，生成对应的可执行文件

用法: 
  python build.py          # 仅 Nuitka 编译
  python build.py --setup  # Nuitka 编译 + Inno Setup 安装包 (Windows)
"""
import subprocess
import sys
import platform
import shutil
from pathlib import Path


def build_nuitka():
    """Nuitka 编译"""
    is_windows = platform.system() == "Windows"
    
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
    
    if is_windows:
        cmd.extend([
            "--windows-console-mode=disable",
            "--output-filename=GeneScreen.exe",
        ])
        icon_path = Path("ui/resources/icons/app.ico")
        if icon_path.exists():
            cmd.append(f"--windows-icon-from-ico={icon_path}")
    else:
        cmd.extend([
            "--static-libpython=no",
            "--output-filename=GeneScreen",
        ])
    
    cmd.append("main.py")
    
    print(f"平台: {platform.system()}")
    print("开始 Nuitka 编译...")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    print("\nNuitka 编译完成！")


def build_installer():
    """使用 Inno Setup 生成安装包 (仅 Windows)"""
    if platform.system() != "Windows":
        print("Inno Setup 仅支持 Windows")
        return False
    
    # 查找 Inno Setup 编译器
    iscc_paths = [
        r"D:\Inno Setup 6\ISCC.exe",
        r"C:\Program Files (x86)\Inno Setup 6\ISCC.exe",
        r"C:\Program Files\Inno Setup 6\ISCC.exe",
        shutil.which("ISCC"),
    ]
    
    iscc = None
    for p in iscc_paths:
        if p and Path(p).exists():
            iscc = p
            break
    
    if not iscc:
        print("未找到 Inno Setup，请安装: https://jrsoftware.org/isdl.php")
        return False
    
    iss_file = Path("installer.iss")
    if not iss_file.exists():
        print(f"未找到 {iss_file}")
        return False
    
    print(f"\n使用 Inno Setup 生成安装包...")
    subprocess.run([iscc, str(iss_file)], check=True)
    print("\n安装包生成完成！输出: dist/GeneScreen_Setup_3.0.exe")
    return True


def main():
    # 1. Nuitka 编译
    build_nuitka()
    
    # 2. 如果指定 --setup，生成安装包
    if "--setup" in sys.argv:
        build_installer()
    else:
        print("\n提示: 运行 'python build.py --setup' 可生成 Windows 安装包")


if __name__ == "__main__":
    main()
