#!/usr/bin/env python
"""
Nuitka 打包脚本 - GeneScreen 1.0
自动识别平台，生成 GitHub Actions 使用的 GUI 构建目录

用法:
  python build.py
"""
import platform
import shutil
import subprocess
import sys
from pathlib import Path


APP_NAME = "GeneScreen"
APP_VERSION = "1.0"


WINDOWS_QT_RUNTIME_DLLS = {
    "concrt140.dll",
    "msvcp140.dll",
    "msvcp140_1.dll",
    "msvcp140_2.dll",
    "msvcp140_codecvt_ids.dll",
    "pyside6.abi3.dll",
    "qt6core.dll",
    "qt6gui.dll",
    "qt6network.dll",
    "qt6pdf.dll",
    "qt6svg.dll",
    "qt6widgets.dll",
    "shiboken6.abi3.dll",
    "vcamp140.dll",
    "vccorlib140.dll",
    "vcomp140.dll",
    "vcruntime140.dll",
    "vcruntime140_1.dll",
}


def _copy_named_dlls(source_dir: Path, target_dir: Path, names: set[str]):
    if not source_dir.exists():
        return

    wanted = {name.lower() for name in names}
    for dll_path in source_dir.glob("*.dll"):
        if dll_path.name.lower() in wanted:
            target_path = target_dir / dll_path.name
            shutil.copy2(dll_path, target_path)
            print(f"Included runtime DLL: {target_path}")


def _include_windows_qt_runtime_dlls():
    """Copy PySide6/Shiboken runtime DLLs that Nuitka may miss on Windows."""
    dist_dirs = sorted(Path("dist").glob("*.dist"))
    if not dist_dirs:
        raise RuntimeError("Nuitka .dist directory not found after build")
    target_dir = dist_dirs[0]

    try:
        import PySide6
        import shiboken6
    except ImportError as exc:
        raise RuntimeError("PySide6/Shiboken is required for Windows packaging") from exc

    _copy_named_dlls(Path(PySide6.__file__).parent, target_dir, WINDOWS_QT_RUNTIME_DLLS)
    _copy_named_dlls(Path(shiboken6.__file__).parent, target_dir, WINDOWS_QT_RUNTIME_DLLS)

    python3_dll = Path(sys.base_prefix) / "python3.dll"
    if python3_dll.exists():
        target_path = target_dir / python3_dll.name
        shutil.copy2(python3_dll, target_path)
        print(f"Included runtime DLL: {target_path}")


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
    if system == "Windows":
        _include_windows_qt_runtime_dlls()
    print("\nNuitka build completed.")


def main():
    build_nuitka()


if __name__ == "__main__":
    main()
