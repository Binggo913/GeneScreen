# GeneScreen 1.0 实现规划

## 目标

为 Windows 用户提供**开箱即用**的 GUI 体验：

- 用户只需安装 BLAST+
- 双击 EXE 即可使用
- 无需安装 Python、Conda 等环境

---

## 技术选型

### GUI 框架：PySide6

| 特性 | 说明 |
|------|------|
| 原生体验 | Qt 框架，响应快速 |
| 跨平台 | Windows/macOS/Linux |
| 组件丰富 | 表格、树形、进度条、文件对话框 |
| 打包友好 | PyInstaller 支持良好 |

### 数据管理：SQLite

| 特性 | 说明 |
|------|------|
| 轻量级 | 单文件数据库，无需安装 |
| 可靠性 | 支持事务、ACID |
| 查询能力 | SQL 查询，支持索引 |
| Python 内置 | 无需额外依赖 |

### 打包方案：PyInstaller

预计打包体积：80-100MB（比 Streamlit 方案小一半）

---

## 数据库设计

### 表结构

```sql
-- 基因组表
CREATE TABLE genomes (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT UNIQUE NOT NULL,           -- 基因组名称（用户指定）
    display_name TEXT,                    -- 显示名称
    fasta_path TEXT NOT NULL,             -- FASTA 文件路径
    annotation_path TEXT,                 -- GFF/GTF 注释文件路径
    source TEXT DEFAULT 'custom',         -- 来源: custom/igv/ensembl
    assembly TEXT,                        -- 组装版本
    species TEXT,                         -- 物种名
    description TEXT,                     -- 描述
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- 分析历史表
CREATE TABLE analysis_history (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    mode TEXT NOT NULL,                   -- gene_id/location/sequence
    ref_genome_id INTEGER,                -- 参考基因组
    qry_genome_id INTEGER,                -- 目标基因组
    input_value TEXT,                     -- 输入值（Gene ID/位置/序列文件）
    identity REAL DEFAULT 90,             -- Identity 阈值
    output_dir TEXT,                      -- 输出目录
    status TEXT DEFAULT 'pending',        -- pending/running/completed/failed
    snp_count INTEGER,                    -- SNP 数量
    indel_count INTEGER,                  -- Indel 数量
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (ref_genome_id) REFERENCES genomes(id),
    FOREIGN KEY (qry_genome_id) REFERENCES genomes(id)
);

-- 索引
CREATE INDEX idx_genomes_name ON genomes(name);
CREATE INDEX idx_history_created ON analysis_history(created_at);
```

### 数据库文件位置

```
~/.genescreen/
├── genescreen.db          # SQLite 数据库
└── genomes/               # 下载的基因组文件
    ├── rice_nipponbare/
    │   ├── genome.fa
    │   └── annotation.gff3
    └── ...
```

---

## 目录结构

```
GeneScreen_1.0/
├── main.py                   # 程序入口
├── core/                     # 核心业务逻辑
│   ├── __init__.py
│   ├── database.py           # SQLite 数据库管理
│   ├── genome_manager.py     # 基因组管理（增删改查）
│   ├── sequence_extractor.py # 序列提取
│   ├── blast_aligner.py      # BLAST 比对
│   ├── variant_caller.py     # 变异检测
│   └── visualizer.py         # 可视化生成
├── ui/                       # PySide6 界面
│   ├── __init__.py
│   ├── main_window.py        # 主窗口
│   ├── widgets/              # 自定义组件
│   │   ├── __init__.py
│   │   ├── genome_selector.py    # 基因组选择器
│   │   ├── analysis_panel.py     # 分析面板
│   │   ├── result_viewer.py      # 结果查看器
│   │   └── genome_manager_dialog.py  # 基因组管理对话框
│   ├── pages/                # 功能页面
│   │   ├── __init__.py
│   │   ├── gene_id_page.py       # Gene ID 模式
│   │   ├── location_page.py      # Location 模式
│   │   ├── sequence_page.py      # Sequence 模式
│   │   └── history_page.py       # 历史记录
│   └── resources/            # 资源文件
│       ├── icons/
│       └── styles.qss        # Qt 样式表
├── utils/                    # 工具函数
│   ├── __init__.py
│   ├── blast_check.py        # BLAST+ 检测
│   ├── file_utils.py         # 文件操作
│   └── download.py           # 下载工具
├── build/                    # 打包相关
│   ├── build.py              # 打包脚本
│   ├── GeneScreen.spec       # PyInstaller 配置
│   └── icon.ico              # 应用图标
├── requirements.txt          # pip 依赖
├── environment.yml           # 开发环境
└── README.md                 # 使用文档
```

### 打包后分发结构

```
GeneScreen_1.0_Windows/
├── GeneScreen.exe            # 主程序
├── README.txt                # 使用说明
└── _internal/                # 依赖文件
```

---

## 界面设计

### 主窗口布局

```
┌─────────────────────────────────────────────────────────────────┐
│  GeneScreen 1.0                                    [─] [□] [×]  │
├─────────────────────────────────────────────────────────────────┤
│  ┌──────────┐  ┌─────────────────────────────────────────────┐  │
│  │          │  │                                             │  │
│  │  导航栏   │  │              内容区域                        │  │
│  │          │  │                                             │  │
│  │ ○ Gene ID │  │  ┌─────────────────────────────────────┐   │  │
│  │ ○ Location│  │  │  参数设置区                          │   │  │
│  │ ○ Sequence│  │  │  [参考基因组 ▼] [目标基因组 ▼]       │   │  │
│  │           │  │  │  [Gene ID: ________________]         │   │  │
│  │ ─────────│  │  │  [Identity: ====●==== 90%]           │   │  │
│  │ ○ 基因组  │  │  │  [输出目录: ________] [浏览]         │   │  │
│  │ ○ 历史    │  │  │  [        开始分析        ]          │   │  │
│  │ ○ 设置    │  │  └─────────────────────────────────────┘   │  │
│  │          │  │                                             │  │
│  │          │  │  ┌─────────────────────────────────────┐   │  │
│  │          │  │  │  结果展示区                          │   │  │
│  │          │  │  │  [比对结果] [变异] [可视化] [报告]    │   │  │
│  │          │  │  │  ┌─────────────────────────────┐    │   │  │
│  │          │  │  │  │                             │    │   │  │
│  │          │  │  │  │     表格 / 图片 / HTML      │    │   │  │
│  │          │  │  │  │                             │    │   │  │
│  │          │  │  │  └─────────────────────────────┘    │   │  │
│  │          │  │  └─────────────────────────────────────┘   │  │
│  └──────────┘  └─────────────────────────────────────────────┘  │
├─────────────────────────────────────────────────────────────────┤
│  状态栏: BLAST+ ✓ | 基因组: 5 个 | 就绪                          │
└─────────────────────────────────────────────────────────────────┘
```

### 基因组管理对话框

```
┌─────────────────────────────────────────────────────────────────┐
│  基因组管理                                              [×]    │
├─────────────────────────────────────────────────────────────────┤
│  ┌─────────────────────────────────────────────────────────┐    │
│  │ 已添加的基因组                                          │    │
│  ├─────────────────────────────────────────────────────────┤    │
│  │ 名称          │ 来源    │ 物种        │ 操作            │    │
│  ├───────────────┼─────────┼─────────────┼─────────────────┤    │
│  │ Nipponbare    │ custom  │ Rice        │ [详情] [删除]   │    │
│  │ ZS97          │ custom  │ Rice        │ [详情] [删除]   │    │
│  │ IRGSP-1.0     │ igv     │ Rice        │ [详情] [删除]   │    │
│  └─────────────────────────────────────────────────────────┘    │
│                                                                 │
│  [+ 添加本地基因组]  [🔍 搜索在线基因组]  [↓ 下载]              │
│                                                                 │
│  ┌─────────────────────────────────────────────────────────┐    │
│  │ 搜索在线基因组: [rice____________] [搜索]               │    │
│  ├─────────────────────────────────────────────────────────┤    │
│  │ ☐ IRGSP-1.0 (Oryza sativa japonica) - IGV              │    │
│  │ ☐ oryza_sativa (Rice) - Ensembl Plants                 │    │
│  │ ☐ oryza_indica (Indica rice) - Ensembl Plants          │    │
│  └─────────────────────────────────────────────────────────┘    │
│                                                                 │
│  下载进度: [████████████░░░░░░░░] 60%  正在下载 FASTA...        │
└─────────────────────────────────────────────────────────────────┘
```

---

## 核心模块设计

### core/database.py

```python
"""SQLite 数据库管理"""
import sqlite3
from pathlib import Path
from contextlib import contextmanager

class Database:
    def __init__(self, db_path=None):
        if db_path is None:
            db_path = Path.home() / ".genescreen" / "genescreen.db"
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_tables()

    @contextmanager
    def connection(self):
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        try:
            yield conn
            conn.commit()
        finally:
            conn.close()

    def _init_tables(self):
        """初始化数据库表"""
        with self.connection() as conn:
            conn.executescript('''
                CREATE TABLE IF NOT EXISTS genomes (...);
                CREATE TABLE IF NOT EXISTS analysis_history (...);
            ''')

    # CRUD 操作
    def add_genome(self, name, fasta_path, annotation_path=None, **kwargs): ...
    def get_genome(self, name): ...
    def list_genomes(self): ...
    def delete_genome(self, name): ...
    def add_history(self, **kwargs): ...
    def get_history(self, limit=50): ...
```

### core/genome_manager.py

```python
"""基因组管理器 - 整合数据库和文件操作"""
from .database import Database

class GenomeManager:
    def __init__(self):
        self.db = Database()
        self.cache_dir = Path.home() / ".genescreen" / "genomes"

    def get(self, name) -> tuple[str, str]:
        """获取基因组路径 (fasta, annotation)"""
        genome = self.db.get_genome(name)
        if genome:
            return genome['fasta_path'], genome['annotation_path']
        return None, None

    def add_custom(self, name, fasta_path, annotation_path=None):
        """添加自定义基因组"""
        # 验证文件存在
        # 创建 FASTA 索引
        # 写入数据库
        pass

    def download_ensembl(self, species_name, progress_callback=None):
        """从 Ensembl Plants 下载"""
        pass

    def download_igv(self, genome_id, progress_callback=None):
        """从 IGV 下载"""
        pass

    def search_online(self, keyword) -> list:
        """搜索在线基因组"""
        pass

    def remove(self, name):
        """删除基因组"""
        pass
```

---

## PySide6 界面实现

### ui/main_window.py

```python
"""主窗口"""
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QStackedWidget, QListWidget, QStatusBar
)
from PySide6.QtCore import Qt

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("GeneScreen 1.0")
        self.setMinimumSize(1200, 800)

        # 中央部件
        central = QWidget()
        self.setCentralWidget(central)
        layout = QHBoxLayout(central)

        # 左侧导航
        self.nav_list = QListWidget()
        self.nav_list.addItems([
            "Gene ID 模式",
            "Location 模式",
            "Sequence 模式",
            "───────────",
            "基因组管理",
            "历史记录",
            "设置"
        ])
        self.nav_list.setMaximumWidth(150)
        self.nav_list.currentRowChanged.connect(self._on_nav_changed)

        # 右侧内容区
        self.stack = QStackedWidget()
        # 添加各个页面...

        layout.addWidget(self.nav_list)
        layout.addWidget(self.stack, 1)

        # 状态栏
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self._update_status()

    def _on_nav_changed(self, index):
        if index != 3:  # 跳过分隔线
            self.stack.setCurrentIndex(index if index < 3 else index - 1)

    def _update_status(self):
        blast_ok = check_blast()
        genome_count = len(genome_manager.list_all())
        self.status_bar.showMessage(
            f"BLAST+ {'✓' if blast_ok else '✗'} | 基因组: {genome_count} 个 | 就绪"
        )
```

### ui/pages/gene_id_page.py

```python
"""Gene ID 模式页面"""
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QComboBox, QLineEdit, QSlider, QPushButton,
    QTabWidget, QTableWidget, QLabel, QTextEdit,
    QFileDialog, QProgressBar
)
from PySide6.QtCore import Qt, QThread, Signal

class AnalysisWorker(QThread):
    """后台分析线程"""
    progress = Signal(str)  # 日志信息
    finished = Signal(dict)  # 结果
    error = Signal(str)

    def __init__(self, ref, qry, gene_id, identity, output_dir):
        super().__init__()
        self.ref = ref
        self.qry = qry
        self.gene_id = gene_id
        self.identity = identity
        self.output_dir = output_dir

    def run(self):
        try:
            # 执行分析...
            self.progress.emit("提取基因序列...")
            # ...
            self.progress.emit("运行 BLAST 比对...")
            # ...
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))


class GeneIDPage(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)

        # 参数区
        form = QFormLayout()
        self.ref_combo = QComboBox()
        self.qry_combo = QComboBox()
        self.gene_id_input = QLineEdit()
        self.identity_slider = QSlider(Qt.Horizontal)
        self.identity_slider.setRange(50, 100)
        self.identity_slider.setValue(90)
        self.output_btn = QPushButton("浏览...")

        form.addRow("参考基因组:", self.ref_combo)
        form.addRow("目标基因组:", self.qry_combo)
        form.addRow("Gene ID:", self.gene_id_input)
        form.addRow("Identity 阈值:", self.identity_slider)
        form.addRow("输出目录:", self.output_btn)

        self.run_btn = QPushButton("开始分析")
        self.run_btn.clicked.connect(self._run_analysis)

        # 日志区
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumHeight(100)

        # 结果区
        self.result_tabs = QTabWidget()
        self.result_tabs.addTab(QTableWidget(), "比对结果")
        self.result_tabs.addTab(QTableWidget(), "变异")
        self.result_tabs.addTab(QLabel(), "可视化")
        self.result_tabs.addTab(QTextEdit(), "报告")

        layout.addLayout(form)
        layout.addWidget(self.run_btn)
        layout.addWidget(self.log_text)
        layout.addWidget(self.result_tabs, 1)

    def _run_analysis(self):
        # 启动后台线程
        self.worker = AnalysisWorker(...)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.start()

    def _on_progress(self, msg):
        self.log_text.append(msg)

    def _on_finished(self, result):
        # 显示结果
        pass
```

---

## 实现步骤

### 阶段 1：项目初始化

1. 创建目录结构
2. 配置 requirements.txt
3. 实现 core/database.py（SQLite 管理）
4. 实现 main.py（程序入口）

### 阶段 2：核心模块

1. 实现 core/genome_manager.py（基因组管理）
2. 实现 core/sequence_extractor.py（序列提取）
3. 实现 core/blast_aligner.py（BLAST 比对）
4. 实现 core/visualizer.py（可视化生成）
5. 添加 BLAST+ 检测

### 阶段 3：主界面框架

1. 实现 ui/main_window.py（主窗口）
2. 实现导航和页面切换
3. 实现状态栏

### 阶段 4：基因组管理

1. 实现 ui/widgets/genome_manager_dialog.py
2. 基因组列表展示（从 SQLite 读取）
3. 添加本地基因组
4. 搜索在线基因组
5. 下载功能（带进度条）
6. 删除功能

### 阶段 5：分析页面

1. 实现 Gene ID 模式页面
2. 实现 Location 模式页面
3. 实现 Sequence 模式页面
4. 实现后台分析线程
5. 实现结果展示

### 阶段 6：历史记录

1. 实现历史记录页面
2. 从 SQLite 读取历史
3. 重新打开结果

### 阶段 7：打包与测试

1. 编写 PyInstaller spec
2. 测试打包
3. 优化体积
4. 编写文档

---

## main.py 入口设计

```python
"""GeneScreen 1.0 主程序入口"""
import sys
from PySide6.QtWidgets import QApplication, QMessageBox
from PySide6.QtGui import QIcon

from utils.blast_check import check_blast
from ui.main_window import MainWindow


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("GeneScreen")
    app.setApplicationVersion("1.0")

    # 检查 BLAST+
    if not check_blast():
        QMessageBox.critical(
            None,
            "BLAST+ 未安装",
            "请先安装 BLAST+：\n\n"
            "下载地址：\n"
            "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n\n"
            "安装后重新启动 GeneScreen。"
        )
        sys.exit(1)

    # 加载样式
    # with open("ui/resources/styles.qss") as f:
    #     app.setStyleSheet(f.read())

    # 显示主窗口
    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
```

---

## 依赖清单

### environment.yml（开发环境）

```yaml
name: genescreen3
channels:
  - conda-forge
dependencies:
  - python=3.10
  - ripgrep>=14        # GFF 文件搜索（~8MB，比纯 Python 快）
  - cairo              # cairosvg 运行时依赖
  - pip:
    - PySide6>=6.5     # GUI 框架
    - pyfaidx          # FASTA 索引/提取
    - biopython        # BLAST 结果解析
    - cairosvg         # SVG 转 PNG
    - pyinstaller      # 打包工具
```

### 安装命令

```bash
micromamba create -f environment.yml
micromamba activate genescreen3
```

---

## PyInstaller 打包配置

### build/GeneScreen.spec

```python
# -*- mode: python ; coding: utf-8 -*-

a = Analysis(
    ['../main.py'],
    pathex=[],
    binaries=[],
    datas=[
        ('../core', 'core'),
        ('../ui', 'ui'),
        ('../utils', 'utils'),
    ],
    hiddenimports=[
        'PySide6.QtCore',
        'PySide6.QtWidgets',
        'PySide6.QtGui',
        'pyfaidx',
        'Bio',
        'Bio.Blast',
        'Bio.Blast.NCBIXML',
        'cairosvg',
        'pandas',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'tkinter',
        'matplotlib',
        'numpy.testing',
    ],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='GeneScreen',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,  # 无控制台窗口
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='icon.ico',
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='GeneScreen',
)
```

### 打包命令

```bash
cd build
pyinstaller GeneScreen.spec
```

打包后文件位于 `build/dist/GeneScreen/`

---

## 用户使用流程

### 前置条件

1. 安装 BLAST+：https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
2. 下载 GeneScreen_1.0_Windows.zip 并解压

### 使用步骤

1. 双击 `GeneScreen.exe`
2. 等待浏览器自动打开（约 3 秒）
3. 在 Web 界面中操作：
   - 首次使用：在「基因组管理」页面添加/下载基因组
   - 选择分析模式（Gene ID / Location / Sequence）
   - 填写参数，点击「开始分析」
   - 查看结果，下载文件
4. 关闭命令行窗口即可退出

---

## 注意事项

1. **BLAST+ 必须手动安装**
   - EXE 中不包含 BLAST+（体积太大且需要 PATH 配置）
   - 启动时自动检测，未安装会弹窗提示

2. **ripgrep 打包进 EXE**
   - 用于快速搜索 GFF 文件（~8MB）
   - PyInstaller 打包时需要包含 rg.exe

3. **基因组存储位置**
   - 默认：`C:\Users\<用户名>\.genescreen\genomes\`
   - 可在设置中修改

4. **大文件处理**
   - 基因组文件使用本地路径，不上传
   - 结果文件提供 ZIP 下载

5. **打包体积优化**
   - 使用 UPX 压缩
   - 排除不必要的依赖
   - 预计最终体积：150-200MB

6. **错误处理**
   - BLAST 未安装：弹窗提示下载链接
   - 基因组不存在：引导到基因组管理页面
   - 比对失败：显示详细日志

---

## 预计工作量

| 阶段 | 内容 | 文件数 |
|------|------|--------|
| 阶段 1 | 项目初始化 | 4 |
| 阶段 2 | 核心模块 | 6 |
| 阶段 3 | 主界面框架 | 2 |
| 阶段 4 | 基因组管理 | 2 |
| 阶段 5 | 分析页面 | 4 |
| 阶段 6 | 历史记录 | 1 |
| 阶段 7 | 打包与测试 | 3 |

**总计约 22 个文件**

---

## 方案优势总结

| 对比项 | Streamlit 方案 | PySide6 + SQLite 方案 |
|--------|---------------|----------------------|
| 打包体积 | ~200MB | ~80-100MB |
| 启动速度 | 3-5秒（启动服务器） | <1秒 |
| 响应速度 | 有网络延迟 | 即时响应 |
| 数据管理 | JSON 文件 | SQLite 数据库 |
| 离线使用 | 支持 | 支持 |
| 用户体验 | 浏览器界面 | 原生桌面应用 |
| 开发复杂度 | 低 | 中 |

**PySide6 + SQLite 方案更适合作为专业的桌面工具分发。**
