# GeneScreen 1.0

[![Release](https://img.shields.io/github/v/release/Binggo913/GeneScreen)](https://github.com/Binggo913/GeneScreen/releases)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

基因组比对与变异分析工具 - 桌面版

![screenshot](docs/screenshot.png)

## 功能特性

- **Gene ID 模式**: 输入基因 ID，从参考基因组提取序列并与目标基因组比对
- **Location 模式**: 输入染色体坐标，提取指定区域并与目标基因组比对
- **Sequence 模式**: 输入序列，与参考基因组进行比对（支持多序列批量分析）
- **变异检测**: 自动检测 SNP 和 Indel 变异
- **可视化报告**: 生成交互式 HTML 分析报告，包含 LINKVIEW 可视化
- **基因组管理**: 支持本地基因组、IGV 公共基因组、Ensembl Plants 植物基因组
- **历史记录**: 自动保存分析历史，支持快速查看和重新分析

## 下载安装

### Windows

1. 从 [Releases](https://github.com/Binggo913/GeneScreen/releases) 下载对应架构的 `GeneScreen_1.0_windows_<arch>_<timestamp>.zip`
2. 解压后运行 `GeneScreen.exe`
3. 安装 BLAST+（见下方说明）

### macOS

1. 从 [Releases](https://github.com/Binggo913/GeneScreen/releases) 下载对应架构的 `GeneScreen_1.0_macos_<arch>_<timestamp>.dmg`
2. 打开 DMG 后将 GeneScreen 拖到 Applications
3. 安装 BLAST+（见下方说明）

## 安装 BLAST+

GeneScreen 依赖 NCBI BLAST+ 进行序列比对，需要单独安装：

**Windows:**
1. 下载 [BLAST+ Windows 安装包](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
2. 运行安装程序，安装时勾选"Add to PATH"
3. 重启电脑或重新打开终端

**macOS:**
```bash
brew install blast
```

**验证安装:**
```bash
blastn -version
```

## 使用方法

### Gene ID 模式

1. 选择参考基因组（需要包含 GFF 注释文件）
2. 选择查询基因组
3. 输入基因 ID（如 `LOC_Os06g10990`），支持实时搜索匹配
4. 设置 Identity 阈值（默认 90%）
5. 点击"开始分析"

### Location 模式

1. 选择参考基因组
2. 选择查询基因组
3. 输入位置，支持多种格式：
   - `Chr1:1000000-1050000`
   - `Chr1:1000000..1050000`
   - `Chr1 1000000 1050000`
4. 点击"开始分析"

### Sequence 模式

1. 选择参考基因组
2. 输入或粘贴序列（FASTA 格式，支持多序列批量分析）
3. 点击"开始分析"

## 基因组管理

点击侧边栏的"基因组管理"按钮：

- **添加本地基因组**: 添加本地 FASTA 和 GFF 文件，自动建立索引
- **搜索在线基因组**: 搜索 IGV 和 Ensembl Plants 数据库
- **下载基因组**: 自动下载、解压并索引基因组文件

## 输出文件

分析完成后，输出目录包含：

| 文件 | 说明 |
|------|------|
| `*.fasta` | 提取的序列文件 |
| `*.coords` | 比对坐标文件 |
| `*.snps` | SNP/Indel 变异文件 |
| `*.blast.xml` | BLAST 原始结果 |
| `*.report.html` | HTML 分析报告（含 LINKVIEW 可视化） |

## 开发

### 环境配置

```bash
# 使用 micromamba/conda
micromamba create -f environment.yml
micromamba activate genescreen

# 或使用 pip
pip install -r requirements.txt
```

### 运行

```bash
python main.py
```

### 打包

```bash
# 本地仅执行 Nuitka 构建；正式四架构产物由 GitHub Actions 生成
python build.py
```

## 系统要求

- Windows 10/11 (x64) 或 macOS 10.15+
- BLAST+ 2.12+
- 4GB+ 内存（推荐 8GB+）

## 许可证

MIT License

## 作者

xbzhang
