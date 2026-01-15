# GeneScreen

[![Release](https://img.shields.io/github/v/release/Binggo913/GeneScreen)](https://github.com/Binggo913/GeneScreen/releases)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Build](https://github.com/Binggo913/GeneScreen/actions/workflows/build.yml/badge.svg)](https://github.com/Binggo913/GeneScreen/actions)

基因组序列比对与变异分析工具，支持 Gene ID / Location / Sequence 三种模式，自动检测 SNP/Indel 变异，生成交互式可视化报告。

## 功能特性

- **多模式比对**：支持 Gene ID、染色体坐标、FASTA 序列三种输入方式
- **变异检测**：自动检测 SNP 和 Indel 变异
- **可视化报告**：生成交互式 HTML 分析报告，包含 LINKVIEW 共线性可视化
- **基因组管理**：支持本地基因组、IGV 公共基因组、Ensembl Plants 植物基因组
- **跨平台**：Windows / macOS / Linux 全平台支持

## 版本选择

| 版本 | 适用场景 | 平台支持 |
|------|----------|----------|
| [GeneScreen_GUI](GeneScreen_GUI/) | 桌面用户，图形界面操作 | Windows, macOS |
| [GeneScreen_CLI](GeneScreen_CLI/) | 命令行用户，批量分析，服务器部署 | Windows, macOS, Linux |

## 快速开始

### GUI 版本（推荐新手使用）

1. 从 [Releases](https://github.com/Binggo913/GeneScreen/releases) 下载安装包
   - Windows: `GeneScreen_Setup_x.x.x.exe`
   - macOS (Apple Silicon): `GeneScreen_macOS_arm64.dmg`
   - macOS (Intel): `GeneScreen_macOS_x64.dmg`
2. 安装 [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
3. 运行 GeneScreen，开始分析

### CLI 版本

```bash

# 1. 克隆仓库
git clone https://github.com/Binggo913/GeneScreen.git
cd GeneScreen/GeneScreen_CLI

# 2. 创建环境
# micromamba（推荐）
micromamba create -f environment.yml
micromamba activate genescreen

# conda（传统）
conda env create -f environment.yml
conda activate genescreen

# 3. 下载基因组
python bin/GenomeManager.py download IRGSP-1.0

# 4. 运行分析
python bin/GeneScreen.py -ref IRGSP-1.0 -qry ZS97 -gid LOC_Os06g10990 -o output/
```

## 项目结构

```
GeneScreen/
├── GeneScreen_GUI/          # 图形界面版本
│   ├── main.py              # GUI 入口
│   ├── core/                # 核心逻辑
│   ├── ui/                  # 界面组件
│   └── build.py             # 打包脚本
│
├── GeneScreen_CLI/          # 命令行版本
│   └── bin/
│       ├── GeneScreen.py          # 主程序
│       ├── GeneScreenVisualizer.py # 可视化模块
│       ├── GenomeManager.py       # 基因组管理
│       └── LINKVIEW.py            # 共线性可视化
│
└── .github/workflows/       # CI/CD 配置
    └── build.yml            # 跨平台自动构建
```

## 依赖

### 环境管理器

| 工具 | 说明 |
|------|------|
| [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) | 轻量级环境管理（推荐） |
| [conda](https://docs.conda.io/en/latest/miniconda.html) | 传统环境管理 |

### 核心依赖

| 工具 | 用途 |
|------|------|
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) | 序列比对 |
| [pyfaidx](https://github.com/mdshw5/pyfaidx) | FASTA 序列提取 |
| [biopython](https://biopython.org/) | BLAST 结果解析 |
| [ripgrep](https://github.com/BurntSushi/ripgrep) | 高速文本搜索 |
| [cairosvg](https://cairosvg.org/) | SVG 转 PNG |

## 输出文件

| 文件 | 说明 |
|------|------|
| `*.fasta` | 提取的序列文件 |
| `*.coords` | 比对坐标文件 |
| `*.snps` | SNP/Indel 变异文件 |
| `*.blast.xml` | BLAST 原始结果 |
| `*.svg` / `*.png` | 共线性可视化图 |
| `report.html` | 交互式 HTML 分析报告 |

## 文档

- [GUI 版本详细文档](GeneScreen_GUI/README.md)
- [CLI 版本详细文档](GeneScreen_CLI/README.md)

## 系统要求

- Windows 10/11 (x64) 或 macOS 10.15+ 或 Linux
- BLAST+ 2.12+
- Python 3.11+（CLI 版本）
- 4GB+ 内存（推荐 8GB+）

## 许可证

MIT License

## 作者

xbzhang
