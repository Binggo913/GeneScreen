# GeneScreen 3.0

基因组比对与变异分析工具 - 跨平台桌面版

## 功能特性

- **Gene ID 模式**: 输入基因 ID，从参考基因组提取序列并与目标基因组比对
- **Location 模式**: 输入染色体坐标，提取指定区域并与目标基因组比对
- **Sequence 模式**: 输入序列，与参考基因组进行比对
- **变异检测**: 自动检测 SNP 和 Indel 变异
- **可视化报告**: 生成 HTML 格式的分析报告
- **基因组管理**: 支持本地基因组、IGV 公共基因组、Ensembl Plants 植物基因组

## 系统要求

- Windows 10/11 (64-bit) 或 Linux
- BLAST+ 2.12+ (需要单独安装)

## 安装

### 1. 安装 BLAST+

从 NCBI 下载并安装 BLAST+:
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

安装后确保 `blastn` 和 `makeblastdb` 在系统 PATH 中。

验证安装:
```bash
blastn -version
```

### 2. 运行 GeneScreen

- **Windows**: 双击 `GeneScreen.exe`
- **Linux**: 运行 `./GeneScreen`

## 使用方法

### Gene ID 模式

1. 选择参考基因组（需要包含 GFF 注释文件）
2. 选择查询基因组
3. 输入基因 ID（如 `LOC_Os06g10990`）
4. 设置 Identity 阈值
5. 选择输出目录
6. 点击"开始分析"

### Location 模式

1. 选择参考基因组
2. 选择查询基因组
3. 输入位置（如 `Chr1:1000000-1050000`）
4. 设置 Identity 阈值
5. 选择输出目录
6. 点击"开始分析"

### Sequence 模式

1. 选择参考基因组
2. 输入或粘贴序列（FASTA 格式或纯序列）
3. 设置 Identity 阈值
4. 选择输出目录
5. 点击"开始分析"

## 基因组管理

点击侧边栏的"基因组管理"按钮，可以:

- **添加自定义基因组**: 添加本地 FASTA 和 GFF 文件
- **搜索在线基因组**: 搜索 IGV 和 Ensembl Plants 数据库
- **下载基因组**: 自动下载并索引基因组文件

## 输出文件

分析完成后，输出目录包含:

- `*.fasta` - 提取的序列文件
- `*.coords` - 比对坐标文件
- `*.snps` - SNP/Indel 变异文件
- `*.blast.xml` - BLAST 原始结果
- `*.report.html` - HTML 分析报告

## 开发

### 环境配置

```bash
# Linux/WSL
micromamba create -f environment.yml
micromamba activate genescreen

# Windows
micromamba create -f environment_windows.yml
micromamba activate genescreen

# 或使用 pip
pip install -r requirements.txt
```

### 运行

```bash
python main.py
```

### 打包

使用 Nuitka 编译成独立可执行文件：

```bash
python build.py
```

脚本会自动识别当前平台（Windows/Linux），生成对应的可执行文件到 `dist/` 目录。

## 依赖

- PySide6 - GUI 框架
- pyfaidx - FASTA 索引和序列提取
- biopython - BLAST 结果解析
- ripgrep - 快速 GFF 搜索（已打包）
- nuitka - 打包工具

## 许可证

MIT License

## 作者

xbzhang
