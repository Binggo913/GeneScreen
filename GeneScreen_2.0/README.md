# GeneScreen 2.0

基因组比对与变异分析工具，跨平台版本（Windows/macOS/Linux）。

使用 BLAST+ + pyfaidx + ripgrep，实现完全跨平台兼容。

## 目录结构

```
GeneScreen_2.0/
├── bin/                        # 脚本目录
│   ├── GeneScreen.py           # 主程序
│   ├── GeneScreenVisualizer.py # 可视化模块
│   ├── GenomeManager.py        # 基因组管理
│   ├── LINKVIEW.py             # LINKVIEW 可视化工具
│   └── interval.py             # LINKVIEW 依赖模块
├── environment.yml             # Linux/macOS 环境
├── environment_windows.yml     # Windows 环境
└── README.md
```

## 依赖说明

| 工具 | 用途 | Windows | Linux/macOS |
|------|------|---------|-------------|
| **BLAST+** | 序列比对 | ✅ 官方支持 | ✅ |
| **ripgrep** | 文本搜索 | ✅ 官方支持 | ✅ |
| **pyfaidx** | 序列提取 | ✅ pip | ✅ |
| **biopython** | BLAST 解析 | ✅ pip | ✅ |
| **cairosvg** | SVG 转 PNG | ✅ pip | ✅ |

## 环境安装

### Linux / macOS

#### 1. 安装 Micromamba（如未安装）

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

> 💡 Tip: 如果想用 `mamba` 命令代替 `micromamba`：
>
> ```bash
> alias mamba=micromamba
> ```

更多信息参考 [Micromamba 官方安装指南](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)

#### 2. 创建环境

```bash
micromamba create -f environment.yml
micromamba activate genescreen
```

---

### Windows

#### 步骤 1：安装 BLAST+

下载 BLAST+ Windows 安装包：
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-win64.exe

运行安装程序，安装时会自动添加到 PATH。

#### 步骤 2：创建 Python 环境

```powershell
# 安装 Micromamba（如未安装）
Invoke-Expression ((Invoke-WebRequest -Uri https://micro.mamba.pm/install.ps1).Content)

# 创建环境
micromamba create -f environment_windows.yml
micromamba activate genescreen
```


## 一、基因组库管理

使用 `bin/GenomeManager.py` 管理基因组，支持自定义添加和从 IGV 下载。

基因组缓存目录：`~/.genescreen/genomes/`

### 1.1 查看已有基因组

```bash
python bin/GenomeManager.py list
```

### 1.2 搜索 IGV 公共基因组

```bash
python bin/GenomeManager.py search rice
python bin/GenomeManager.py search human
```

### 1.3 下载 IGV 基因组

```bash
python bin/GenomeManager.py download IRGSP-1.0
python bin/GenomeManager.py download hg38
```

### 1.4 添加自定义基因组

```bash
# 只添加基因组序列
python bin/GenomeManager.py add Nippon /path/to/Nipponbare.fa

# 添加基因组 + 注释文件
python bin/GenomeManager.py add Nippon /path/to/Nipponbare.fa -a /path/to/Nipponbare.gff3
```

### 1.5 删除基因组

```bash
python bin/GenomeManager.py remove IRGSP-1.0
```

---

## 二、使用方法

### 模式 1：Gene ID

从参考基因组中提取指定基因，与目标基因组比对。

```bash
# 单个基因
python bin/GeneScreen.py -ref Nippon -qry ZS97 -gid LOC_Os06g10990 -o output/

# 多个基因
python bin/GeneScreen.py -ref Nippon -qry ZS97 -gidl gene_list.txt -o output/

# 提取基因 + 上下游区域（上游 1000bp，下游 500bp）
python bin/GeneScreen.py -ref Nippon -qry ZS97 -gid LOC_Os06g10990 -u 1000 -d 500 -o output/
```

> 💡 **上下游延伸说明**：
> - `-u` / `--upstream`：上游延伸长度（默认 0）
> - `-d` / `--downstream`：下游延伸长度（默认 0）
> - 正链基因：上游在 5' 端（左侧），下游在 3' 端（右侧）
> - 负链基因：上游在 3' 端（右侧），下游在 5' 端（左侧）
> - 可视化图中会用大括号标注基因区域位置

### 模式 2：Location

从参考基因组中提取指定区域，与目标基因组比对。

```bash
# 使用位置文件
python bin/GeneScreen.py -ref Nippon -qry ZS97 -loc positions.txt -o output/

# 直接指定区域
python bin/GeneScreen.py -ref Nippon -qry ZS97 -loc Chr1:1000000-1050000 -o output/
```

位置文件格式（Tab 分隔）：

```
Chr1    1000000    1050000
Chr2    2000000    2100000    my_region
```

### 模式 3：Sequence

用户提供 FASTA 序列，与参考基因组比对。

```bash
python bin/GeneScreen.py -ref Nippon -seq query.fasta -o output/
```

---

## 三、比对参数

BLAST 使用 identity 阈值过滤比对结果：

```bash
# 默认 90% identity
python bin/GeneScreen.py -ref Nippon -qry ZS97 -gid LOC_Os06g10990 -o output/

# 自定义 identity 阈值
python bin/GeneScreen.py -ref Nippon -qry ZS97 -gid LOC_Os06g10990 -o output/ --identity 95
```

---

## 四、输出文件

| 文件 | 说明 |
|------|------|
| `{id}.fasta` | 提取的序列 |
| `{id}.blast.xml` | BLAST 原始结果 |
| `{id}.coords` | 比对坐标 |
| `{id}.snps` | SNP/Indel 变异 |
| `{id}.linkview.gff3` | 相对坐标注释（Gene ID 模式） |
| `{id}.svg` / `{id}.png` | 可视化图 |
| `report.html` | HTML 报告 |

---

## 五、参数速查

| 参数 | 说明 | 必需 |
|------|------|------|
| `-ref` | 参考基因组 | ✓ |
| `-qry` | 目标基因组 | 模式1/2 |
| `-ra` | 注释文件（覆盖默认） | 可选 |
| `-gid` | 单个 Gene ID | 模式1 |
| `-gidl` | Gene ID 列表 | 模式1 |
| `-u` / `--upstream` | 上游延伸长度（bp） | 模式1 可选 |
| `-d` / `--downstream` | 下游延伸长度（bp） | 模式1 可选 |
| `-loc` | 位置文件或区域字符串 | 模式2 |
| `-seq` | 序列文件 | 模式3 |
| `--identity` | 最小 identity 阈值 | 可选 |
| `-o` | 输出目录 | ✓ |
