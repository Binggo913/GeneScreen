# GeneScreen 3.0 开发任务清单

## 阶段 1：项目初始化

- [x] 1.1 创建目录结构
- [x] 1.2 创建 requirements.txt / environment.yml
- [x] 1.3 实现 core/database.py（SQLite 数据库管理）
- [x] 1.4 实现 main.py（程序入口，待 UI 完成后可运行）

## 阶段 2：核心模块（从 2.0 整合）

> 核心业务逻辑从 2.0 复用，主要工作是适配 SQLite 和整理代码结构

- [x] 2.1 实现 utils/blast_check.py（BLAST+ 检测）
- [x] 2.2 整合 core/genome_manager.py（从 2.0 复用，改用 SQLite）
- [x] 2.3 整合 core/analysis.py（整合 2.0 的 SequenceExtractor + BlastAligner + 变异检测）
- [x] 2.4 整合 core/visualizer.py（从 2.0 复用 GeneScreenVisualizer + LINKVIEW）

## 阶段 3：主界面框架

- [x] 3.1 实现 ui/main_window.py（主窗口）
- [x] 3.2 实现 ui/resources/styles.qss（样式表）

## 阶段 4：基因组管理

- [x] 4.1 实现 ui/widgets/genome_selector.py（基因组选择器）
- [x] 4.2 实现 ui/widgets/genome_manager_dialog.py（基因组管理对话框）

## 阶段 5：分析页面

- [x] 5.1 实现 ui/pages/gene_id_page.py（Gene ID 模式）
- [x] 5.2 实现 ui/pages/location_page.py（Location 模式）
- [x] 5.3 实现 ui/pages/sequence_page.py（Sequence 模式）
- [x] 5.4 实现 ui/widgets/result_viewer.py（结果查看器）

## 阶段 6：历史记录

- [x] 6.1 实现 ui/pages/history_page.py（历史记录页面）

## 阶段 7：打包与测试

- [x] 7.1 创建 build/GeneScreen.spec（PyInstaller 配置，包含 ripgrep）
- [x] 7.2 测试打包 EXE（已由用户手动完成）
- [x] 7.3 编写 README.md（用户文档）

---

## 进度统计

- 总任务数：18
- 已完成：18
- 进度：100%

## 备注

- 核心模块从 GeneScreen_2.0/bin/ 复用，不需要重写
- 依赖：保留 ripgrep（~8MB），搜索 GFF 比纯 Python 快
- 打包时需要把 rg.exe 一起打包
- Task 7.2 (测试打包 EXE) 需要手动执行验证
