# Sanger 16S rRNA 测序分析工具包

<div align="center">

![Python Version](https://img.shields.io/badge/python-3.7%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey)
![Status](https://img.shields.io/badge/status-stable-brightgreen)

**一个智能、易用的 Sanger 测序数据分析工具，专为 16S rRNA 基因测序设计**

[功能特性](#-主要特性) • [快速开始](#-快速开始) • [使用文档](#-详细使用指南) • [常见问题](#-常见问题) • [贡献指南](#-贡献与反馈)

</div>

---

## ✨ 主要特性

- 🎯 **智能文件名识别** - 自动识别多种命名规则，无需手动重命名
- 🔄 **完整分析流程** - 从原始测序文件到 BLAST 物种鉴定一键完成
- 💾 **断点续传** - BLAST 分析支持中断后继续，节省时间
- 📊 **专业报告** - 自动生成 Excel 和文本格式的分析报告
- ⚙️ **灵活配置** - 通过配置文件自定义所有参数
- 🖥️ **友好界面** - 支持交互式和批处理两种模式

## 📋 功能模块

### 1. 智能文件名解析器 (`filename_parser.py`)

自动识别以下命名格式：

| 格式类型 | 示例 | 说明 |
|---------|------|------|
| 标准格式 | `0001_31525103103397_(3sm1)_[16SF].ab1` | 编号_订单号_(样本ID)_[方向] |
| 下划线分隔 | `sample1_F.ab1` | 样本ID_方向 |
| 横杠分隔 | `sample1-F.ab1` | 样本ID-方向 |
| 前缀格式 | `F_sample1.ab1` | 方向_样本ID |
| 后缀格式 | `sample1.F.ab1` | 样本ID.方向 |
| 16S标记 | `sample1_16SF.ab1` | 样本ID_16S方向 |

**不支持的格式？** 只需在 `config.yaml` 中添加自定义规则！

### 2. 序列组装模块 (`graph_based_assembly.py`)

- 基于图论的序列组装算法
- 自动配对正反向测序文件
- 质量控制和末端修剪
- 生成高质量的一致性序列

### 3. BLAST 分析模块 (`ncbi_blast_analysis.py`)

- 在线 NCBI BLAST 物种鉴定
- 支持多种数据库（nt, RefSeq, 16S rRNA）
- 断点续传功能
- 详细的比对结果和统计信息

### 4. 完整流程管理器 (`sanger_16s_pipeline.py`)

- 一键完成全流程分析
- 交互式向导模式
- 批处理自动化模式

## 🚀 快速开始

### 系统要求

- Python 3.7 或更高版本
- 操作系统：Windows / Linux / macOS

### 安装

#### 方法一：从 GitHub 克隆（推荐）

```bash
# 克隆仓库
git clone https://github.com/YOUR_USERNAME/sanger-16s-analysis.git
cd sanger-16s-analysis

# 安装依赖
pip install -r requirements.txt
```

#### 方法二：直接下载

1. 下载最新版本的 [ZIP 文件](https://github.com/YOUR_USERNAME/sanger-16s-analysis/archive/refs/heads/main.zip)
2. 解压到本地目录
3. 安装依赖：

```bash
pip install biopython pandas openpyxl pyyaml
```

### 基本使用

#### 方式一：交互式模式（推荐新手）

```bash
python sanger_16s_pipeline.py
```

程序会引导您完成以下步骤：
1. 选择输入目录
2. 自动扫描和识别文件
3. 确认样本信息
4. 序列组装
5. BLAST 分析（可选）

#### 方式二：批处理模式

```bash
# 完整分析（组装 + BLAST）
python sanger_16s_pipeline.py -i ./data --batch

# 只进行序列组装
python sanger_16s_pipeline.py -i ./data --batch --skip-blast
```

#### 方式三：单独运行各模块

```bash
# 1. 扫描文件
python filename_parser.py ./data

# 2. 序列组装
python graph_based_assembly.py

# 3. BLAST 分析
python ncbi_blast_analysis.py -i assembled_sequences.fasta
```

## 📖 详细使用指南

### 查看支持的文件命名规则

```bash
python sanger_16s_pipeline.py --show-patterns
```

### 只扫描文件（不进行分析）

```bash
python sanger_16s_pipeline.py -i ./data --scan-only
```

### 自定义配置

编辑 `config.yaml` 文件：

```yaml
# 添加自定义文件名规则
filename_patterns:
  - name: "我的自定义格式"
    pattern: '^Sample_([^_]+)_([FR])\.ab1$'
    sample_group: 1
    direction_group: 2
    description: "格式: Sample_样本ID_方向"

# 调整组装参数
assembly:
  min_overlap: 20        # 最小重叠长度
  min_quality: 20        # 最小质量分数
  trim_ends: true        # 是否修剪低质量末端

# 调整 BLAST 参数
blast:
  database: "nt"         # 数据库: nt, refseq_rna, 16S_ribosomal_RNA
  max_hits: 10           # 每个序列返回的最大匹配数
  delay: 10              # 每次查询间隔（秒）
```

### BLAST 分析选项

```bash
# 使用 RefSeq 数据库（更快）
python ncbi_blast_analysis.py -d refseq_rna

# 调整延迟时间
python ncbi_blast_analysis.py --delay 15

# 从头开始（不使用断点续传）
python ncbi_blast_analysis.py --no-resume

# 查看帮助
python ncbi_blast_analysis.py -h
```

## 📂 输出文件说明

### 序列组装结果

位于 `sanger_assembly_results/` 目录：

- `assembled_sequences.fasta` - 组装后的序列（FASTA 格式）
- `assembly_report.xlsx` - 详细的组装报告（Excel 格式）
- `assembly_report.txt` - 文本格式报告
- `quality_plots/` - 质量分数可视化图表

### BLAST 分析结果

位于 `blast_results/` 目录：

- `blast_results_YYYYMMDD_HHMMSS.xlsx` - BLAST 结果（Excel 格式）
- `blast_report_YYYYMMDD_HHMMSS.txt` - 文本格式报告
- `blast_progress.json` - 进度保存文件（用于断点续传）

### 日志文件

- `analysis.log` - 完整的分析日志
- `blast_analysis.log` - BLAST 分析日志

## 🔧 常见问题

### Q1: 我的文件命名格式不被识别怎么办？

**A:** 有两种解决方案：

1. **添加自定义规则**（推荐）：
   - 编辑 `config.yaml` 文件
   - 在 `filename_patterns` 部分添加您的命名规则
   - 使用正则表达式匹配文件名

2. **重命名文件**：
   - 使用工具批量重命名为支持的格式
   - 例如：`sample1_F.ab1`, `sample1_R.ab1`

### Q2: BLAST 分析太慢怎么办？

**A:** 几个建议：

1. 使用更快的数据库：`-d refseq_rna` 或 `-d 16S_ribosomal_RNA`
2. 减少返回的匹配数：`-n 5`
3. 利用断点续传：中断后重新运行会自动继续
4. 考虑使用本地 BLAST（需要下载数据库）

### Q3: 如何处理单端测序数据？

**A:** 工具会自动识别：
- 配对的文件会进行组装
- 单个文件会直接使用（不组装）
- 在报告中会明确标注

### Q4: 组装质量不理想怎么办？

**A:** 调整 `config.yaml` 中的参数：

```yaml
assembly:
  min_overlap: 30        # 增加最小重叠长度
  min_quality: 25        # 提高质量阈值
  trim_threshold: 0.1    # 调整修剪阈值
```

### Q5: BLAST 分析中断了怎么办？

**A:** 不用担心！
- 进度会自动保存到 `blast_progress.json`
- 重新运行相同命令即可继续
- 使用 `--no-resume` 参数可以从头开始

## 🎯 使用场景示例

### 场景 1：测序公司返回的标准格式数据

```bash
# 文件格式: 0001_订单号_(样本ID)_[16SF].ab1
python sanger_16s_pipeline.py -i ./sequencing_results --batch
```

### 场景 2：实验室自己测序的简单命名

```bash
# 文件格式: sample1_F.ab1, sample1_R.ab1
python sanger_16s_pipeline.py -i ./lab_data --batch
```

### 场景 3：只想快速查看文件是否能被识别

```bash
python sanger_16s_pipeline.py -i ./data --scan-only
```

### 场景 4：已有组装序列，只做 BLAST

```bash
python ncbi_blast_analysis.py -i my_sequences.fasta -d refseq_rna
```

### 场景 5：大批量样本，分批处理

```bash
# 第一批
python sanger_16s_pipeline.py -i ./batch1 --batch

# 第二批
python sanger_16s_pipeline.py -i ./batch2 --batch

# 合并结果（手动或使用脚本）
```

## 📊 输出示例

### Excel 报告示例

| 样本ID | 物种名称 | 相似度 | 覆盖度 | E-value | 序列长度 |
|--------|---------|--------|--------|---------|---------|
| 3sm1 | Escherichia coli | 99.5% | 100% | 0.0 | 1450 |
| 3sm2 | Bacillus subtilis | 98.2% | 99% | 0.0 | 1420 |

### 文本报告示例

```
样本: 3sm1
最佳匹配: Escherichia coli strain K-12
相似度: 99.5%
覆盖度: 100%
E-value: 0.0
序列长度: 1450 bp
```

## 🔬 技术细节

### 序列组装算法

1. **质量控制**：过滤低质量碱基
2. **末端修剪**：移除低质量区域
3. **重叠检测**：基于图论的最优重叠搜索
4. **一致性序列**：生成高质量共有序列

### BLAST 参数

- **程序**：blastn（核苷酸比对）
- **默认数据库**：nt（所有核苷酸）
- **E-value 阈值**：0.001
- **最大目标序列数**：10

## 🤝 贡献与反馈

我们欢迎任何形式的贡献！如果您：
- 🐛 发现了 bug
- 💡 有功能建议
- 📝 需要支持新的文件格式
- 🔧 想要贡献代码

请查看 [贡献指南](CONTRIBUTING.md) 并提交 Issue 或 Pull Request！

### 联系方式

- 📧 提交 [Issue](https://github.com/YOUR_USERNAME/sanger-16s-analysis/issues)
- 💬 参与 [Discussions](https://github.com/YOUR_USERNAME/sanger-16s-analysis/discussions)

## 📝 更新日志

查看完整的 [CHANGELOG.md](CHANGELOG.md)

### v1.0.0 (2025-11-01)
- ✨ 智能文件名识别系统
- ✨ 完整的分析流程
- ✨ BLAST 断点续传
- ✨ 交互式和批处理模式
- ✨ 灵活的配置系统

## 📄 许可证

本项目采用 [MIT License](LICENSE) 开源许可证。

这意味着您可以自由地：
- ✅ 使用本软件进行商业或非商业用途
- ✅ 修改源代码
- ✅ 分发本软件
- ✅ 私人使用

唯一的要求是保留版权声明和许可证声明。

## 👨‍🔬 致谢

本工具使用以下优秀的开源库：
- [Biopython](https://biopython.org/) - 生物信息学工具
- [Pandas](https://pandas.pydata.org/) - 数据处理
- [PyYAML](https://pyyaml.org/) - 配置文件解析

---

<div align="center">

**祝您科研顺利！** 🎉

如有问题，请查看文档或提交 [Issue](https://github.com/YOUR_USERNAME/sanger-16s-analysis/issues)

⭐ 如果这个项目对您有帮助，欢迎给个 Star！

Made with ❤️ for the scientific community

</div>

