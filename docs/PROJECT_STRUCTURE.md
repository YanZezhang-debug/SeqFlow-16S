# 项目结构说明

本文档详细说明了 Sanger 16S rRNA 测序分析工具包的项目结构。

## 📁 目录结构

```
sanger-16s-analysis/
│
├── 📄 README.md                    # 项目主文档
├── 📄 LICENSE                      # MIT 开源许可证
├── 📄 CHANGELOG.md                 # 版本更新日志
├── 📄 CONTRIBUTING.md              # 贡献指南
├── 📄 requirements.txt             # Python 依赖包列表
├── 📄 config.yaml                  # 默认配置文件
├── 📄 .gitignore                   # Git 忽略文件列表
│
├── 📂 docs/                        # 文档目录
│   ├── PROJECT_STRUCTURE.md        # 项目结构说明（本文件）
│   ├── API.md                      # API 文档（待添加）
│   └── TROUBLESHOOTING.md          # 故障排除指南（待添加）
│
├── 📂 examples/                    # 示例和教程
│   ├── README.md                   # 示例说明
│   ├── custom_config.yaml          # 自定义配置示例
│   └── sample_data/                # 示例数据（可选）
│
├── 📂 src/ 或 核心 Python 脚本      # 源代码
│   ├── 📄 sanger_16s_pipeline.py   # 主流程管理器
│   ├── 📄 filename_parser.py       # 文件名解析模块
│   ├── 📄 graph_based_assembly.py  # 序列组装模块
│   └── 📄 ncbi_blast_analysis.py   # BLAST 分析模块
│
├── 📂 tests/                       # 单元测试（待添加）
│   ├── test_filename_parser.py
│   ├── test_assembly.py
│   └── test_blast.py
│
├── 📂 sanger_assembly_results/     # 序列组装输出（运行后生成）
│   ├── assembled_sequences.fasta
│   ├── assembly_report.xlsx
│   ├── assembly_report.txt
│   ├── assembly_summary.csv
│   └── quality_report.txt
│
├── 📂 blast_results/               # BLAST 分析输出（运行后生成）
│   ├── blast_results_*.xlsx
│   ├── blast_report_*.txt
│   └── blast_progress.json
│
└── 📂 原始数据文件/                 # 用户的测序数据
    ├── *.ab1                       # ABI 测序文件
    ├── *.fasta                     # FASTA 序列文件
    ├── *.seq                       # 序列文件
    └── *.pdf                       # 测序质量报告
```

## 📄 核心文件说明

### 主要脚本

#### 1. `sanger_16s_pipeline.py`
**完整流程管理器**

- **功能**：协调整个分析流程
- **模式**：
  - 交互式模式：引导用户逐步完成分析
  - 批处理模式：自动化执行全流程
- **主要功能**：
  - 文件扫描和识别
  - 调用序列组装
  - 调用 BLAST 分析
  - 生成综合报告
- **命令行参数**：
  ```bash
  python sanger_16s_pipeline.py [选项]
  
  -i, --input DIR          输入目录
  -o, --output DIR         输出目录
  -c, --config FILE        配置文件
  --batch                  批处理模式
  --skip-blast             跳过 BLAST 分析
  --scan-only              仅扫描文件
  --show-patterns          显示支持的命名格式
  ```

#### 2. `filename_parser.py`
**智能文件名解析器**

- **功能**：识别和解析测序文件名
- **支持格式**：
  - 标准格式：`0001_订单号_(样本ID)_[16SF].ab1`
  - 简单格式：`sample1_F.ab1`
  - 横杠格式：`sample1-F.ab1`
  - 前缀格式：`F_sample1.ab1`
  - 后缀格式：`sample1.F.ab1`
  - 16S标记：`sample1_16SF.ab1`
- **扩展性**：通过配置文件添加自定义规则
- **主要类**：
  - `FileNameParser`: 文件名解析主类
  - `SampleInfo`: 样本信息数据类

#### 3. `graph_based_assembly.py`
**序列组装引擎**

- **功能**：组装正反向测序序列
- **算法**：基于图论的序列组装
- **主要步骤**：
  1. 质量控制和过滤
  2. 末端修剪
  3. 重叠区域检测
  4. 一致性序列生成
- **输出**：
  - FASTA 格式的组装序列
  - 详细的质量报告
  - Excel 格式的统计表
- **主要类**：
  - `GraphBasedAssembler`: 组装器主类
  - `SequenceNode`: 序列节点
  - `OverlapGraph`: 重叠图

#### 4. `ncbi_blast_analysis.py`
**BLAST 物种鉴定模块**

- **功能**：在线 NCBI BLAST 查询
- **特性**：
  - 支持多种数据库
  - 断点续传功能
  - 自动重试机制
  - 进度保存和恢复
- **输出**：
  - Excel 格式的 BLAST 结果
  - 文本格式的详细报告
  - JSON 格式的进度文件
- **主要类**：
  - `NCBIBlastAnalyzer`: BLAST 分析器
  - `BlastResult`: 结果数据类

### 配置文件

#### `config.yaml`
**主配置文件**

包含所有可配置参数：
- 文件名识别规则
- 序列组装参数
- BLAST 分析参数
- 输出设置
- 日志设置

结构：
```yaml
filename_patterns:
  - name: "格式名称"
    pattern: "正则表达式"
    sample_group: 1
    direction_group: 2

assembly:
  min_overlap: 20
  min_quality: 20
  trim_ends: true

blast:
  database: "nt"
  max_hits: 10
  delay: 10
```

### 文档文件

#### `README.md`
- 项目概述
- 功能特性
- 安装说明
- 使用指南
- 常见问题

#### `CHANGELOG.md`
- 版本历史
- 功能更新
- Bug 修复记录

#### `CONTRIBUTING.md`
- 贡献指南
- 代码规范
- 提交流程
- 开发环境设置

#### `LICENSE`
- MIT 开源许可证

## 📂 输出目录说明

### `sanger_assembly_results/`
序列组装结果目录

**文件列表：**

1. **`assembled_sequences.fasta`**
   - 格式：FASTA
   - 内容：所有组装后的序列
   - 用途：用于后续 BLAST 分析

2. **`assembly_report.xlsx`**
   - 格式：Excel
   - 内容：详细的组装报告
   - 包含：样本ID、序列长度、质量分数、组装状态等

3. **`assembly_report.txt`**
   - 格式：纯文本
   - 内容：与 Excel 相同，便于查看

4. **`assembly_summary.csv`**
   - 格式：CSV
   - 内容：组装统计摘要
   - 用途：便于导入其他分析工具

5. **`quality_report.txt`**
   - 格式：纯文本
   - 内容：质量控制详细信息

### `blast_results/`
BLAST 分析结果目录

**文件列表：**

1. **`blast_results_YYYYMMDD_HHMMSS.xlsx`**
   - 格式：Excel
   - 内容：BLAST 比对结果
   - 包含：物种名称、相似度、覆盖度、E-value 等

2. **`blast_report_YYYYMMDD_HHMMSS.txt`**
   - 格式：纯文本
   - 内容：详细的 BLAST 报告

3. **`blast_progress.json`**
   - 格式：JSON
   - 内容：分析进度
   - 用途：断点续传

## 🔧 数据流程

```
原始测序文件 (.ab1, .fasta, .seq)
    ↓
[filename_parser.py]
识别文件名，提取样本信息
    ↓
[graph_based_assembly.py]
质量控制 → 序列组装 → 生成一致性序列
    ↓
assembled_sequences.fasta
    ↓
[ncbi_blast_analysis.py]
在线 BLAST 查询 → 物种鉴定
    ↓
BLAST 结果报告 (.xlsx, .txt)
```

## 📊 文件格式说明

### 输入文件格式

#### ABI 测序文件 (.ab1)
- 二进制格式
- 包含原始测序数据
- 包含质量分数
- 使用 Biopython 解析

#### FASTA 文件 (.fasta, .fa)
- 纯文本格式
- 序列描述行：以 `>` 开头
- 序列数据：多行核苷酸序列

#### SEQ 文件 (.seq)
- 纯文本格式
- 通常只包含序列数据
- 可能不包含质量信息

### 输出文件格式

#### FASTA 输出
```
>sample1
ATCGATCGATCG...
>sample2
GCTAGCTAGCTA...
```

#### Excel 输出
- 使用 openpyxl 库生成
- 包含多个工作表
- 格式化的表格数据

#### CSV 输出
- 逗号分隔值
- UTF-8 编码
- 便于导入其他工具

## 🔍 代码组织

### 模块化设计

每个脚本都可以独立运行：

```python
# 独立使用文件名解析器
from filename_parser import FileNameParser
parser = FileNameParser()
result = parser.parse("sample1_F.ab1")

# 独立使用组装器
from graph_based_assembly import GraphBasedAssembler
assembler = GraphBasedAssembler()
assembled = assembler.assemble(forward_seq, reverse_seq)

# 独立使用 BLAST 分析
from ncbi_blast_analysis import NCBIBlastAnalyzer
analyzer = NCBIBlastAnalyzer()
results = analyzer.blast("sequences.fasta")
```

### 配置管理

使用 YAML 配置文件：
- 集中管理所有参数
- 易于修改和扩展
- 支持多配置文件

### 日志系统

统一的日志记录：
- 文件日志：保存到 `.log` 文件
- 控制台日志：实时显示进度
- 多级别：DEBUG, INFO, WARNING, ERROR

## 🚀 扩展开发

### 添加新功能

1. **添加新的文件格式支持**
   - 修改 `filename_parser.py`
   - 在 `config.yaml` 中添加规则

2. **添加新的组装算法**
   - 在 `graph_based_assembly.py` 中添加新类
   - 保持接口一致

3. **支持本地 BLAST**
   - 修改 `ncbi_blast_analysis.py`
   - 添加本地数据库路径配置

4. **添加 GUI 界面**
   - 使用 tkinter 或 PyQt
   - 调用现有的核心模块

### 测试开发

创建单元测试：
```python
# tests/test_filename_parser.py
import pytest
from filename_parser import FileNameParser

def test_parse_standard_format():
    parser = FileNameParser()
    result = parser.parse("0001_12345_(sample1)_[16SF].ab1")
    assert result['sample_id'] == 'sample1'
    assert result['direction'] == 'F'
```

## 📝 维护建议

### 版本控制

- 使用 Git 管理代码
- 遵循语义化版本号
- 及时更新 CHANGELOG.md

### 代码质量

- 遵循 PEP 8 规范
- 添加类型注解
- 编写文档字符串
- 定期代码审查

### 文档更新

- 功能变更时更新 README
- 添加新功能时更新示例
- 保持文档与代码同步

## 🔗 相关资源

- [Biopython 文档](https://biopython.org/wiki/Documentation)
- [NCBI BLAST 文档](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [Python 正则表达式](https://docs.python.org/3/library/re.html)
- [YAML 语法](https://yaml.org/spec/1.2/spec.html)

---

如有问题或建议，请提交 [Issue](https://github.com/YOUR_USERNAME/sanger-16s-analysis/issues)。

