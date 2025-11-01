# 使用示例

本目录包含了使用 Sanger 16S rRNA 测序分析工具包的各种示例。

## 📁 目录结构

```
examples/
├── README.md                    # 本文件
├── quick_start.sh              # 快速开始脚本（Linux/macOS）
├── quick_start.bat             # 快速开始脚本（Windows）
├── custom_config.yaml          # 自定义配置示例
└── sample_data/                # 示例数据（如果包含）
    └── README.md               # 数据说明
```

## 🚀 快速开始示例

### 示例 1：基本使用（交互式模式）

最简单的使用方式，程序会引导您完成所有步骤：

```bash
# 进入项目根目录
cd ..

# 运行交互式流程
python sanger_16s_pipeline.py
```

程序会提示您：
1. 选择输入目录
2. 确认识别的文件
3. 选择是否进行 BLAST 分析
4. 查看分析结果

### 示例 2：批处理模式（完整分析）

适合自动化处理，包含序列组装和 BLAST 分析：

```bash
# Linux/macOS
python sanger_16s_pipeline.py -i ./your_data_directory --batch

# Windows
python sanger_16s_pipeline.py -i .\your_data_directory --batch
```

### 示例 3：只进行序列组装（跳过 BLAST）

如果只需要组装序列，不需要物种鉴定：

```bash
python sanger_16s_pipeline.py -i ./your_data_directory --batch --skip-blast
```

### 示例 4：查看支持的文件命名格式

在处理数据前，先检查您的文件命名是否被支持：

```bash
python sanger_16s_pipeline.py --show-patterns
```

### 示例 5：扫描文件（不进行分析）

只想看看哪些文件会被识别：

```bash
python sanger_16s_pipeline.py -i ./your_data_directory --scan-only
```

## 🔧 高级使用示例

### 示例 6：使用自定义配置

创建自定义配置文件并使用：

```bash
# 复制示例配置
cp examples/custom_config.yaml my_config.yaml

# 编辑配置文件
# 修改 my_config.yaml 中的参数

# 使用自定义配置运行
python sanger_16s_pipeline.py -i ./data -c my_config.yaml --batch
```

### 示例 7：单独运行 BLAST 分析

如果您已经有组装好的序列文件：

```bash
# 使用默认数据库（nt）
python ncbi_blast_analysis.py -i assembled_sequences.fasta

# 使用 RefSeq 数据库（更快）
python ncbi_blast_analysis.py -i assembled_sequences.fasta -d refseq_rna

# 使用 16S 专用数据库
python ncbi_blast_analysis.py -i assembled_sequences.fasta -d 16S_ribosomal_RNA

# 调整查询间隔（避免被限制）
python ncbi_blast_analysis.py -i assembled_sequences.fasta --delay 15

# 返回更多匹配结果
python ncbi_blast_analysis.py -i assembled_sequences.fasta -n 20
```

### 示例 8：BLAST 断点续传

如果 BLAST 分析中断了：

```bash
# 直接重新运行相同的命令，会自动继续
python ncbi_blast_analysis.py -i assembled_sequences.fasta

# 如果想从头开始（不使用断点续传）
python ncbi_blast_analysis.py -i assembled_sequences.fasta --no-resume
```

## 📊 实际使用场景

### 场景 1：测序公司返回的数据

```bash
# 文件格式: 0001_订单号_(样本ID)_[16SF].ab1
# 目录结构:
# sequencing_results/
#   ├── 0001_31525103103397_(3sm1)_[16SF].ab1
#   ├── 0002_31525103103397_(3sm1)_[16SR].ab1
#   ├── 0003_31525103103398_(3sm2)_[16SF].ab1
#   └── 0004_31525103103398_(3sm2)_[16SR].ab1

# 运行完整分析
python sanger_16s_pipeline.py -i ./sequencing_results --batch
```

### 场景 2：实验室简单命名

```bash
# 文件格式: sample1_F.ab1, sample1_R.ab1
# 目录结构:
# lab_data/
#   ├── sample1_F.ab1
#   ├── sample1_R.ab1
#   ├── sample2_F.ab1
#   └── sample2_R.ab1

# 运行完整分析
python sanger_16s_pipeline.py -i ./lab_data --batch
```

### 场景 3：混合格式文件

```bash
# 如果文件命名格式不统一，先扫描看看能识别多少
python sanger_16s_pipeline.py -i ./mixed_data --scan-only

# 查看支持的格式
python sanger_16s_pipeline.py --show-patterns

# 根据扫描结果，可能需要：
# 1. 重命名部分文件
# 2. 或在 config.yaml 中添加自定义规则
```

### 场景 4：大批量样本分批处理

```bash
# 第一批（样本 1-10）
python sanger_16s_pipeline.py -i ./batch1 --batch -o ./results/batch1

# 第二批（样本 11-20）
python sanger_16s_pipeline.py -i ./batch2 --batch -o ./results/batch2

# 第三批（样本 21-30）
python sanger_16s_pipeline.py -i ./batch3 --batch -o ./results/batch3

# 最后合并结果（手动或使用 Excel）
```

## 🎯 自定义配置示例

### 添加自定义文件名规则

编辑 `config.yaml` 或创建新的配置文件：

```yaml
filename_patterns:
  # 添加您的自定义格式
  - name: "实验室格式A"
    pattern: '^Lab_([A-Z0-9]+)_([FR])\.ab1$'
    sample_group: 1
    direction_group: 2
    description: "格式: Lab_样本编号_方向.ab1"
  
  - name: "实验室格式B"
    pattern: '^(\d{4})_([FR])_16S\.ab1$'
    sample_group: 1
    direction_group: 2
    description: "格式: 编号_方向_16S.ab1"
```

### 调整组装参数

```yaml
assembly:
  min_overlap: 30          # 增加最小重叠长度（更严格）
  min_quality: 25          # 提高质量阈值
  trim_ends: true          # 启用末端修剪
  trim_threshold: 0.1      # 修剪阈值
  window_size: 20          # 质量评估窗口大小
```

### 调整 BLAST 参数

```yaml
blast:
  database: "refseq_rna"   # 使用 RefSeq（更快）
  max_hits: 5              # 减少返回结果数量
  evalue: 0.001            # E-value 阈值
  delay: 12                # 查询间隔（秒）
  timeout: 300             # 超时时间（秒）
```

## 📝 输出文件说明

### 序列组装结果

```
sanger_assembly_results/
├── assembled_sequences.fasta    # 组装后的序列（FASTA 格式）
├── assembly_report.xlsx         # 详细报告（Excel）
├── assembly_report.txt          # 文本报告
├── assembly_summary.csv         # 摘要（CSV）
└── quality_report.txt           # 质量评估报告
```

### BLAST 分析结果

```
blast_results/
├── blast_results_20251101_143022.xlsx    # BLAST 结果（Excel）
├── blast_report_20251101_143022.txt      # 文本报告
└── blast_progress.json                    # 进度文件（断点续传用）
```

### 日志文件

```
analysis.log           # 完整分析日志
blast_analysis.log     # BLAST 分析日志
```

## ❓ 常见问题

### Q1: 如何处理单端测序数据？

**A:** 工具会自动识别单端数据，不会尝试组装，直接使用原始序列。

### Q2: 文件名不被识别怎么办？

**A:** 
1. 使用 `--show-patterns` 查看支持的格式
2. 在 `config.yaml` 中添加自定义规则
3. 或重命名文件为支持的格式

### Q3: BLAST 分析太慢？

**A:**
1. 使用更快的数据库：`-d refseq_rna`
2. 减少返回结果：`-n 5`
3. 利用断点续传功能

### Q4: 如何批量处理多个项目？

**A:** 编写简单的批处理脚本：

**Linux/macOS (bash):**
```bash
#!/bin/bash
for dir in project1 project2 project3; do
    echo "Processing $dir..."
    python sanger_16s_pipeline.py -i "./$dir" --batch
done
```

**Windows (PowerShell):**
```powershell
$projects = @("project1", "project2", "project3")
foreach ($dir in $projects) {
    Write-Host "Processing $dir..."
    python sanger_16s_pipeline.py -i ".\$dir" --batch
}
```

## 💡 提示和技巧

1. **首次使用**：建议先用交互式模式熟悉流程
2. **大批量数据**：使用批处理模式并启用断点续传
3. **测试配置**：先用 `--scan-only` 测试文件识别
4. **保存日志**：日志文件对于问题排查很有帮助
5. **定期备份**：重要的分析结果记得备份

## 📚 更多资源

- [主文档](../README.md)
- [配置文件说明](../config.yaml)
- [贡献指南](../CONTRIBUTING.md)
- [更新日志](../CHANGELOG.md)

## 🆘 获取帮助

如有问题：
- 查看 [常见问题](../README.md#-常见问题)
- 提交 [Issue](https://github.com/YOUR_USERNAME/sanger-16s-analysis/issues)
- 参与 [Discussions](https://github.com/YOUR_USERNAME/sanger-16s-analysis/discussions)

---

祝您使用愉快！🎉

