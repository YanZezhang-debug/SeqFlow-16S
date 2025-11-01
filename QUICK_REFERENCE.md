# SeqFlow-16S 快速参考卡片

## 📋 项目信息

| 项目 | 信息 |
|------|------|
| **项目名称** | SeqFlow-16S |
| **版本** | v1.0.0 |
| **作者** | Zhang Yanze |
| **许可证** | MIT |
| **GitHub** | https://github.com/YanZezhang-debug/SeqFlow-16S |
| **标语** | From Sequencing Files to Species Identification - In One Flow |

---

## 🚀 快速开始

### 安装
```bash
git clone https://github.com/YanZezhang-debug/SeqFlow-16S.git
cd SeqFlow-16S
pip install -r requirements.txt
```

### 运行
```bash
# 交互模式
python sanger_16s_pipeline.py

# 批处理模式
python sanger_16s_pipeline.py -i ./data --batch

# 跳过 BLAST
python sanger_16s_pipeline.py -i ./data --batch --skip-blast
```

---

## 📁 核心文件

| 文件 | 说明 |
|------|------|
| `sanger_16s_pipeline.py` | 主程序入口 |
| `filename_parser.py` | 文件名解析 |
| `graph_based_assembly.py` | 序列组装 |
| `ncbi_blast_analysis.py` | BLAST 分析 |
| `config.yaml` | 配置文件 |
| `requirements.txt` | 依赖列表 |

---

## 🎯 核心功能

- ✅ 智能文件名识别（6+ 格式）
- ✅ 自动序列组装
- ✅ NCBI BLAST 鉴定
- ✅ 多格式报告生成
- ✅ 断点续传支持
- ✅ 灵活配置管理

---

## 📖 文档导航

| 文档 | 用途 |
|------|------|
| [README.md](README.md) | 中文主文档 |
| [README_EN.md](README_EN.md) | 英文主文档 |
| [CHANGELOG.md](CHANGELOG.md) | 版本历史 |
| [CONTRIBUTING.md](CONTRIBUTING.md) | 贡献指南 |
| [docs/PROJECT_STRUCTURE.md](docs/PROJECT_STRUCTURE.md) | 项目结构 |
| [docs/GITHUB_GUIDE.md](docs/GITHUB_GUIDE.md) | GitHub 指南 |
| [examples/README.md](examples/README.md) | 使用示例 |
| [PROJECT_COMPLETION.md](PROJECT_COMPLETION.md) | 项目总结 |

---

## ⚙️ 配置要点

### 文件名模式
```yaml
filename_patterns:
  - name: "Standard Format"
    pattern: "^\\d+_[^_]+_\\(([^)]+)\\)_\\[16S([FR])\\]\\.ab1$"
    sample_group: 1
    direction_group: 2
```

### 组装参数
```yaml
assembly:
  min_overlap: 20
  min_quality: 20
  trim_ends: true
  trim_length: 20
```

### BLAST 参数
```yaml
blast:
  database: "nt"
  max_hits: 10
  delay: 10
```

---

## 📊 输出文件

### 组装结果
- `assembled_sequences.fasta` - 组装序列
- `assembly_report.xlsx` - Excel 报告
- `assembly_report.txt` - 文本报告
- `assembly_summary.csv` - 汇总数据

### BLAST 结果
- `blast_results_*.xlsx` - BLAST Excel 报告
- `blast_report_*.txt` - BLAST 文本报告
- `blast_progress.json` - 进度文件

---

## 🔧 命令行选项

```bash
python sanger_16s_pipeline.py [OPTIONS]

选项:
  -i, --input DIR       输入目录
  -o, --output DIR      输出目录
  -c, --config FILE     配置文件
  --batch               批处理模式
  --skip-blast          跳过 BLAST
  --scan-only           仅扫描文件
  --show-patterns       显示支持的模式
  -h, --help            帮助信息
```

---

## 🐛 常见问题

### 1. 文件无法识别
- 检查文件命名格式
- 使用 `--show-patterns` 查看支持的格式
- 在 config.yaml 中添加自定义模式

### 2. BLAST 查询失败
- 检查网络连接
- 等待后重试（自动重试 3 次）
- 使用 `--skip-blast` 跳过 BLAST

### 3. 组装质量低
- 检查原始序列质量
- 调整 `min_quality` 参数
- 增加 `trim_length` 值

---

## 📞 获取帮助

- **文档**: 查看 docs/ 目录
- **示例**: 查看 examples/ 目录
- **Issues**: https://github.com/YanZezhang-debug/SeqFlow-16S/issues
- **Discussions**: https://github.com/YanZezhang-debug/SeqFlow-16S/discussions

---

## 🤝 贡献

欢迎贡献！请查看 [CONTRIBUTING.md](CONTRIBUTING.md)

1. Fork 项目
2. 创建分支
3. 提交更改
4. 推送分支
5. 创建 PR

---

## 📜 许可证

MIT License - 详见 [LICENSE](LICENSE)

---

## ⭐ 支持项目

如果觉得有用，请给项目点个 Star！⭐

https://github.com/YanZezhang-debug/SeqFlow-16S

---

<div align="center">

**Made with ❤️ for the Bioinformatics Community**

</div>

