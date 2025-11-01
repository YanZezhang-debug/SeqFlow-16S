# SeqFlow-16S: Intelligent Sanger 16S rRNA Analysis Toolkit

<div align="center">

![Python Version](https://img.shields.io/badge/python-3.7%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey)
![Status](https://img.shields.io/badge/status-stable-brightgreen)
![Maintenance](https://img.shields.io/badge/maintained-yes-brightgreen)

**From Sequencing Files to Species Identification - In One Flow**

[中文文档](README.md) | [Documentation](docs/) | [Examples](examples/) | [Contributing](CONTRIBUTING.md)

</div>

---

## 🌟 Overview

**SeqFlow-16S** is a comprehensive, intelligent toolkit designed for automated analysis of Sanger 16S rRNA sequencing data. It streamlines the entire workflow from raw sequencing files to species identification, eliminating the need for manual file sorting, sequence assembly, and BLAST queries.

### 🎯 Key Features

- **🎯 Intelligent Filename Recognition** - Automatically identifies 6+ naming conventions with customizable pattern support
- **🔄 Complete Analysis Pipeline** - One-click workflow from raw files to species identification
- **💾 Resume Capability** - BLAST analysis supports interruption and continuation
- **📊 Professional Reports** - Auto-generates Excel and text format reports
- **⚙️ Flexible Configuration** - YAML-based configuration with all parameters customizable
- **🖥️ Dual Modes** - Interactive wizard and batch automation modes
- **🌐 Multi-platform** - Works seamlessly on Windows, Linux, and macOS
- **📚 Rich Documentation** - Comprehensive guides, examples, and API documentation

---

## 📋 Table of Contents

- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Features](#-features)
- [Usage](#-usage)
- [Configuration](#-configuration)
- [Output](#-output)
- [Examples](#-examples)
- [Documentation](#-documentation)
- [Contributing](#-contributing)
- [License](#-license)
- [Citation](#-citation)

---

## 🚀 Installation

### Prerequisites

- Python 3.7 or higher
- Internet connection (for NCBI BLAST queries)

### Install Dependencies

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/seqflow-16s.git
cd seqflow-16s

# Install required packages
pip install -r requirements.txt
```

### Required Packages

- `biopython` - Sequence file parsing and BLAST
- `openpyxl` - Excel report generation
- `pyyaml` - Configuration file parsing
- `requests` - HTTP requests for NCBI API

---

## ⚡ Quick Start

### Interactive Mode (Recommended for First-Time Users)

```bash
python sanger_16s_pipeline.py
```

Follow the interactive prompts to:
1. Specify input directory containing sequencing files
2. Configure analysis parameters
3. Run sequence assembly
4. Perform BLAST species identification
5. Generate comprehensive reports

### Batch Mode (For Automation)

```bash
# Basic usage
python sanger_16s_pipeline.py -i ./data --batch

# With custom configuration
python sanger_16s_pipeline.py -i ./data -c custom_config.yaml --batch

# Skip BLAST analysis
python sanger_16s_pipeline.py -i ./data --batch --skip-blast

# Scan files only (no analysis)
python sanger_16s_pipeline.py -i ./data --scan-only
```

---

## ✨ Features

### 1. 🎯 Intelligent Filename Recognition

Automatically recognizes multiple naming conventions:

| Format | Example | Description |
|--------|---------|-------------|
| Standard | `0001_12345_(sample1)_[16SF].ab1` | Order number with sample ID |
| Simple | `sample1_F.ab1` | Basic sample + direction |
| Dash | `sample1-F.ab1` | Dash-separated format |
| Prefix | `F_sample1.ab1` | Direction prefix |
| Suffix | `sample1.F.ab1` | Dot-separated suffix |
| 16S Marker | `sample1_16SF.ab1` | With 16S identifier |

**Extensible**: Add custom patterns via configuration file!

### 2. 🧬 Graph-Based Sequence Assembly

Advanced assembly algorithm with:
- **Quality Control**: Automatic filtering of low-quality bases
- **End Trimming**: Removes unreliable terminal sequences
- **Overlap Detection**: Finds optimal alignment between forward and reverse reads
- **Consensus Generation**: Creates high-quality consensus sequences
- **Quality Metrics**: Comprehensive quality assessment

### 3. 🔍 NCBI BLAST Species Identification

Robust BLAST analysis featuring:
- **Multiple Databases**: Support for nt, 16S_ribosomal_RNA, and more
- **Resume Capability**: Automatic progress saving and recovery
- **Smart Retry**: Handles network errors gracefully
- **Rate Limiting**: Respects NCBI usage policies
- **Detailed Results**: Top hits with alignment details, E-values, and identity scores

### 4. 📊 Comprehensive Reporting

Professional output formats:
- **Excel Reports**: Multi-sheet workbooks with formatted tables
- **Text Reports**: Human-readable summaries
- **FASTA Files**: Standard sequence format for downstream analysis
- **CSV Files**: Easy import into other tools
- **Quality Reports**: Detailed QC metrics

### 5. ⚙️ Flexible Configuration

YAML-based configuration system:
```yaml
# Customize filename patterns
filename_patterns:
  - name: "Custom Format"
    pattern: "^([A-Z0-9]+)_(F|R)\\.ab1$"
    sample_group: 1
    direction_group: 2

# Adjust assembly parameters
assembly:
  min_overlap: 20
  min_quality: 20
  trim_ends: true
  trim_length: 20

# Configure BLAST settings
blast:
  database: "nt"
  max_hits: 10
  delay: 10
```

---

## 📖 Usage

### Command-Line Interface

```bash
python sanger_16s_pipeline.py [OPTIONS]

Options:
  -i, --input DIR          Input directory containing sequencing files
  -o, --output DIR         Output directory for results (default: current dir)
  -c, --config FILE        Custom configuration file (default: config.yaml)
  --batch                  Run in batch mode (non-interactive)
  --skip-blast             Skip BLAST analysis
  --scan-only              Scan files only, no analysis
  --show-patterns          Display supported filename patterns
  -h, --help               Show help message
```

### Individual Modules

Each module can be used independently:

#### Filename Parser

```python
from filename_parser import FileNameParser

parser = FileNameParser()
result = parser.parse("sample1_F.ab1")
print(f"Sample: {result['sample_id']}, Direction: {result['direction']}")
```

#### Sequence Assembly

```python
from graph_based_assembly import GraphBasedAssembler

assembler = GraphBasedAssembler()
assembled = assembler.assemble_from_files("sample1_F.ab1", "sample1_R.ab1")
print(f"Assembled length: {len(assembled.seq)}")
```

#### BLAST Analysis

```python
from ncbi_blast_analysis import NCBIBlastAnalyzer

analyzer = NCBIBlastAnalyzer()
results = analyzer.blast_fasta("sequences.fasta", database="nt")
```

---

## ⚙️ Configuration

### Configuration File Structure

The `config.yaml` file contains all configurable parameters:

```yaml
# Filename recognition patterns
filename_patterns:
  - name: "Standard Format"
    pattern: "^\\d+_[^_]+_\\(([^)]+)\\)_\\[16S([FR])\\]\\.ab1$"
    sample_group: 1
    direction_group: 2
  # Add more patterns...

# Assembly parameters
assembly:
  min_overlap: 20          # Minimum overlap length
  min_quality: 20          # Minimum quality score
  trim_ends: true          # Enable end trimming
  trim_length: 20          # Bases to trim from each end
  match_score: 2           # Score for matching bases
  mismatch_penalty: -1     # Penalty for mismatches
  gap_penalty: -2          # Penalty for gaps

# BLAST parameters
blast:
  database: "nt"           # NCBI database
  max_hits: 10             # Maximum hits to return
  delay: 10                # Delay between queries (seconds)
  timeout: 300             # Query timeout (seconds)
  max_retries: 3           # Maximum retry attempts

# Output settings
output:
  save_excel: true         # Generate Excel reports
  save_text: true          # Generate text reports
  save_fasta: true         # Save assembled sequences
  save_csv: true           # Generate CSV files

# Logging
logging:
  level: "INFO"            # DEBUG, INFO, WARNING, ERROR
  file: "seqflow.log"      # Log file name
```

---

## 📁 Output

### Directory Structure

```
project_directory/
├── sanger_assembly_results/
│   ├── assembled_sequences.fasta      # Assembled sequences
│   ├── assembly_report.xlsx           # Detailed assembly report
│   ├── assembly_report.txt            # Text format report
│   ├── assembly_summary.csv           # Summary statistics
│   └── quality_report.txt             # Quality control metrics
│
└── blast_results/
    ├── blast_results_YYYYMMDD_HHMMSS.xlsx   # BLAST results
    ├── blast_report_YYYYMMDD_HHMMSS.txt     # Detailed BLAST report
    └── blast_progress.json                   # Progress file (for resume)
```

### Report Contents

#### Assembly Report
- Sample ID
- Sequence length
- Average quality score
- GC content
- Assembly status
- Quality metrics

#### BLAST Report
- Query sequence ID
- Top hit species
- Percentage identity
- Query coverage
- E-value
- Accession number
- Alignment details

---

## 📚 Examples

### Example 1: Basic Analysis

```bash
# Analyze all files in a directory
python sanger_16s_pipeline.py -i ./my_sequences --batch
```

### Example 2: Custom Configuration

```bash
# Use custom settings
python sanger_16s_pipeline.py -i ./data -c my_config.yaml --batch
```

### Example 3: Assembly Only

```bash
# Skip BLAST analysis
python sanger_16s_pipeline.py -i ./data --batch --skip-blast
```

### Example 4: File Scanning

```bash
# Check which files will be recognized
python sanger_16s_pipeline.py -i ./data --scan-only
```

### Example 5: Python API

```python
from sanger_16s_pipeline import SangerPipeline

# Create pipeline instance
pipeline = SangerPipeline(
    input_dir="./data",
    config_file="config.yaml"
)

# Run analysis
pipeline.run(skip_blast=False)

# Access results
assembly_results = pipeline.get_assembly_results()
blast_results = pipeline.get_blast_results()
```

For more examples, see the [examples](examples/) directory.

---

## 📖 Documentation

- **[Project Structure](docs/PROJECT_STRUCTURE.md)** - Detailed project organization
- **[GitHub Guide](docs/GITHUB_GUIDE.md)** - Publishing and maintenance guide
- **[Quick Publish Guide](PUBLISH_TO_GITHUB.md)** - Quick start for GitHub
- **[Examples](examples/README.md)** - Usage examples and tutorials
- **[Changelog](CHANGELOG.md)** - Version history and updates
- **[Contributing](CONTRIBUTING.md)** - How to contribute

---

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Ways to Contribute

- 🐛 Report bugs
- 💡 Suggest new features
- 📝 Improve documentation
- 🔧 Submit pull requests
- ⭐ Star the project

### Development Setup

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/seqflow-16s.git
cd seqflow-16s

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run tests (when available)
python -m pytest tests/
```

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2025 Zhang Yanze

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction...
```

---

## 📊 Project Statistics

![GitHub stars](https://img.shields.io/github/stars/YOUR_USERNAME/seqflow-16s?style=social)
![GitHub forks](https://img.shields.io/github/forks/YOUR_USERNAME/seqflow-16s?style=social)
![GitHub watchers](https://img.shields.io/github/watchers/YOUR_USERNAME/seqflow-16s?style=social)

![GitHub issues](https://img.shields.io/github/issues/YOUR_USERNAME/seqflow-16s)
![GitHub pull requests](https://img.shields.io/github/issues-pr/YOUR_USERNAME/seqflow-16s)
![GitHub last commit](https://img.shields.io/github/last-commit/YOUR_USERNAME/seqflow-16s)
![GitHub contributors](https://img.shields.io/github/contributors/YOUR_USERNAME/seqflow-16s)

---

## 📮 Contact & Support

- **Issues**: [GitHub Issues](https://github.com/YOUR_USERNAME/seqflow-16s/issues)
- **Discussions**: [GitHub Discussions](https://github.com/YOUR_USERNAME/seqflow-16s/discussions)
- **Email**: your.email@example.com

---

## 🙏 Acknowledgments

- **Biopython** - For excellent bioinformatics tools
- **NCBI** - For providing BLAST services
- **Community** - For feedback and contributions

---

## 📚 Citation

If you use SeqFlow-16S in your research, please cite:

```bibtex
@software{seqflow16s2025,
  author = {Zhang, Yanze},
  title = {SeqFlow-16S: Intelligent Sanger 16S rRNA Analysis Toolkit},
  year = {2025},
  url = {https://github.com/YOUR_USERNAME/seqflow-16s},
  version = {1.0.0}
}
```

---

## 🗺️ Roadmap

### Version 1.x
- [x] Core functionality
- [x] Interactive and batch modes
- [x] BLAST resume capability
- [x] Comprehensive documentation

### Version 2.0 (Planned)
- [ ] GUI interface
- [ ] Local BLAST support
- [ ] Multi-threading for faster processing
- [ ] Advanced quality filtering
- [ ] Phylogenetic tree generation
- [ ] Docker container support

---

## ⭐ Star History

If you find this project helpful, please consider giving it a star! ⭐

[![Star History Chart](https://api.star-history.com/svg?repos=YOUR_USERNAME/seqflow-16s&type=Date)](https://star-history.com/#YOUR_USERNAME/seqflow-16s&Date)

---

## 📈 Performance

- **Processing Speed**: ~10-20 samples per minute (assembly)
- **BLAST Queries**: Respects NCBI rate limits (3 queries per second)
- **Memory Usage**: Minimal (<100MB for typical datasets)
- **File Support**: .ab1, .fasta, .seq formats

---

## 🔒 Security

- No sensitive data is stored or transmitted
- All NCBI queries use HTTPS
- Configuration files are local only
- No external dependencies with known vulnerabilities

---

## 🌍 Supported Platforms

| Platform | Status | Notes |
|----------|--------|-------|
| Windows 10/11 | ✅ Fully Supported | Tested on Windows 10/11 |
| Linux | ✅ Fully Supported | Tested on Ubuntu 20.04+ |
| macOS | ✅ Fully Supported | Tested on macOS 11+ |

---

## 💡 Tips & Best Practices

1. **Organize Files**: Keep forward and reverse reads in the same directory
2. **Naming Convention**: Use consistent filename patterns
3. **Quality Check**: Review assembly quality reports before BLAST
4. **Network**: Ensure stable internet connection for BLAST queries
5. **Configuration**: Customize settings for your specific needs
6. **Backup**: Keep original sequencing files safe

---

<div align="center">

**Made with ❤️ for the Bioinformatics Community**

[⬆ Back to Top](#seqflow-16s-intelligent-sanger-16s-rrna-analysis-toolkit)

</div>

