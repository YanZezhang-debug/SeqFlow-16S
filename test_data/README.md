# Test Data for SeqFlow-16S

## 📁 Overview

This directory contains sample Sanger sequencing files (`.ab1` format) for testing and demonstration purposes.

## 📊 Sample Information

| Sample ID | Description | Forward File | Reverse File |
|-----------|-------------|--------------|--------------|
| 3sm1 | Sample 1 | 0001_31525103103397_(3sm1)_[16SF].ab1 | 0002_31525103103397_(3sm1)_[16SR].ab1 |
| 3sm2 | Sample 2 | 0003_31525103103398_(3sm2)_[16SF].ab1 | 0004_31525103103398_(3sm2)_[16SR].ab1 |
| 3sm3 | Sample 3 | 0005_31525103103399_(3sm3)_[16SF].ab1 | 0006_31525103103399_(3sm3)_[16SR].ab1 |
| 3sm4 | Sample 4 | 0007_31525103103400_(3sm4)_[16SF].ab1 | 0008_31525103103400_(3sm4)_[16SR].ab1 |
| 3sm5 | Sample 5 | 0009_31525103103401_(3sm5)_[16SF].ab1 | 0010_31525103103401_(3sm5)_[16SR].ab1 |
| 3sm6 | Sample 6 | 0011_31525103103402_(3sm6)_[16SF].ab1 | 0012_31525103103402_(3sm6)_[16SR].ab1 |
| 3sm7 | Sample 7 | 0013_31525103103403_(3sm7)_[16SF].ab1 | 0014_31525103103403_(3sm7)_[16SR].ab1 |
| 1sm6 | Sample 8 | 0015_31525103103404_(1sm6)_[16SF].ab1 | 0016_31525103103404_(1sm6)_[16SR].ab1 |

## 🧬 File Format

These files follow the **Standard Format** naming convention:

```
<序号>_<日期>_(<样本编号>)_[16S<方向>].ab1
```

Where:
- **序号**: Sequential number (e.g., 0001, 0002)
- **日期**: Date identifier (e.g., 31525103103397)
- **样本编号**: Sample ID (e.g., 3sm1, 3sm2)
- **方向**: Direction - F (Forward) or R (Reverse)

## 🚀 Quick Test

### Test Assembly Only

```bash
# From project root
python sanger_16s_pipeline.py -i test_data -o test_output --skip-blast
```

This will:
- Parse all 16 ab1 files
- Assemble 8 sample sequences
- Generate assembly reports
- Skip BLAST queries (for quick testing)

### Full Pipeline Test

```bash
# From project root
python sanger_16s_pipeline.py -i test_data -o test_output --batch
```

This will:
- Parse and assemble sequences
- Perform NCBI BLAST queries
- Generate complete identification reports

**Note**: Full pipeline requires internet connection and takes ~10-15 minutes.

## 📊 Expected Output

After running the pipeline, you should see:

```
test_output/
├── assembled_sequences.fasta          # 8 assembled sequences
├── assembly_report.xlsx               # Excel report
├── assembly_report.txt                # Text report
├── assembly_summary.csv               # CSV summary
├── blast_results_<timestamp>.xlsx     # BLAST results (if not skipped)
└── blast_report_<timestamp>.txt       # BLAST text report (if not skipped)
```

## 🔍 Assembly Quality

These test samples should produce high-quality assemblies:
- Overlap length: 400-600 bp
- Consensus length: 1200-1400 bp
- Quality score: 95-100%

## 📝 Notes

- These are real Sanger sequencing data from bacterial 16S rRNA gene analysis
- Files are in ABI format (.ab1), containing both sequence and quality scores
- Each sample has forward (F) and reverse (R) reads for bidirectional sequencing
- Data is suitable for testing all features of SeqFlow-16S

## 🤝 Usage in Examples

These test files are used in:
- [Basic Usage Example](../examples/README.md#basic-usage)
- [Batch Processing Example](../examples/README.md#batch-processing)
- [Custom Configuration Example](../examples/README.md#custom-configuration)

## 📄 License

These test data files are provided for demonstration and testing purposes only.

