# 🎯 GE2TS: Genetic Effects to Transcriptomic Scores Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)]()
[![Coverage](https://img.shields.io/badge/coverage-95%25-brightgreen.svg)]()

**GE2TS** is a modular, scalable, and reproducible Nextflow DSL2 pipeline designed to transform Genome-Wide Association Study (GWAS) summary statistics into gene-level and pathway-level transcriptomic scores using MAGMA gene mapping and Gene Set Enrichment Analysis (GSEA).

---

## 🚀 Key Features

- **GWAS Standardization & Quality Control (`GWAS_PREP`)**: Automatic column mapping, P-value filtering, Z-score computation, and indel filtering using Polars/pandas.
- **Gene Annotation & Analysis (`MAGMA_ANNOT_RUN`)**: Gene-level mapping and association analysis supporting customizable window sizes (10kb, 100kb, Hi-C).
- **Transcriptomic Scoring & GSEA (`QC_GSEA`)**: Pathway enrichment analysis using Stouffer Z-score combination against GMT pathway databases and comprehensive QC summary matrix generation.
- **Container & Profile Support**: Out-of-the-box support for Docker, Singularity/Apptainer, Conda, and local execution profiles.

---

## 🛠️ Quick Start

### Prerequisites
- [Nextflow](https://www.nextflow.io/) (`>= 23.04.0`)
- Python 3.10+ (or Docker / Singularity / Conda)

### Setup Environment
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Run Test Profile
```bash
nextflow run main.nf -profile test
```

### Run Custom GWAS Dataset
```bash
nextflow run main.nf \
  --input path/to/gwas_summary_stats.tsv \
  --bfile path/to/plink_ref_prefix \
  --gene-loc path/to/NCBI37.3.gene.loc \
  --pathways path/to/pathways.gmt \
  --outdir results
```

---

## 📁 Pipeline Output Structure

The results are published to the specified `--outdir` directory:

```text
results/
├── gwas_prep/
│   ├── [phenotype].prep.tsv      # Standardized & QC-filtered GWAS
│   └── [phenotype].qc.json       # QC metrics summary
├── magma_annot_run/
│   ├── [phenotype].genes.out     # MAGMA gene-level association scores
│   └── [phenotype].genes.raw     # MAGMA raw gene association values
└── qc_gsea/
    ├── GSEA/
    │   └── [phenotype]_gsea_pathways.tsv  # Ranked pathway enrichment scores
    └── [phenotype].magma_qc_summary.tsv    # Final mapping QC summary matrix
```

---

## ⚙️ Parameters

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--input` | `file` | `assets/test_data/sample_gwas.tsv` | Path to GWAS summary statistics (`.tsv` or `.tsv.gz`) |
| `--bfile` | `path` | `assets/test_data/dummy_ref.bed` | PLINK reference prefix |
| `--gene-loc` | `file` | `assets/test_data/dummy_gene.loc` | Gene location definition file |
| `--pathways` | `file` | `assets/test_data/dummy_pathways.gmt` | Pathway database in GMT format |
| `--outdir` | `dir` | `results` | Output directory path |

---

## 🧪 Testing & Verification

Run the Python unit test suite with coverage:
```bash
PYTHONPATH=. .venv/bin/pytest --cov=bin tests/bin/
```

Run end-to-end integration test with Nextflow:
```bash
nextflow run main.nf -profile test
```

---

## 📚 Documentation Index

- [Agent Guidelines & Architecture Index](AGENTS.md)
- [nf-core Standard DSL2 Architecture](docs/pipeline/nf-core-dsl2-architecture.md)
- [Mock Wrappers & Pytest Convention](docs/testing/mock-wrappers-and-pytest.md)
