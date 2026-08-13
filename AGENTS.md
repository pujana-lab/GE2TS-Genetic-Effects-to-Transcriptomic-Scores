# 🤖 Agent Guidelines & Index

Welcome to the **GE2TS** repository. This project implements a modular Nextflow DSL2 pipeline to convert GWAS summary statistics into transcriptomic scores with pathway analysis (GSEA).

## 🛠️ Verification & Testing Commands

- **Run Python Unit Tests**: `PYTHONPATH=. .venv/bin/pytest tests/bin/`
- **Run Unit Tests with Coverage**: `PYTHONPATH=. .venv/bin/pytest --cov=bin tests/bin/`
- **Run End-to-End Nextflow Test Workflow**: `nextflow run main.nf -profile test`

## 📚 Project Documentation Index

- [nf-core Standard DSL2 Architecture](docs/pipeline/nf-core-dsl2-architecture.md)
- [Mock Wrappers and Pytest Convention](docs/testing/mock-wrappers-and-pytest.md)
