---
name: "adapt_ge2ts_nf_core_pipeline"
description: "Adapt GE2TS pipeline to nf-core standard DSL2 architecture using GWAS_MAGMA_PREP design patterns and dummy test datasets."
created_at: "2026-08-13T00:00:00Z"

created_by:
  tool: "opencode"
  model:
    name: "Gemini"
    version: "3.6-flash"
    reasoning_effort: "high"

implemented_by:
  tool: "opencode"
  model:
    name: "Gemini"
    version: "3.6-flash"
    reasoning_effort: "high"

last_implementation_at: "2026-08-13T00:00:00Z"
has_completed_all_phases: false
---

# Goal

Adapt and implement the `GE2TS` pipeline following standard nf-core Nextflow DSL2 architectural patterns. The goal is to build a reproducible, modular pipeline for converting GWAS summary statistics into transcriptomic scores with pathway analysis (GSEA), drawing reference patterns from `gwas-magma-pipeline`.

# Context

- [gets.md](gets.md): High-level specification of the GE2TS pipeline goals, modules, and expected outputs.
- Reference Pipeline (`gwas-magma-pipeline` at `~/Documentos/projects/tools/gwas-magma-pipeline`): Source of logic for GWAS standardization, quality control, indel filtering, and MAGMA annotations.
- Target Architecture: nf-core standard Nextflow DSL2 structure with modular subfolders (`bin/`, `modules/local/`, `workflows/`, `assets/`, `conf/`, `assets/test_data/`).

# Public Contracts

- Pipeline Input/Output Parameters: Defined in `nextflow.config` and validated via `nextflow_schema.json`.
- Nextflow DSL2 Module Channels:
  - `GWAS_PREP`: Accepts `tuple val(meta), path(gwas)` -> Emits `tuple val(meta), path("*.prep.tsv")`, `path("*.qc.log")`.
  - `MAGMA_ANNOT_RUN`: Accepts `tuple val(meta), path(prep_gwas)` -> Emits `tuple val(meta), path("*.genes.out")`, `path("*.genes.raw")`.
  - `QC_GSEA`: Accepts `tuple val(meta), path(genes_out)` -> Emits `tuple val(meta), path("GSEA/")`, `path("magma_mapping_qc_summary.tsv")`.
- Python Scripts CLI Contracts:
  - `bin/gwas_prep.py`: `--input`, `--genome`, `--out-tsv`, `--qc-log`.
  - `bin/qc_gsea.py`: `--genes-out`, `--pathway-db`, `--out-dir`, `--summary-tsv`.
- Test Suites: Pytest unit tests under `tests/bin/` and Nextflow workflow execution tests using dummy test datasets in `assets/test_data/`.

# Phases

## Phase 1: Project Scaffolding, Dummy Test Data & `GWAS_Prep` Module

### Description
Establish the nf-core directory structure, create small synthetic/dummy test datasets (dummy GWAS, dummy reference, dummy pathway DB) for rapid pipeline verification, and implement the end-to-end `GWAS_PREP` vertical slice (Python script `bin/gwas_prep.py`, Nextflow module `modules/local/gwas_prep.nf`, and unit tests).

### To-do Actions List
- [x] Scaffold standard nf-core directory structure (`bin/`, `modules/local/`, `workflows/`, `assets/test_data/`, `conf/`, `tests/`).
- [x] Create synthetic test datasets in `assets/test_data/` (`sample_gwas.tsv`, `dummy_ref.bed`, `dummy_pathways.gmt`).
- [x] Implement `bin/gwas_prep.py` for column standardization, SNP QC, allele alignment, and indel filtering based on `gwas-magma-pipeline`.
- [x] Create unit tests in `tests/bin/test_gwas_prep.py` to verify `gwas_prep.py` CLI and filtering logic on dummy data.
- [x] Implement Nextflow DSL2 module `modules/local/gwas_prep.nf` wrapping `bin/gwas_prep.py`.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [x] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

## Phase 2: `MAGMA_Annot_Run` Module

### Description
Implement the `MAGMA_ANNOT_RUN` vertical slice to handle gene annotation (10kb, 100kb, or Hi-C windows) and gene-based association analysis using MAGMA, including container/tool execution and error handling.

### To-do Actions List
- [x] Implement Nextflow DSL2 module `modules/local/magma_annot_run.nf` handling annotation and gene analysis steps.
- [x] Add configuration handling in `conf/modules.config` for MAGMA window sizes (10kb, 100kb, Hi-C) and binary arguments.
- [x] Create test case validating `MAGMA_ANNOT_RUN` module execution using dummy/mocked MAGMA outputs.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [x] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

## Phase 3: `QC_GSEA` Module & Transcriptomic Scoring

### Description
Implement the `QC_GSEA` vertical slice to process MAGMA gene output files, perform Gene Set Enrichment Analysis (GSEA) against pathway databases, and compile the final mapping QC summary matrix.

### To-do Actions List
- [ ] Implement `bin/qc_gsea.py` to run pathway enrichment on gene scores and format output matrices.
- [ ] Create unit tests in `tests/bin/test_qc_gsea.py` verifying pathway scoring logic with dummy pathways GMT.
- [ ] Implement Nextflow DSL2 module `modules/local/qc_gsea.nf` wrapping `bin/qc_gsea.py`.
- [ ] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

## Phase 4: Workflow Integration, Configuration & End-to-End Pipeline

### Description
Assemble the complete end-to-end Nextflow workflow in `workflows/ge2ts.nf` and `main.nf`, configure `nextflow.config`, `nextflow_schema.json`, and profile options (docker, singularity, conda, test), and execute an end-to-end test using the dummy datasets.

### To-do Actions List
- [ ] Build the main workflow `workflows/ge2ts.nf` linking `GWAS_PREP`, `MAGMA_ANNOT_RUN`, and `QC_GSEA`.
- [ ] Create entry point `main.nf` with parameter validation via schema.
- [ ] Configure `nextflow.config` with container definitions, default parameters, and execution profiles (`test`, `docker`, `singularity`, `conda`).
- [ ] Run full end-to-end integration test with `nextflow run main.nf -profile test`.
- [ ] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

# Next step

The next step is to implement Phase 3 to create the `QC_GSEA` module for pathway enrichment and transcriptomic score calculation.

MAGMA gene analysis leveled up with 🧬 📊 by 🐢 💨 (Turbotuga™, [Codely](https://codely.com)’s mascot)
