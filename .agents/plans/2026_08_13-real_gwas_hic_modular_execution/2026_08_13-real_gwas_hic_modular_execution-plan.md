---
name: "real_gwas_hic_modular_execution"
description: "Adapt GE2TS pipeline for real GWAS datasets with CSV/delimiter auto-detection, flexible column mapping, Hi-C BED annotation, and modular step-by-step execution."
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
has_completed_all_phases: true
---

# Goal

Enhance the GE2TS pipeline to seamlessly process real-world GWAS datasets (such as CIMBA BRCA1 with CSV comma separation and alternate column names), support Hi-C BED region annotation for MAGMA, and provide modular entry points and stop flags for step-by-step execution (GWAS Prep, MAGMA, GSEA).

# Context

- [gets.md](gets.md): High-level specification of pipeline goals.
- Reference Pipeline (`gwas-magma-pipeline`): Source of logic for CSV summary statistics, Hi-C annotations, and step-by-step execution.
- Target Architecture: nf-core standard Nextflow DSL2 structure with modular subfolders (`bin/`, `modules/local/`, `workflows/`, `assets/`, `conf/`).

# Public Contracts

- Pipeline Input/Output Parameters (`nextflow.config`, `nextflow_schema.json`):
  - `--input`: Raw GWAS input file.
  - `--input_prep`: Pre-prepared GWAS input file to skip Phase 1.
  - `--input_genes_out`: Pre-computed MAGMA genes output to skip Phases 1 and 2.
  - `--hic_bed`: Optional Hi-C mappable regions BED file.
  - `--stop_at_prep` / `--skip_magma`: Stop after GWAS preparation.
  - `--skip_gsea`: Stop after MAGMA gene analysis.
- `bin/gwas_prep.py` CLI:
  - Automatic delimiter detection (comma vs tab).
  - Synonym column mapping (`MarkerName` -> `SNP`, `position` -> `BP`, `Allele1`/`Allele2` -> `A1`/`A2`, `Effect` -> `BETA`, `StdErr` -> `SE`, `P-value` -> `P`).
- `bin/magma_wrapper.py` & `modules/local/magma_annot_run.nf`:
  - Supports `--hic-bed` parameter for integration of distal Hi-C interactions.
  - Handles PLINK binary reference tuple (`.bed`, `.bim`, `.fam`).
- Test Suites (`tests/bin/`):
  - Unit tests for CSV parsing, column renaming, Hi-C integration, and modular workflow branches.

# Phases

## Phase 1: Real GWAS CSV Auto-Delimiter & Column Mapping

### Description
Enhance `bin/gwas_prep.py` to automatically detect CSV vs TSV delimiters and map non-standard column names (e.g. `MarkerName`, `position`, `Allele1`, `Allele2`, `Effect`, `StdErr`, `P-value`) to standard GE2TS names, adding comprehensive unit tests.

### To-do Actions List
- [x] Update `bin/gwas_prep.py` with separator auto-detection and synonym column renaming dictionary.
- [x] Add unit test in `tests/bin/test_gwas_prep.py` covering comma-separated GWAS with alternate column headers.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [x] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

## Phase 2: Hi-C BED Annotation & PLINK Multi-file Handling

### Description
Extend `bin/magma_wrapper.py` and `modules/local/magma_annot_run.nf` to support optional Hi-C region annotation via `--hic-bed` and robust handling of PLINK reference file bundles.

### To-do Actions List
- [x] Update `bin/magma_wrapper.py` to accept `--hic-bed` and integrate Hi-C regions into annotation or mock generation.
- [x] Update `modules/local/magma_annot_run.nf` to handle PLINK reference sets and optional Hi-C BED files.
- [x] Add unit test in `tests/bin/test_magma_wrapper.py` verifying Hi-C BED integration in mock mode.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [x] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

## Phase 3: Modular Pipeline Entry Points & Step Control Flags

### Description
Configure `main.nf`, `workflows/ge2ts.nf`, and `nextflow.config` to support modular execution entry points (`--input_prep`, `--input_genes_out`) and stop flags (`--stop_at_prep`, `--skip_gsea`), validating with end-to-end tests.

### To-do Actions List
- [x] Update `workflows/ge2ts.nf` and `main.nf` to branch channels based on `--input`, `--input_prep`, or `--input_genes_out`.
- [x] Implement stop/skip conditional logic (`--stop_at_prep`, `--skip_gsea`) in `main.nf` and `workflows/ge2ts.nf`.
- [x] Update `nextflow_schema.json` and `nextflow.config` with new parameters.
- [x] Run end-to-end test with real CIMBA/dummy sample inputs and verify all workflows.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [x] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

# Next step

All phases of the real GWAS, Hi-C, and modular execution adaptation plan have been fully completed and verified!

Pipeline ready for real GWAS and Hi-C analysis with 🚀 🧬 by 🐢 💨 (Turbotuga™, [Codely](https://codely.com)’s mascot)
