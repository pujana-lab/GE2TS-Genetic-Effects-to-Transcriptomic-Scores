---
name: "add_cli_window_parameter"
description: "Add configurable --window parameter to command line interface, nextflow.config, nextflow_schema.json, and MAGMA wrapper/module."
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

Expose the MAGMA annotation window size as a command-line parameter (`--window`, default `"10,10"`) in Nextflow, passing it down to the `MAGMA_ANNOT_RUN` module and `bin/magma_wrapper.py`.

# Context

- `nextflow.config`: Contains pipeline parameters.
- `nextflow_schema.json`: Parameter validation schema.
- `main.nf`: CLI help message and parameter logging.
- `conf/modules.config`: Module ext.window settings.
- `modules/local/magma_annot_run.nf`: Nextflow module wrapping MAGMA.
- `bin/magma_wrapper.py`: Python wrapper script handling window annotation arguments.

# Public Contracts

- CLI Parameter: `--window` (string, default `"10,10"`).
- Schema definition in `nextflow_schema.json`.

# Phases

## Phase 1: Pipeline Configuration & Schema Update

### Description
Add `params.window = '10,10'` to `nextflow.config`, include it in `nextflow_schema.json`, document it in `main.nf` help message, and update `conf/modules.config` to reference `params.window`.

### To-do Actions List
- [x] Add `window = '10,10'` to `params` in `nextflow.config`.
- [x] Add `window` definition in `nextflow_schema.json`.
- [x] Update `main.nf` help message and execution log to display `--window`.
- [x] Update `conf/modules.config` to set `ext.window = params.window`.
- [x] Verify the changes in terms of typechecking, linting and tests using pytest.
- [x] STOP. Present the changes to the user for review and suggest commit messages.

## Phase 2: Module & Wrapper Integration

### Description
Ensure `modules/local/magma_annot_run.nf` and `bin/magma_wrapper.py` correctly forward and utilize the configurable window argument.

### To-do Actions List
- [x] Verify `modules/local/magma_annot_run.nf` uses `params.window` or `task.ext.window`.
- [x] Add unit test verifying custom window argument in `bin/magma_wrapper.py`.
- [x] Run end-to-end Nextflow test with custom `--window "100,100"`.
- [x] Verify the changes in terms of typechecking, linting and tests using pytest.
- [x] STOP. Present the changes to the user for review and suggest commit messages.

# Next step

All phases of the `--window` parameter addition have been successfully completed and verified!

Window configuration ready with 🪟 ⚙️ by 🐢 💨 (Turbotuga™, [Codely](https://codely.com)’s mascot)
