# 🎯 Mock Wrappers and Pytest Convention

## 💡 Convention

All CLI helper tools in `bin/` must support standalone mock execution flags or fallbacks when external bioinformatic binaries (e.g. MAGMA executable) are absent. Unit tests must be placed in `tests/bin/` and executed using `pytest` to maintain high code coverage (>80%).

## 🏆 Benefits

- Enables continuous integration (CI) and local testing without requiring heavy external bioinformatic binaries.
- Ensures fast feedback loop for python script changes before running full Nextflow executions.
- Guarantees high test coverage across data parsing, filtering, and pathway scoring logic.

## 👀 Examples

### ✅ Good: Mock mode fallback in CLI script and corresponding pytest suite

```python
if has_magma and not mock:
    subprocess.run([magma_bin, ...])
else:
    # Generate mock outputs for testing
    write_mock_outputs(out_prefix)
```

```python
def test_run_magma_mock():
    run_magma(..., mock=True)
    assert os.path.exists("test.genes.out")
```

### ❌ Bad: Hardcoding external binary execution without fallback

```python
# Fails in CI when binary is missing
subprocess.run(["/usr/local/bin/magma", ...], check=True)
```

## 🧐 Real world examples

- [bin/magma_wrapper.py](../../bin/magma_wrapper.py)
- [tests/bin/test_magma_wrapper.py](../../tests/bin/test_magma_wrapper.py)
- [tests/bin/test_gwas_prep.py](../../tests/bin/test_gwas_prep.py)
- [tests/bin/test_qc_gsea.py](../../tests/bin/test_qc_gsea.py)

## 🔗 Related agreements

- [AGENTS.md](../../AGENTS.md)

Doc created by 🐢 💨 (Turbotuga™, [Codely](https://codely.com)’s mascot)
