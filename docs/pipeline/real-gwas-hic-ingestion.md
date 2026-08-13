# 🎯 Real GWAS Ingestion, Delimiter Auto-Detection, and Hi-C Integration Convention

## 💡 Convention

GWAS ingestion scripts (`bin/gwas_prep.py`) must automatically detect file delimiters (CSV vs TSV) and map non-standard column headers (e.g. `MarkerName`, `position`, `Allele1`, `Allele2`, `Effect`, `StdErr`, `P-value`, `onco_icogs_*`) to standard GE2TS names (`SNP`, `BP`, `A1`, `A2`, `BETA`, `SE`, `P`). Furthermore, MAGMA execution modules must support optional Hi-C mappable regions BED files (`--hic_bed`) and configurable annotation window sizes (`--window`).

## 🏆 Benefits

- Enables seamless ingestion of diverse public GWAS summary statistics formats without manual file formatting.
- Increases biological mapping accuracy by incorporating high-resolution chromatin conformation (Hi-C) interaction data.
- Gives users full command-line flexibility over annotation window sizes (e.g., `10,10` vs `100,100`).

## 👀 Examples

### ✅ Good: Auto-delimiter detection and synonym column mapping in Python

```python
sep = detect_delimiter(input_path)
df = pd.read_csv(input_path, sep=sep)
df = df.rename(columns=SYNONYM_MAP)
```

### ❌ Bad: Hardcoding tab separation and strict fixed column names

```python
# Fails when ingesting comma-separated CSV files or alternate column headers
df = pd.read_csv(input_path, sep='\t')
assert "CHR" in df.columns
```

## 🧐 Real world examples

- [bin/gwas_prep.py](../../bin/gwas_prep.py)
- [bin/magma_wrapper.py](../../bin/magma_wrapper.py)
- [modules/local/magma_annot_run.nf](../../modules/local/magma_annot_run.nf)

## 🔗 Related agreements

- [nf-core Standard DSL2 Architecture](nf-core-dsl2-architecture.md)
- [Modular Step-by-Step Pipeline Execution](modular-step-by-step-pipeline-execution.md)
- [AGENTS.md](../../AGENTS.md)

Doc created by 🐢 💨 (Turbotuga™, [Codely](https://codely.com)’s mascot)
