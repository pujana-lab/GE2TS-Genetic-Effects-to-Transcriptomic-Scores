# 🎯 Modular Step-by-Step Pipeline Execution Convention

## 💡 Convention

Pipelines in GE2TS must support modular execution points (`--input`, `--input_prep`, `--input_genes_out`) and inter-step control flags (`--stop_at_prep`, `--skip_magma`, `--skip_gsea`), alongside automatic delimiter detection and synonym column mapping for real-world GWAS datasets.

## 🏆 Benefits

- Allows users and researchers to resume, debug, or start at any intermediate step without re-running upstream compute-heavy processes.
- Provides flexibility to ingest pre-prepared GWAS summary statistics or pre-computed MAGMA gene scores directly.
- Handles diverse real-world file formats (CSV vs TSV) and alternate header nomenclatures transparently.

## 👀 Examples

### ✅ Good: Modular channel branching and step control in Nextflow workflows

```nextflow
if (!params.input_prep && !params.input_genes_out) {
    GWAS_PREP(ch_gwas)
    ch_current_prep = GWAS_PREP.out.gwas_prep
} else {
    ch_current_prep = ch_prep
}
```

### ❌ Bad: Monolithic pipelines requiring full re-execution from raw inputs on every change

```nextflow
// Forcing full pipeline execution from raw GWAS every single time
GWAS_PREP(ch_gwas)
MAGMA_ANNOT_RUN(GWAS_PREP.out.gwas_prep, ...)
QC_GSEA(MAGMA_ANNOT_RUN.out.genes_out, ...)
```

## 🧐 Real world examples

- [main.nf](../../main.nf)
- [workflows/ge2ts.nf](../../workflows/ge2ts.nf)
- [bin/gwas_prep.py](../../bin/gwas_prep.py)
- [bin/magma_wrapper.py](../../bin/magma_wrapper.py)

## 🔗 Related agreements

- [nf-core Standard DSL2 Architecture](nf-core-dsl2-architecture.md)
- [AGENTS.md](../../AGENTS.md)

Doc created by 🐢 💨 (Turbotuga™, [Codely](https://codely.com)’s mascot)
