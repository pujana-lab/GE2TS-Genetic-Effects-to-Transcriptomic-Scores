# 🎯 nf-core Standard DSL2 Architecture Convention

## 💡 Convention

Bioinformatic pipelines in GE2TS must follow the standard **nf-core Nextflow DSL2** architecture rather than heavy DDD / Hexagonal layers. All process logic is encapsulated in modular DSL2 process files under `modules/local/`, orchestrated by `workflows/ge2ts.nf`, and executed via `main.nf`. Helper processing logic belongs in executable Python scripts inside `bin/`.

## 🏆 Benefits

- High community adoption, readability, and maintainability for bioinformaticians.
- Native Nextflow caching, parallelism, and multi-environment container execution (Docker, Apptainer, Conda).
- Clear separation between workflow orchestration (`modules/local/*.nf`) and data transformation logic (`bin/*.py`).

## 👀 Examples

### ✅ Good: Standard DSL2 process wrapping a clean Python CLI script

```nextflow
process GWAS_PREP {
    tag "$meta.id"
    label '$meta.id'

    input:
    tuple val(meta), path(gwas)

    output:
    tuple val(meta), path("*.prep.tsv"), emit: gwas_prep

    script:
    """
    gwas_prep.py --input ${gwas} --output-tsv ${meta.id}.prep.tsv
    """
}
```

### ❌ Bad: Heavy DDD entities/adapters abstraction inside bioinformatic scripts

```python
class GWASRecordEntity:
    def __init__(self, rsid: str, pvalue: float):
        self._rsid = rsid
        self._pvalue = pvalue
```

## 🧐 Real world examples

- [main.nf](../../main.nf)
- [workflows/ge2ts.nf](../../workflows/ge2ts.nf)
- [modules/local/gwas_prep.nf](../../modules/local/gwas_prep.nf)
- [modules/local/magma_annot_run.nf](../../modules/local/magma_annot_run.nf)
- [modules/local/qc_gsea.nf](../../modules/local/qc_gsea.nf)

## 🔗 Related agreements

- [AGENTS.md](../../AGENTS.md)

Doc created by 🐢 💨 (Turbotuga™, [Codely](https://codely.com)’s mascot)
