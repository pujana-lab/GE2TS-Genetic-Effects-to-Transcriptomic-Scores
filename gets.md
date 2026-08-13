# 🎯 GE2TS: Genetic Effects to Transcriptomic Scores Pipeline

## 1. Objectiu
Pipeline modular en Nextflow per transformar GWAS en puntuacions d'efectes genètics amb anàlisi de vies (GSEA).

## 2. Estructura de Mòduls
- `GWAS_Prep`: Standardize, QC, Indel Filtering, Match.
- `MAGMA_Annot_Run`: 10kb, 100kb, Hi-C Annotations & MAGMA Analysis.
- `QC_GSEA`: GSEA Pathway enrichment & Final QC Matrix.

## 3. Flux de Treball (Pipeline)
1. **INPUT**: GWAS brut + Referències.
2. **PROCESS**: Execució modular.
3. **OUTPUT**:
   - `[phenotype].magma.txt` (Ready)
   - `[phenotype]_[mode].genes.out` (Analysis)
   - `[phenotype]_[mode]_GSEA/` (Pathways)
   - `magma_mapping_qc_summary.tsv` (Resum)

## 4. Configuració
- Adaptat per ser asèptic a les dades (Data-agnostic).
- Ús de containers (Apptainer/Docker) o Conda.
- Documentació completa per a la beca inclosa.
