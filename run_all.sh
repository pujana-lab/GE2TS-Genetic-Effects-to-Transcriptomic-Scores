#!/bin/bash
# GE2TS Batch Execution Script for CIMBA (BRCA1, BRCA2, TNBC)
# Place and run this script from /home/luis/Escritorio/gwas_to_gsea/

GE2TS_HOME="/home/luis/Documentos/projects/tools/GE2TS-Genetic-Effects-to-Transcriptomic-Scores"
WORK_DIR="/home/luis/Escritorio/gwas_to_gsea"

PATHWAYS="${GE2TS_HOME}/assets/gmt_data/h.all.v2026.1.Hs.symbols.gmt"
BFILE="${WORK_DIR}/reference/g1000_eur"
GENE_LOC="${WORK_DIR}/reference/NCBI37.3/NCBI37.3.gene.loc"
HIC_BED="${WORK_DIR}/reference/normal_breast_consensus_hic_SNP_mappable_regions_hg19.bed"
MAGMA_BIN="/home/luis/Documentos/projects/tools/gwas-magma-pipeline/bin/magma_v1.10/magma"

echo "=== Starting GE2TS Batch Processing (9 Runs with Real MAGMA) ==="

# 1. BRCA1
for WINDOW in "10,10" "100,100"; do
    echo "Running BRCA1 window ${WINDOW}..."
    nextflow run "${GE2TS_HOME}/main.nf" \
      --input "${WORK_DIR}/input_data/cimba_onco_icogs_brca1_combined_results.txt" \
      --bfile "${BFILE}" \
      --gene_loc "${GENE_LOC}" \
      --magma_bin "${MAGMA_BIN}" \
      --pathways "${PATHWAYS}" \
      --window "${WINDOW}" \
      --outdir "${WORK_DIR}/results_BRCA1_${WINDOW//,/}" \
      -resume
done

echo "Running BRCA1 with 10kb + Hi-C..."
nextflow run "${GE2TS_HOME}/main.nf" \
  --input "${WORK_DIR}/input_data/cimba_onco_icogs_brca1_combined_results.txt" \
  --bfile "${BFILE}" \
  --gene_loc "${GENE_LOC}" \
  --magma_bin "${MAGMA_BIN}" \
  --hic_bed "${HIC_BED}" \
  --pathways "${PATHWAYS}" \
  --window "10,10" \
  --outdir "${WORK_DIR}/results_BRCA1_10kb_hic" \
  -resume


# 2. BRCA2
for WINDOW in "10,10" "100,100"; do
    echo "Running BRCA2 window ${WINDOW}..."
    nextflow run "${GE2TS_HOME}/main.nf" \
      --input "${WORK_DIR}/input_data/cimba_onco_icogs_brca2_combined_results.txt" \
      --bfile "${BFILE}" \
      --gene_loc "${GENE_LOC}" \
      --magma_bin "${MAGMA_BIN}" \
      --pathways "${PATHWAYS}" \
      --window "${WINDOW}" \
      --outdir "${WORK_DIR}/results_BRCA2_${WINDOW//,/}" \
      -resume
done

echo "Running BRCA2 with 10kb + Hi-C..."
nextflow run "${GE2TS_HOME}/main.nf" \
  --input "${WORK_DIR}/input_data/cimba_onco_icogs_brca2_combined_results.txt" \
  --bfile "${BFILE}" \
  --gene_loc "${GENE_LOC}" \
  --magma_bin "${MAGMA_BIN}" \
  --hic_bed "${HIC_BED}" \
  --pathways "${PATHWAYS}" \
  --window "10,10" \
  --outdir "${WORK_DIR}/results_BRCA2_10kb_hic" \
  -resume


# 3. TNBC META
for WINDOW in "10,10" "100,100"; do
    echo "Running TNBC window ${WINDOW}..."
    nextflow run "${GE2TS_HOME}/main.nf" \
      --input "${WORK_DIR}/input_data/CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt" \
      --bfile "${BFILE}" \
      --gene_loc "${GENE_LOC}" \
      --magma_bin "${MAGMA_BIN}" \
      --pathways "${PATHWAYS}" \
      --window "${WINDOW}" \
      --outdir "${WORK_DIR}/results_TNBC_${WINDOW//,/}" \
      -resume
done

echo "Running TNBC with 10kb + Hi-C..."
nextflow run "${GE2TS_HOME}/main.nf" \
  --input "${WORK_DIR}/input_data/CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt" \
  --bfile "${BFILE}" \
  --gene_loc "${GENE_LOC}" \
  --magma_bin "${MAGMA_BIN}" \
  --hic_bed "${HIC_BED}" \
  --pathways "${PATHWAYS}" \
  --window "10,10" \
  --outdir "${WORK_DIR}/results_TNBC_10kb_hic" \
  -resume

echo "=== All 9 GE2TS Pipelines Completed Successfully ==="
