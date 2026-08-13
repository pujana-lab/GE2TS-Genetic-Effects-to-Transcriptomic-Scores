process QC_GSEA {
    tag "$meta.id"
    label 'process_low'

    container 'community.wave.seqera.io/library/pandas_scipy:latest'

    input:
    tuple val(meta), path(genes_out)
    path pathway_gmt

    output:
    tuple val(meta), path("GSEA/*.tsv")             , emit: gsea_results
    tuple val(meta), path("*.magma_qc_summary.tsv"), emit: qc_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    qc_gsea.py \\
        --genes-out ${genes_out} \\
        --pathway-gmt ${pathway_gmt} \\
        --out-dir GSEA \\
        --qc-summary-tsv ${prefix}.magma_qc_summary.tsv \\
        ${args}
    """
}
