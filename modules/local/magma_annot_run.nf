process MAGMA_ANNOT_RUN {
    tag "$meta.id"
    label 'process_medium'

    container 'community.wave.seqera.io/library/magma_pandas:latest'

    input:
    tuple val(meta), path(prep_gwas)
    path bfile
    path gene_loc

    output:
    tuple val(meta), path("*.genes.out"), emit: genes_out
    tuple val(meta), path("*.genes.raw"), emit: genes_raw

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def window = task.ext.window ?: '10,10'
    """
    magma_wrapper.py \\
        --gwas ${prep_gwas} \\
        --bfile ${bfile} \\
        --gene-loc ${gene_loc} \\
        --window ${window} \\
        --out-prefix ${prefix} \\
        ${args}
    """
}
