process MAGMA_ANNOT_RUN {
    tag "$meta.id"
    label 'process_medium'

    container 'community.wave.seqera.io/library/magma_pandas:latest'

    input:
    tuple val(meta), path(prep_gwas)
    path bfile
    path gene_loc
    path hic_bed, stageAs: 'hic_regions.bed'

    output:
    tuple val(meta), path("*.genes.out"), emit: genes_out
    tuple val(meta), path("*.genes.raw"), emit: genes_raw

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def window = task.ext.window ?: '10,10'
    def bfile_prefix = bfile instanceof List ? bfile[0].baseName : file(bfile).baseName
    def hic_arg = hic_bed && hic_bed.name != 'NO_FILE' ? "--hic-bed hic_regions.bed" : ''
    """
    magma_wrapper.py \\
        --gwas ${prep_gwas} \\
        --bfile ${bfile_prefix} \\
        --gene-loc ${gene_loc} \\
        --window ${window} \\
        --out-prefix ${prefix} \\
        ${hic_arg} \\
        ${args}
    """
}
