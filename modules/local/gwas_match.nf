process GWAS_MATCH {
    tag "$meta.id"
    label 'process_medium'

    container 'community.wave.seqera.io/library/polars_pandas:latest'

    input:
    tuple val(meta), path(gwas)
    path bfile_files

    output:
    tuple val(meta), path("*.matched.tsv"), emit: matched_gwas

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bim_file = bfile_files instanceof List ? bfile_files.find { it.name.endsWith('.bim') } : bfile_files
    """
    gwas_match.py \\
        --gwas ${gwas} \\
        --bim ${bim_file} \\
        --output-tsv ${prefix}.matched.tsv
    """
}
