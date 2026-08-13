process GWAS_PREP {
    tag "$meta.id"
    label 'process_low'

    container 'community.wave.seqera.io/library/polars_pandas:latest'

    input:
    tuple val(meta), path(gwas)

    output:
    tuple val(meta), path("*.prep.tsv"), emit: gwas_prep
    tuple val(meta), path("*.qc.json")  , emit: qc_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gwas_prep.py \\
        --input ${gwas} \\
        --output-tsv ${prefix}.prep.tsv \\
        --qc-summary ${prefix}.qc.json \\
        ${args}
    """
}
