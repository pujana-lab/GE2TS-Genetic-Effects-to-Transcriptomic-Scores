include { GWAS_PREP       } from '../modules/local/gwas_prep'
include { MAGMA_ANNOT_RUN } from '../modules/local/magma_annot_run'

workflow GE2TS {
    take:
    ch_gwas     // channel: [ val(meta), path(gwas) ]
    ch_bfile    // channel: path(bfile)
    ch_gene_loc // channel: path(gene_loc)

    main:
    GWAS_PREP(ch_gwas)

    MAGMA_ANNOT_RUN(
        GWAS_PREP.out.gwas_prep,
        ch_bfile,
        ch_gene_loc
    )

    emit:
    gwas_prep  = GWAS_PREP.out.gwas_prep
    qc_summary = GWAS_PREP.out.qc_summary
    genes_out  = MAGMA_ANNOT_RUN.out.genes_out
    genes_raw  = MAGMA_ANNOT_RUN.out.genes_raw
}
