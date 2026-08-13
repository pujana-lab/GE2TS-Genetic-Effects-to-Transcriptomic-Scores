include { GWAS_PREP } from '../modules/local/gwas_prep'

workflow GE2TS {
    take:
    ch_gwas // channel: [ val(meta), path(gwas) ]

    main:
    GWAS_PREP(ch_gwas)

    emit:
    gwas_prep  = GWAS_PREP.out.gwas_prep
    qc_summary = GWAS_PREP.out.qc_summary
}
