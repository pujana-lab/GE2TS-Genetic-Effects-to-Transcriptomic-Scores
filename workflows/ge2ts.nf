include { GWAS_PREP       } from '../modules/local/gwas_prep'
include { GWAS_MATCH      } from '../modules/local/gwas_match'
include { MAGMA_ANNOT_RUN } from '../modules/local/magma_annot_run'
include { QC_GSEA         } from '../modules/local/qc_gsea'

workflow GE2TS {
    take:
    ch_gwas         // channel: [ val(meta), path(gwas) ]
    ch_prep         // channel: [ val(meta), path(prep_gwas) ]
    ch_genes_out    // channel: [ val(meta), path(genes_out) ]
    ch_bfile        // channel: path(bfile)
    ch_gene_loc     // channel: path(gene_loc)
    ch_hic_bed      // channel: path(hic_bed)
    ch_pathways     // channel: path(pathway_gmt)

    main:
    // 1. GWAS Prep (Phase 1)
    if (!params.input_prep && !params.input_genes_out) {
        GWAS_PREP(ch_gwas)
        ch_current_prep = GWAS_PREP.out.gwas_prep
        ch_qc_summary   = GWAS_PREP.out.qc_summary
    } else {
        ch_current_prep = ch_prep
        ch_qc_summary   = Channel.empty()
    }

    // 1.5. GWAS Match against Reference BIM IDs
    if (!params.input_genes_out && !params.stop_at_prep && !params.skip_magma) {
        GWAS_MATCH(
            ch_current_prep,
            ch_bfile
        )
        ch_matched_gwas = GWAS_MATCH.out.matched_gwas
    } else {
        ch_matched_gwas = ch_current_prep
    }

    // 2. MAGMA Annotation & Analysis (Phase 2)
    if (params.stop_at_prep || params.skip_magma) {
        ch_current_genes = Channel.empty()
        ch_genes_raw     = Channel.empty()
    } else {
        ch_magma_input = params.input_genes_out ? ch_genes_out : ch_matched_gwas
        if (!params.input_genes_out) {
            MAGMA_ANNOT_RUN(
                ch_magma_input,
                ch_bfile,
                ch_gene_loc,
                ch_hic_bed
            )
            ch_current_genes = MAGMA_ANNOT_RUN.out.genes_out
            ch_genes_raw     = MAGMA_ANNOT_RUN.out.genes_raw
        } else {
            ch_current_genes = ch_genes_out
            ch_genes_raw     = Channel.empty()
        }
    }

    // 3. GSEA Pathway Enrichment & QC (Phase 3)
    if (!params.stop_at_prep && !params.skip_magma && !params.skip_gsea) {
        QC_GSEA(
            ch_current_genes,
            ch_pathways,
            ch_gene_loc
        )
        ch_gsea_results = QC_GSEA.out.gsea_results
        ch_magma_qc     = QC_GSEA.out.qc_summary
    } else {
        ch_gsea_results = Channel.empty()
        ch_magma_qc     = Channel.empty()
    }

    emit:
    gwas_prep    = ch_current_prep
    genes_out    = ch_current_genes
    gsea_results = ch_gsea_results
}
