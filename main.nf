#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    GE2TS: Genetic Effects to Transcriptomic Scores Pipeline
========================================================================================
*/

include { GE2TS } from './workflows/ge2ts'

workflow {
    // 1. Create input channels
    if (!params.input) {
        error "Please specify a GWAS summary statistics file with --input <path>"
    }

    ch_input = Channel.fromPath(params.input, checkIfExists: true)
        .map { file ->
            def meta = [ id: file.getSimpleName() ]
            return [ meta, file ]
        }

    ch_bfile    = params.bfile ? Channel.fromPath(params.bfile, checkIfExists: true) : Channel.empty()
    ch_gene_loc = params.gene_loc ? Channel.fromPath(params.gene_loc, checkIfExists: true) : Channel.empty()

    // 2. Run workflow
    GE2TS(
        ch_input,
        ch_bfile,
        ch_gene_loc
    )
}
