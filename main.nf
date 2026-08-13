#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    GE2TS: Genetic Effects to Transcriptomic Scores Pipeline
========================================================================================
*/

include { GE2TS } from './workflows/ge2ts'

workflow {
    // 1. Create input channel
    if (!params.input) {
        error "Please specify a GWAS summary statistics file with --input <path>"
    }

    ch_input = Channel.fromPath(params.input, checkIfExists: true)
        .map { file ->
            def meta = [ id: file.getSimpleName() ]
            return [ meta, file ]
        }

    // 2. Run workflow
    GE2TS(ch_input)
}
