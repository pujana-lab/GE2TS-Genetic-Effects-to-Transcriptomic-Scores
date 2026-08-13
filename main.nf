#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    GE2TS: Genetic Effects to Transcriptomic Scores Pipeline
========================================================================================
*/

include { GE2TS } from './workflows/ge2ts'

def helpMessage() {
    log.info"""
    ========================================================================================
    GE2TS: Genetic Effects to Transcriptomic Scores Pipeline v${workflow.manifest.version}
    ========================================================================================
    Usage:
      nextflow run main.nf --input <gwas_file> --outdir <output_dir> [options]

    Mandatory arguments:
      --input [file]         Path to GWAS summary statistics file (.tsv or .tsv.gz)
      --outdir [directory]   Path to output directory

    Optional arguments:
      --bfile [path]         PLINK binary reference dataset prefix
      --gene-loc [file]      Gene location file
      --pathways [file]      Pathway database file (.gmt)
      --help                 Display this help message
    ========================================================================================
    """.stripIndent()
}

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    log.info """
    ========================================================================================
    GE2TS Pipeline Execution
    ========================================================================================
    GWAS Input  : ${params.input}
    Reference   : ${params.bfile}
    Gene Location: ${params.gene_loc}
    Pathways GMT: ${params.pathways}
    Output Dir  : ${params.outdir}
    ========================================================================================
    """.stripIndent()

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
    ch_pathways = params.pathways ? Channel.fromPath(params.pathways, checkIfExists: true) : Channel.empty()

    GE2TS(
        ch_input,
        ch_bfile,
        ch_gene_loc,
        ch_pathways
    )
}
