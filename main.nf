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
      --outdir [directory]   Path to output directory

    Input options (choose one):
      --input [file]         Path to raw GWAS summary statistics file (.tsv / .csv / .gz)
      --input_prep [file]    Path to pre-prepared GWAS file (skips GWAS_PREP step)
      --input_genes_out [file] Path to pre-computed MAGMA genes output (skips steps 1 & 2)

    Reference options:
      --bfile [path]         PLINK binary reference dataset prefix (.bed/.bim/.fam)
      --gene-loc [file]      Gene location file
      --hic_bed [file]       Optional Hi-C mappable regions BED file
      --pathways [file]      Pathway database file (.gmt)

    Execution control flags:
      --stop_at_prep         Stop pipeline after GWAS preparation step
      --skip_magma           Alias for stop_at_prep
      --skip_gsea            Stop pipeline after MAGMA gene analysis step
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
    GE2TS Pipeline Execution v${workflow.manifest.version}
    ========================================================================================
    Input GWAS  : ${params.input ?: params.input_prep ?: params.input_genes_out}
    Reference   : ${params.bfile}
    Gene Loc    : ${params.gene_loc}
    Hi-C BED    : ${params.hic_bed ?: 'None'}
    Pathways GMT: ${params.pathways}
    Output Dir  : ${params.outdir}
    ========================================================================================
    """.stripIndent()

    if (!params.input && !params.input_prep && !params.input_genes_out) {
        error "Please specify an input file with --input <path>, --input_prep <path>, or --input_genes_out <path>"
    }

    ch_gwas         = params.input ? Channel.fromPath(params.input, checkIfExists: true).map { file -> [ [id: file.getSimpleName()], file ] } : Channel.empty()
    ch_prep         = params.input_prep ? Channel.fromPath(params.input_prep, checkIfExists: true).map { file -> [ [id: file.getSimpleName().replaceAll(/\.prep$/, '')], file ] } : Channel.empty()
    ch_genes_out    = params.input_genes_out ? Channel.fromPath(params.input_genes_out, checkIfExists: true).map { file -> [ [id: file.getSimpleName().replaceAll(/\.genes$/, '')], file ] } : Channel.empty()

    ch_bfile        = params.bfile ? Channel.fromPath(params.bfile, checkIfExists: true) : Channel.empty()
    ch_gene_loc     = params.gene_loc ? Channel.fromPath(params.gene_loc, checkIfExists: true) : Channel.empty()
    ch_hic_bed      = params.hic_bed ? Channel.fromPath(params.hic_bed, checkIfExists: true) : Channel.value(file('NO_FILE'))
    ch_pathways     = params.pathways ? Channel.fromPath(params.pathways, checkIfExists: true) : Channel.empty()

    GE2TS(
        ch_gwas,
        ch_prep,
        ch_genes_out,
        ch_bfile,
        ch_gene_loc,
        ch_hic_bed,
        ch_pathways
    )
}
