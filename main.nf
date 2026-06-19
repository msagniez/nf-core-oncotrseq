#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oncotrseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/oncotrseq
    Website: https://nf-co.re/oncotrseq
    Slack  : https://nfcore.slack.com/channels/oncotrseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oncotrseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oncotrseq_pipeline'
include { ONCOTRSEQ               } from './workflows/oncotrseq'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_ONCOTRSEQ {

    take:
    samplesheet_fastq  // channel: fastq samples - needs alignment/quantification before classification
    samplesheet_matrix // channel: already-quantified matrix samples - go straight to classification
    genome_fasta        // channel: resolved (and indexed) reference genome - fastq samples only
    ref_classif          // channel: resolved reference classification per sample
    biotype               // channel: biotype per sample

    main:

    //
    // WORKFLOW: Run pipeline
    //
    ONCOTRSEQ (
        samplesheet_fastq,
        samplesheet_matrix,
        genome_fasta,
        ref_classif,
        biotype
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    ch_samplesheet_fastq  = PIPELINE_INITIALISATION.out.fastq
    ch_samplesheet_matrix = PIPELINE_INITIALISATION.out.matrix
    ch_genome_fasta        = PIPELINE_INITIALISATION.out.genome_fasta

    ch_samplesheet_fastq.view  { "Fastq samplesheet channel:  ${it}" }
    ch_samplesheet_matrix.view { "Matrix samplesheet channel: ${it}" }
    ch_genome_fasta.view       { "Genome fasta channel: ${it}" }
    PIPELINE_INITIALISATION.out.ref_classif.view { "Reference classification channel: ${it}" }
    PIPELINE_INITIALISATION.out.biotype.view { "Biotype channel: ${it}" }

    ONCOTRSEQ (
        ch_samplesheet_fastq,
        ch_samplesheet_matrix,
        ch_genome_fasta,
        PIPELINE_INITIALISATION.out.ref_classif,
        PIPELINE_INITIALISATION.out.biotype
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/