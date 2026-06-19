/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_oncotrseq_pipeline'
include { MAPPING                   } from '../subworkflows/local/mapping/mapping.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONCOTRSEQ {

    take:
    samplesheet_fastq
    samplesheet_matrix
    genome_fasta
    ref_classif
    biotype

    main:

    ch_versions = Channel.empty()
    ch_mapping = Channel.empty()

    samplesheet_fastq.view { "Fastq samplesheet channel: ${it}" }
    samplesheet_matrix.view { "Matrix samplesheet channel: ${it}" }
    genome_fasta.view { "Genome fasta channel: ${it}" }
    ref_classif.view { "Reference classification channel: ${it}" }
    biotype.view { "Biotype channel: ${it}" }

    //MAPPING (
    //    ch_samplesheet,
    //    ch_ref
    //)


    emit:
    multiqc_report = Channel.empty()
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
