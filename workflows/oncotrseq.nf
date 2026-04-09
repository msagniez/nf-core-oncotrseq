/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Basecalling subworkflows
include { BASECALL_SIMPLEX   } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'

// Core analysis subworkflows
include { MAPPING            } from '../subworkflows/local/mapping/mapping'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONCOTRSEQ {

    take:
    samplesheet             // channel: samplesheet read in from --input
    skip_mapping             // channel: skip_mapping read in from --skip_mapping
    skip_basecalling        // channel: skip_basecalling read in from --skip_basecalling
    demux                   // channel: demux read in from --demux
    demux_samplesheet       // channel: demux_samplesheet read in from --demux_samplesheet
    ref_transcriptome       // channel: ref read in from --ref

    main:

    ch_versions = Channel.empty()
    ch_mapping = Channel.empty()

     if (params.skip_mapping) {

        ch_mapping = ch_samplesheet.map { meta, fastqFiles, gref, tref -> 
        tuple(meta, fastqFiles, tref)
        }
        MAPPING (
            ch_mapping
        )

        ch_seqkit = MAPPING.out.seqkit

    } else if (params.skip_basecalling) {

        ch_mapping = ch_samplesheet.map { meta, fastqFiles, gref, tref -> 
        tuple(meta, fastqFiles, tref)
        }
        MAPPING (
            ch_mapping
        )

        ch_seqkit = MAPPING.out.seqkit
    } else {

        if (params.demux) {

            BASECALL_MULTIPLEX (
                samplesheet,
                demux_samplesheet
            )

            MAPPING (
                BASECALL_MULTIPLEX.out.fastq,
                ref_transcriptome
            )

            ch_fastq = BASECALL_MULTIPLEX.out.fastq

            ch_seqkit   = BASECALL_MULTIPLEX.out.stats_pass
            ch_versions = BASECALL_MULTIPLEX.out.versions

        } else {

            BASECALL_SIMPLEX (
                samplesheet
            )

            MAPPING (
                BASECALL_SIMPLEX.out.fastq,
                ref_transcriptome
            )

            ch_fastq = BASECALL_SIMPLEX.out.fastq

            ch_seqkit   = BASECALL_SIMPLEX.out.stats_pass
            ch_versions = BASECALL_SIMPLEX.out.versions
        }


    emit:
    multiqc_report = Channel.empty()
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
