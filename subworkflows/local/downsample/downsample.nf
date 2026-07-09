/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oncotrseq: Downsampling Subworkflow
    - Handles downsampling of reads and mapping to reference using ontime
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ONTIME_RANGE_FILTER_FASTQ              } from '../../../modules/local/ontime/main.nf'                     // ontime downsampling
include { paramsSummaryMap                       } from 'plugin/nf-schema'                                          // Parameter summary
include { softwareVersionsToYAML                 } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'       // Version reporting

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN DOWNSAMPLE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOWNSAMPLE {
    take:
    in_ch

    main:

    // Fork based on downsamp_hours ; FULL inputs are not processed through ontime.
    in_ch
        .branch { meta, biotype, input, zero, downsamp_hours ->
            full: downsamp_hours == 'FULL'
            downsamp: downsamp_hours != 'FULL'
        }
        .set { ch_branched }

    ch_ontime_in = ch_branched.downsamp
        .map { meta, biotype, input, zero, downsamp_hours ->
            tuple(meta + [biotype: biotype], input, zero, downsamp_hours)
        }

    ch_full_out = ch_branched.full
        .map { meta, biotype, input, zero, downsamp_hours ->
            tuple(meta, biotype, input)
        }

    // Downsample reads using ontime
    ONTIME_RANGE_FILTER_FASTQ(ch_ontime_in)

    //Add biotype back to the downsamp channel
    downsampled_ch = ONTIME_RANGE_FILTER_FASTQ.out.fastq
        .map { meta, input ->
            def biotype = meta.biotype
            def clean_meta = meta.clone()
            clean_meta.remove('biotype')

            tuple(clean_meta, biotype, input)
        }


    //Add the downsampled reads to the full channel
    out_channel = downsampled_ch
        .mix(ch_full_out)

    // Collect versions from all modules
    ch_versions = ONTIME_RANGE_FILTER_FASTQ.out.versions


    emit:
    reads         = out_channel                       // Final downsampled reads channel
    versions      = ch_versions                       // All tool versions

}