/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_oncotrseq_pipeline'
include { modifyMetaId              } from '../subworkflows/local/utils_nfcore_oncotrseq_pipeline'
include { MAPPING                   } from '../subworkflows/local/mapping/mapping.nf'
include { DOWNSAMPLE                } from '../subworkflows/local/downsample/downsample.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Parse a --downsamp value such as "1h" or "1h,2h,72h,FULL" into a list of
// [label, hours] pairs, e.g. [ ["1h", 1], ["2h", 2], ["72h", 72], ["FULL", "FULL"] ]
// The special token 'FULL' (case-insensitive) means: keep the complete,
// unfiltered set of reads for that sample alongside the timed subsets.
def parseDownsampParam(String downsampStr) {
    if (!downsampStr?.trim()) {
        return []
    }

    return downsampStr.split(',').collect { token ->
        token = token.trim()

        if (token.toUpperCase() == 'FULL') {
            return ['FULL', 'FULL']
        }
        
        def matcher = token =~ /^(\d+)h$/
        if (!matcher.matches()) {
            error "Invalid --downsamp value '${token}'. Only integer hours are supported (e.g. '1h,2h,72h,FULL')."
        }

        def hours = matcher[0][1].toInteger()
        [token, hours]
    }
}

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


    ch_samplesheet = samplesheet_fastq.map { meta, biotype, tumor_type, input, ref_genome, ref_rna, gtf ->
        tuple(meta, biotype, input)
    }

    // Optional time-based downsampling: --downsamp 1h  or  --downsamp 1h,2h,72h,FULL
    // If --downsamp is not specified, the original samplesheet is used as-is.
    // Each sample is expanded into one copy per requested timepoint, tagged
    // with a suffix on the meta id (e.g. Sample_1h, Sample_2h, Sample_72h, Sample_FULL).,
    // so downstream results for every timepoint are kept separate.
    if (params.downsamp) {
        ch_downsamp_times = Channel.fromList(parseDownsampParam(params.downsamp))

        ch_downsamp_in = ch_samplesheet
            .combine(ch_downsamp_times)
            .map { meta, biotype, input, downsamp_label, downsamp_hours ->
                def meta_ds = modifyMetaId(meta, 'add_suffix', '', '', "_${downsamp_label}")
                meta_ds.original_id = meta.id


                tuple(meta_ds, biotype, input, 0, downsamp_hours)
            }

        DOWNSAMPLE(ch_downsamp_in)
        


        ch_mapping_input = DOWNSAMPLE.out.reads
    } else {
        ch_mapping_input = ch_samplesheet
    }

    MAPPING (
        ch_mapping_input,
        genome_fasta
    )


    emit:
    multiqc_report = Channel.empty()
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/