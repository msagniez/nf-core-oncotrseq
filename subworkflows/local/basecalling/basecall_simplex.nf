/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DORADO_BASECALL                           } from '../../../modules/local/dorado/main.nf'
include { SAMTOOLS_QSFILTER                         } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_PASS } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_FAIL } from '../../../modules/local/samtools/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BASECALL_SIMPLEX {

    //TODO Add reports for read stats figure and tables

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    ch_input = ch_samplesheet
        .map { meta, model, modif, input, _qc_type, ubam ->
             tuple(meta, input, ubam, model, modif)}

    DORADO_BASECALL(ch_input)

    SAMTOOLS_QSFILTER(DORADO_BASECALL.out.ubam)

    // Add pass and fail to meta in tuples for output naming

    ch_ubam_pass = SAMTOOLS_QSFILTER.out.ubam_pass
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')
            def new_meta = meta.id + '_' + meta_suffix
            tuple(project:meta.project, id:new_meta, ubam)
            }

    ch_ubam_fail = SAMTOOLS_QSFILTER.out.ubam_fail
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')
            def new_meta = meta.id + '_' + meta_suffix
            tuple(project:meta.project, id:new_meta, ubam)
            }

    SAMTOOLS_TOFASTQ_PASS(ch_ubam_pass)
    SAMTOOLS_TOFASTQ_FAIL(ch_ubam_fail)

    ch_fastq = SAMTOOLS_TOFASTQ_PASS.out.fq
        .map { meta, fastq ->
            def meta_restore = meta.id.replace('_pass', '')
            tuple(project:meta.project, id:meta_restore, fastq) }
        .join(ch_samplesheet)
        .map { meta, fastq, model, _modif, _input, qc_type, _ubam ->
        tuple(meta, fastq, model, qc_type)}



    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_versions = DORADO_BASECALL.out.versions
        .mix(SAMTOOLS_QSFILTER.out.versions)
        .mix(SAMTOOLS_TOFASTQ_PASS.out.versions)
        .mix(SAMTOOLS_TOFASTQ_FAIL.out.versions)



    emit:
    fastq          = ch_fastq
    versions       = ch_versions

}