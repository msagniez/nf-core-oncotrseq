/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DORADO_BASECALL                           } from '../../../modules/local/dorado/main.nf'
include { DORADO_DEMULTIPLEX                        } from '../../../modules/local/dorado/main.nf'
include { SAMTOOLS_QSFILTER                         } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_PASS } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_FAIL } from '../../../modules/local/samtools/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BASECALL_MULTIPLEX {

    //TODO Add reports for read stats figure and tables

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    ch_samplesheet_to_basecall = ch_samplesheet
        .map { meta, _barcode, model, modif, input, _qc_type, _species, _kit, ubam ->
            tuple(id:meta.project,input,ubam,model,modif) }             // Meta project as meta id to basecall only once per sequencing project
        .unique()

    ch_demux_corrected_barcode = ch_samplesheet
        .map { meta, barcode, _model, _modif, _input, _qc_type, _species, kit, _ubam ->
            def new_barcode = kit + '_' + barcode
            tuple(id:meta.project, meta.id, new_barcode) }

    DORADO_BASECALL(ch_samplesheet_to_basecall)

    ch_to_dmux = ch_samplesheet
        .map { meta, _barcode, _model, _modif, _input, _qc_type, _species, kit, _ubam ->
            tuple(id:meta.project,kit) }
        .join(DORADO_BASECALL.out.ubam)

    DORADO_DEMULTIPLEX(ch_to_dmux)

     split_bams_ch = DORADO_DEMULTIPLEX.out.demux_ubam
        .flatMap { meta, bams ->
            bams.collect { bam -> tuple(meta, bam) }                    // Split output of Dorado multiplex in a multiple tuples of one bam per barcode
        }
        .map { meta, bam ->
            def barcode = bam.baseName.replaceAll(/^[^_]*_/, '')
            def new_meta = [ id: "${meta.id}_${barcode}" ]
            tuple(new_meta, bam)
        }

    SAMTOOLS_QSFILTER(split_bams_ch)

     // Get the sample_ids from the demux_samplesheet that specify each barcode to each sample_id
    ch_new_sample_ids_pass = ch_demux_corrected_barcode
        .map { meta, sample_id, barcode ->
            def new_meta = meta.id + '_' + barcode
            tuple(new_meta, sample_id)    // flatten key to string
        }
        .combine(
            SAMTOOLS_QSFILTER.out.ubam_pass.map { meta, bam ->
                tuple(meta.id, bam)      // also flatten key
            },
            by: 0
        )
        .map { _meta, sampleid, ubam ->                         // Get rid of barcodes here and use real sample_id as meta
            tuple(id: sampleid, ubam) }
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')       // Add pass to meta in tuples for output naming
            def new_meta = meta.id + '_' + meta_suffix
            tuple(id:new_meta, ubam)
        }

    ch_new_sample_ids_fail = ch_demux_corrected_barcode
        .map { meta, sample_id, barcode ->
            def new_meta = meta.id + '_' + barcode
            tuple(new_meta, sample_id)    // flatten key to string
        }
        .combine(
            SAMTOOLS_QSFILTER.out.ubam_fail.map { meta, bam ->
                tuple(meta.id, bam)      // also flatten key
            },
            by: 0
        )
        .map { _meta, sampleid, ubam ->                         // Get rid of barcodes here and use real sample_id as meta
            tuple(id: sampleid, ubam) }
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')       // Add fail to meta in tuples for output naming
            def new_meta = meta.id + '_' + meta_suffix
            tuple(id:new_meta, ubam)
            }

    SAMTOOLS_TOFASTQ_PASS(ch_new_sample_ids_pass)
    SAMTOOLS_TOFASTQ_FAIL(ch_new_sample_ids_fail)

    ch_fastq_pass = SAMTOOLS_TOFASTQ_PASS.out.fq
        .map { meta, fq ->
         def meta_restore = meta.id.replace('_pass', '')
         tuple(id:meta_restore, fq)}

    // Add project in meta of resulting fastq files :

    ch_samplesheet_to_rejoin = ch_samplesheet
        .map { meta, _barcode, model, _modif, _input, qc_type, species, _kit, _ubam ->
        tuple(meta, model, qc_type, species)}

    ch_fastq = ch_samplesheet_to_rejoin
        .map { meta, model, qc_type, species ->
            tuple(id:meta.id, meta.project, model, qc_type, species) }
        .join(ch_fastq_pass, by:0)
        .map { meta_sample, project, model, qc_type, species, fq ->
            tuple(project:project, id:meta_sample.id, fq, model, qc_type, species) }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_versions = DORADO_BASECALL.out.versions
        .mix(DORADO_DEMULTIPLEX.out.versions)
        .mix(SAMTOOLS_QSFILTER.out.versions)
        .mix(SAMTOOLS_TOFASTQ_PASS.out.versions)
        .mix(SAMTOOLS_TOFASTQ_FAIL.out.versions)


    emit:
    fastq          = ch_fastq
    versions       = ch_versions              // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/