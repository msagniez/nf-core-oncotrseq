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

include { ONCOTRSEQ  } from './workflows/oncotrseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oncotrseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oncotrseq_pipeline'
include { DORADO_DOWNLOAD_LIST    } from './modules/local/dorado/main.nf'
include { DORADO_DOWNLOAD_MODEL   } from './modules/local/dorado/main.nf'
include { selectLatestModel       } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { selectLatestModif       } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { selectModelDownload     } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'


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
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    tumor_type  // channel: samplesheet read in from --input, contains only tumor type
    basecall_model  // channel : basecalling model used with dorado

    main:

    //
    // WORKFLOW: Run pipeline
    //
    ONCOTRSEQ (
        samplesheet,
        demux,
        ref,
        basecall_model,
        tumor_type
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    if (!params.reference_cache_dir) {
        throw new IllegalArgumentException(
            'Please provide --reference_cache_dir to stage reference assets ' +
                'and Dorado models.'
        )
    }

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.ubam_samplesheet,
        params.demux_samplesheet,
    )

    // Load model channels from parameters:
    if (params.skip_basecalling || params.skip_mapping ) {

        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .map { meta, input, _ubam ->
            tuple(meta, input) }

        ch_model = channel.of(params.basecall_model)

    } else {

        ch_model_dir = channel.fromPath("${params.reference_cache_dir}")
        def model_resolved = selectLatestModel(params.basecall_model,file("${params.reference_cache_dir}/dorado_models"))
        def modif_resolved = params.m_bases && model_resolved ?
            selectLatestModif(file("${params.reference_cache_dir}/dorado_models"), model_resolved, params.m_bases) :
            null
        ch_modif = channel.fromPath("${projectDir}/assets/NOMOD")

        if (model_resolved) {
            ch_model = channel.fromPath(model_resolved)
        }

        if (modif_resolved) {
            ch_modif = channel.fromPath(modif_resolved)
        }

        // download both
        if (!model_resolved && !modif_resolved && params.m_bases) {

            DORADO_DOWNLOAD_LIST()

            ch_model_to_download_resolved = selectBaseModelDownload(
                DORADO_DOWNLOAD_LIST.out.list,
                params.basecall_model
            )

            ch_model_to_download = ch_model_to_download_resolved
                .combine(ch_model_dir)
                .map { type, model -> tuple("base", type, model)}

            ch_modif_to_download = selectModifiedModelDownload(
                DORADO_DOWNLOAD_LIST.out.list,
                ch_model_to_download_resolved,
                params.m_bases
            )
                .combine(ch_model_dir)
                .map { type, model -> tuple("modif", type, model)}

            ch_in_download = ch_model_to_download
                .mix(ch_modif_to_download)

            DORADO_DOWNLOAD_MODEL(ch_in_download)

            ch_model = DORADO_DOWNLOAD_MODEL.out.model
                .filter { type, model -> type == "base" }
                .map { type, model -> model }

            ch_modif = DORADO_DOWNLOAD_MODEL.out.model
                .filter { type, model -> type == "modif" }
                .map { type, model -> model }

        // download only model
        } else if (!model_resolved && !params.m_bases) {

            DORADO_DOWNLOAD_LIST()

            ch_model_to_download = selectBaseModelDownload(
                DORADO_DOWNLOAD_LIST.out.list,
                params.basecall_model
            )
            .combine(ch_model_dir)
            .map { type, model -> tuple("base", type, model)}

            ch_in_download = ch_model_to_download

            DORADO_DOWNLOAD_MODEL(ch_model_to_download)

            ch_model = DORADO_DOWNLOAD_MODEL.out.model
                .map { type, model -> model }

        // download only modif
        } else if (model_resolved && !modif_resolved && params.m_bases) {

            DORADO_DOWNLOAD_LIST()

            ch_model_match = ch_model.map { it -> it.name }

            ch_modif_to_download = selectModifiedModelDownload(
                DORADO_DOWNLOAD_LIST.out.list,
                ch_model_match,
                params.m_bases
            )
                .combine(ch_model_dir)
                .map { type, model -> tuple("modif", type, model)}

            DORADO_DOWNLOAD_MODEL(ch_modif_to_download)

            ch_modif = DORADO_DOWNLOAD_MODEL.out.model
                .map { type, model -> model }
        }

        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .combine(ch_model)
            .combine(ch_modif)
    }

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_ONCOTRSEQ (
        ch_input,
        PIPELINE_INITIALISATION.out.demux_sheet,
        PIPELINE_INITIALISATION.out.ref_ch,
        ch_model,
        PIPELINE_INITIALISATION.out.bed_sheet,
        PIPELINE_INITIALISATION.out.tumor_type
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
