/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oncotrseq: Quantifying Subworkflow
    - Handles quantification of reads to transcriptomic reference and quantification matrix handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { OARFISH_QUANTIFY_GENOMEALIGNMENTPROJECTION  } from '../../../modules/local/oarfish/main.nf'                    // oarfish quantification from genomic alignments
include { QUANT_CLEANTRANSCRIPTIDS                    } from '../../../modules/local/utils/quant_cleantranscriptids.nf'  // Clean transcript IDs in quantification matrix
include { QUARTO_TABLE                                } from '../../../modules/local/quarto/main.nf'                     // Reporting (optional)
include { paramsSummaryMap                            } from 'plugin/nf-schema'                                          // Parameter summary
include { paramsSummaryMultiqc                        } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'       // MultiQC summary
include { softwareVersionsToYAML                      } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'       // Version reporting
include { methodsDescriptionText                      } from '../../../subworkflows/local/utils_nfcore_oncotrseq_pipeline' // Methods for MultiQC

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN QUANTIFYING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QUANTIFYING {
    take:
    in_ch

    main:
    ch_versions = Channel.empty()

    OARFISH_QUANTIFY_GENOMEALIGNMENTPROJECTION(in_ch)

    QUANT_CLEANTRANSCRIPTIDS(OARFISH_QUANTIFY_GENOMEALIGNMENTPROJECTION.out.quant)

    ch_quant = QUANT_CLEANTRANSCRIPTIDS.out.quant
        .map { meta, quant -> [ meta + [ quantifier: 'oarfish' ], quant ] }


    // Collect versions from all modules
    ch_versions = OARFISH_QUANTIFY_GENOMEALIGNMENTPROJECTION.out.versions
                    .mix(QUANT_CLEANTRANSCRIPTIDS.out.versions)


    emit:
    quant          = ch_quant
    versions       = ch_versions                        // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/