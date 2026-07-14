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
include { OARFISH_QUANTIFY_GENOMEALIGNMENTPROJECTION  } from '../../../modules/local/oarfish/main.nf'                   // oarfish quantification from genomic alignments
include { modifyMetaId                                } from '../utils_nfcore_oncotrseq_pipeline'
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


    // Collect versions from all modules
    ch_versions = OARFISH_QUANTIFY_GENOMEALIGNMENTPROJECTION.out.versions


    emit:
    quant          = OARFISH_QUANTIFY_GENOMEALIGNMENTPROJECTION.out.quant
    versions       = ch_versions                        // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/