/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oncotrseq: Classifying Subworkflow
    - Handles classification for ALL/AML/Other tumor types from transcript-level quantification matrices as well as a cascade mode based on MnM results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GENE_SUMMARIZE as GENE_SUMMARIZE_MNM        } from '../../../modules/local/tximport/main.nf'                   // summarize transcript-level quantification to gene-level
include { REMOVE_DUPLICATE_GENES                     } from '../../../modules/local/utils/remove_duplicate_genes.nf'     // Remove duplicate genes in quantification matrix
include { MNM_PREDICT                                 } from '../../../modules/local/MnM/main.nf'                        // MnM classification
//include { ALLCATCHR_PREDICT                           } from '../../../modules/local/allcatchr/main.nf'                  // ALLcatchR classification
//include { MDALL_PREDICT                               } from '../../../modules/local/mdall/main.nf'                        // MD-ALL classification
include { QUARTO_TABLE                                } from '../../../modules/local/quarto/main.nf'                     // Reporting (optional)
include { paramsSummaryMap                            } from 'plugin/nf-schema'                                          // Parameter summary
include { paramsSummaryMultiqc                        } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'       // MultiQC summary
include { softwareVersionsToYAML                      } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'       // Version reporting
include { methodsDescriptionText                      } from '../../../subworkflows/local/utils_nfcore_oncotrseq_pipeline' // Methods for MultiQC

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN CLASSIFYING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLASSIFYING {
    take:
    classif_in

    main:
    ch_versions = Channel.empty()
    classif_out = Channel.empty()

    // Fork based on tumor_type to run different classification strategies.
    ch_classif_forked = classif_in
        .branch { meta, biotype, tumor_type, quant, ref_rna, gtf, quantifier ->
            all: tumor_type == 'all'
            aml: tumor_type == 'aml'
            cascade: tumor_type == 'cascade'
            other: true
        }
        .set { ch_classif }

    

    // MnM classification for any tumor_type (other, all, aml, cascade)
    //Collapse to gene-level with tximport --gene-names ; use TPM file
    GENE_SUMMARIZE_MNM(
        classif_in.map { meta, biotype, tumor_type, quant, ref_rna, gtf, quantifier -> tuple(meta, quant) },
        classif_in.map { meta, biotype, tumor_type, quant, ref_rna, gtf, quantifier -> gtf },
        true, //use --gene-names option
        false, //do not use --ignore-tx-version option
        classif_in.map { meta, biotype, tumor_type, quant, ref_rna, gtf, quantifier -> quantifier }
    )

    REMOVE_DUPLICATE_GENES(GENE_SUMMARIZE_MNM.out.tpm)

    //Apply MnM
    MNM_PREDICT(REMOVE_DUPLICATE_GENES.out.quant)

    // Collect type and subtype + scores for cascade branch
    

    
    // ALL classification for all tumor_type only
    // Collapse to gene-level with --ignore-tx-version ; use COUNTS file (SIGNATURE - ALLcatchR - MD-ALL)

    //Apply ALLcatchR
    //Apply MD-ALL
    //Apply SIGNATURE

    // Collect type and subtype + scores for cascade branch



    // AML classification for aml tumor_type only
    // Collapse to gene-level ; use TPM file (Maylis) ; use COUNTS file (AMLmapR)
    // Collapse to gene-level with --gene-names ; use TPM file (AttentionAML)

    //Apply Maylis (reformat + apply + rearrange results)
    //Apply AMLmapR (Transpose + reformat + apply)
    //Apply AttentionAML (Transpose + apply)

    // Collect type and subtype + scores for cascade branch

    // Collect versions from all modules
    ch_versions = ch_versions.mix(GENE_SUMMARIZE_MNM.out.versions.first())
                    .mix(REMOVE_DUPLICATE_GENES.out.versions.first())


    emit:
    classif        = classif_out
    versions       = ch_versions                        // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

