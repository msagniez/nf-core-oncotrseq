/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oncotrseq: Mapping Subworkflow
    - Handles mapping of reads to reference using minimap2 and downstream processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MINIMAP2_ALIGN                         } from '../../../modules/local/minimap2/main.nf'                   // minimap2 alignment
include { SAMTOOLS_TOBAM                         } from '../../../modules/local/samtools/main.nf'                   // Convert SAM to BAM
include { SAMTOOLS_SORT as SAMTOOLS_SORT         } from '../../../modules/local/samtools/main.nf'                   // Sort BAM
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_CHUNK } from '../../../modules/local/samtools/main.nf'                   // Merge BAMs
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_FINAL } from '../../../modules/local/samtools/main.nf'                   // Merge BAMs
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FULL  } from '../../../modules/local/samtools/main.nf'                   // Index BAM
include { SAMTOOLS_TOFASTQ                       } from '../../../modules/local/samtools/main.nf'
include { CRAMINO_STATS                          } from '../../../modules/local/cramino/main.nf'
include { SEQKIT_STATS                           } from '../../../modules/local/seqkit/main.nf'                     // Coverage stats
include { modifyMetaId                           } from '../utils_nfcore_oncotrseq_pipeline'
include { QUARTO_TABLE                           } from '../../../modules/local/quarto/main.nf'                     // Reporting (optional)
include { paramsSummaryMap                       } from 'plugin/nf-schema'                                          // Parameter summary
include { paramsSummaryMultiqc                   } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'       // MultiQC summary
include { softwareVersionsToYAML                 } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'       // Version reporting
include { methodsDescriptionText                 } from '../../../subworkflows/local/utils_nfcore_oncotrseq_pipeline' // Methods for MultiQC

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN MAPPING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAPPING {
    take:
    in_ch
    ref


    main:
    ch_versions = Channel.empty()

    ch_ref = ref
        .map { meta, ref, ref_fasta, _ref_fai ->
            tuple(meta, ref, ref_fasta) }

    in_ch
        .map { meta, reads ->
            // Ensure 'reads' is a list and flatten it
            def dir_list = reads instanceof List ? reads.flatten() : [reads]
            def dir = file(dir_list[0])

            if (dir.isDirectory()) {
                // Collect all FASTQ files with common extensions
                def files = dir.listFiles().findAll { f ->
                    f.name ==~ /.*\.(fastq|fq)(\.gz)?$/
                }
                return tuple(meta, files)
            } else {
                return tuple(meta, [dir])
            }
        }
        .set { in_ch }

    ch_mapping_in = in_ch
        .join(ch_ref)
        .map { meta, fastq, ref, ref_fasta ->
            def meta_ref = modifyMetaId(meta, 'add_suffix', '', '', "_${ref}")
            tuple(meta_ref, fastq, ref_fasta)
            }

    // Run minimap2 alignment
    MINIMAP2_ALIGN(ch_mapping_in)

    // Convert SAM to BAM
    SAMTOOLS_TOBAM(MINIMAP2_ALIGN.out.sam)
    // Sort and index BAM
    SAMTOOLS_SORT(SAMTOOLS_TOBAM.out.bamfile)
    SAMTOOLS_INDEX_FULL(SAMTOOLS_SORT.out.sortedbam)

    // Restore meta ID by removing ref id

    ch_ref_id = ch_ref
        .map { meta, ref, _ref_fasta ->
            def meta_ref = modifyMetaId(meta, 'add_suffix', '', '', "_${ref}")
            tuple(meta_ref, meta.id) }

    bam_ch = SAMTOOLS_INDEX_FULL.out.bamfile_index
        .join(ch_ref_id)
        .map { meta, bam, bai, meta_restore ->
            tuple(id:meta_restore, bam, bai) }

    // Compute coverage stats
    CRAMINO_STATS(bam_ch)

    // Collect versions from all modules
    ch_versions = MINIMAP2_ALIGN.out.versions
        .mix(SAMTOOLS_TOBAM.out.versions)
        .mix(SAMTOOLS_SORT.out.versions)
        .mix(SAMTOOLS_INDEX_FULL.out.versions)
        .mix(CRAMINO_STATS.out.versions)

    
    emit:
    bam      = bam_ch                                 // Final sorted BAM with index
    coverage = CRAMINO_STATS.out.stats                // Coverage stats
    //seqkit   = ch_seqkit_out                          // Input fastq stats
    versions = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/