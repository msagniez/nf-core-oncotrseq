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
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_CDNA  } from '../../../modules/local/minimap2/main.nf'                   // minimap2 alignment for cDNA
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_DRNA  } from '../../../modules/local/minimap2/main.nf'                   // minimap2 alignment for dRNA
include { SAMTOOLS_TOBAM                         } from '../../../modules/local/samtools/main.nf'                   // Convert SAM to BAM
include { SAMTOOLS_SORT                          } from '../../../modules/local/samtools/main.nf'                   // Sort BAM
include { SAMTOOLS_COLLATE                       } from '../../../subworkflows/local/utils/samtools_collate.nf'     // Collate BAM
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_CHUNK } from '../../../modules/local/samtools/main.nf'                   // Merge BAMs
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_FINAL } from '../../../modules/local/samtools/main.nf'                   // Merge BAMs
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FULL  } from '../../../modules/local/samtools/main.nf'                   // Index BAM
include { SAMTOOLS_TOFASTQ                       } from '../../../modules/local/samtools/main.nf'
include { CRAMINO_STATS                          } from '../../../modules/local/cramino/main.nf'
include { SEQKIT_STATS                           } from '../../../modules/local/seqkit/main.nf'                     // Coverage stats
include { QUARTO_TABLE_COLNAMES as QUARTO_TABLE1 } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION                         } from '../../../modules/local/quarto/main.nf'
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
        .map { meta, ref_name, ref_fasta, _ref_fai ->
            tuple(meta, ref_name, ref_fasta) }

    SEQKIT_STATS(in_ch.map { meta, biotype, reads -> tuple(meta, reads) })
    ch_seqkit_out = SEQKIT_STATS.out.stats

    in_ch
        .map { meta, biotype, reads ->
            // Ensure 'reads' is a list and flatten it
            def dir_list = reads instanceof List ? reads.flatten() : [reads]
            def dir = file(dir_list[0])

            if (dir.isDirectory()) {
                // Collect all FASTQ files with common extensions
                def files = dir.listFiles().findAll { f ->
                    f.name ==~ /.*\.(fastq|fq)(\.gz)?$/
                }
                return tuple(meta, biotype, files)
            } else {
                return tuple(meta, biotype, [dir])
            }
        }
        .set { in_ch }

    in_ch_keyed  = in_ch.map { meta, biotype, fastq -> tuple(meta.id, meta, biotype, fastq) }
    ch_ref_keyed = ch_ref.map { meta, ref_name, ref_fasta -> tuple(meta.id, ref_name, ref_fasta) }

    ch_mapping_in = in_ch_keyed
        .combine(ch_ref_keyed, by: 0)
        .map { sample_id, meta, biotype, fastq, ref_name, ref_fasta ->
            def meta_ref = meta + [ref_name: ref_name]
            tuple(meta_ref, biotype, fastq, ref_fasta)
    }

    // Fork based on biotype to run different mapping strategies.
    ch_mapping_in
        .branch { meta, biotype, fastq, ref_fasta ->
            drna: biotype == 'drna'
            cdna: true
        }
        .set { ch_branched }

    ch_mapping_drna = ch_branched.drna
    ch_mapping_cdna = ch_branched.cdna

    // Run minimap2 alignment
    MINIMAP2_ALIGN_CDNA(ch_mapping_cdna.map { meta, biotype, fastq, ref_fasta -> tuple(meta, fastq, ref_fasta) })
    MINIMAP2_ALIGN_DRNA(ch_mapping_drna.map { meta, biotype, fastq, ref_fasta -> tuple(meta, fastq, ref_fasta) })

    //Merge channels back together
    ch_mapping_out = Channel.empty()
    ch_mapping_out = MINIMAP2_ALIGN_CDNA.out.sam
        .mix(MINIMAP2_ALIGN_DRNA.out.sam)

    // Convert SAM to BAM
    SAMTOOLS_TOBAM(ch_mapping_out)
    // Sort and index BAM
    SAMTOOLS_SORT(SAMTOOLS_TOBAM.out.bamfile)
    SAMTOOLS_INDEX_FULL(SAMTOOLS_SORT.out.sortedbam)

    bam_ch = SAMTOOLS_INDEX_FULL.out.bamfile_index

    // Compute coverage stats
    CRAMINO_STATS(bam_ch)

    // Collate bam for downstream oarfish processing
    SAMTOOLS_COLLATE(SAMTOOLS_TOBAM.out.bamfile)


    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TABLE INPUT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Get only necessary lines from cramino summary

    def desired_stats = [
    'Number of alignments',
    '% from total alignments',
    'Yield [Gb]',
    'Mean coverage',
    'N50',
    'Mean length',
    'Mean identity'
    ]

    ch_cramino_parsed = CRAMINO_STATS.out.stats
        .map { meta, table ->
            def stats = table.readLines()
                .findAll { it.contains('\t') }
                .collectEntries { line ->
                    def parts = line.split('\t', 2)
                    [(parts[0]): parts[1]]
                }
            tuple(meta, stats)
        }

    ch_mean_coverage = ch_cramino_parsed
        .map { meta, stats -> tuple(meta, stats['Mean coverage'] ?: 'NA') }

    ch_cramino_summary = ch_cramino_parsed
    .map { meta, stats ->
        def project = file(params.outdir).name
        tuple(project, meta, stats)
    }
    .collectFile(sort: true) { project, meta, stats ->
        def sample  = meta.timepoint ? "${meta.id}_${meta.timepoint}" : meta.id
        def refname = meta.ref_name.toString().replaceAll(/\.(fa|fasta)(\.gz)?$/, '')
        def row = ([sample, refname] + desired_stats.collect { stats[it] ?: 'NA' }).join('\t')
        return [ "${project}_mapping_stats.tsv", "${row}\n" ]
    }

    ch_quarto_table = ch_cramino_summary
        .map { table ->
            def meta_project = table.name.replace('_mapping_stats.tsv', '')
            tuple(id: meta_project, table)
        }
        .map { meta, table ->
            def caption   = "Summary mapping stats for ${meta.id} on BAM files"
            def col_names = "Sample, Ref, # Alignments, % from total alignments, Yield [Gb], Mean Coverage, N50, Mean length, Mean identity"
            def section   = "Mapping_QC"
            def process   = "mapping-qc-${meta.id}"
            tuple(meta, table, caption, col_names, section, process)
        }
    

    QUARTO_TABLE1(
        ch_quarto_table
        )

    ch_section_inputs = QUARTO_TABLE1.out.quarto_table


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QUARTO_SECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_section_description = channel.of("Summary Statistics for BAM files")

    ch_section_inputs = ch_section_inputs
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }.combine(ch_section_description)

    QUARTO_SECTION(
        ch_section_inputs
    )

    ch_section = QUARTO_SECTION.out.quarto_section

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT VERSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_versions = MINIMAP2_ALIGN_CDNA.out.versions
        .mix(MINIMAP2_ALIGN_DRNA.out.versions)
        .mix(SAMTOOLS_TOBAM.out.versions)
        .mix(SAMTOOLS_SORT.out.versions)
        .mix(SAMTOOLS_COLLATE.out.versions)
        .mix(SAMTOOLS_INDEX_FULL.out.versions)
        .mix(CRAMINO_STATS.out.versions)
        .mix(QUARTO_TABLE1.out.versions)
        .mix(QUARTO_SECTION.out.versions)


    emit:
    sortedbam      = bam_ch                             // Final sorted BAM with index
    collatedbam    = SAMTOOLS_COLLATE.out.collatedbam   // Collated BAM
    coverage       = CRAMINO_STATS.out.stats            // Coverage stats
    seqkit         = ch_seqkit_out                      // Input fastq stats
    section        = QUARTO_SECTION.out.quarto_section  // Quarto section for mapping stats
    versions       = ch_versions                        // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/