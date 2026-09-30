/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_oncotrseq_pipeline'
include { CAT_FASTQ                 } from '../modules/local/utils/cat_fastq.nf'
include { DOWNSAMPLE                } from '../subworkflows/local/downsample/downsample.nf'
include { MAPPING                   } from '../subworkflows/local/mapping/mapping.nf'
include { QUANTIFYING               } from '../subworkflows/local/quantifying/quantifying.nf'
include { QUARTO_TEXT               } from '../modules/local/quarto/main.nf'
include { QUARTO_SECTION            } from '../modules/local/quarto/main.nf'
include { QUARTO_REPORT             } from '../modules/local/quarto/main.nf'
include { QUARTO_TABLE              } from '../modules/local/quarto/main.nf'

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

    if (params.skip_quantification) {
        log.info "Skipping quantification step as requested by user."
        classif_in = samplesheet_matrix

    } else {
        
    

        ch_samplesheet = samplesheet_fastq.map { meta, biotype, tumor_type, input, ref_genome, ref_rna, gtf ->
            tuple(meta, biotype, input)
        }

        // Optional time-based downsampling: --downsamp 1h  or  --downsamp 1h,2h,72h,FULL
        // If --downsamp is not specified, the original samplesheet is used as-is.
        // Each sample is expanded into one copy per requested timepoint, tagged with
        // an explicit meta.timepoint field (e.g. "1h", "2h", "72h", "FULL") so that
        // results for every timepoint stay distinguishable downstream - meta.id
        // itself is left untouched.
        if (params.downsamp) {

            CAT_FASTQ(ch_samplesheet)
            ch_samplesheet_cat = CAT_FASTQ.out.reads

            ch_downsamp_times = Channel.fromList(parseDownsampParam(params.downsamp))

            ch_downsamp_in = ch_samplesheet_cat
                .combine(ch_downsamp_times)
                .map { meta, biotype, input, downsamp_label, downsamp_hours ->
                    def meta_ds = meta + [timepoint: downsamp_label]
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

        quant_ref = genome_fasta
            .map { meta, ref_name, ref_fasta, ref_fai -> tuple(meta.id, ref_fasta) }
            .join(ref_classif.map { meta, ref_fa, ref_gtf, ref_fai -> tuple(meta.id, ref_gtf) })
            .join(biotype.map { meta, biotype -> tuple(meta.id, biotype) })

        quant_in = MAPPING.out.collatedbam
            .map { meta, collated_bam -> tuple(meta.id, meta, collated_bam) }
            .combine(quant_ref, by: 0)
            .map { sample_id, meta, collated_bam, ref_fasta, ref_gtf, biotype ->
                def gtf_label = ref_gtf.name
                    .replaceAll(/\.gtf(\.gz)?$/, '')
                    .replace('.', '')
                    .replace('_', '')
                def meta_q = meta + [ref_gtf: gtf_label]
                tuple(meta_q, biotype, collated_bam, ref_fasta, ref_gtf)
            }

        QUANTIFYING (
            quant_in
        )

        classif_in = QUANTIFYING.out.quant
            .map { meta, quant -> tuple(meta.id, meta, quant) }
            .combine(samplesheet_fastq.map { meta, biotype, tumor_type, input, ref_genome, ref_rna, gtf -> tuple(meta.id, biotype, tumor_type, ref_genome, ref_rna, gtf) }, by: 0)
            .map { sample_id, meta, quant, biotype, tumor_type, ref_genome, ref_rna, gtf ->
                def meta_c = [id: "${meta.id}_${meta.timepoint}_${meta.ref_name}_${meta.ref_gtf}"]
                tuple(meta_c, biotype, tumor_type, quant, ref_rna, gtf)
            }
    }

    // Fork based on tumor_type to run different classification strategies.
    ch_classif_forked = classif_in
        .branch { meta, biotype, tumor_type, quant, ref_rna, gtf ->
            all: tumor_type == 'all'
            aml: tumor_type == 'aml'
            cascade: tumor_type == 'cascade'
            other: true
        }
        .set { ch_classif }

    ch_classif.all.view()
    

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Collect sections for report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_sections = MAPPING.out.section
                    
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOFTWARE VERSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    // Collect versions
    ch_versions = ch_versions
        .mix(MAPPING.out.versions)

    ch_samples = MAPPING.out.seqkit
    .map { meta, stats ->
        def meta_reduced = file(params.outdir).name
        tuple(id: meta_reduced)
    }
    .unique()

    // Extract all versions into a single channel of values
    versions = softwareVersionsToYAML(ch_versions)
    // Collapse the channel of versions into a single value
    versions = versions.collect().map { it.join('\n\n') }
    versions = ch_samples
        .combine(versions)

    // Give it an ID of versions
    versions = versions
        .map {
            versions_out ->
            def section = "Versions"
            def process = "versions"

            [versions_out[0], versions_out[1]] + [section, process]
            }

    QUARTO_TEXT(
        versions
        )

    ch_section_description = channel.of("Software Versions")

    ch_section_inputs = QUARTO_TEXT.out.quarto_text
        .combine(ch_section_description)

    QUARTO_SECTION(
        ch_section_inputs
    )
    // Add the versions to the channel of sections for every report

    ch_sections = ch_sections.mix(QUARTO_SECTION.out.quarto_section)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPORT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_template = channel.fromPath(params.report_template)
    ch_subtitle = channel.of('OncoTRseq: Nanopore RNA-Seq for Fusion Detection and Expression Profiling')
    ch_title    = channel.of('MPGI OncoTRseq Report')

    //Fix sections order in final report
    def sectionOrder = [
        'Mapping_QC',
        'Software Versions'
    ]

    ch_report_sections = ch_sections
        .groupTuple()
        .map { meta, sections, filePaths, reports ->
            def zipped = [sections, filePaths, reports].transpose()
            zipped = zipped.sort(false) { a, b ->
                def ia = sectionOrder.indexOf(a[0])
                def ib = sectionOrder.indexOf(b[0])
                if (ia < 0) ia = sectionOrder.size()
                if (ib < 0) ib = sectionOrder.size()
                ia <=> ib
            }

            def sortedSections  = zipped.collect { it[0] }
            def sortedFilePaths = zipped.collect { it[1] }
            def sortedReports   = zipped.collect { it[2] }

            [meta, sortedSections, sortedFilePaths, sortedReports]
        }
        .combine(ch_title)
        .combine(ch_subtitle)
        .combine(ch_template)

    QUARTO_REPORT(
        ch_report_sections
    )

    ch_report = QUARTO_REPORT.out.report


    emit:
    multiqc_report = Channel.empty()
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/