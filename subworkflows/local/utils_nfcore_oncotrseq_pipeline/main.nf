//
// Subworkflow with functionality specific to the nf-core/oncotrseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME  } from '../../../modules/local/samtools/main.nf'    // indexes the alignment genome (--genome / ref_genome column)
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CLASSIF } from '../../../modules/local/samtools/main.nf'    // indexes the classification reference (ref_classif column)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// File extensions used to decide whether an input row is raw reads (fastq)
// or an already-quantified matrix. No basecalling happens in this pipeline,
// so this is purely a routing decision for downstream classification.
def FASTQ_EXTENSIONS_REGEX  = /.*\.(fastq|fq)(\.gz)?$/
def MATRIX_EXTENSIONS_REGEX = /.*\.(csv|tsv)(\.gz)?$/

// Accepted biotypes for classification
def VALID_BIOTYPES = ['cdna', 'drna']

// Genome ID aliases -> canonical genome key used to name files in --ref_dir / on UCSC.
// Only used for the alignment genome (ref_genome / --genome), never for ref_classif/gtf.
def GENOME_ALIAS_MAP = [
    'hg38'   : 'hg38',
    'GRCh38' : 'hg38',
    'hg19'   : 'hg19',
    'GRCh37' : 'hg19',
    'hs1'    : 'hs1',
    'CHM13'  : 'hs1'
]
def DEFAULT_GENOME        = 'hg38'
def UCSC_ALLOWED_GENOMES  = ['hg38', 'hg19', 'hs1']

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Resolve a samplesheet value (ref_classif or gtf column) to an existing file.
    // The column may contain either a full path, or just a filename to be looked
    // up inside --ref_dir. Always required - no silent fallback.
    //
    def resolveRefFile = { value, label, meta ->
        if (!value) {
            throw new IllegalArgumentException(
                "${label} is required for sample '${meta.id}'"
            )
        }

        def direct = file(value)
        if (direct.exists()) {
            return direct
        }

        if (!params.ref_dir) {
            throw new IllegalArgumentException(
                "${label} '${value}' is not an existing path for sample '${meta.id}', " +
                "and --ref_dir is not set to look it up by filename"
            )
        }

        def in_ref_dir = file("${params.ref_dir}/${value}")
        if (!in_ref_dir.exists()) {
            throw new IllegalArgumentException(
                "${label} '${value}' was not found as a path, nor as a file named '${value}' " +
                "in --ref_dir (${params.ref_dir}) for sample '${meta.id}'"
            )
        }

        return in_ref_dir
    }

    //
    // Parse samplesheet once and validate every row:
    //   - the input file (fastq or matrix) must exist
    //   - biotype must be 'cdna' or 'drna'
    //   - the input file extension must match either a fastq or a matrix pattern
    //   - ref_classif and gtf are mandatory and resolved to existing files
    //     (full path as given, or filename looked up in --ref_dir)
    // Expected columns: meta, biotype, tumor_type, input_file, ref_genome, ref_classif, gtf
    // `input_file` is either a fastq (raw reads) or a quantification matrix.
    //
    ch_samplesheet_raw = channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, biotype, tumor_type, input_file, ref_genome, ref_classif, gtf ->

            if (!input_file.exists()) {
                throw new IllegalArgumentException(
                    "Input file not found for sample '${meta.id}': ${input_file}"
                )
            }

            if (!(biotype in VALID_BIOTYPES)) {
                throw new IllegalArgumentException(
                    "Invalid biotype '${biotype}' for sample '${meta.id}': must be one of ${VALID_BIOTYPES}"
                )
            }

            def is_fastq  = input_file.name ==~ FASTQ_EXTENSIONS_REGEX
            def is_matrix = input_file.name ==~ MATRIX_EXTENSIONS_REGEX
            if (!is_fastq && !is_matrix) {
                throw new IllegalArgumentException(
                    "Unsupported input file '${input_file.name}' for sample '${meta.id}': " +
                    "expected a fastq (.fastq/.fq[.gz]) or a matrix (.csv/.tsv[.gz])"
                )
            }

            // Classification reference (fasta) and its annotation (gtf) are a separate
            // reference from the alignment genome (ref_genome / genome_fasta below),
            // and are mandatory for every sample regardless of fastq/matrix.
            def resolved_ref_classif = resolveRefFile(ref_classif, 'ref_classif', meta)
            def resolved_gtf         = resolveRefFile(gtf, 'gtf', meta)

            return [ meta, biotype, tumor_type, input_file, ref_genome, resolved_ref_classif, resolved_gtf ]
        }

    //
    // Branch into fastq (needs alignment/quantification before classification)
    // vs matrix (already quantified, goes straight to classification). Extension
    // was already validated above, so this only routes. ref_classif/gtf are kept
    // in both branches; only ref_genome (alignment-only) is dropped for matrix.
    //
    ch_samplesheet_raw
        .branch { meta, biotype, tumor_type, input_file, ref_genome, ref_classif, gtf ->
            fastq:  input_file.name ==~ FASTQ_EXTENSIONS_REGEX
                return [ meta, biotype, tumor_type, input_file, ref_genome, ref_classif, gtf ]
            matrix: true
                return [ meta, biotype, tumor_type, input_file, ref_classif, gtf ]
        }
        .set { ch_samplesheet_branched }

    ch_samplesheet_fastq  = ch_samplesheet_branched.fastq
    ch_samplesheet_matrix = ch_samplesheet_branched.matrix

    ch_samplesheet_fastq.view  { "Fastq input:  ${it}" }
    ch_samplesheet_matrix.view { "Matrix input: ${it}" }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Alignment genome resolution - fastq samples only (matrix samples are already
    quantified and never need an alignment reference). This is independent from
    ref_classif/gtf above.

    Priority per sample:
      1. ref_genome column is an existing FASTA file path     -> used directly
      2. ref_genome column (or --genome) is a genome ID
         and --ref_dir is set                                 -> looked up by alias in --ref_dir
      3. neither set                                          -> fall back to --genome
         (default 'hg38'); looked up in --ref_dir if set,
         otherwise downloaded from the UCSC FTP mirror

    Missing .fai indexes are generated below with SAMTOOLS_FAIDX_GENOME.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def faExts  = ['.fa', '.fasta', '.fa.gz', '.fasta.gz']
    def faiExts = ['.fa.fai', '.fasta.fai', '.fa.gz.fai', '.fasta.gz.fai']

    def findFa = { dir, aliases ->
        dir?.listFiles()?.find { f ->
            aliases.any { a -> f.name.toLowerCase().contains(a.toLowerCase()) } &&
            faExts.any { ext -> f.name.toLowerCase().endsWith(ext) }
        }
    }

    def findFai = { dir, aliases ->
        dir?.listFiles()?.find { f ->
            aliases.any { a -> f.name.toLowerCase().contains(a.toLowerCase()) } &&
            faiExts.any { ext -> f.name.toLowerCase().endsWith(ext) }
        }
    }

    ch_samplesheet_fastq
        .map { meta, _biotype, _tumor_type, _input_file, ref_genome, _ref_classif, _gtf ->

            def genome_input = ref_genome ?: params.genome
            def fa  = null
            def fai = null

            // CASE 1: explicit FASTA path
            if (genome_input && file(genome_input).exists()) {
                fa = file(genome_input)
                def fai_candidate = file("${fa}.fai")
                fai = fai_candidate.exists() ? fai_candidate : null
            }
            // CASE 2: genome ID -> lookup in --ref_dir
            else if (genome_input) {
                if (!params.ref_dir) {
                    throw new IllegalArgumentException(
                        "Genome '${genome_input}' provided for sample '${meta.id}' but --ref_dir is not set"
                    )
                }

                def canonical = GENOME_ALIAS_MAP[genome_input] ?: genome_input
                def aliases   = [canonical, genome_input]
                def ref_dir   = file(params.ref_dir)

                fa  = findFa(ref_dir, aliases)
                fai = findFai(ref_dir, aliases)

                if (!fa) {
                    throw new IllegalArgumentException(
                        "Could not find FASTA for genome '${genome_input}' in ${params.ref_dir} (sample: ${meta.id})"
                    )
                }
            }
            // CASE 3: no genome specified -> default genome, --ref_dir lookup, else UCSC download
            else {
                def fallback  = params.genome ?: DEFAULT_GENOME
                def canonical = GENOME_ALIAS_MAP[fallback] ?: fallback

                if (params.ref_dir) {
                    def ref_dir = file(params.ref_dir)
                    fa  = findFa(ref_dir, [canonical])
                    fai = findFai(ref_dir, [canonical])
                }

                if (!fa) {
                    log.warn("No matching genome found in --ref_dir -- staging default genome '${canonical}' from UCSC for sample '${meta.id}'")
                    if (!UCSC_ALLOWED_GENOMES.contains(canonical)) {
                        log.warn("Genome '${canonical}' may not exist on the UCSC FTP mirror")
                    }
                    fa  = file("ftp://hgdownload.cse.ucsc.edu/goldenPath/${canonical}/bigZips/${canonical}.fa.gz")
                    fai = null
                }
            }

            return tuple(meta, fa, fai)
        }
        .set { ch_genome_ref }

    ch_genome_ref
        .branch { meta, fa, fai ->
            needs_index: fai == null
                return [ meta, fa ]
            has_index: true
                return [ meta, fa, fai ]
        }
        .set { ch_genome_split }

    ch_genome_needs_index = ch_genome_split.needs_index
    ch_genome_with_index  = ch_genome_split.has_index

    ch_genome_needs_index.view { "Genome indexing required: ${it}" }

    SAMTOOLS_FAIDX_GENOME (
        ch_genome_needs_index
    )

    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions.first())

    // The module copies/renames the fasta alongside the new .fai and emits both
    // together as (meta, fasta, fai) - use that directly, no join needed.
    ch_genome_indexed = SAMTOOLS_FAIDX_GENOME.out.fasta_index

    ch_genome_final = ch_genome_indexed.mix(ch_genome_with_index)

    ch_genome_final.view { "Genome ready: ${it}" }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Classification reference indexing - both fastq and matrix samples.
    ref_classif and gtf were already resolved to existing files above; only the
    .fai index may still be missing, and is generated here if so.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_classif_raw = ch_samplesheet_fastq
        .map { meta, _biotype, _tumor_type, _input_file, _ref_genome, ref_classif, gtf -> [ meta, ref_classif, gtf ] }
        .mix(
            ch_samplesheet_matrix
                .map { meta, _biotype, _tumor_type, _input_file, ref_classif, gtf -> [ meta, ref_classif, gtf ] }
        )

    ch_classif_raw
        .branch { meta, ref_classif, gtf ->
            def fai_candidate = file("${ref_classif}.fai")
            needs_index: !fai_candidate.exists()
                return [ meta, ref_classif, gtf ]
            has_index: true
                return [ meta, ref_classif, gtf, fai_candidate ]
        }
        .set { ch_classif_split }

    ch_classif_needs_index = ch_classif_split.needs_index
    ch_classif_with_index  = ch_classif_split.has_index

    ch_classif_needs_index.view { "Classification reference indexing required: ${it}" }

    SAMTOOLS_FAIDX_CLASSIF (
        ch_classif_needs_index.map { meta, ref_classif, _gtf -> tuple(meta, ref_classif) }
    )

    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_CLASSIF.out.versions.first())

    // The module emits (meta, fasta, fai) but doesn't know about gtf - rejoin it
    // from the pre-indexing channel, keyed by meta, then reorder to (meta, fasta, gtf, fai).
    ch_classif_indexed = ch_classif_needs_index
        .map { meta, _ref_classif, gtf -> tuple(meta, gtf) }
        .join(SAMTOOLS_FAIDX_CLASSIF.out.fasta_index)
        .map { meta, gtf, fasta, fai -> tuple(meta, fasta, gtf, fai) }

    ch_classif_final = ch_classif_indexed.mix(ch_classif_with_index)

    ch_classif_final.view { "Classification reference ready: ${it}" }

    //
    // Extract biotype per sample (both branches)
    //
    ch_biotype = ch_samplesheet_fastq
        .map { meta, biotype, _tumor_type, _input_file, _ref_genome, _ref_classif, _gtf -> [ meta, biotype ] }
        .mix(
            ch_samplesheet_matrix
                .map { meta, biotype, _tumor_type, _input_file, _ref_classif, _gtf -> [ meta, biotype ] }
        )

    ch_biotype.view { "Biotype: ${it}" }

    emit:
    fastq                        = ch_samplesheet_fastq
    matrix                       = ch_samplesheet_matrix
    genome_fasta                 = ch_genome_final          // (meta, fasta, fai)        - alignment genome, fastq only
    ref_classif                  = ch_classif_final          // (meta, fasta, gtf, fai)   - classification reference, all samples
    biotype                      = ch_biotype
    versions                     = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Function to modify metadata id field in a flexible way
// This is a generalized function that can handle various metadata transformations
//
def modifyMetaId(Map meta, String operation, String search_string = '', String replace_string = '', def suffix = '') {
    // Clone metadata and normalize to strings
    def new_meta = meta.collectEntries { k, v -> [k, v?.toString()] }

    switch(operation) {
        case 'remove_suffix':
            if (new_meta.id && suffix) {
                new_meta.id = new_meta.id.replaceFirst(suffix.toString(), '')
            }
            break

        case 'add_suffix':
            if (new_meta.id && suffix) {
                new_meta.id = new_meta.id + suffix
            }
            break

        case 'replace':
            if (new_meta.id && search_string) {
                new_meta.id = new_meta.id.replace(search_string, replace_string ?: '')
            }
            break

        case 'prefix':
            if (new_meta.id && suffix) {
                new_meta.id = suffix + new_meta.id
            }
            break
    }

    return new_meta
}