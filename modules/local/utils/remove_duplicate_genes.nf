process REMOVE_DUPLICATE_GENES {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(quant)

    output:
    tuple val(meta), path("*.quant"), emit: quant
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    // must differ from the staged input name (upstream emits "<id>.clean.quant")
    def prefix = task.ext.prefix ?: "${meta.id}.dedup"
    """
    # Keep the header and the first occurrence of each ID in column 1, preserving order
    awk $args 'BEGIN {FS=OFS="\\t"} NR==1 || !seen[\$1]++' ${quant} > ${prefix}.quant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -Wversion 2>&1 | head -n1 | sed 's/,.*//' || echo unknown)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.dedup"
    """
    touch ${prefix}.quant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: stub
    END_VERSIONS
    """
}