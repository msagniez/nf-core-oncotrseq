process SAMTOOLS_COLLATE {
    // TODO : SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/chusj-pigu/samtools:b195aca24376fa3482000f5bcdc804ac36d9da0b'

    label "process_medium_cpu"              // Label for mpgi drac memory alloc
    label "process_medium_memory"           // Label for mpgi drac memory alloc
    label "process_medium_low_time"         // Label for mpgi drac time alloc

    tag "$meta.id"

    input:
    tuple val(meta),
        path(in_bam)

    output:
    tuple val(meta),
        path("*.collated.bam"),
        emit: collatedbam
    path "versions.yml",
        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args
        if (task.ext.args) {
            args = task.ext.args
        } else if (meta.id.toString().contains('oarfish')) {      // Sort by read ID required by oarfish
            args = '-n'
        } else {
            args = ''
        }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus
    """
    samtools \\
        collate \\
        -@ ${threads} \\
        ${args} \\
        ${in_bam} \\
        -o ${prefix}.collated.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}