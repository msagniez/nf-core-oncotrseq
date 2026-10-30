process QUANT_CLEANTRANSCRIPTIDS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
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
    // output must differ from the staged input name
    def prefix = task.ext.prefix ?: "${meta.id}.clean"
    """
    # Check whether transcript IDs (column 1) contain gene info, e.g. ENST...|ENSG...|...
    PIPE_COUNT=\$(head -n 10 ${quant} | awk '{print \$1}' | grep -c '|' || true)

    if [ "\$PIPE_COUNT" -eq 0 ]; then
        echo "transcript ids of ${quant} passed the check; format is correct"
        cp -L ${quant} ${prefix}.quant
    else
        echo "reformatting ${quant} transcript ids"
        awk $args 'BEGIN {FS=OFS="\\t"}
            {
                split(\$1, arr, "|")
                \$1 = arr[1]
                print
            }' ${quant} > ${prefix}.quant
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -Wversion 2>&1 | head -n1 | sed 's/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.clean"
    """
    touch ${prefix}.quant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: stub
    END_VERSIONS
    """
}