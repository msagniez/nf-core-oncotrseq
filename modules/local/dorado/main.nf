nextflow.enable.dsl=2

def DORADO_CONTAINER =
    'ghcr.io/chusj-pigu/dorado:851b0a37ecbff6b3b952c811b31dad48ca6bf995'

process DORADO_BASECALL {
    container DORADO_CONTAINER

    // nf-core resource label
    label "process_high"

    // MPGI DRAC resource labels
    label "process_medium_low_cpu"
    label "process_higher_memory"
    label "process_medium_high_time"
    label "process_gpu"

    tag "$meta.id"

    input:
    tuple val(meta),
        path(pod5),
        path(ubam),
        path(model),
        path(model_mh)

    output:
    tuple val(meta),
        path("*.bam"),
        emit: ubam
    path "versions.yml",
        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def device = params.device != null ? "-x $params.device" : ""
    def mod = model_mh.name != "NOMOD"
        ? "--modified-bases-models ${model_mh}"
        : ""
    def multi = params.demux ? "--no-trim" : ""
    def resume = ubam.name != 'NOFILE'
        ? "--resume-from $ubam > ${prefix}_unaligned_final.bam"
        : "> ${prefix}_unaligned.bam"
    """
    dorado basecaller \\
        $args \\
        $device \\
        $model \\
        $mod \\
        $pod5 \\
        $multi \\
        $resume

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | \
            sed 's/^.*dorado //; s/Using.*\$//')
    END_VERSIONS
    """
}

process DORADO_DEMULTIPLEX {
    container DORADO_CONTAINER

    // nf-core resource label
    label "process_high"

    // MPGI DRAC resource labels
    label "process_high_cpu"
    label "process_medium_low_memory"
    label "process_medium_low_time"

    tag "$meta.id"

    input:
    tuple val(meta),
        val(kit),
        path(bam)

    output:
    tuple val(meta),
        path("${meta.id}"),
        emit: demux
    path "versions.yml",
        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dorado \\
        demux \\
        $args \\
        --kit-name ${kit} \\
        --output-dir $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | \
            sed 's/^.*dorado //; s/Using.*\$//')
    END_VERSIONS
    """
}

process DORADO_DOWNLOAD_LIST {
    container DORADO_CONTAINER

    // nf-core resource label
    label "process_low"

    label "local"

    // MPGI DRAC resource labels
    label "process_single_cpu"
    label "process_very_low_memory"
    label "process_very_low_time"

    output:
    path("*.json"),
        emit: list
    path "versions.yml",
        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    dorado download --list-structured > model_list.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | \
            sed 's/^.*dorado //; s/Using.*\$//')
    END_VERSIONS
    """
}

process DORADO_DOWNLOAD_MODEL {
    container DORADO_CONTAINER

    // nf-core resource label
    label "process_low"

    label "local"

    // MPGI DRAC resource labels
    label "process_single_cpu"
    label "process_very_low_memory"
    label "process_very_low_time"

    tag "$model"

    input:
    tuple val(type),
        val(model),
        path(dir)

    output:
    tuple val(type),
        path("${dir}/dorado_models/${model}"),
        emit:model
    path "versions.yml",
        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p ${dir}/dorado_models
    dorado download --model ${model} --models-directory ${dir}/dorado_models

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | \
            sed 's/^.*dorado //; s/Using.*\$//')
    END_VERSIONS
    """
}

process DORADO_TRIM {
    container 'ghcr.io/chusj-pigu/dorado:851b0a37ecbff6b3b952c811b31dad48ca6bf995'
    label "process_high"                    // nf-core label
    label "process_high_cpu"                 // Label for mpgi drac cpu alloc
    label "process_medium_low_memory"        // Label for mpgi drac memory alloc
    label "process_medium_low_time"          // Label for mpgi drac time alloc

    tag "$meta.id"

    input:
    tuple val(meta),
        val(kit),
        path(bam)

    output:
    tuple val(meta),
        path("*.bam"),
        emit: bam
    path "versions.yml",
        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus
    """
    dorado \\
        trim \\
        ${args} \\
        ${bam} \\
        -t ${threads} \\
        -k ${kit} > ${prefix}_trimmed.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed 's/^.*dorado //; s/Using.*\$//')
    END_VERSIONS
    """
}
