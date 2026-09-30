process CAT_FASTQ {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), val(biotype), path(reads)

    output:
    tuple val(meta), val(biotype), path("${meta.id}.merged.fastq.gz"), emit: reads

    script:
    """
    find -L ${reads} -type f \\( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \\) \\
        | sort \\
        | xargs zcat -f \\
        | gzip -c > ${meta.id}.merged.fastq.gz

    # fail early if nothing was found
    if [ ! -s ${meta.id}.merged.fastq.gz ]; then
        echo "No fastq reads found for ${meta.id}" >&2
        exit 1
    fi
    """

    stub:
    """
    echo "" | gzip > ${meta.id}.merged.fastq.gz
    """
}