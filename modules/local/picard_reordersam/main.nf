process PICARD_REORDERSAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0' :
        'biocontainers/picard:3.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(meta), path("*.reordered.bam"), path("*.reordered.bam.bai"), emit: bam
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[PICARD ReorderSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        ReorderSam \\
        I=${bam} \\
        O=${prefix}.reordered.bam \\
        SD=${dict} \\
        ALLOW_INCOMPLETE_DICT_CONCORDANCE=true \\
        CREATE_INDEX=true \\
        $args

    # Rename index to match expected naming convention
    mv ${prefix}.reordered.bai ${prefix}.reordered.bam.bai || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard ReorderSam --version 2>&1 | grep -o 'Version:.*' | cut -f2 -d: | tr -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reordered.bam
    touch ${prefix}.reordered.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard ReorderSam --version 2>&1 | grep -o 'Version:.*' | cut -f2 -d: | tr -d ' ')
    END_VERSIONS
    """
}
