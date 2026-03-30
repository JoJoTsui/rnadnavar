process SALMON_QUANT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::salmon=1.11.4 conda-forge::libstdcxx-ng>=14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.11.4--h43eeafb_0' :
        'quay.io/biocontainers/salmon:1.11.4--h43eeafb_0' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${prefix}/quant.sf")              , emit: quant_sf
    tuple val(meta), path("${prefix}/lib_format_counts.json"), emit: lib_format_counts
    tuple val(meta), path("${prefix}/cmd_info.json")         , emit: cmd_info
    tuple val(meta), path("${prefix}/aux_info/")             , emit: aux_info
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"

    def read_args
    if (reads.size() == 2) {
        // Paired-end: reads may be single files or lists (multi-lane)
        def r1 = reads[0] instanceof List ? reads[0].join(' ') : reads[0]
        def r2 = reads[1] instanceof List ? reads[1].join(' ') : reads[1]
        read_args = "-1 ${r1} -2 ${r2}"
    } else {
        // Single-end
        def r = reads[0] instanceof List ? reads[0].join(' ') : reads[0]
        read_args = "-r ${r}"
    }

    """
    salmon quant \\
        --index ${index} \\
        ${args} \\
        ${read_args} \\
        -p ${task.cpus} \\
        --validateMappings \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version 2>&1 | sed 's/salmon //g')
    END_VERSIONS
    """
}
