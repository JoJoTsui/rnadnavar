process VT_DECOMPOSE {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vt:2015.11.10--h5ef6573_4':
        'biocontainers/vt:2015.11.10--h5ef6573_4' }"

    input:
    tuple val(meta), path(vcf), path(intervals)

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: vcf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$vcf" == "${prefix}.vcf" || "$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    def bed = intervals ? "-i ${intervals}" : ""
    def VERSION = "2015.11.10" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    vt decompose \\
        -o ${prefix}.vcf \\
        ${bed} \\
        ${args} \\
        ${vcf}

    bgzip ${args2} --threads ${task.cpus} ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vt: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$vcf" == "${prefix}.vcf" || "$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    def VERSION = "2015.11.10" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    # Keep the stub graph semantically valid: downstream consensus parses the
    # decomposed VCF, so an empty placeholder is not a usable test artifact.
    # The synthetic run does not need decomposition itself; preserving the
    # input VCF gives later processes a valid, deterministic record set.
    if [[ "${vcf}" == *.gz ]]; then
        cp "${vcf}" "${prefix}.vcf.gz"
    else
        gzip -c "${vcf}" > "${prefix}.vcf.gz"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vt: ${VERSION}
    END_VERSIONS
    """
}
