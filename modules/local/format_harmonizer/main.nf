process FORMAT_HARMONIZER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.harmonized.vcf.gz"), path("*.harmonized.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    harmonize_vcf_format.py \\
        --vcf ${vcf} \\
        --caller ${meta.variantcaller} \\
        --output ${prefix}.harmonized.vcf.gz \\
        ${args}

    tabix -p vcf ${prefix}.harmonized.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        cyvcf2: \$(python -c "import cyvcf2; print(cyvcf2.__version__)")
    END_VERSIONS
    """
}
