process VCF_FILTERING {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(fasta)
    path whitelist
    path blacklist

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.filtered.vcf.stripped.vcf.gz"), path("*.filtered.vcf.stripped.vcf.gz.tbi"), emit: vcf_stripped
    path "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def whitelist_arg = whitelist ? "--whitelist ${whitelist}" : ""
    def blacklist_arg = blacklist ? "--blacklist ${blacklist}" : ""
    def ref_arg = fasta ? "--ref ${fasta}" : ""

    """
    filter_vcf.py \\
        -i ${vcf} \\
        -o ${prefix}.filtered.vcf.gz \\
        --output_stripped ${prefix}.filtered.vcf.stripped.vcf.gz \\
        --filter_multiallelic \\
        ${whitelist_arg} \\
        ${blacklist_arg} \\
        ${ref_arg} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        cyvcf2: \$(python -c "import cyvcf2; print(cyvcf2.__version__)")
    END_VERSIONS
    """
}
