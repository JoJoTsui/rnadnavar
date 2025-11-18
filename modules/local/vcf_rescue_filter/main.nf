process VCF_RESCUE_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dna_thr = task.ext.dna_thr ?: 2
    def rna_thr = task.ext.rna_thr ?: 2

    """
    filter_rescue_vcf.py \\
        --input ${vcf} \\
        --output ${prefix}.filtered.vcf.gz \\
        --dna_threshold ${dna_thr} \\
        --rna_threshold ${rna_thr}
    
    tabix -p vcf ${prefix}.filtered.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        cyvcf2: \$(python -c "import cyvcf2; print(cyvcf2.__version__)")
    END_VERSIONS
    """
}
