process RNA_EDITING_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path rediportal_vcf
    path rediportal_tbi
    val min_rna_support

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_support = min_rna_support ?: 2
    """
    annotate_rna_editing.py \\
        -i ${vcf} \\
        -r ${rediportal_vcf} \\
        -o ${prefix}.rna_annotated.vcf.gz \\
        --min-rna-support ${min_support} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotate_rna_editing: 1.0.0
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pysam: \$(python -c "import pysam; print(pysam.__version__)" 2>/dev/null || echo "unknown")
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//' || echo "unknown")
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^.*tabix //; s/ .*\$//' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rna_annotated.vcf.gz
    touch ${prefix}.rna_annotated.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotate_rna_editing: 1.0.0
        pysam: 0.22.0
        bcftools: 1.21
        tabix: 1.21
    END_VERSIONS
    """
}