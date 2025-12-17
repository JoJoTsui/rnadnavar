process COSMIC_GNOMAD_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path cosmic_vcf
    path cosmic_tbi
    path gnomad_dir
    val verbose_logging

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    path "*_stats.json", emit: stats, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def verbose = verbose_logging ? '--verbose' : ''
    def cosmic_arg = cosmic_vcf && cosmic_vcf.name != 'NO_FILE' ? "--cosmic ${cosmic_vcf}" : ''
    def gnomad_arg = gnomad_dir && gnomad_dir.name != 'NO_FILE' ? "--gnomad ${gnomad_dir}" : ''
    def stats_output = verbose_logging ? "--stats-output ${prefix}_stats.json" : ''

    """
    annotate_cosmic_gnomad.py \\
        --input ${vcf} \\
        ${cosmic_arg} \\
        ${gnomad_arg} \\
        --output ${prefix}.vcf.gz \\
        --workers ${task.cpus} \\
        ${verbose} \\
        ${stats_output} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pysam: \$(python -c "try: import pysam; print(pysam.__version__);" 2>/dev/null || echo "not available")
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//' || echo "not available")
        htslib: \$(tabix --version 2>&1 | head -n1 | sed 's/^.*tabix (htslib) //; s/ .*\$//' || echo "not available")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stats_file = verbose_logging ? "${prefix}_stats.json" : ''
    """
    touch ${prefix}.cosmic_gnomad_annotated.vcf.gz
    touch ${prefix}.cosmic_gnomad_annotated.vcf.gz.tbi
    ${stats_file ? "touch ${stats_file}" : ''}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
        pysam: 0.22.1
        bcftools: 1.21
        htslib: 1.21
    END_VERSIONS
    """
}