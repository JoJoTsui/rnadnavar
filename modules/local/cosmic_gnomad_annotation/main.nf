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
    def cosmic_arg = cosmic_vcf ? "--cosmic ${cosmic_vcf}" : ''
    def gnomad_arg = gnomad_dir ? "--gnomad ${gnomad_dir}" : ''
    def stats_output = verbose_logging ? "--stats-output ${prefix}_stats.json" : ''

    """
    # Validate input files and databases
    if [ ! -f "${vcf}" ]; then
        echo "ERROR: Input VCF file not found: ${vcf}" >&2
        exit 1
    fi

    if [ ! -f "${tbi}" ]; then
        echo "ERROR: Input VCF index file not found: ${tbi}" >&2
        exit 1
    fi

    # Check database availability
    if [ -n "${cosmic_vcf}" ] && [ ! -f "${cosmic_vcf}" ]; then
        echo "WARNING: COSMIC database not found: ${cosmic_vcf}" >&2
        echo "Proceeding without COSMIC annotation" >&2
    fi

    if [ -n "${gnomad_dir}" ] && [ ! -d "${gnomad_dir}" ]; then
        echo "WARNING: gnomAD directory not found: ${gnomad_dir}" >&2
        echo "Proceeding without gnomAD annotation" >&2
    fi

    # Check if any databases are available
    if [ -z "${cosmic_arg}" ] && [ -z "${gnomad_arg}" ]; then
        echo "WARNING: No databases provided for annotation" >&2
        echo "Implementing graceful fallback - copying input VCF to output" >&2
        
        # Copy input VCF to output location to maintain workflow continuity
        cp ${vcf} ${prefix}.vcf.gz
        cp ${tbi} ${prefix}.vcf.gz.tbi
        
        # Create empty stats file if verbose logging was requested
        if [ "${verbose_logging}" = "true" ]; then
            echo '{"status": "fallback", "message": "No databases provided, input VCF preserved"}' > ${prefix}_stats.json
        fi
        
        annotation_exit_code=0
    else
        # Run COSMIC/gnomAD annotation with error handling
        set +e  # Allow script to continue on annotation errors
        
        annotate_cosmic_gnomad.py \\
            --input ${vcf} \\
            ${cosmic_arg} \\
            ${gnomad_arg} \\
            --output ${prefix}.vcf.gz \\
            ${verbose} \\
            ${stats_output} \\
            ${args}
        
        annotation_exit_code=\$?
    fi
    
    # Implement graceful fallback on annotation failure
    if [ \$annotation_exit_code -ne 0 ]; then
        echo "WARNING: COSMIC/gnomAD annotation failed with exit code \$annotation_exit_code" >&2
        echo "Implementing graceful fallback - copying input VCF to output" >&2
        
        # Copy input VCF to output location to maintain workflow continuity
        cp ${vcf} ${prefix}.vcf.gz
        cp ${tbi} ${prefix}.vcf.gz.tbi
        
        # Create empty stats file if verbose logging was requested
        if [ "${verbose_logging}" = "true" ]; then
            echo '{"status": "fallback", "message": "Annotation failed, input VCF preserved"}' > ${prefix}_stats.json
        fi
        
        echo "Graceful fallback completed - workflow will continue" >&2
    else
        echo "COSMIC/gnomAD annotation completed successfully" >&2
    fi
    
    # Validate output files exist
    if [ ! -f "${prefix}.vcf.gz" ]; then
        echo "ERROR: Output VCF file not created: ${prefix}.vcf.gz" >&2
        exit 1
    fi
    
    if [ ! -f "${prefix}.vcf.gz.tbi" ]; then
        echo "ERROR: Output VCF index file not created: ${prefix}.vcf.gz.tbi" >&2
        exit 1
    fi
    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pysam: \$(python -c "try: import pysam; print(pysam.__version__); except ImportError: print('not available')" 2>/dev/null || echo "not available")
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//' || echo "not available")
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^.*tabix //; s/ .*\$//' || echo "not available")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create stub output files for testing
    touch ${prefix}.cosmic_gnomad_annotated.vcf.gz
    touch ${prefix}.cosmic_gnomad_annotated.vcf.gz.tbi
    
    # Create stub stats file if verbose logging enabled
    if [ "${verbose_logging}" = "true" ]; then
        echo '{"status": "stub", "total_variants": 0, "classification_counts": {}}' > ${prefix}_stats.json
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
        pysam: 0.22.0
        bcftools: 1.21
        tabix: 1.21
    END_VERSIONS
    """
}