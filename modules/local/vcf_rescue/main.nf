process VCF_RESCUE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta),
          path(dna_consensus_vcf), path(dna_consensus_tbi),
          path(rna_consensus_vcf), path(rna_consensus_tbi),
          path(dna_caller_vcfs), path(dna_caller_tbis), val(dna_callers),
          path(rna_caller_vcfs), path(rna_caller_tbis), val(rna_callers)

    output:
    tuple val(meta), path("*.rescued.vcf.gz"), path("*.rescued.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def snv_thr = task.ext.snv_thr ?: 2
    def indel_thr = task.ext.indel_thr ?: 2
    
    // Build DNA caller VCF arguments
    def dna_vcf_args = ""
    if (dna_caller_vcfs && dna_caller_vcfs.size() > 0) {
        dna_vcf_args = dna_caller_vcfs.collect { "--dna_vcf ${it}" }.join(' ')
    }
    
    // Build RNA caller VCF arguments
    def rna_vcf_args = ""
    if (rna_caller_vcfs && rna_caller_vcfs.size() > 0) {
        rna_vcf_args = rna_caller_vcfs.collect { "--rna_vcf ${it}" }.join(' ')
    }

    """
    run_rescue_vcf.py \\
        --dna_consensus ${dna_consensus_vcf} \\
        --rna_consensus ${rna_consensus_vcf} \\
        ${dna_vcf_args} \\
        ${rna_vcf_args} \\
        --out_prefix ${prefix}.rescued \\
        --snv_thr ${snv_thr} \\
        --indel_thr ${indel_thr}
    
    tabix -p vcf ${prefix}.rescued.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        cyvcf2: \$(python -c "import cyvcf2; print(cyvcf2.__version__)")
    END_VERSIONS
    """
}
