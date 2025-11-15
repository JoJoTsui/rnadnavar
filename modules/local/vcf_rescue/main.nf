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
          path(dna_vcfs), path(dna_tbis),
          path(rna_vcfs), path(rna_tbis)

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

    """
    mkdir -p dna_vcfs rna_vcfs
    
    # Link DNA caller VCFs
    for vcf in ${dna_vcfs}; do
        ln -s \$(readlink -f \$vcf) dna_vcfs/
    done
    for tbi in ${dna_tbis}; do
        ln -s \$(readlink -f \$tbi) dna_vcfs/
    done
    
    # Link RNA caller VCFs
    for vcf in ${rna_vcfs}; do
        ln -s \$(readlink -f \$vcf) rna_vcfs/
    done
    for tbi in ${rna_tbis}; do
        ln -s \$(readlink -f \$tbi) rna_vcfs/
    done
    
    run_rescue_vcf.py \\
        --dna_consensus ${dna_consensus_vcf} \\
        --rna_consensus ${rna_consensus_vcf} \\
        --dna_vcfs dna_vcfs/ \\
        --rna_vcfs rna_vcfs/ \\
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
