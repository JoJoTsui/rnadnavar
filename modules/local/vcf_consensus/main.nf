process VCF_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis), val(callers)

    output:
    tuple val(meta), path("*.consensus.vcf.gz"), path("*.consensus.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def snv_thr = task.ext.snv_thr ?: 2
    def indel_thr = task.ext.indel_thr ?: 2

    """
    mkdir -p inputs
    
    # Link VCF files
    for vcf in ${vcfs}; do
        ln -s \$(readlink -f \$vcf) inputs/
    done
    
    # Link index files
    for tbi in ${tbis}; do
        tbi_basename=\$(basename \$tbi)
        if [ ! -e "inputs/\${tbi_basename}" ]; then
            ln -s \$(readlink -f \$tbi) inputs/
        fi
    done
    
    run_consensus_vcf.py \\
        --input_dir inputs/ \\
        --out_prefix ${prefix}.consensus \\
        --snv_thr ${snv_thr} \\
        --indel_thr ${indel_thr}
    
    # Compress and index
    # bgzip -c ${prefix}.consensus.vcf > ${prefix}.consensus.vcf.gz
    tabix -p vcf ${prefix}.consensus.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        cyvcf2: \$(python -c "import cyvcf2; print(cyvcf2.__version__)")
    END_VERSIONS
    """
}
