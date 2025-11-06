process DEEPSOMATIC {
    tag "$meta.id"
    label 'process_high'

    // Use local installation instead of container
    // conda "${moduleDir}/environment.yml"
    
    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor)
    tuple val(meta2), path(intervals)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")      ,  emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi")  ,  emit: vcf_tbi
    tuple val(meta), path("${prefix}.g.vcf.gz")    ,  emit: gvcf
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"),  emit: gvcf_tbi
    path "versions.yml"                            ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--regions=${intervals}" : ""
    def VERSION = '1.9.0'
    
    // Configure DeepSomatic binary path and model type
    def deepsomatic_bin = params.deepsomatic_bin_path ?: '/opt/deepvariant/bin/deepsomatic/run_deepsomatic'
    def model_type = params.deepsomatic_model_type ?: (params.wes ? 'WES' : 'WGS')
    // def intermediate_dir = params.deepsomatic_intermediate_results_dir ?: 'tmp'
    
    // Handle sample names from metadata or use defaults
    def sample_name_tumor = meta.tumor_id ?: meta.id + "_tumor"
    def sample_name_normal = meta.normal_id ?: meta.id + "_normal"

    """
    ${deepsomatic_bin} \\
        --model_type=${model_type} \\
        --ref=${fasta} \\
        --reads_normal=${input_normal} \\
        --reads_tumor=${input_tumor} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        --sample_name_tumor="${sample_name_tumor}" \\
        --sample_name_normal="${sample_name_normal}" \\
        --num_shards=${task.cpus} \\
        ${regions} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepsomatic: $VERSION
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.9.0'
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepsomatic: $VERSION
    END_VERSIONS
    """
}
