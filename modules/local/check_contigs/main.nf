process CHECK_CONTIGS {
    tag "$meta.id"
    label 'process_alignment_dictionary'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(alignment), path(index)
    path(fasta)
    path(fai)
    path(dict)
    val(policy)
    val(stage)

    output:
    tuple val(meta), path(alignment), path(index), env(NEEDS_NORMALIZATION), emit: alignment_with_status
    tuple val(meta), path("*.dictionary_audit.json")                         , emit: audit
    path "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${stage}"
    """
    samtools view -H ${alignment} > ${prefix}.alignment_header.sam
    samtools idxstats ${alignment} > ${prefix}.idxstats.tsv

    NEEDS_NORMALIZATION=\$(check_alignment_dictionary.py \\
        --header ${prefix}.alignment_header.sam \\
        --idxstats ${prefix}.idxstats.tsv \\
        --reference-dict ${dict} \\
        --alignment ${alignment} \\
        --sample-id '${meta.id}' \\
        --stage ${stage} \\
        --policy ${policy} \\
        --output ${prefix}.dictionary_audit.json)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${stage}"
    """
    NEEDS_NORMALIZATION="false"
    echo '{"schema_version":1,"sample_id":"${meta.id}","stage":"${stage}","policy":"${policy}","decision":"compatible"}' > ${prefix}.dictionary_audit.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
