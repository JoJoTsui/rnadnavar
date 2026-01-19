process CHECK_CONTIGS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(meta), path(bam), path(bai), env(NEEDS_REORDER), emit: bam_with_status
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract contigs from BAM header
    samtools view -H ${bam} | grep "^@SQ" | cut -f2 | sed 's/SN://' | sort > bam_contigs.txt

    # Extract contigs from reference dict
    grep "^@SQ" ${dict} | cut -f2 | sed 's/SN://' | sort > ref_contigs.txt

    # Find contigs in BAM that are not in reference
    EXTRA_CONTIGS=\$(comm -23 bam_contigs.txt ref_contigs.txt | wc -l)

    if [ "\$EXTRA_CONTIGS" -gt 0 ]; then
        echo "Found \$EXTRA_CONTIGS contigs in BAM not present in reference"
        comm -23 bam_contigs.txt ref_contigs.txt | head -20
        NEEDS_REORDER="true"
    else
        echo "All BAM contigs are present in reference"
        NEEDS_REORDER="false"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    NEEDS_REORDER="false"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
