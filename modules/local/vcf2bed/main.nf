process VCF2BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path('*.bed'), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools view \\
        -O v \\
        -o - ${vcf} | \\
    awk 'BEGIN{OFS="\\t"} !/^#/ {
        start = \$2 - 1
        end = \$2 + length(\$4) - 1
        print \$1, start, end
    }' > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
