process SORT_MERGE_BED_AUDITED {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path('*.bed'), emit: bed
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n ${input} \
      | awk 'BEGIN{OFS="\\t"} { if (NF < 3) next; if (seen && \$1 == chrom && \$2 <= end) { if (\$3 > end) end=\$3 } else { if (seen) print chrom, start, end; chrom=\$1; start=\$2; end=\$3; seen=1 } } END { if (seen) print chrom, start, end }' \
      > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        sort: coreutils
    END_VERSIONS
    """
}
