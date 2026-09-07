process FILTER_HISAT_SPLICESITES {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(splicesites), path(fai)

    output:
    tuple val(meta), path('*.splice_sites.txt'), emit: splicesites
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    awk 'BEGIN{FS=OFS="\\t"}
      FNR==NR { contig_len[\$1]=\$2; next }
      NF < 4 { bad++; next }
      !(\$1 in contig_len) { excluded++; next }
      \$2 !~ /^[0-9]+\$/ || \$3 !~ /^[0-9]+\$/ || \$2 < 0 || \$2 >= \$3 || \$3 > contig_len[\$1] { bad++; next }
      { print }
      END { if (bad) { printf("malformed or out-of-range splice rows: %d\\n", bad) > "/dev/stderr"; exit 2 } }' \\
      ${fai} ${splicesites} > ${prefix}.splice_sites.txt
    test -s ${prefix}.splice_sites.txt || { echo 'no splice sites remain after reference audit' >&2; exit 2; }

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        awk: coreutils
    END_VERSIONS
    """
}
