process FILTER_HISAT_SPLICESITES {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(splicesites), path(fai)

    output:
    tuple val(meta), path('*.splice_sites.txt'), emit: splicesites
    tuple val(meta), path('*.splice_sites.audit.tsv'), emit: audit
    tuple val(meta), path('*.splice_sites.excluded.tsv'), emit: excluded
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    awk -v excluded_file=${prefix}.splice_sites.excluded.tsv 'BEGIN{FS=OFS="\\t"; print "reason\\trow" > excluded_file}
      FNR==NR { contig_len[\$1]=\$2; next }
      NF < 4 { print "malformed", \$0 > excluded_file; bad++; next }
      \$2 !~ /^[0-9]+\$/ || \$3 !~ /^[0-9]+\$/ || \$2 < 0 || \$2 >= \$3 { print "invalid_coordinates", \$0 > excluded_file; bad++; next }
      !(\$1 in contig_len) { print "unknown_contig", \$0 > excluded_file; excluded++; next }
      \$3 > contig_len[\$1] { print "out_of_range", \$0 > excluded_file; bad++; next }
      { print }
      END { if (bad) { printf("malformed or out-of-range splice rows: %d\\n", bad) > "/dev/stderr"; exit 2 } }' \\
      ${fai} ${splicesites} > ${prefix}.splice_sites.txt
    test -s ${prefix}.splice_sites.txt || { echo 'no splice sites remain after reference audit' >&2; exit 2; }
    printf 'source\\tderivative\\texcluded_unknown_contig\\n' > ${prefix}.splice_sites.audit.tsv
    printf '%s\\t%s\\t' "\$(sha256sum ${splicesites} | cut -d' ' -f1)" "\$(sha256sum ${prefix}.splice_sites.txt | cut -d' ' -f1)" >> ${prefix}.splice_sites.audit.tsv
    awk 'BEGIN{FS="\\t"} FNR==NR { contig[\$1]=1; next } NF >= 1 && !(\$1 in contig) { n++ } END { print n+0 }' ${fai} ${splicesites} >> ${prefix}.splice_sites.audit.tsv

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        awk: coreutils
    END_VERSIONS
    """
}
