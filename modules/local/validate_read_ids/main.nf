process VALIDATE_READ_IDS_AUDITED {
    tag "$meta.id"
    label 'process_single'
    input:
    tuple val(meta), path(read_ids)

    output:
    tuple val(meta), path('*.txt'), emit: read_ids
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    test -s ${read_ids} || { echo "nonempty realignment candidates yielded zero usable paired read IDs for ${meta.id}" >&2; exit 2; }
    cp ${read_ids} ${prefix}.read_ids.txt
    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        coreutils: bundled
    END_VERSIONS
    """
}
