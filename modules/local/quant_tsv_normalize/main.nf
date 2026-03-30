process QUANT_TSV_NORMALIZE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta), path(quant_sf)

    output:
    tuple val(meta), path("quant.tsv"), emit: quant_tsv
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 <<-EOF
    import sys

    with open('${quant_sf}') as f:
        lines = f.readlines()

    if len(lines) < 2:
        print("ERROR: quant.sf is empty or has fewer than 2 lines", file=sys.stderr)
        sys.exit(1)

    with open('quant.tsv', 'w') as out:
        for line in lines:
            parts = line.rstrip('\\n').split('\\t')
            if parts[0] != 'Name':  # skip header
                parts[0] = parts[0].split('|')[0]
            out.write('\\t'.join(parts) + '\\n')
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
