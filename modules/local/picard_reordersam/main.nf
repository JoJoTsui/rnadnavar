process PICARD_REORDERSAM {
    tag "$meta.id"
    label 'process_reordersam'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.4.0--hdfd78af_0' :
        'biocontainers/picard:3.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(meta), path("*.reordered.bam"), path("*.reordered.bam.bai"), emit: bam
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Use 80% of available memory for JVM heap
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[PICARD ReorderSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    // Calculate MAX_RECORDS_IN_RAM based on available memory
    // ~1.5KB per record, use ~60% of heap for records buffer
    def max_records = ((avail_mem * 0.6 * 1024) / 1.5).intValue()
    // Parallel GC threads based on available CPUs
    def gc_threads = Math.min(task.cpus ?: 4, 8)
    """
    picard \\
        -Xmx${avail_mem}M \\
        -XX:ParallelGCThreads=${gc_threads} \\
        ReorderSam \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.reordered.bam \\
        --SEQUENCE_DICTIONARY ${dict} \\
        --ALLOW_INCOMPLETE_DICT_CONCORDANCE true \\
        --MAX_RECORDS_IN_RAM ${max_records} \\
        --COMPRESSION_LEVEL 1 \\
        --VALIDATION_STRINGENCY SILENT \\
        --CREATE_INDEX true \\
        $args

    # Rename index to match expected naming convention
    mv ${prefix}.reordered.bai ${prefix}.reordered.bam.bai || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard ReorderSam --version 2>&1 | grep -o 'Version:.*' | cut -f2 -d: | tr -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reordered.bam
    touch ${prefix}.reordered.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard ReorderSam --version 2>&1 | grep -o 'Version:.*' | cut -f2 -d: | tr -d ' ')
    END_VERSIONS
    """
}

/*
 * OPTIMIZATION NOTES (Picard 3.4.0 Verified):
 * 
 * 1. --ALLOW_INCOMPLETE_DICT_CONCORDANCE true
 *    Enables filtering of reads from contigs not in reference (instead of crashing)
 * 
 * 2. --MAX_RECORDS_IN_RAM (dynamically calculated)
 *    Default is only 500,000. Increased based on available memory for speed.
 * 
 * 3. --COMPRESSION_LEVEL 1
 *    Default is 5. Lower compression = faster writes (larger output file)
 * 
 * 4. --VALIDATION_STRINGENCY SILENT
 *    Per manual: "Setting stringency to SILENT can improve performance"
 * 
 * 5. -XX:ParallelGCThreads
 *    Optimizes JVM garbage collection for multi-core systems
 * 
 * WARNING: Do NOT use --ALLOW_CONTIG_LENGTH_DISCORDANCE (default: false)
 *          Per manual: "Highly dangerous" - only for mismatched assemblies
 */
