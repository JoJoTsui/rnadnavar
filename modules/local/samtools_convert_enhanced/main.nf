process SAMTOOLS_CONVERT_ENHANCED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.bam")  , emit: bam ,   optional: true
    tuple val(meta), path("*.cram") , emit: cram,   optional: true
    tuple val(meta), path("*.bai")  , emit: bai ,   optional: true
    tuple val(meta), path("*.crai") , emit: crai,   optional: true
    path  "versions.yml"            , emit: versions

    when:
    // Enhanced when condition handling with null safety
    task.ext.when == null || task.ext.when == true

    script:
    // Input validation before processing
    if (!meta) {
        error "SAMTOOLS_CONVERT_ENHANCED: meta is null or empty"
    }
    if (!meta.id) {
        error "SAMTOOLS_CONVERT_ENHANCED: meta.id is required but missing"
    }
    if (!input) {
        error "SAMTOOLS_CONVERT_ENHANCED: input file is null"
    }
    if (!input.exists()) {
        error "SAMTOOLS_CONVERT_ENHANCED: input file does not exist: ${input}"
    }
    if (!index) {
        error "SAMTOOLS_CONVERT_ENHANCED: index file is null"
    }
    if (!index.exists()) {
        error "SAMTOOLS_CONVERT_ENHANCED: index file does not exist: ${index}"
    }
    if (!fasta) {
        error "SAMTOOLS_CONVERT_ENHANCED: fasta reference is null"
    }
    if (!fasta.exists()) {
        error "SAMTOOLS_CONVERT_ENHANCED: fasta reference does not exist: ${fasta}"
    }
    if (!fai) {
        error "SAMTOOLS_CONVERT_ENHANCED: fasta index is null"
    }
    if (!fai.exists()) {
        error "SAMTOOLS_CONVERT_ENHANCED: fasta index does not exist: ${fai}"
    }

    // Validate input file format
    def input_extension = input.getExtension().toLowerCase()
    if (!(input_extension in ['bam', 'cram'])) {
        error "SAMTOOLS_CONVERT_ENHANCED: unsupported input format: ${input_extension}. Expected: bam or cram"
    }

    // Validate index file matches input file
    def expected_index_ext = input_extension == 'bam' ? 'bai' : 'crai'
    def actual_index_ext = index.getExtension().toLowerCase()
    if (actual_index_ext != expected_index_ext) {
        error "SAMTOOLS_CONVERT_ENHANCED: index file extension mismatch. Expected: ${expected_index_ext}, got: ${actual_index_ext}"
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_extension = input_extension == "bam" ? "cram" : "bam"

    // Log conversion details for debugging
    log.info "SAMTOOLS_CONVERT_ENHANCED: Converting ${input.getSimpleName()} (${input_extension}) to ${output_extension} for sample ${meta.id}"

    """
    # Validate samtools is available
    if ! command -v samtools &> /dev/null; then
        echo "ERROR: samtools command not found" >&2
        exit 1
    fi

    # Validate input files are readable
    if [ ! -r "${input}" ]; then
        echo "ERROR: Cannot read input file: ${input}" >&2
        exit 1
    fi

    if [ ! -r "${index}" ]; then
        echo "ERROR: Cannot read index file: ${index}" >&2
        exit 1
    fi

    if [ ! -r "${fasta}" ]; then
        echo "ERROR: Cannot read reference fasta: ${fasta}" >&2
        exit 1
    fi

    if [ ! -r "${fai}" ]; then
        echo "ERROR: Cannot read fasta index: ${fai}" >&2
        exit 1
    fi

    # Check available disk space (warn if less than 2GB)
    available_space=\$(df . | tail -1 | awk '{print \$4}')
    if [ "\$available_space" -lt 2097152 ]; then
        echo "WARNING: Low disk space available: \${available_space}KB" >&2
    fi

    # Perform conversion with error handling
    echo "Starting conversion: ${input} -> ${prefix}.${output_extension}"
    
    samtools view \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        $args \\
        $input \\
        -o ${prefix}.${output_extension}
    
    conversion_exit_code=\$?
    if [ \$conversion_exit_code -ne 0 ]; then
        echo "ERROR: samtools view failed with exit code \$conversion_exit_code" >&2
        exit \$conversion_exit_code
    fi

    # Verify output file was created and is not empty
    if [ ! -f "${prefix}.${output_extension}" ]; then
        echo "ERROR: Output file was not created: ${prefix}.${output_extension}" >&2
        exit 1
    fi

    if [ ! -s "${prefix}.${output_extension}" ]; then
        echo "ERROR: Output file is empty: ${prefix}.${output_extension}" >&2
        exit 1
    fi

    echo "Conversion completed successfully. Creating index..."

    # Create index with error handling
    samtools index -@${task.cpus} ${prefix}.${output_extension}
    
    index_exit_code=\$?
    if [ \$index_exit_code -ne 0 ]; then
        echo "ERROR: samtools index failed with exit code \$index_exit_code" >&2
        exit \$index_exit_code
    fi

    # Verify index file was created
    index_extension=\$(if [ "${output_extension}" = "bam" ]; then echo "bai"; else echo "crai"; fi)
    if [ ! -f "${prefix}.${output_extension}.\$index_extension" ]; then
        echo "ERROR: Index file was not created: ${prefix}.${output_extension}.\$index_extension" >&2
        exit 1
    fi

    echo "Indexing completed successfully."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_extension = input.getExtension().toLowerCase()
    def output_extension = input_extension == "bam" ? "cram" : "bam"
    def index_extension = output_extension == "bam" ? "bai" : "crai"

    """
    touch ${prefix}.${output_extension}
    touch ${prefix}.${output_extension}.${index_extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}