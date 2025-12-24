process SAMTOOLS_CONVERT_ENHANCED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(fasta_meta), path(fasta)
    tuple val(fai_meta), path(fai)

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.sam"), emit: sam, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.csi"), emit: csi, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input.getExtension()
    def reference = fasta && fai ? "--reference ${fasta}" : ""
    
    // Enhanced error handling and validation
    def validation_log = "validation_${prefix}.log"
    
    // Determine output format and file extension
    def output_fmt = ""
    def output_ext = ""
    def index_ext = ""
    
    if (file_type == "bam") {
        output_fmt = "CRAM"
        output_ext = "cram"
        index_ext = "crai"
    } else if (file_type == "cram") {
        output_fmt = "BAM"
        output_ext = "bam"
        index_ext = "bai"
    } else if (file_type == "sam") {
        output_fmt = "BAM"
        output_ext = "bam"
        index_ext = "bai"
    } else {
        error "Unsupported input file type: ${file_type}"
    }
    
    """
    # Enhanced SAMTOOLS_CONVERT with comprehensive validation and error handling
    
    echo "=== SAMTOOLS_CONVERT_ENHANCED: Starting conversion ===" > ${validation_log}
    echo "Input file: ${input}" >> ${validation_log}
    echo "Input type: ${file_type}" >> ${validation_log}
    echo "Output format: ${output_fmt}" >> ${validation_log}
    echo "Reference: ${fasta ?: 'None'}" >> ${validation_log}
    echo "Sample ID: ${meta.id}" >> ${validation_log}
    echo "Patient ID: ${meta.patient}" >> ${validation_log}
    echo "Status: ${meta.status}" >> ${validation_log}
    echo "Timestamp: \$(date)" >> ${validation_log}
    
    # Pre-execution validation
    echo "=== Pre-execution validation ===" >> ${validation_log}
    
    # Validate input file exists and is readable
    if [[ ! -f "${input}" ]]; then
        echo "ERROR: Input file does not exist: ${input}" >> ${validation_log}
        exit 1
    fi
    
    if [[ ! -r "${input}" ]]; then
        echo "ERROR: Input file is not readable: ${input}" >> ${validation_log}
        exit 1
    fi
    
    # Check input file size
    input_size=\$(stat -c%s "${input}")
    echo "Input file size: \${input_size} bytes" >> ${validation_log}
    
    if [[ \${input_size} -lt 1000 ]]; then
        echo "WARNING: Input file is very small (\${input_size} bytes)" >> ${validation_log}
    fi
    
    # Validate index file if provided
    if [[ -n "${index}" && "${index}" != "null" ]]; then
        if [[ ! -f "${index}" ]]; then
            echo "ERROR: Index file does not exist: ${index}" >> ${validation_log}
            exit 1
        fi
        
        if [[ ! -r "${index}" ]]; then
            echo "ERROR: Index file is not readable: ${index}" >> ${validation_log}
            exit 1
        fi
        
        index_size=\$(stat -c%s "${index}")
        echo "Index file size: \${index_size} bytes" >> ${validation_log}
    fi
    
    # Validate reference files if provided
    if [[ -n "${fasta}" && "${fasta}" != "null" ]]; then
        if [[ ! -f "${fasta}" ]]; then
            echo "ERROR: Reference FASTA does not exist: ${fasta}" >> ${validation_log}
            exit 1
        fi
        
        if [[ ! -r "${fasta}" ]]; then
            echo "ERROR: Reference FASTA is not readable: ${fasta}" >> ${validation_log}
            exit 1
        fi
        
        fasta_size=\$(stat -c%s "${fasta}")
        echo "Reference FASTA size: \${fasta_size} bytes" >> ${validation_log}
    fi
    
    if [[ -n "${fai}" && "${fai}" != "null" ]]; then
        if [[ ! -f "${fai}" ]]; then
            echo "ERROR: Reference FASTA index does not exist: ${fai}" >> ${validation_log}
            exit 1
        fi
        
        if [[ ! -r "${fai}" ]]; then
            echo "ERROR: Reference FASTA index is not readable: ${fai}" >> ${validation_log}
            exit 1
        fi
        
        fai_size=\$(stat -c%s "${fai}")
        echo "Reference FASTA index size: \${fai_size} bytes" >> ${validation_log}
    fi
    
    # Check available disk space
    available_space=\$(df -B1 . | tail -1 | awk '{print \$4}')
    echo "Available disk space: \${available_space} bytes" >> ${validation_log}
    
    # Estimate required space (conservative estimate: 2x input size)
    required_space=\$((input_size * 2))
    echo "Estimated required space: \${required_space} bytes" >> ${validation_log}
    
    if [[ \${available_space} -lt \${required_space} ]]; then
        echo "WARNING: Available disk space may be insufficient" >> ${validation_log}
        echo "Available: \${available_space}, Required: \${required_space}" >> ${validation_log}
    fi
    
    # Check memory usage before starting
    echo "Memory usage before conversion:" >> ${validation_log}
    free -h >> ${validation_log} 2>/dev/null || echo "Memory info not available" >> ${validation_log}
    
    echo "=== Starting samtools conversion ===" >> ${validation_log}
    
    # Construct and validate samtools command
    samtools_cmd="samtools view ${args} ${reference} -o ${prefix}.${output_ext} ${input}"
    echo "Command: \${samtools_cmd}" >> ${validation_log}
    
    # Execute samtools conversion with error handling
    set +e  # Don't exit on error immediately
    
    # Run the conversion
    eval "\${samtools_cmd}" 2>&1 | tee -a ${validation_log}
    conversion_exit_code=\${PIPESTATUS[0]}
    
    echo "Conversion exit code: \${conversion_exit_code}" >> ${validation_log}
    
    if [[ \${conversion_exit_code} -ne 0 ]]; then
        echo "ERROR: samtools conversion failed with exit code \${conversion_exit_code}" >> ${validation_log}
        echo "=== Error diagnostics ===" >> ${validation_log}
        
        # Capture additional diagnostic information
        echo "Current directory contents:" >> ${validation_log}
        ls -la >> ${validation_log}
        
        echo "Memory usage after failed conversion:" >> ${validation_log}
        free -h >> ${validation_log} 2>/dev/null || echo "Memory info not available" >> ${validation_log}
        
        echo "Disk usage:" >> ${validation_log}
        df -h . >> ${validation_log}
        
        # Check if partial output was created
        if [[ -f "${prefix}.${output_ext}" ]]; then
            partial_size=\$(stat -c%s "${prefix}.${output_ext}")
            echo "Partial output file created: \${partial_size} bytes" >> ${validation_log}
            rm -f "${prefix}.${output_ext}"  # Clean up partial file
        fi
        
        exit \${conversion_exit_code}
    fi
    
    set -e  # Re-enable exit on error
    
    # Post-conversion validation
    echo "=== Post-conversion validation ===" >> ${validation_log}
    
    # Validate output file was created
    if [[ ! -f "${prefix}.${output_ext}" ]]; then
        echo "ERROR: Output file was not created: ${prefix}.${output_ext}" >> ${validation_log}
        exit 1
    fi
    
    # Check output file size
    output_size=\$(stat -c%s "${prefix}.${output_ext}")
    echo "Output file size: \${output_size} bytes" >> ${validation_log}
    
    if [[ \${output_size} -lt 1000 ]]; then
        echo "WARNING: Output file is very small (\${output_size} bytes)" >> ${validation_log}
    fi
    
    # Validate output file format using samtools
    echo "Validating output file format..." >> ${validation_log}
    samtools quickcheck "${prefix}.${output_ext}" 2>&1 | tee -a ${validation_log}
    quickcheck_exit_code=\${PIPESTATUS[0]}
    
    if [[ \${quickcheck_exit_code} -ne 0 ]]; then
        echo "ERROR: Output file failed format validation" >> ${validation_log}
        exit 1
    fi
    
    echo "Output file format validation passed" >> ${validation_log}
    
    # Create index if output is BAM or CRAM
    if [[ "${output_ext}" == "bam" || "${output_ext}" == "cram" ]]; then
        echo "Creating index for ${output_ext} file..." >> ${validation_log}
        
        index_cmd="samtools index ${prefix}.${output_ext}"
        echo "Index command: \${index_cmd}" >> ${validation_log}
        
        eval "\${index_cmd}" 2>&1 | tee -a ${validation_log}
        index_exit_code=\${PIPESTATUS[0]}
        
        if [[ \${index_exit_code} -ne 0 ]]; then
            echo "ERROR: Index creation failed with exit code \${index_exit_code}" >> ${validation_log}
            exit \${index_exit_code}
        fi
        
        # Validate index file was created
        if [[ ! -f "${prefix}.${output_ext}.${index_ext}" ]]; then
            echo "ERROR: Index file was not created: ${prefix}.${output_ext}.${index_ext}" >> ${validation_log}
            exit 1
        fi
        
        index_size=\$(stat -c%s "${prefix}.${output_ext}.${index_ext}")
        echo "Index file size: \${index_size} bytes" >> ${validation_log}
        
        echo "Index creation completed successfully" >> ${validation_log}
    fi
    
    # Final memory and disk usage check
    echo "=== Final resource usage ===" >> ${validation_log}
    echo "Memory usage after conversion:" >> ${validation_log}
    free -h >> ${validation_log} 2>/dev/null || echo "Memory info not available" >> ${validation_log}
    
    echo "Disk usage after conversion:" >> ${validation_log}
    df -h . >> ${validation_log}
    
    echo "Final directory contents:" >> ${validation_log}
    ls -la >> ${validation_log}
    
    echo "=== SAMTOOLS_CONVERT_ENHANCED: Conversion completed successfully ===" >> ${validation_log}
    echo "Conversion completed at: \$(date)" >> ${validation_log}
    
    # Display validation log for debugging
    echo "=== Validation Log Summary ==="
    tail -20 ${validation_log}
    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input.getExtension()
    
    def output_fmt = ""
    def output_ext = ""
    
    if (file_type == "bam") {
        output_fmt = "CRAM"
        output_ext = "cram"
    } else if (file_type == "cram") {
        output_fmt = "BAM"
        output_ext = "bam"
    } else if (file_type == "sam") {
        output_fmt = "BAM"
        output_ext = "bam"
    }
    
    """
    touch ${prefix}.${output_ext}
    touch ${prefix}.${output_ext}.bai
    touch ${prefix}.${output_ext}.crai
    touch ${prefix}.${output_ext}.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' || echo "1.19.2")
    END_VERSIONS
    """
}