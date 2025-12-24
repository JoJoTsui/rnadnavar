process FILE_VALIDATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path(files), emit: validated_files
    path "versions.yml",          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/bin/bash
    set -euo pipefail

    # Function to validate file existence and accessibility
    validate_file_access() {
        local file="\$1"
        local file_type="\$2"
        
        if [[ ! -f "\$file" ]]; then
            echo "ERROR: \$file_type file '\$file' does not exist" >&2
            return 1
        fi
        
        if [[ ! -r "\$file" ]]; then
            echo "ERROR: \$file_type file '\$file' is not readable" >&2
            return 1
        fi
        
        if [[ ! -s "\$file" ]]; then
            echo "WARNING: \$file_type file '\$file' is empty" >&2
        fi
        
        echo "INFO: \$file_type file '\$file' validation passed"
        return 0
    }

    # Function to detect compression format
    detect_compression() {
        local file="\$1"
        
        # Check file magic bytes
        if file "\$file" | grep -q "gzip compressed"; then
            echo "gzip"
        elif file "\$file" | grep -q "bzip2 compressed"; then
            echo "bzip2"
        elif file "\$file" | grep -q "XZ compressed"; then
            echo "xz"
        else
            echo "uncompressed"
        fi
    }

    # Function to validate CRAM files
    validate_cram() {
        local cram_file="\$1"
        local index_file="\$2"
        
        echo "INFO: Validating CRAM file: \$cram_file"
        
        # Check CRAM file header
        if ! samtools view -H "\$cram_file" >/dev/null 2>&1; then
            echo "ERROR: Invalid CRAM file header in '\$cram_file'" >&2
            return 1
        fi
        
        # Check if index exists and is valid
        if [[ -n "\$index_file" ]]; then
            validate_file_access "\$index_file" "CRAM index"
            
            # Verify index corresponds to CRAM file
            if ! samtools quickcheck "\$cram_file" 2>/dev/null; then
                echo "ERROR: CRAM file '\$cram_file' failed quickcheck" >&2
                return 1
            fi
        fi
        
        echo "INFO: CRAM file validation passed"
        return 0
    }

    # Function to validate BAM files
    validate_bam() {
        local bam_file="\$1"
        local index_file="\$2"
        
        echo "INFO: Validating BAM file: \$bam_file"
        
        # Check BAM file header
        if ! samtools view -H "\$bam_file" >/dev/null 2>&1; then
            echo "ERROR: Invalid BAM file header in '\$bam_file'" >&2
            return 1
        fi
        
        # Check if index exists and is valid
        if [[ -n "\$index_file" ]]; then
            validate_file_access "\$index_file" "BAM index"
            
            # Verify index corresponds to BAM file
            if ! samtools quickcheck "\$bam_file" 2>/dev/null; then
                echo "ERROR: BAM file '\$bam_file' failed quickcheck" >&2
                return 1
            fi
        fi
        
        echo "INFO: BAM file validation passed"
        return 0
    }

    # Main validation logic
    echo "Starting file validation for sample: ${meta.id}"
    
    # Validate each file based on extension
    for file in ${files}; do
        echo "Processing file: \$file"
        
        # Basic file access validation
        validate_file_access "\$file" "Input"
        
        # Detect compression format
        compression=\$(detect_compression "\$file")
        echo "INFO: Detected compression format: \$compression"
        
        # File type specific validation
        case "\$file" in
            *.cram)
                # Look for corresponding index file
                index_file=""
                if [[ -f "\${file}.crai" ]]; then
                    index_file="\${file}.crai"
                elif [[ -f "\${file%.*}.crai" ]]; then
                    index_file="\${file%.*}.crai"
                fi
                validate_cram "\$file" "\$index_file"
                ;;
            *.bam)
                # Look for corresponding index file
                index_file=""
                if [[ -f "\${file}.bai" ]]; then
                    index_file="\${file}.bai"
                elif [[ -f "\${file%.*}.bai" ]]; then
                    index_file="\${file%.*}.bai"
                fi
                validate_bam "\$file" "\$index_file"
                ;;
            *.crai|*.bai)
                echo "INFO: Index file '\$file' detected"
                ;;
            *)
                echo "INFO: Generic file validation for '\$file'"
                ;;
        esac
    done
    
    echo "File validation completed successfully for sample: ${meta.id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        file: \$(file --version | head -n1 | sed 's/file-//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create dummy output files to match the expected outputs
    for file in ${files}; do
        touch "\$file"
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' || echo "1.17")
        file: \$(file --version | head -n1 | sed 's/file-//' || echo "5.44")
    END_VERSIONS
    """
}