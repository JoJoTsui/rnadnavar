process VCF_FORMAT_VALIDATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0' :
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(hisat2_index, stageAs: "hisat2_index/*")

    output:
    tuple val(meta), path(vcf), path(tbi), emit: validated_vcf
    path "versions.yml",                   emit: versions
    path "*_validation_report.txt",        emit: validation_report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def hisat2_idx = hisat2_index ? "--hisat2-index ${hisat2_index}" : ""
    
    """
    #!/bin/bash
    set -euo pipefail

    # Initialize validation report
    REPORT="${prefix}_validation_report.txt"
    echo "VCF Format Validation Report for ${meta.id}" > \$REPORT
    echo "Generated on: \$(date)" >> \$REPORT
    echo "=========================================" >> \$REPORT
    echo "" >> \$REPORT

    # Function to log validation results
    log_result() {
        local test_name="\$1"
        local status="\$2"
        local message="\$3"
        
        echo "[\$status] \$test_name: \$message" | tee -a \$REPORT
    }

    # Function to validate VCF header
    validate_vcf_header() {
        local vcf_file="\$1"
        
        echo "INFO: Validating VCF header structure"
        
        # Check if file has proper VCF header
        if ! bcftools view -h "\$vcf_file" | head -1 | grep -q "^##fileformat=VCF"; then
            log_result "VCF Header" "FAIL" "Missing or invalid VCF fileformat declaration"
            return 1
        fi
        
        # Check for required header lines
        local header_lines=\$(bcftools view -h "\$vcf_file")
        
        # Check for contig information
        if ! echo "\$header_lines" | grep -q "^##contig="; then
            log_result "VCF Header" "WARN" "No contig information found in header"
        else
            local contig_count=\$(echo "\$header_lines" | grep -c "^##contig=")
            log_result "VCF Header" "PASS" "Found \$contig_count contig declarations"
        fi
        
        # Check for column header line
        if ! echo "\$header_lines" | tail -1 | grep -q "^#CHROM"; then
            log_result "VCF Header" "FAIL" "Missing column header line (#CHROM...)"
            return 1
        fi
        
        log_result "VCF Header" "PASS" "VCF header structure is valid"
        return 0
    }

    # Function to validate coordinate system
    validate_coordinate_system() {
        local vcf_file="\$1"
        
        echo "INFO: Validating coordinate system"
        
        # Check if VCF has any records
        local record_count=\$(bcftools view -H "\$vcf_file" | wc -l)
        if [[ \$record_count -eq 0 ]]; then
            log_result "Coordinate System" "WARN" "VCF file contains no variant records"
            return 0
        fi
        
        # Check coordinate consistency (positions should be positive integers)
        local invalid_coords=\$(bcftools view -H "\$vcf_file" | awk '\$2 !~ /^[1-9][0-9]*\$/ {print NR}' | wc -l)
        if [[ \$invalid_coords -gt 0 ]]; then
            log_result "Coordinate System" "FAIL" "Found \$invalid_coords records with invalid coordinates"
            return 1
        fi
        
        # Check if coordinates are sorted (basic check on first 1000 records)
        local unsorted_count=\$(bcftools view -H "\$vcf_file" | head -1000 | awk 'prev_chr==\$1 && prev_pos>\$2 {count++} {prev_chr=\$1; prev_pos=\$2} END {print count+0}')
        if [[ \$unsorted_count -gt 0 ]]; then
            log_result "Coordinate System" "WARN" "VCF may not be properly sorted (found \$unsorted_count unsorted positions in first 1000 records)"
        else
            log_result "Coordinate System" "PASS" "Coordinate system validation passed"
        fi
        
        return 0
    }

    # Function to detect and validate compression
    validate_compression() {
        local vcf_file="\$1"
        local tbi_file="\$2"
        
        echo "INFO: Validating compression format"
        
        # Detect compression format
        local compression=\$(file "\$vcf_file" | grep -o "gzip\\|bzip2\\|XZ" || echo "uncompressed")
        log_result "Compression" "INFO" "Detected format: \$compression"
        
        # If compressed, check for index
        if [[ "\$compression" == "gzip" ]]; then
            if [[ -f "\$tbi_file" ]]; then
                # Validate tabix index
                if tabix -l "\$vcf_file" >/dev/null 2>&1; then
                    log_result "Compression" "PASS" "Valid gzip compression with tabix index"
                else
                    log_result "Compression" "FAIL" "Invalid or corrupted tabix index"
                    return 1
                fi
            else
                log_result "Compression" "WARN" "Compressed VCF without tabix index"
            fi
        elif [[ "\$compression" == "uncompressed" ]]; then
            log_result "Compression" "PASS" "Uncompressed VCF file"
        else
            log_result "Compression" "WARN" "Unsupported compression format: \$compression"
        fi
        
        return 0
    }

    # Function to validate HISAT2 index format
    validate_hisat2_index() {
        local index_dir="\$1"
        
        if [[ -z "\$index_dir" ]]; then
            log_result "HISAT2 Index" "SKIP" "No HISAT2 index provided"
            return 0
        fi
        
        echo "INFO: Validating HISAT2 index format"
        
        # Look for HISAT2 index files (*.ht2 or *.ht2l)
        local ht2_files=\$(find "\$index_dir" -name "*.ht2" -o -name "*.ht2l" | wc -l)
        
        if [[ \$ht2_files -eq 0 ]]; then
            log_result "HISAT2 Index" "FAIL" "No HISAT2 index files (*.ht2 or *.ht2l) found"
            return 1
        fi
        
        # Check for minimum required index files (usually 8 files for standard index)
        if [[ \$ht2_files -lt 6 ]]; then
            log_result "HISAT2 Index" "WARN" "Found only \$ht2_files index files, may be incomplete"
        else
            log_result "HISAT2 Index" "PASS" "Found \$ht2_files HISAT2 index files"
        fi
        
        # List found index files
        echo "Found HISAT2 index files:" >> \$REPORT
        find "\$index_dir" -name "*.ht2" -o -name "*.ht2l" | sed 's/^/  /' >> \$REPORT
        echo "" >> \$REPORT
        
        return 0
    }

    # Function to validate VCF content structure
    validate_vcf_content() {
        local vcf_file="\$1"
        
        echo "INFO: Validating VCF content structure"
        
        # Check for required VCF columns
        local header_cols=\$(bcftools view -h "\$vcf_file" | tail -1 | tr '\t' '\n' | wc -l)
        if [[ \$header_cols -lt 8 ]]; then
            log_result "VCF Content" "FAIL" "VCF has only \$header_cols columns, minimum 8 required"
            return 1
        fi
        
        # Check for sample columns
        local sample_cols=\$((header_cols - 9))
        if [[ \$sample_cols -gt 0 ]]; then
            log_result "VCF Content" "PASS" "VCF contains \$sample_cols sample column(s)"
        else
            log_result "VCF Content" "PASS" "VCF contains no sample columns (sites-only VCF)"
        fi
        
        # Validate a sample of records
        local total_records=\$(bcftools view -H "\$vcf_file" | wc -l)
        if [[ \$total_records -gt 0 ]]; then
            local sample_size=\$((\$total_records < 100 ? \$total_records : 100))
            local invalid_records=\$(bcftools view -H "\$vcf_file" | head -\$sample_size | awk 'NF<8 {count++} END {print count+0}')
            
            if [[ \$invalid_records -gt 0 ]]; then
                log_result "VCF Content" "FAIL" "Found \$invalid_records invalid records in sample of \$sample_size"
                return 1
            else
                log_result "VCF Content" "PASS" "All sampled records (\$sample_size/\$total_records) have valid structure"
            fi
        fi
        
        return 0
    }

    # Main validation workflow
    echo "Starting VCF format validation for: ${vcf}"
    echo ""

    # Run all validation checks
    validation_failed=0
    
    validate_vcf_header "${vcf}" || validation_failed=1
    validate_coordinate_system "${vcf}" || validation_failed=1
    validate_compression "${vcf}" "${tbi}" || validation_failed=1
    validate_vcf_content "${vcf}" || validation_failed=1
    
    # Validate HISAT2 index if provided
    if [[ -n "${hisat2_index}" ]]; then
        validate_hisat2_index "${hisat2_index}" || validation_failed=1
    fi

    # Final validation summary
    echo "" >> \$REPORT
    echo "=========================================" >> \$REPORT
    if [[ \$validation_failed -eq 0 ]]; then
        echo "OVERALL RESULT: VALIDATION PASSED" >> \$REPORT
        log_result "Overall" "PASS" "All validation checks completed successfully"
    else
        echo "OVERALL RESULT: VALIDATION FAILED" >> \$REPORT
        log_result "Overall" "FAIL" "One or more validation checks failed"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^.*tabix //; s/ .*\$//')
        file: \$(file --version | head -n1 | sed 's/file-//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create dummy output files to match the expected outputs
    touch "${vcf}"
    if [[ -n "${tbi}" && "${tbi}" != "[]" ]]; then
        touch "${tbi}"
    fi
    echo "STUB: VCF validation would be performed here" > ${prefix}_validation_report.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' || echo "1.17")
        tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*tabix //; s/ .*\$//' || echo "1.17")
        file: \$(file --version | head -n1 | sed 's/file-//' || echo "5.44")
    END_VERSIONS
    """
}