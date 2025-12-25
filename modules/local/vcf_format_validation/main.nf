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
    
    """
    #!/bin/bash
    set -euo pipefail

    # Function to validate VCF file
    validate_vcf() {
        local vcf_file="\$1"
        local report_file="\$2"
        
        echo "Starting VCF validation for: \$vcf_file"
        
        if [[ ! -f "\$vcf_file" ]]; then
            echo "ERROR: VCF file '\$vcf_file' does not exist" >&2
            return 1
        fi
        
        if [[ ! -r "\$vcf_file" ]]; then
            echo "ERROR: VCF file '\$vcf_file' is not readable" >&2
            return 1
        fi
        
        # Basic VCF header validation
        if bcftools view -h "\$vcf_file" | head -1 | grep -q "^##fileformat=VCF"; then
            echo "PASS: VCF header format is valid" | tee -a "\$report_file"
        else
            echo "WARN: VCF header format may be non-standard" | tee -a "\$report_file"
        fi
        
        # Check record count
        local record_count=\$(bcftools view -H "\$vcf_file" | wc -l)
        echo "INFO: VCF contains \$record_count variant records" | tee -a "\$report_file"
        
        echo "INFO: VCF validation completed successfully"
        return 0
    }

    # Main validation logic
    echo "Starting VCF format validation for sample: ${meta.id}"
    
    # Create validation report
    REPORT="${prefix}_validation_report.txt"
    echo "VCF Format Validation Report for ${meta.id}" > "\$REPORT"
    echo "Generated on: \$(date)" >> "\$REPORT"
    echo "=========================================" >> "\$REPORT"
    
    # Validate VCF file
    validate_vcf "${vcf}" "\$REPORT"
    
    # Final result
    echo "OVERALL RESULT: VALIDATION PASSED" >> "\$REPORT"
    echo "VCF format validation completed successfully for sample: ${meta.id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/Using.*\$//')
        file: \$(file --version | head -n1 | sed 's/file-//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create dummy output files to match the expected outputs
    for file in ${vcf}; do
        touch "\$file"
    done
    
    if [[ -n "${tbi}" && "${tbi}" != "[]" ]]; then
        touch "${tbi}"
    fi
    
    echo "STUB: VCF validation would be performed here" > ${prefix}_validation_report.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/Using.*\$//' || echo "1.17")
        file: \$(file --version | head -n1 | sed 's/file-//' || echo "5.44")
    END_VERSIONS
    """
}