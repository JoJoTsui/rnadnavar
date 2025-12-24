//
// Comprehensive input validation subworkflow for VCF realignment optimization
//

include { FILE_VALIDATION }      from '../../../modules/local/file_validation/main'
include { VCF_FORMAT_VALIDATION } from '../../../modules/local/vcf_format_validation/main'

workflow INPUT_VALIDATION {
    
    take:
    vcf_channel     // channel: [meta, vcf, tbi]
    cram_channel    // channel: [meta, cram, crai] 
    hisat2_index    // path: hisat2 index directory (optional)
    
    main:
    
    ch_versions = Channel.empty()
    
    //
    // Validate CRAM files and their indices
    //
    cram_files_for_validation = cram_channel
        .map { meta, cram, crai ->
            def files = [cram]
            if (crai) {
                files.add(crai)
            }
            [meta, files]
        }
    
    FILE_VALIDATION (
        cram_files_for_validation
    )
    ch_versions = ch_versions.mix(FILE_VALIDATION.out.versions)
    
    // Reconstruct CRAM channel after validation
    validated_cram = FILE_VALIDATION.out.validated_files
        .map { meta, files ->
            def cram = files.find { it.name.endsWith('.cram') }
            def crai = files.find { it.name.endsWith('.crai') }
            [meta, cram, crai]
        }
    
    //
    // Validate VCF files with comprehensive format checking
    //
    VCF_FORMAT_VALIDATION (
        vcf_channel,
        hisat2_index ? hisat2_index.map { meta, path -> path } : []
    )
    ch_versions = ch_versions.mix(VCF_FORMAT_VALIDATION.out.versions)
    
    // Use the validated VCF from VCF_FORMAT_VALIDATION
    validated_vcf = VCF_FORMAT_VALIDATION.out.validated_vcf
    
    emit:
    vcf              = validated_vcf                           // channel: [meta, vcf, tbi]
    cram             = validated_cram                          // channel: [meta, cram, crai]
    validation_reports = VCF_FORMAT_VALIDATION.out.validation_report  // channel: validation_report.txt
    versions         = ch_versions                             // channel: versions.yml
}