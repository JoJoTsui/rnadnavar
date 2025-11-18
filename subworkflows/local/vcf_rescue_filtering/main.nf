// 
// VCF Rescue Filtering Workflow - Modality-aware variant classification
//
include { VCF_RESCUE_FILTER } from '../../../modules/local/vcf_rescue_filter/main'

workflow VCF_RESCUE_FILTERING {
    take:
    vcf_rescue      // channel: [ [meta], vcf, tbi ]
    
    main:
    versions = Channel.empty()
    filtered_vcf = Channel.empty()
    
    // Apply rescue-specific filtering
    VCF_RESCUE_FILTER(vcf_rescue)
    filtered_vcf = VCF_RESCUE_FILTER.out.vcf
    versions = versions.mix(VCF_RESCUE_FILTER.out.versions)
    
    emit:
    vcf      = filtered_vcf
    versions = versions
}
