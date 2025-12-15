//
// VCF Rescue Post-Processing Workflow - Modular pipeline for extensible post-rescue processing
//
include { RNA_EDITING_ANNOTATION } from '../../../modules/local/rna_editing_annotation/main'
include { VCF_RESCUE_FILTERING   } from '../vcf_rescue_filtering/main'

workflow VCF_RESCUE_POST_PROCESSING {
    take:
    vcf_rescue              // channel: [meta, vcf, tbi]
    rediportal_vcf          // channel: path
    rediportal_tbi          // channel: path
    min_rna_support         // val: integer
    enable_rna_annotation   // val: boolean

    main:
    versions = Channel.empty()
    
    // RNA editing annotation (conditional)
    if (enable_rna_annotation) {
        RNA_EDITING_ANNOTATION(
            vcf_rescue,
            rediportal_vcf,
            rediportal_tbi,
            min_rna_support
        )
        vcf_for_filtering = RNA_EDITING_ANNOTATION.out.vcf
        versions = versions.mix(RNA_EDITING_ANNOTATION.out.versions)
    } else {
        vcf_for_filtering = vcf_rescue
    }
    
    // Rescue VCF filtering
    VCF_RESCUE_FILTERING(vcf_for_filtering)
    
    versions = versions.mix(VCF_RESCUE_FILTERING.out.versions)

    emit:
    vcf                 = VCF_RESCUE_FILTERING.out.vcf
    vcf_stripped        = VCF_RESCUE_FILTERING.out.vcf_stripped
    versions            = versions
}