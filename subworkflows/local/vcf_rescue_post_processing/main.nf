//
// VCF Rescue Post-Processing Workflow - Modular pipeline for extensible post-rescue processing
// REFACTORED: Uses runtime channel handling instead of parse-time conditionals
//
include { COSMIC_GNOMAD_ANNOTATION } from '../../../modules/local/cosmic_gnomad_annotation/main'
include { RNA_EDITING_ANNOTATION   } from '../../../modules/local/rna_editing_annotation/main'
include { VCF_RESCUE_FILTERING     } from '../vcf_rescue_filtering/main'

workflow VCF_RESCUE_POST_PROCESSING {
    take:
    vcf_rescue                      // channel: [meta, vcf, tbi]
    rediportal_vcf                  // channel: path
    rediportal_tbi                  // channel: path
    min_rna_support                 // val: integer
    enable_rna_annotation           // val: boolean
    cosmic_vcf                      // channel: path
    cosmic_tbi                      // channel: path
    gnomad_dir                      // channel: path
    enable_cosmic_gnomad_annotation // val: boolean
    cosmic_gnomad_verbose           // val: boolean

    main:
    versions = Channel.empty()
    
    // Validate min_rna_support
    def validated_min_rna_support = (min_rna_support < 1) ? 2 : min_rna_support
    
    // Debug logging only when enabled
    if (params.debug_verbose) {
        log.info "[DEBUG] VCF_RESCUE_POST_PROCESSING: cosmic_gnomad=${enable_cosmic_gnomad_annotation}, rna_annotation=${enable_rna_annotation}"
    }
    
    // ============================================================================
    // COSMIC/gnomAD Annotation - Uses .ifEmpty() to handle missing databases
    // This allows proper runtime channel handling instead of parse-time checks
    // ============================================================================
    
    // Create safe channels with NO_FILE placeholder for empty database channels
    cosmic_vcf_safe = cosmic_vcf.ifEmpty { file('NO_FILE') }
    cosmic_tbi_safe = cosmic_tbi.ifEmpty { file('NO_FILE') }
    gnomad_dir_safe = gnomad_dir.ifEmpty { file('NO_FILE') }
    
    // Branch VCF based on whether COSMIC/gnomAD annotation is enabled
    // This is evaluated at runtime, not parse time
    vcf_for_cosmic = vcf_rescue.branch {
        run_cosmic: enable_cosmic_gnomad_annotation
        skip_cosmic: true
    }
    
    // Run COSMIC/gnomAD annotation only when enabled
    // The process has internal handling for NO_FILE placeholders
    COSMIC_GNOMAD_ANNOTATION(
        vcf_for_cosmic.run_cosmic,
        cosmic_vcf_safe,
        cosmic_tbi_safe,
        gnomad_dir_safe,
        cosmic_gnomad_verbose
    )
    
    // Merge annotated and skipped VCFs
    vcf_after_cosmic_gnomad = COSMIC_GNOMAD_ANNOTATION.out.vcf.mix(vcf_for_cosmic.skip_cosmic)
    versions = versions.mix(COSMIC_GNOMAD_ANNOTATION.out.versions)
    
    // ============================================================================
    // RNA Editing Annotation - Uses .ifEmpty() to handle missing REDIportal files
    // ============================================================================
    
    // Create safe channels for REDIportal files
    rediportal_vcf_safe = rediportal_vcf.ifEmpty { file('NO_FILE') }
    rediportal_tbi_safe = rediportal_tbi.ifEmpty { file('NO_FILE') }
    
    // Branch VCF based on whether RNA annotation is enabled
    vcf_for_rna = vcf_after_cosmic_gnomad.branch {
        run_rna: enable_rna_annotation
        skip_rna: true
    }
    
    // Run RNA editing annotation only when enabled
    RNA_EDITING_ANNOTATION(
        vcf_for_rna.run_rna,
        rediportal_vcf_safe,
        rediportal_tbi_safe,
        validated_min_rna_support
    )
    
    // Merge annotated and skipped VCFs
    vcf_for_filtering = RNA_EDITING_ANNOTATION.out.vcf.mix(vcf_for_rna.skip_rna)
    versions = versions.mix(RNA_EDITING_ANNOTATION.out.versions)
    
    // Debug output (only when enabled)
    if (params.debug_verbose) {
        vcf_for_filtering.view { meta, vcf, tbi -> "[DEBUG] VCF for filtering: ${meta.id}" }
    }
    
    // Rescue VCF filtering (always executed)
    VCF_RESCUE_FILTERING(vcf_for_filtering)
    
    versions = versions.mix(VCF_RESCUE_FILTERING.out.versions)

    emit:
    vcf                 = VCF_RESCUE_FILTERING.out.vcf
    vcf_stripped        = VCF_RESCUE_FILTERING.out.vcf_stripped
    versions            = versions
}