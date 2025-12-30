//
// VCF Rescue Post-Processing Workflow - Modular pipeline for extensible post-rescue processing
// Enhanced with comprehensive error handling and resource management
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
    
    // Log workflow parameters for debugging
    log.info "=== VCF Rescue Post-Processing Workflow ==="
    log.info "COSMIC/gnomAD annotation enabled: ${enable_cosmic_gnomad_annotation}"
    log.info "COSMIC database provided: ${cosmic_vcf ? 'Yes' : 'No'}"
    log.info "gnomAD database provided: ${gnomad_dir ? 'Yes' : 'No'}"
    log.info "COSMIC/gnomAD verbose logging: ${cosmic_gnomad_verbose}"
    log.info "RNA annotation enabled: ${enable_rna_annotation}"
    log.info "Min RNA support threshold: ${min_rna_support}"
    log.info "REDIportal VCF provided: ${rediportal_vcf ? 'Yes' : 'No'}"
    log.info "REDIportal VCF channel: ${rediportal_vcf}"
    log.info "REDIportal TBI channel: ${rediportal_tbi}"
    
    // Determine actual execution flags based on input validation
    // Note: Channels need explicit check - empty channels are truthy in Groovy
    def cosmic_vcf_provided = cosmic_vcf != null && !cosmic_vcf.toString().contains('DataflowBroadcast')
    def gnomad_dir_provided = gnomad_dir != null && !gnomad_dir.toString().contains('DataflowBroadcast')
    def rediportal_provided = rediportal_vcf != null && !rediportal_vcf.toString().contains('DataflowBroadcast')
    
    def run_cosmic_gnomad = enable_cosmic_gnomad_annotation && (cosmic_vcf_provided || gnomad_dir_provided)
    def run_rna_annotation = enable_rna_annotation && rediportal_provided
    def validated_min_rna_support = (min_rna_support < 1) ? 2 : min_rna_support
    
    // Log validation results with detailed channel inspection
    log.info "Conditional evaluation: run_cosmic_gnomad = ${run_cosmic_gnomad}, run_rna_annotation = ${run_rna_annotation}"
    if (enable_cosmic_gnomad_annotation && !run_cosmic_gnomad) {
        log.warn "COSMIC/gnomAD annotation enabled but no databases provided - skipping"
    }
    if (enable_rna_annotation && !run_rna_annotation) {
        log.warn "RNA annotation enabled (${enable_rna_annotation}) but REDIportal VCF not provided or empty - skipping"
        log.warn "  rediportal_vcf channel content: ${rediportal_vcf}"
    }
    if (min_rna_support < 1) {
        log.warn "Invalid min_rna_support value: ${min_rna_support}. Using default value of 2"
    }
    
    // COSMIC/gnomAD annotation (conditional with error handling)
    if (run_cosmic_gnomad) {
        // log.info "COSMIC/gnomAD annotation will be performed"
        // if (cosmic_vcf) log.info "  - COSMIC database: ${cosmic_vcf}"
        // if (gnomad_dir) log.info "  - gnomAD database: ${gnomad_dir}"
        
        COSMIC_GNOMAD_ANNOTATION(
            vcf_rescue,
            cosmic_vcf,
            cosmic_tbi,
            gnomad_dir,
            cosmic_gnomad_verbose
        )
        
        vcf_after_cosmic_gnomad = COSMIC_GNOMAD_ANNOTATION.out.vcf
        versions = versions.mix(COSMIC_GNOMAD_ANNOTATION.out.versions)
    } else {
        vcf_after_cosmic_gnomad = vcf_rescue
    }

    // RNA editing annotation (conditional with error handling)
    if (run_rna_annotation) {
        // log.info "RNA editing annotation will be performed using REDIportal database"
        
        RNA_EDITING_ANNOTATION(
            vcf_after_cosmic_gnomad,
            rediportal_vcf,
            rediportal_tbi,
            validated_min_rna_support
        )
        
        vcf_for_filtering = RNA_EDITING_ANNOTATION.out.vcf
        versions = versions.mix(RNA_EDITING_ANNOTATION.out.versions)
    } else {
        vcf_for_filtering = vcf_after_cosmic_gnomad
    }
    
    // Rescue VCF filtering (always executed)
    VCF_RESCUE_FILTERING(vcf_for_filtering)
    
    versions = versions.mix(VCF_RESCUE_FILTERING.out.versions)

    emit:
    vcf                 = VCF_RESCUE_FILTERING.out.vcf
    vcf_stripped        = VCF_RESCUE_FILTERING.out.vcf_stripped
    versions            = versions
}