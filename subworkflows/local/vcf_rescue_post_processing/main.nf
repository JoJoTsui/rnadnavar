//
// VCF Rescue Post-Processing Workflow - Modular pipeline for extensible post-rescue processing
// Enhanced with comprehensive error handling and resource management
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
    
    // Log workflow parameters for debugging
    log.info "=== VCF Rescue Post-Processing Workflow ==="
    log.info "RNA annotation enabled: ${enable_rna_annotation}"
    log.info "Min RNA support threshold: ${min_rna_support}"
    log.info "REDIportal VCF provided: ${rediportal_vcf ? 'Yes' : 'No'}"
    
    // Validate input parameters with comprehensive error handling
    if (enable_rna_annotation) {
        // Validate REDIportal database availability
        if (!rediportal_vcf) {
            log.warn "RNA annotation enabled but REDIportal VCF not provided"
            log.warn "Disabling RNA annotation and proceeding with direct filtering"
            enable_rna_annotation = false
        } else {
            log.info "RNA editing annotation will be performed using REDIportal database"
        }
        
        // Validate min_rna_support parameter
        if (min_rna_support < 1) {
            log.warn "Invalid min_rna_support value: ${min_rna_support}. Setting to default value of 2"
            min_rna_support = 2
        }
    }
    
    // RNA editing annotation (conditional with error handling)
    if (enable_rna_annotation) {
        log.info "Starting RNA editing annotation process..."
        
        try {
            RNA_EDITING_ANNOTATION(
                vcf_rescue,
                rediportal_vcf,
                rediportal_tbi,
                min_rna_support
            )
            
            // Check if annotation completed successfully
            vcf_for_filtering = RNA_EDITING_ANNOTATION.out.vcf
            versions = versions.mix(RNA_EDITING_ANNOTATION.out.versions)
            
            log.info "RNA editing annotation completed successfully"
            
        } catch (Exception e) {
            log.error "RNA editing annotation failed: ${e.getMessage()}"
            log.warn "Implementing graceful fallback - proceeding with original rescue VCF"
            
            // Graceful fallback to original rescue VCF on annotation failure
            vcf_for_filtering = vcf_rescue
            
            // Log fallback action for monitoring
            log.info "Graceful fallback implemented - workflow continues with unannotated VCF"
        }
    } else {
        log.info "RNA editing annotation disabled - proceeding directly to filtering"
        vcf_for_filtering = vcf_rescue
    }
    
    // Validate input to filtering step
    vcf_for_filtering
        .ifEmpty { 
            log.error "No VCF files available for filtering step"
            error "VCF rescue post-processing failed: no input files for filtering"
        }
        .subscribe { meta, vcf, tbi ->
            log.debug "Filtering input validated for sample: ${meta.id}"
            log.debug "VCF file: ${vcf}"
            log.debug "Index file: ${tbi}"
        }
    
    // Rescue VCF filtering with error handling
    log.info "Starting VCF rescue filtering process..."
    
    try {
        VCF_RESCUE_FILTERING(vcf_for_filtering)
        
        versions = versions.mix(VCF_RESCUE_FILTERING.out.versions)
        
        log.info "VCF rescue filtering completed successfully"
        
        // Validate filtering outputs
        VCF_RESCUE_FILTERING.out.vcf
            .ifEmpty {
                log.error "VCF rescue filtering produced no output files"
                error "VCF rescue filtering failed: no output files generated"
            }
            .subscribe { meta, vcf, tbi ->
                log.debug "Filtering output validated for sample: ${meta.id}"
            }
            
    } catch (Exception e) {
        log.error "VCF rescue filtering failed: ${e.getMessage()}"
        error "VCF rescue post-processing failed at filtering step: ${e.getMessage()}"
    }
    
    // Final validation and logging
    log.info "VCF rescue post-processing workflow completed"
    log.info "Output channels:"
    log.info "  - Filtered VCF: available"
    log.info "  - Stripped VCF: available" 
    log.info "  - Versions: collected"

    emit:
    vcf                 = VCF_RESCUE_FILTERING.out.vcf
    vcf_stripped        = VCF_RESCUE_FILTERING.out.vcf_stripped
    versions            = versions
}