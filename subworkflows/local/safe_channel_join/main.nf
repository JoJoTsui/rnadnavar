//
// SAFE_CHANNEL_JOIN: Implement robust channel joining with validation
//

workflow SAFE_CHANNEL_JOIN {
    take:
        vcf_channel     // [meta, vcf, tbi]
        cram_channel    // [meta, cram, crai]
    
    main:
        versions = Channel.empty()
        
        // Prepare VCF channel for joining with patient ID validation
        vcf_keyed = vcf_channel
            .filter { meta, vcf, tbi -> 
                if (!meta.patient) {
                    log.error "VCF channel missing patient ID: ${meta}"
                    return false
                }
                if (!vcf.exists()) {
                    log.error "VCF file does not exist: ${vcf}"
                    return false
                }
                if (!tbi.exists()) {
                    log.error "VCF index file does not exist: ${tbi}"
                    return false
                }
                return true
            }
            .map { meta, vcf, tbi -> 
                // Create join key using patient ID
                [meta.patient, meta, vcf, tbi] 
            }
        
        // Prepare CRAM channel for joining with patient ID validation
        cram_keyed = cram_channel
            .filter { meta, cram, crai -> 
                if (!meta.patient) {
                    log.error "CRAM channel missing patient ID: ${meta}"
                    return false
                }
                if (!cram.exists()) {
                    log.error "CRAM file does not exist: ${cram}"
                    return false
                }
                if (!crai.exists()) {
                    log.error "CRAM index file does not exist: ${crai}"
                    return false
                }
                return true
            }
            .map { meta, cram, crai -> 
                // Create join key using patient ID
                [meta.patient, meta, cram, crai] 
            }
        
        // Perform safe join with validation and error handling
        // Join by patient ID (first element) with failOnMismatch: false to handle missing matches gracefully
        joined_raw = cram_keyed.join(vcf_keyed, by: 0, failOnMismatch: false)
        
        // Validate and process joined results
        joined = joined_raw
            .filter { patient, cram_meta, cram, crai, vcf_meta, vcf, tbi -> 
                // Validate all components exist after join
                if (!patient) {
                    log.error "Join resulted in null patient ID"
                    return false
                }
                if (!cram_meta || !vcf_meta) {
                    log.error "Join resulted in null metadata for patient: ${patient}"
                    return false
                }
                if (!cram || !vcf) {
                    log.error "Join resulted in null files for patient: ${patient}"
                    return false
                }
                
                // Validate patient IDs match between CRAM and VCF
                if (cram_meta.patient != vcf_meta.patient) {
                    log.error "Patient ID mismatch in join: CRAM=${cram_meta.patient}, VCF=${vcf_meta.patient}"
                    return false
                }
                
                // Validate sample status compatibility
                if (cram_meta.status != vcf_meta.status) {
                    log.warn "Status mismatch for patient ${patient}: CRAM=${cram_meta.status}, VCF=${vcf_meta.status}"
                }
                
                return true
            }
            .map { patient, cram_meta, cram, crai, vcf_meta, vcf, tbi ->
                // Merge metadata safely, preserving essential fields from CRAM meta
                // and adding VCF file references without storing them in meta to prevent circular references
                def merged_meta = [
                    patient: cram_meta.patient,
                    sample: cram_meta.sample,
                    status: cram_meta.status,
                    id: cram_meta.id
                ]
                
                // Add safe additional fields
                if (cram_meta.single_end != null) {
                    merged_meta.single_end = cram_meta.single_end
                }
                if (cram_meta.data_type) {
                    merged_meta.data_type = cram_meta.data_type
                }
                
                // Add VCF information as separate fields to avoid circular references
                merged_meta.vcf_sample = vcf_meta.sample
                merged_meta.vcf_id = vcf_meta.id
                
                // Return joined data with merged metadata
                [merged_meta, cram, crai, vcf, tbi]
            }
        
        // Log successful joins for monitoring
        joined_logged = joined.map { meta, cram, crai, vcf, tbi ->
            log.info "Successfully joined patient ${meta.patient}: CRAM=${cram.getSimpleName()}, VCF=${vcf.getSimpleName()}"
            [meta, cram, crai, vcf, tbi]
        }
    
    emit:
        joined = joined_logged
        versions = versions
}