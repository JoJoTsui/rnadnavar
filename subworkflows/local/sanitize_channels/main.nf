//
// SANITIZE_CHANNELS: Clean and validate channel data to prevent StackOverflowError
//

workflow SANITIZE_CHANNELS {
    take:
        vcf_channel     // [meta, vcf, tbi]
        cram_channel    // [meta, cram, crai]
    
    main:
        versions = Channel.empty()
        
        // Validate and clean VCF channel
        clean_vcf = vcf_channel
            .filter { meta, vcf, tbi ->
                // Validate essential fields exist
                if (!meta.patient) {
                    log.warn "Missing patient ID in VCF meta: ${meta}"
                    return false
                }
                if (!meta.id) {
                    log.warn "Missing sample ID in VCF meta: ${meta}"
                    return false
                }
                if (meta.status == null) {
                    log.warn "Missing status in VCF meta: ${meta}"
                    return false
                }
                return true
            }
            .map { meta, vcf, tbi ->
                // Remove potential circular references and file objects from meta
                // Keep only essential fields to prevent StackOverflowError
                def clean_meta = [
                    patient: meta.patient,
                    sample: meta.id + "_realign", 
                    status: meta.status,
                    id: meta.id + "_realign"
                ]
                
                // Add any additional safe fields that don't contain file references
                if (meta.single_end != null) {
                    clean_meta.single_end = meta.single_end
                }
                if (meta.data_type) {
                    clean_meta.data_type = meta.data_type
                }
                
                [clean_meta, vcf, tbi]
            }
        
        // Validate and clean CRAM channel
        clean_cram = cram_channel
            .filter { meta, cram, crai ->
                // Validate essential fields exist
                if (!meta.patient) {
                    log.warn "Missing patient ID in CRAM meta: ${meta}"
                    return false
                }
                if (!meta.id) {
                    log.warn "Missing sample ID in CRAM meta: ${meta}"
                    return false
                }
                if (meta.status == null) {
                    log.warn "Missing status in CRAM meta: ${meta}"
                    return false
                }
                
                // Validate file accessibility
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
                // Remove potential circular references and file objects from meta
                // Keep only essential fields to prevent StackOverflowError
                def clean_meta = [
                    patient: meta.patient,
                    sample: meta.id + "_realign",
                    status: meta.status, 
                    id: meta.id + "_realign"
                ]
                
                // Add any additional safe fields that don't contain file references
                if (meta.single_end != null) {
                    clean_meta.single_end = meta.single_end
                }
                if (meta.data_type) {
                    clean_meta.data_type = meta.data_type
                }
                
                [clean_meta, cram, crai]
            }
    
    emit:
        vcf = clean_vcf
        cram = clean_cram
        versions = versions
}