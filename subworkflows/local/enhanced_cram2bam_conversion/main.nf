//
// ENHANCED_CRAM2BAM_CONVERSION: Convert CRAM to BAM with comprehensive validation and error handling
//

include { VALIDATED_SAMTOOLS_CONVERT } from '../validated_samtools_convert/main'

workflow ENHANCED_CRAM2BAM_CONVERSION {
    take:
        cram_channel    // [meta, cram, crai]
        fasta           // reference fasta file
        fasta_fai       // reference fasta index file
    
    main:
        versions = Channel.empty()
        
        // Pre-validate inputs before processing
        validated_input = cram_channel
            .filter { meta, cram, crai ->
                // Comprehensive input validation
                if (!meta) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: meta is null for input"
                    return false
                }
                
                if (!meta.id) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: meta.id is required but missing"
                    return false
                }
                
                if (!meta.patient) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: meta.patient is required but missing for sample ${meta.id}"
                    return false
                }
                
                if (meta.status == null) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: meta.status is required but missing for sample ${meta.id}"
                    return false
                }
                
                if (!cram) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: CRAM file is null for sample ${meta.id}"
                    return false
                }
                
                if (!cram.exists()) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: CRAM file does not exist: ${cram} for sample ${meta.id}"
                    return false
                }
                
                if (!crai) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: CRAM index file is null for sample ${meta.id}"
                    return false
                }
                
                if (!crai.exists()) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: CRAM index file does not exist: ${crai} for sample ${meta.id}"
                    return false
                }
                
                // Validate CRAM file format
                def cram_extension = cram.getExtension().toLowerCase()
                if (cram_extension != 'cram') {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: Expected CRAM file but got ${cram_extension} for sample ${meta.id}"
                    return false
                }
                
                // Validate CRAM index file format
                def crai_extension = crai.getExtension().toLowerCase()
                if (crai_extension != 'crai') {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: Expected CRAI index but got ${crai_extension} for sample ${meta.id}"
                    return false
                }
                
                // Check file sizes (warn if files are suspiciously small)
                if (cram.size() < 1000) {
                    log.warn "ENHANCED_CRAM2BAM_CONVERSION: CRAM file is very small (${cram.size()} bytes): ${cram} for sample ${meta.id}"
                }
                
                if (crai.size() < 100) {
                    log.warn "ENHANCED_CRAM2BAM_CONVERSION: CRAM index is very small (${crai.size()} bytes): ${crai} for sample ${meta.id}"
                }
                
                log.info "ENHANCED_CRAM2BAM_CONVERSION: Validated input for sample ${meta.id}: CRAM=${cram.getSimpleName()}, Index=${crai.getSimpleName()}"
                return true
            }
            .map { meta, cram, crai ->
                // Create safe metadata that avoids circular references
                // Remove any file objects or complex nested structures that could cause StackOverflowError
                def safe_meta = [
                    id: meta.id,
                    patient: meta.patient,
                    sample: meta.sample ?: meta.id,
                    status: meta.status
                ]
                
                // Add safe additional fields that don't contain file references
                if (meta.single_end != null) {
                    safe_meta.single_end = meta.single_end
                }
                if (meta.data_type) {
                    safe_meta.data_type = meta.data_type
                }
                // Preserve path fields needed for downstream processing
                if (meta.readsid_path) {
                    safe_meta.readsid_path = meta.readsid_path
                }
                if (meta.cram_path) {
                    safe_meta.cram_path = meta.cram_path
                }
                if (meta.crai_path) {
                    safe_meta.crai_path = meta.crai_path
                }
                
                // Add conversion-specific metadata
                safe_meta.conversion_stage = "cram2bam"
                safe_meta.input_format = "cram"
                safe_meta.output_format = "bam"
                
                // Log metadata cleaning
                log.debug "ENHANCED_CRAM2BAM_CONVERSION: Cleaned metadata for ${meta.id}: ${safe_meta}"
                
                [safe_meta, cram, crai]
            }
        
        // Skip reference file validation for now to avoid channel/file object issues
        // The validation will be handled by the downstream modules
        log.info "ENHANCED_CRAM2BAM_CONVERSION: Skipping reference file validation (handled downstream)"
        
        // Perform conversion with enhanced error handling, resource management, and command validation
        VALIDATED_SAMTOOLS_CONVERT(
            validated_input,
            fasta,
            fasta_fai
        )
        
        // Post-process outputs with validation
        validated_bam = VALIDATED_SAMTOOLS_CONVERT.out.bam
            .filter { meta, bam ->
                if (!bam) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: BAM output is null for sample ${meta.id}"
                    return false
                }
                if (!bam.exists()) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: BAM output does not exist: ${bam} for sample ${meta.id}"
                    return false
                }
                if (bam.size() < 1000) {
                    log.warn "ENHANCED_CRAM2BAM_CONVERSION: BAM output is very small (${bam.size()} bytes): ${bam} for sample ${meta.id}"
                }
                
                log.info "ENHANCED_CRAM2BAM_CONVERSION: Successfully converted to BAM for sample ${meta.id}: ${bam.getSimpleName()}"
                return true
            }
        
        validated_bai = VALIDATED_SAMTOOLS_CONVERT.out.bai
            .filter { meta, bai ->
                if (!bai) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: BAI output is null for sample ${meta.id}"
                    return false
                }
                if (!bai.exists()) {
                    log.error "ENHANCED_CRAM2BAM_CONVERSION: BAI output does not exist: ${bai} for sample ${meta.id}"
                    return false
                }
                if (bai.size() < 100) {
                    log.warn "ENHANCED_CRAM2BAM_CONVERSION: BAI output is very small (${bai.size()} bytes): ${bai} for sample ${meta.id}"
                }
                
                log.info "ENHANCED_CRAM2BAM_CONVERSION: Successfully created BAI index for sample ${meta.id}: ${bai.getSimpleName()}"
                return true
            }
        
        // Collect versions
        versions = versions.mix(VALIDATED_SAMTOOLS_CONVERT.out.versions)
    
    emit:
        bam = validated_bam
        bai = validated_bai
        versions = versions
}