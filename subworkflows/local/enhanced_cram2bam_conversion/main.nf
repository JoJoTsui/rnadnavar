//
// ENHANCED_CRAM2BAM_CONVERSION: Convert CRAM to BAM with comprehensive validation and error handling
//

include { SAMTOOLS_CONVERT_ENHANCED } from '../../../modules/local/samtools_convert_enhanced/main'

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
                
                // Add conversion-specific metadata
                safe_meta.conversion_stage = "cram2bam"
                safe_meta.input_format = "cram"
                safe_meta.output_format = "bam"
                
                // Log metadata cleaning
                log.debug "ENHANCED_CRAM2BAM_CONVERSION: Cleaned metadata for ${meta.id}: ${safe_meta}"
                
                [safe_meta, cram, crai]
            }
        
        // Validate reference files (these are single files, not channels)
        // Perform validation directly without channel operations
        if (!fasta) {
            error "ENHANCED_CRAM2BAM_CONVERSION: Reference FASTA is null"
        }
        if (!fasta.exists()) {
            error "ENHANCED_CRAM2BAM_CONVERSION: Reference FASTA does not exist: ${fasta}"
        }
        
        def fasta_extension = fasta.getExtension().toLowerCase()
        if (!(fasta_extension in ['fa', 'fasta'])) {
            error "ENHANCED_CRAM2BAM_CONVERSION: Invalid FASTA extension: ${fasta_extension}"
        }
        
        if (fasta.size() < 1000) {
            log.warn "ENHANCED_CRAM2BAM_CONVERSION: FASTA file is very small (${fasta.size()} bytes): ${fasta}"
        }
        
        if (!fasta_fai) {
            error "ENHANCED_CRAM2BAM_CONVERSION: Reference FASTA index is null"
        }
        if (!fasta_fai.exists()) {
            error "ENHANCED_CRAM2BAM_CONVERSION: Reference FASTA index does not exist: ${fasta_fai}"
        }
        
        def fai_extension = fasta_fai.getExtension().toLowerCase()
        if (fai_extension != 'fai') {
            error "ENHANCED_CRAM2BAM_CONVERSION: Invalid FASTA index extension: ${fai_extension}"
        }
        
        if (fasta_fai.size() < 100) {
            log.warn "ENHANCED_CRAM2BAM_CONVERSION: FASTA index is very small (${fasta_fai.size()} bytes): ${fasta_fai}"
        }
        
        log.info "ENHANCED_CRAM2BAM_CONVERSION: Validated reference files: FASTA=${fasta.getSimpleName()}, FAI=${fasta_fai.getSimpleName()}"
        
        // Create channels for the validated reference files
        validated_fasta = Channel.of([[id: "fasta"], fasta])
        validated_fasta_fai = Channel.of([[id: "fasta_fai"], fasta_fai])
        
        // Perform conversion with enhanced error handling and resource management
        SAMTOOLS_CONVERT_ENHANCED(
            validated_input,
            validated_fasta,
            validated_fasta_fai
        )
        
        // Post-process outputs with validation
        validated_bam = SAMTOOLS_CONVERT_ENHANCED.out.bam
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
        
        validated_bai = SAMTOOLS_CONVERT_ENHANCED.out.bai
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
        versions = versions.mix(SAMTOOLS_CONVERT_ENHANCED.out.versions)
    
    emit:
        bam = validated_bam
        bai = validated_bai
        versions = versions
}