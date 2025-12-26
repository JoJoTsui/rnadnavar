//
// PREPARE VCF-BASED REALIGNMENT: extract reads from candidate regions for re-alignment (RNA and DNA normal only)
// OPTIMIZED VERSION: Integrates sanitization, validation, and error handling components
//
include { VCF2BED                                        } from '../../../modules/local/vcf2bed/main'
// Extract read ids for selected regions
include { SAMTOOLS_EXTRACT_READ_IDS                    } from '../../../modules/local/extract_reads_id/main'
// Filter bam for selected regions
include { PICARD_FILTERSAMREADS                        } from '../../../modules/nf-core/picard/filtersamreads/main'
// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT  } from '../bam_convert_samtools/main'
// Realignment with HISAT2
include { FASTQ_ALIGN_HISAT2                            } from '../../nf-core/fastq_align_hisat2/main'

// OPTIMIZATION COMPONENTS: Sanitization, validation, and error handling
include { SANITIZE_CHANNELS                             } from '../sanitize_channels/main'
include { SAFE_CHANNEL_JOIN                             } from '../safe_channel_join/main'
include { ENHANCED_CRAM2BAM_CONVERSION                  } from '../enhanced_cram2bam_conversion/main'
include { PROCESS_MONITORING                            } from '../process_monitoring/main'
include { ERROR_SAFE_PROCESS                            } from '../error_safe_process/main'
include { INPUT_VALIDATION                              } from '../input_validation/main'


workflow BAM_EXTRACT_READS_HISAT2_ALIGN_VCF {
    take:
        input_sample
        vcf_with_candidates              // VCF with candidate regions to extract [meta, vcf, tbi]
        reads_to_realign                 // CRAM/BAM to extract reads from [meta, cram, crai]
        fasta
        fasta_fai
        dict
        hisat2_index
        splicesites
        dna_consensus_maf
        dna_varcall_mafs

    main:
        versions   = Channel.empty()
        bam_mapped = Channel.empty()
        bed        = Channel.empty()

        if (params.step in ['mapping', 'markduplicates', 'splitncigar',
        'prepare_recalibration', 'recalibrate', 'variant_calling', 'norm', 'consensus',
        'realignment'] && !(params.skip_tools && params.skip_tools.split(",").contains("realignment"))) {

            // === STEP 1: INPUT VALIDATION ===
            // Validate all input files and formats before processing (if enabled)
            if (params.enable_input_validation && !params.disable_realignment_optimization) {
                INPUT_VALIDATION(
                    vcf_with_candidates,
                    reads_to_realign,
                    hisat2_index
                )
                versions = versions.mix(INPUT_VALIDATION.out.versions)

                // Use validated inputs for further processing
                validated_vcf = INPUT_VALIDATION.out.vcf
                validated_cram = INPUT_VALIDATION.out.cram
            } else {
                // Skip validation - use inputs directly
                validated_vcf = vcf_with_candidates
                validated_cram = reads_to_realign
            }

            // === STEP 2: SAMPLE FILTERING WITH MONITORING ===
            // Filter samples: RNA (status=2) and DNA normal (status=0) only
            // Exclude DNA tumor (status=1)
            
            // Monitor VCF filtering (if enabled)
            if (params.enable_process_monitoring && !params.disable_realignment_optimization) {
                PROCESS_MONITORING(
                    "VCF_SAMPLE_FILTERING",
                    validated_vcf
                )
                monitored_vcf = PROCESS_MONITORING.out.output
                versions = versions.mix(PROCESS_MONITORING.out.versions)
            } else {
                monitored_vcf = validated_vcf
            }

            vcf_to_realign = monitored_vcf.branch{
                                        norealign: it[0].status == 1
                                        realign:   it[0].status == 2 || it[0].status == 0
            }

            // Monitor CRAM filtering (if enabled)
            if (params.enable_process_monitoring && !params.disable_realignment_optimization) {
                PROCESS_MONITORING(
                    "CRAM_SAMPLE_FILTERING", 
                    validated_cram
                )
                monitored_cram = PROCESS_MONITORING.out.output
                versions = versions.mix(PROCESS_MONITORING.out.versions)
            } else {
                monitored_cram = validated_cram
            }

            reads_to_realign_branch = monitored_cram.branch{
                                        norealign: it[0].status == 1
                                        realign:   it[0].status == 2 || it[0].status == 0
            }

            // === STEP 3: CHANNEL SANITIZATION ===
            // Clean and validate channel data to prevent StackOverflowError (if enabled)
            if (params.enable_channel_sanitization && !params.disable_realignment_optimization) {
                SANITIZE_CHANNELS(
                    vcf_to_realign.realign,
                    reads_to_realign_branch.realign
                )
                versions = versions.mix(SANITIZE_CHANNELS.out.versions)
                sanitized_vcf = SANITIZE_CHANNELS.out.vcf
                sanitized_cram = SANITIZE_CHANNELS.out.cram
            } else {
                // Skip sanitization - use channels directly
                sanitized_vcf = vcf_to_realign.realign
                sanitized_cram = reads_to_realign_branch.realign
            }

            // === STEP 4: SAFE CHANNEL JOIN ===
            // Perform robust channel joining with comprehensive validation (if enabled)
            if (params.enable_safe_channel_join && !params.disable_realignment_optimization) {
                SAFE_CHANNEL_JOIN(
                    sanitized_vcf,
                    sanitized_cram
                )
                versions = versions.mix(SAFE_CHANNEL_JOIN.out.versions)
                joined_data = SAFE_CHANNEL_JOIN.out.joined
            } else {
                // Use basic channel join - both VCF and CRAM are now filtered to RNA samples (status=2)
                joined_data = sanitized_cram
                    .map { meta, cram, crai -> [meta.patient, meta, cram, crai] }
                    .join(
                        sanitized_vcf.map { meta, vcf, tbi -> [meta.patient, meta, vcf, tbi] },
                        by: 0
                    )
                    .map { patient, cram_meta, cram, crai, vcf_meta, vcf, tbi ->
                        [cram_meta, cram, crai, vcf, tbi]
                    }
            }

            // === STEP 5: VCF TO BED CONVERSION WITH ERROR HANDLING ===
            // Convert VCF to BED format for region extraction
            vcf_to_bed_input = joined_data.map { meta, cram, crai, vcf, tbi -> 
                // Store file paths as strings to avoid serialization issues
                def enhanced_meta = meta + [
                    cram_path: cram.toString(),
                    crai_path: crai.toString()
                ]
                [enhanced_meta, vcf, tbi]
            }

            // Monitor VCF2BED conversion (if enabled)
            if (params.enable_process_monitoring && !params.disable_realignment_optimization) {
                PROCESS_MONITORING(
                    "VCF2BED_CONVERSION",
                    vcf_to_bed_input
                )
                monitored_vcf2bed = PROCESS_MONITORING.out.output
                versions = versions.mix(PROCESS_MONITORING.out.versions)
            } else {
                monitored_vcf2bed = vcf_to_bed_input
            }

            VCF2BED(monitored_vcf2bed)
            versions = versions.mix(VCF2BED.out.versions)
            bed = VCF2BED.out.bed

            // === STEP 6: READ ID EXTRACTION WITH MONITORING ===
            // Extract read IDs using BED from VCF2BED
            cram_to_extract = VCF2BED.out.bed.map { meta, bed -> 
                // Reconstruct file objects from stored paths with null checks
                if (!meta.cram_path) {
                    throw new IllegalArgumentException("Missing cram_path in metadata for sample ${meta.id}")
                }
                if (!meta.crai_path) {
                    throw new IllegalArgumentException("Missing crai_path in metadata for sample ${meta.id}")
                }
                def cram_file = file(meta.cram_path)
                def crai_file = file(meta.crai_path)
                [meta, cram_file, crai_file, bed]
            }

            if (params.enable_process_monitoring && !params.disable_realignment_optimization) {
                PROCESS_MONITORING(
                    "READ_ID_EXTRACTION",
                    cram_to_extract
                )
                monitored_extract = PROCESS_MONITORING.out.output
                versions = versions.mix(PROCESS_MONITORING.out.versions)
            } else {
                monitored_extract = cram_to_extract
            }

            SAMTOOLS_EXTRACT_READ_IDS(monitored_extract, fasta)
            versions = versions.mix(SAMTOOLS_EXTRACT_READ_IDS.out.versions)

            // === STEP 7: ENHANCED CRAM TO BAM CONVERSION ===
            // Use enhanced conversion with comprehensive error handling (if enabled)
            cram_to_convert = SAMTOOLS_EXTRACT_READ_IDS.out.read_ids.map { meta, readsid -> 
                // Store read IDs path as string and reconstruct CRAM files with null checks
                if (!meta.cram_path) {
                    throw new IllegalArgumentException("Missing cram_path in metadata for sample ${meta.id}")
                }
                if (!meta.crai_path) {
                    throw new IllegalArgumentException("Missing crai_path in metadata for sample ${meta.id}")
                }
                def cram_file = file(meta.cram_path)
                def crai_file = file(meta.crai_path)
                def enhanced_meta = meta + [readsid_path: readsid.toString()]
                [enhanced_meta, cram_file, crai_file]
            }

            if (params.enable_enhanced_cram_conversion && !params.disable_realignment_optimization) {
                ENHANCED_CRAM2BAM_CONVERSION(
                    cram_to_convert,
                    fasta,
                    fasta_fai
                )
                versions = versions.mix(ENHANCED_CRAM2BAM_CONVERSION.out.versions)
                converted_bam = ENHANCED_CRAM2BAM_CONVERSION.out.bam
            } else {
                // Use basic SAMTOOLS_CONVERT (would need to be implemented or imported)
                // For now, use the enhanced version but without extra monitoring
                ENHANCED_CRAM2BAM_CONVERSION(
                    cram_to_convert,
                    fasta,
                    fasta_fai
                )
                versions = versions.mix(ENHANCED_CRAM2BAM_CONVERSION.out.versions)
                converted_bam = ENHANCED_CRAM2BAM_CONVERSION.out.bam
            }

            // === STEP 8: READ FILTERING WITH ERROR HANDLING ===
            // Apply picard filtersamreads with enhanced error handling
            bam_to_filter = converted_bam.map { meta, bam -> 
                // Reconstruct read IDs file from stored path with null check
                if (!meta.readsid_path) {
                    throw new IllegalArgumentException("Missing readsid_path in metadata for sample ${meta.id}")
                }
                def readsid_file = file(meta.readsid_path)
                [meta, bam, readsid_file]
            }

            if (params.enable_process_monitoring && !params.disable_realignment_optimization) {
                PROCESS_MONITORING(
                    "READ_FILTERING",
                    bam_to_filter
                )
                monitored_filter = PROCESS_MONITORING.out.output
                versions = versions.mix(PROCESS_MONITORING.out.versions)
            } else {
                monitored_filter = bam_to_filter
            }

            PICARD_FILTERSAMREADS(monitored_filter, fasta, 'includeReadList')
            versions = versions.mix(PICARD_FILTERSAMREADS.out.versions)

            // === STEP 9: BAM TO FASTQ CONVERSION WITH MONITORING ===
            // Convert filtered BAM to FASTQ for realignment
            bam_to_fq = PICARD_FILTERSAMREADS.out.bam.join(PICARD_FILTERSAMREADS.out.bai)

            if (params.enable_process_monitoring && !params.disable_realignment_optimization) {
                PROCESS_MONITORING(
                    "BAM_TO_FASTQ_CONVERSION",
                    bam_to_fq
                )
                monitored_bam2fq = PROCESS_MONITORING.out.output
                versions = versions.mix(PROCESS_MONITORING.out.versions)
            } else {
                monitored_bam2fq = bam_to_fq
            }

            interleave_input = false  // Currently don't allow interleaved input
            CONVERT_FASTQ_INPUT(
                                monitored_bam2fq,
                                fasta,
                                fasta_fai.map{it -> [ [ id:"fasta_fai" ], it ]},
                                interleave_input
            )
            versions = versions.mix(CONVERT_FASTQ_INPUT.out.versions)

            // === STEP 10: HISAT2 REALIGNMENT WITH ERROR HANDLING ===
            // Perform realignment with HISAT2
            reads_for_realignment = CONVERT_FASTQ_INPUT.out.reads

            if (params.enable_process_monitoring && !params.disable_realignment_optimization) {
                PROCESS_MONITORING(
                    "HISAT2_REALIGNMENT",
                    reads_for_realignment
                )
                monitored_realign = PROCESS_MONITORING.out.output
                versions = versions.mix(PROCESS_MONITORING.out.versions)
            } else {
                monitored_realign = reads_for_realignment
            }

            // Note: single_end in meta always false for this subworkflow
            FASTQ_ALIGN_HISAT2(
                                monitored_realign.map{meta, reads -> [meta + [single_end:false], reads]},
                                hisat2_index,
                                splicesites,
                                fasta
            )
            versions = versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

            // === STEP 11: OUTPUT PROCESSING WITH VALIDATION ===
            // Process final outputs with validation and monitoring
            raw_bam_mapped = FASTQ_ALIGN_HISAT2.out.bam.join(FASTQ_ALIGN_HISAT2.out.bai)
                .map{meta,bam,bai -> [meta + [ id:meta.sample, data_type:"bam"], bam, bai]}

            if (params.enable_process_monitoring && !params.disable_realignment_optimization) {
                PROCESS_MONITORING(
                    "FINAL_OUTPUT_VALIDATION",
                    raw_bam_mapped
                )
                validated_output = PROCESS_MONITORING.out.output
                versions = versions.mix(PROCESS_MONITORING.out.versions)
            } else {
                validated_output = raw_bam_mapped
            }

            // Clean up metadata for final output
            bam_mapped = validated_output.map{meta, bam, bai -> 
                // Remove temporary fields added during processing
                def clean_meta = meta - meta.subMap(['single_end', 'cram_path', 'crai_path', 'readsid_path'])
                [clean_meta, bam]
            }
        }

    emit:
        bam_mapped         = bam_mapped
        bed                = bed
        dna_consensus_maf  = dna_consensus_maf
        dna_varcall_mafs   = dna_varcall_mafs
        versions           = versions
}
