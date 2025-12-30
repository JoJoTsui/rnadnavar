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
include { ENHANCED_CRAM2BAM_CONVERSION                  } from '../enhanced_cram2bam_conversion/main'
include { validateMeta                                  } from '../utils_nfcore_rnadnavar_pipeline/main'


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
            // Validate input files and formats using lightweight validation
            validated_vcf = vcf_with_candidates.filter { meta, vcf, tbi -> validateMeta(meta) }
            validated_cram = reads_to_realign.filter { meta, cram, crai -> validateMeta(meta) }

            // === STEP 2: SAMPLE FILTERING ===
            // Filter samples: RNA (status=2) and DNA normal (status=0) only
            // Exclude DNA tumor (status=1)
            
            vcf_to_realign = validated_vcf.branch{
                                        norealign: it[0].status == 1
                                        realign:   it[0].status == 2 || it[0].status == 0
            }

            reads_to_realign_branch = validated_cram.branch{
                                        norealign: it[0].status == 1
                                        realign:   it[0].status == 2 || it[0].status == 0
            }

            // === STEP 3: CHANNEL JOIN ===
            // Use basic channel join - both VCF and CRAM are now filtered to RNA samples (status=2)
            // Ensure we join by patient ID
            
            // Prepare VCF for join: [patient, meta, vcf, tbi]
            vcf_keyed = vcf_to_realign.realign
                .map { meta, vcf, tbi -> [meta.patient, meta, vcf, tbi] }
            
            // Prepare CRAM for join: [patient, meta, cram, crai]
            cram_keyed = reads_to_realign_branch.realign
                .map { meta, cram, crai -> [meta.patient, meta, cram, crai] }
            
            // Join
            joined_data = cram_keyed.join(vcf_keyed, by: 0)
                .map { patient, cram_meta, cram, crai, vcf_meta, vcf, tbi ->
                    // Merge metadata safely
                    def merged_meta = cram_meta + [
                        vcf_id: vcf_meta.id,
                        vcf_sample: vcf_meta.sample
                    ]
                    [merged_meta, cram, crai, vcf, tbi]
                }

            // === STEP 4: VCF TO BED CONVERSION ===
            // Convert VCF to BED format for region extraction
            vcf_to_bed_input = joined_data.map { meta, cram, crai, vcf, tbi -> 
                // Store file paths as strings to avoid serialization issues if needed, 
                // but with lightweight meta we might not need to be so aggressive.
                // Keeping it safe though.
                def enhanced_meta = meta + [
                    cram_path: cram.toString(),
                    crai_path: crai.toString()
                ]
                [enhanced_meta, vcf, tbi]
            }

            VCF2BED(vcf_to_bed_input)
            versions = versions.mix(VCF2BED.out.versions)
            bed = VCF2BED.out.bed

            // === STEP 5: READ ID EXTRACTION ===
            // Extract read IDs using BED from VCF2BED
            cram_to_extract = VCF2BED.out.bed.map { meta, bed -> 
                // Reconstruct file objects
                def cram_file = file(meta.cram_path)
                def crai_file = file(meta.crai_path)
                [meta, cram_file, crai_file, bed]
            }

            SAMTOOLS_EXTRACT_READ_IDS(cram_to_extract, fasta)
            versions = versions.mix(SAMTOOLS_EXTRACT_READ_IDS.out.versions)

            // === STEP 6: ENHANCED CRAM TO BAM CONVERSION ===
            cram_to_convert = SAMTOOLS_EXTRACT_READ_IDS.out.read_ids.map { meta, readsid -> 
                def cram_file = file(meta.cram_path)
                def crai_file = file(meta.crai_path)
                def enhanced_meta = meta + [readsid_path: readsid.toString()]
                [enhanced_meta, cram_file, crai_file]
            }

            ENHANCED_CRAM2BAM_CONVERSION(
                cram_to_convert,
                fasta,
                fasta_fai
            )
            versions = versions.mix(ENHANCED_CRAM2BAM_CONVERSION.out.versions)
            converted_bam = ENHANCED_CRAM2BAM_CONVERSION.out.bam

            // === STEP 7: READ FILTERING ===
            // Apply picard filtersamreads
            bam_to_filter = converted_bam.map { meta, bam -> 
                def readsid_file = file(meta.readsid_path)
                [meta, bam, readsid_file]
            }

            PICARD_FILTERSAMREADS(bam_to_filter, fasta, 'includeReadList')
            versions = versions.mix(PICARD_FILTERSAMREADS.out.versions)

            // === STEP 8: BAM TO FASTQ CONVERSION ===
            // Convert filtered BAM to FASTQ for realignment
            bam_to_fq = PICARD_FILTERSAMREADS.out.bam.join(PICARD_FILTERSAMREADS.out.bai)

            interleave_input = false  // Currently don't allow interleaved input
            CONVERT_FASTQ_INPUT(
                                bam_to_fq,
                                fasta,
                                fasta_fai.map{it -> [ [ id:"fasta_fai" ], it ]},
                                interleave_input
            )
            versions = versions.mix(CONVERT_FASTQ_INPUT.out.versions)

            // === STEP 9: HISAT2 REALIGNMENT ===
            // Perform realignment with HISAT2
            reads_for_realignment = CONVERT_FASTQ_INPUT.out.reads

            // Note: single_end in meta always false for this subworkflow
            FASTQ_ALIGN_HISAT2(
                                reads_for_realignment.map{meta, reads -> [meta + [single_end:false], reads]},
                                hisat2_index,
                                splicesites,
                                fasta
            )
            versions = versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

            // === STEP 10: OUTPUT PROCESSING ===
            // Process final outputs
            raw_bam_mapped = FASTQ_ALIGN_HISAT2.out.bam.join(FASTQ_ALIGN_HISAT2.out.bai)
                .map{meta,bam,bai -> [meta + [ id:meta.sample, data_type:"bam"], bam, bai]}

            // Clean up metadata for final output
            bam_mapped = raw_bam_mapped.map{meta, bam, bai -> 
                // Remove temporary fields added during processing
                def clean_meta = meta - meta.subMap(['single_end', 'cram_path', 'crai_path', 'readsid_path', 'vcf_id', 'vcf_sample'])
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
