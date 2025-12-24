//
// PREPARE VCF-BASED REALIGNMENT: extract reads from candidate regions for re-alignment (RNA and DNA normal only)
//
include { VCF2BED                                        } from '../../../modules/local/vcf2bed/main'
// Extract read ids for selected regions
include { SAMTOOLS_EXTRACT_READ_IDS                    } from '../../../modules/local/extract_reads_id/main'
// Filter bam for selected regions
include { PICARD_FILTERSAMREADS                        } from '../../../modules/nf-core/picard/filtersamreads/main'
// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT  } from '../bam_convert_samtools/main'
// Conver CRAM to BAM (picard/filtersamreads can't work with cram)
include { SAMTOOLS_CONVERT as CONVERT_CRAM2BAM         } from '../../../modules/nf-core/samtools/convert/main'
// Realignment with HISAT2
include { FASTQ_ALIGN_HISAT2                            } from '../../nf-core/fastq_align_hisat2/main'


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

        if (params.step in ['mapping', 'markduplicates', 'splitncigar',
        'prepare_recalibration', 'recalibrate', 'variant_calling', 'norm', 'consensus',
        'realignment'] && !(params.skip_tools && params.skip_tools.split(",").contains("realignment"))) {

            // Filter samples: RNA (status=2) and DNA normal (status=0) only
            // Exclude DNA tumor (status=1)
            vcf_to_realign = vcf_with_candidates.branch{
                                        norealign: it[0].status == 1
                                        realign:   it[0].status == 2 || it[0].status == 0
            }

            reads_to_realign_branch = reads_to_realign.branch{
                                        norealign: it[0].status == 1
                                        realign:   it[0].status == 2 || it[0].status == 0
            }

            // Join VCF with CRAM files
            // Map CRAM to add _realign suffix
            reads_to_realign_and_join = reads_to_realign_branch.realign.map{meta, cram, crai ->
                                        [[patient: meta.patient,
                                        sample:  meta.id + "_realign",
                                        status:  meta.status,
                                        id:      meta.id + "_realign"], cram, crai]
            }

            vcf_to_realign_and_join = vcf_to_realign.realign.map{meta, vcf, tbi ->
                                        [[patient: meta.patient,
                                        sample:  meta.id + "_realign",
                                        status:  meta.status,
                                        id:      meta.id + "_realign"], vcf, tbi]
            }

            // Join VCF with CRAM by patient
            cram_vcf_joined = reads_to_realign_and_join.join(vcf_to_realign_and_join)

            // VCF to BED conversion
            vcf_to_bed = cram_vcf_joined.map{meta, cram, crai, vcf, tbi -> [meta + [cram_file:cram, crai_file:crai], vcf, tbi]}
            VCF2BED(vcf_to_bed)

            // Extract read IDs using BED from VCF2BED
            cram_to_extract = VCF2BED.out.bed.map{meta, bed -> [meta, meta.cram_file, meta.crai_file, bed]}
            SAMTOOLS_EXTRACT_READ_IDS(cram_to_extract, fasta)

            // Extract reads
            cram_to_convert = SAMTOOLS_EXTRACT_READ_IDS.out.read_ids.map{meta, readsid -> [meta + [readsid_file:readsid], meta.cram_file, meta.crai_file]}

            // 1) Convert cram 2 bam
            CONVERT_CRAM2BAM(cram_to_convert, fasta, fasta_fai.map{fai -> [[id:"fai"], fai]})
            bam_to_filter = CONVERT_CRAM2BAM.out.bam.map{meta, bam -> [meta, bam, meta.readsid_file]}

            // 2) Apply picard filtersamreads
            PICARD_FILTERSAMREADS(bam_to_filter, fasta, 'includeReadList')  // bam -> filtered_bam

            // Convert to FQ
            bam_to_fq = PICARD_FILTERSAMREADS.out.bam.join(PICARD_FILTERSAMREADS.out.bai)
            interleave_input = false  // Currently don't allow interleaved input
            CONVERT_FASTQ_INPUT(
                                bam_to_fq,
                                fasta,
                                fasta_fai.map{it -> [ [ id:"fasta_fai" ], it ]},
                                interleave_input
            )

            // Align with HISAT2
            reads_for_realignment = CONVERT_FASTQ_INPUT.out.reads
            // Note: single_end in meta always false for this subworkflow TODO: add to samplesheet in future?
            FASTQ_ALIGN_HISAT2(
                                reads_for_realignment.map{meta, reads -> [meta + [single_end:false], reads]},
                                hisat2_index,
                                splicesites,
                                fasta
            )

            // Mix with index add data type and change id to sample
            bam_mapped = FASTQ_ALIGN_HISAT2.out.bam.join(FASTQ_ALIGN_HISAT2.out.bai).map{meta,bam,bai -> [meta + [ id:meta.sample, data_type:"bam"], bam, bai]}
        }

        bam_mapped = bam_mapped.map{meta, bam, bai -> [meta - meta.subMap('single_end'), bam]}

    emit:
        bam_mapped         = bam_mapped
        dna_consensus_maf  = dna_consensus_maf
        dna_varcall_mafs   = dna_varcall_mafs
        versions           = versions
}
