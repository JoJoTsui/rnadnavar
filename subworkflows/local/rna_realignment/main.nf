//
// RNA_REALIGNMENT_WORKFLOW - RNA variant calling on realigned BAMs
//
// This workflow handles RNA variant calling on BAMs that have been realigned with HISAT2.
// It is RNA-specific and does not include DNA steps or rescue (rescue happens separately).
//
include { BAM_GATK_PREPROCESSING    } from '../bam_gatk_preprocessing/main'
include { BAM_VARIANT_CALLING       } from '../bam_variant_calling/main'
include { VCF_NORMALIZE             } from '../vcf_normalize/main'
include { VCF_CONSENSUS_WORKFLOW    } from '../vcf_consensus_workflow/main'
include { VCF_FILTERING             } from '../vcf_filtering/main'
include { MAF_FILTERING              } from '../maf_filtering/main'
include { validateMeta               } from '../utils_nfcore_rnadnavar_pipeline/main'


workflow RNA_REALIGNMENT_WORKFLOW {
    take:
        input_sample
        bam_mapped                  // Realigned BAM from HISAT2 (RNA tumor only)
        dna_normal_cram             // DNA normal CRAM from first pass (for tumor-normal pairing)
        fasta
        fasta_fai
        dict
        dbsnp
        dbsnp_tbi
        pon
        pon_tbi
        known_sites_indels
        known_sites_indels_tbi
        germline_resource
        germline_resource_tbi
        intervals
        intervals_for_preprocessing
        intervals_bed_gz_tbi
        intervals_bed_combined
        intervals_and_num_intervals
        intervals_bed_gz_tbi_combined
        realignment                 // true (indicates this is realignment)

    main:
        reports   = Channel.empty()
        versions  = Channel.empty()

        // Validate inputs
        validated_bam = bam_mapped.filter { meta, bam -> validateMeta(meta) }
        validated_dna_cram = dna_normal_cram.filter { meta, cram, crai -> validateMeta(meta) }

        // GATK preprocessing on realigned BAM (RNA samples only)
        BAM_GATK_PREPROCESSING(
            input_sample,
            validated_bam,
            Channel.empty(),  // No separate cram_mapped for realignment
            fasta,
            fasta_fai,
            dict,
            known_sites_indels,
            known_sites_indels_tbi,
            intervals_for_preprocessing,
            intervals_and_num_intervals,
            realignment
        )

        // Combine realigned RNA tumor CRAM with DNA normal CRAM for tumor-normal pairing
        // RNA tumor samples (status=2) from realignment + DNA normal (status=0) from first pass
        cram_for_variant_calling = BAM_GATK_PREPROCESSING.out.cram_variant_calling
            .mix(validated_dna_cram)

        // Variant calling on combined tumor-normal samples
        BAM_VARIANT_CALLING(
            params.tools,
            cram_for_variant_calling,
            fasta,
            fasta_fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            intervals,
            intervals_bed_gz_tbi,
            intervals_bed_combined,
            intervals_bed_gz_tbi_combined,
            pon,
            pon_tbi,
            input_sample,
            realignment,
            false // params.no_intervals -> false, as we now pass merged realignment intervals
        )

        versions = versions.mix(BAM_VARIANT_CALLING.out.versions)
        reports  = reports.mix(BAM_VARIANT_CALLING.out.reports)

        // Normalize VCFs
        VCF_NORMALIZE(
            BAM_VARIANT_CALLING.out.vcf_to_normalise,
            fasta,
            fasta_fai.map{fai -> [[id:'fai'], fai]},
            input_sample,
            realignment
        )
        versions = versions.mix(VCF_NORMALIZE.out.versions)
        vcf_normalized = VCF_NORMALIZE.out.vcf

        // Note: VEP annotation removed - annotated VCF was not used downstream.
        // Only normalized VCF is used for consensus generation.

        // RNA Consensus (on realigned BAM)
        VCF_CONSENSUS_WORKFLOW(vcf_normalized, input_sample, realignment)
        rna_consensus_vcf = VCF_CONSENSUS_WORKFLOW.out.vcf
        versions = versions.mix(VCF_CONSENSUS_WORKFLOW.out.versions)

        // Note: NO rescue here - rescue happens in separate SECOND_RESCUE_WORKFLOW
        // that combines DNA (first round) with RNA (realigned)

        // Filter consensus VCF
        VCF_FILTERING(rna_consensus_vcf, fasta, input_sample, realignment)
        filtered_rna_vcf = VCF_FILTERING.out.vcf
        versions = versions.mix(VCF_FILTERING.out.versions)

        // Note: MAF filtering can be added here if MAF output is needed
        // For now, focus on VCF outputs

    emit:
        rna_consensus_vcf    = rna_consensus_vcf      // RNA consensus VCF from realigned BAM
        filtered_rna_vcf     = filtered_rna_vcf       // Filtered RNA VCF from realigned BAM
        versions             = versions
        reports              = reports
}
