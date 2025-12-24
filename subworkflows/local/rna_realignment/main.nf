//
// RNA_REALIGNMENT_WORKFLOW - RNA variant calling on realigned BAMs
//
// This workflow handles RNA variant calling on BAMs that have been realigned with HISAT2.
// It is RNA-specific and does not include DNA steps or rescue (rescue happens separately).
//
include { BAM_GATK_PREPROCESSING    } from '../bam_gatk_preprocessing/main'
include { BAM_VARIANT_CALLING       } from '../bam_variant_calling/main'
include { VCF_NORMALIZE             } from '../vcf_normalize/main'
include { VCF_ANNOTATE              } from '../vcf_annotate/main'
include { VCF_CONSENSUS_WORKFLOW    } from '../vcf_consensus_workflow/main'
include { VCF_FILTERING             } from '../vcf_filtering/main'
include { MAF_FILTERING              } from '../maf_filtering/main'


workflow RNA_REALIGNMENT_WORKFLOW {
    take:
        input_sample
        bam_mapped                  // Realigned BAM from HISAT2
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
        vep_cache
        realignment                 // true (indicates this is realignment)

    main:
        reports   = Channel.empty()
        versions  = Channel.empty()

        // GATK preprocessing on realigned BAM (RNA samples)
        BAM_GATK_PREPROCESSING(
            input_sample,
            bam_mapped,
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

        // Variant calling on realigned BAM
        BAM_VARIANT_CALLING(
            params.tools,
            BAM_GATK_PREPROCESSING.out.cram_variant_calling,
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
            Channel.empty()  // No intervals for realignment (targeted regions)
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

        // Annotate VCFs
        VCF_ANNOTATE(
            vcf_normalized.map{ meta, vcf, tbi -> [ meta + [ file_name: vcf.baseName ], vcf, [tbi] ] },
            fasta,
            input_sample,
            realignment,
            vep_cache
        )
        vcf_annotated = VCF_ANNOTATE.out.vcf_ann
        versions = versions.mix(VCF_ANNOTATE.out.versions)
        reports  = reports.mix(VCF_ANNOTATE.out.reports)

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
