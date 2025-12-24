//
// SECOND_RESCUE_WORKFLOW - Cross-modality rescue using realigned RNA VCFs
//
// This workflow performs cross-modality rescue between DNA consensus (from first round)
// and realigned RNA consensus VCFs. It applies full annotation and filtering to the rescued VCF.
//
include { VCF_RESCUE_WORKFLOW              } from '../vcf_rescue_workflow/main'
include { VCF_RESCUE_POST_PROCESSING      } from '../vcf_rescue_post_processing/main'


workflow SECOND_RESCUE_WORKFLOW {
    take:
        first_round_dna_consensus_vcf   // DNA consensus from first round
        realigned_rna_consensus_vcf     // RNA consensus from realigned BAM
        first_round_dna_caller_vcfs     // Optional: Individual DNA caller VCFs
        realigned_rna_caller_vcfs       // Optional: Individual RNA caller VCFs
        fasta
        rediportal_vcf
        rediportal_tbi
        min_rna_support
        enable_rna_annotation
        cosmic_vcf
        cosmic_tbi
        gnomad_dir
        enable_cosmic_gnomad_annotation
        cosmic_gnomad_verbose

    main:
        versions           = Channel.empty()
        second_rescued_vcf = Channel.empty()

        // Run rescue workflow
        VCF_RESCUE_WORKFLOW(
            first_round_dna_consensus_vcf,
            realigned_rna_consensus_vcf,
            first_round_dna_caller_vcfs ?: Channel.empty(),
            realigned_rna_caller_vcfs ?: Channel.empty()
        )
        versions = versions.mix(VCF_RESCUE_WORKFLOW.out.versions)
        rescued_vcf = VCF_RESCUE_WORKFLOW.out.vcf

        // Apply post-rescue processing (annotation, filtering)
        VCF_RESCUE_POST_PROCESSING(
            rescued_vcf,
            rediportal_vcf,
            rediportal_tbi,
            min_rna_support,
            enable_rna_annotation,
            cosmic_vcf,
            cosmic_tbi,
            gnomad_dir,
            enable_cosmic_gnomad_annotation,
            cosmic_gnomad_verbose
        )
        versions = versions.mix(VCF_RESCUE_POST_PROCESSING.out.versions)
        second_rescued_vcf = VCF_RESCUE_POST_PROCESSING.out.vcf

    emit:
        second_rescued_vcf    = second_rescued_vcf     // Final rescued VCF with full annotation
        versions              = versions
}
