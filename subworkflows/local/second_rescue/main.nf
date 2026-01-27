//
// SECOND_RESCUE_WORKFLOW - Cross-modality rescue using realigned RNA VCFs
//
// This workflow performs cross-modality rescue between DNA consensus (from first round)
// and realigned RNA consensus VCFs. It applies full annotation and filtering to the rescued VCF.
//
// IMPORTANT: Individual caller VCFs must be passed for correct N_RNA_CALLERS_SUPPORT calculation.
// Without individual caller VCFs, all variants will have N_RNA_CALLERS_SUPPORT=0 and
// RNA editing annotation will not work correctly.
//
include { VCF_RESCUE_WORKFLOW              } from '../vcf_rescue_workflow/main'
include { VCF_RESCUE_POST_PROCESSING      } from '../vcf_rescue_post_processing/main'
include { VCF_ANNOTATE                     } from '../vcf_annotate/main'


workflow SECOND_RESCUE_WORKFLOW {
    take:
        first_round_dna_consensus_vcf   // DNA consensus from first round
        realigned_rna_consensus_vcf     // RNA consensus from realigned BAM
        first_round_dna_vcf_normalized  // Individual DNA caller VCFs from first round (vcf_normalized)
        realigned_rna_vcf_normalized    // Individual RNA caller VCFs from realignment (vcf_normalized)
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
        input_sample                    // Required for VEP annotation
        vep_cache                       // Required for VEP annotation

    main:
        versions           = Channel.empty()
        reports            = Channel.empty()
        second_rescued_vcf = Channel.empty()

        // Prepare DNA caller VCFs from first round vcf_normalized
        // Filter to DNA samples (status <= 1) and format for VCF_RESCUE_WORKFLOW
        dna_caller_vcfs = first_round_dna_vcf_normalized
            .filter { meta, vcf, tbi -> meta.status <= 1 }
            .map { meta, vcf, tbi ->
                def caller = meta.variantcaller ?: 'unknown'
                [meta, vcf, tbi, caller]
            }

        // Prepare RNA caller VCFs from realignment vcf_normalized
        // Filter to RNA samples (status == 2) and format for VCF_RESCUE_WORKFLOW
        rna_caller_vcfs = realigned_rna_vcf_normalized
            .filter { meta, vcf, tbi -> meta.status == 2 }
            .map { meta, vcf, tbi ->
                def caller = meta.variantcaller ?: 'unknown'
                [meta, vcf, tbi, caller]
            }

        // Run rescue workflow with individual caller VCFs for correct N_RNA_CALLERS_SUPPORT
        VCF_RESCUE_WORKFLOW(
            first_round_dna_consensus_vcf,
            realigned_rna_consensus_vcf,
            dna_caller_vcfs,
            rna_caller_vcfs
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
        second_rescued_vcf_stripped = VCF_RESCUE_POST_PROCESSING.out.vcf_stripped

        // VEP annotation on second rescue filtered stripped VCF
        VCF_ANNOTATE(
            second_rescued_vcf_stripped.map{ meta, vcf, tbi -> [ meta + [ file_name: vcf.baseName ], vcf, [tbi] ] },
            fasta,
            input_sample,
            true,  // realignment = true
            vep_cache
        )
        versions = versions.mix(VCF_ANNOTATE.out.versions)
        reports = reports.mix(VCF_ANNOTATE.out.reports)
        second_rescued_vcf_vep = VCF_ANNOTATE.out.vcf_ann

    emit:
        second_rescued_vcf    = second_rescued_vcf       // Final rescued VCF with annotation (before VEP)
        second_rescued_vcf_stripped = second_rescued_vcf_stripped  // Stripped VCF for VEP input
        second_rescued_vcf_vep = second_rescued_vcf_vep  // VEP-annotated rescue VCF [meta, vcf, tbi]
        versions              = versions
        reports               = reports
}
