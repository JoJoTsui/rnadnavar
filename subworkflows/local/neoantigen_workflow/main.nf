//
// Neoantigen Workflow — Salmon quantification + neoantigen VCF selection
//
// Runs salmon quant on RNA tumor FASTQs (status=2) in quasi-mapping mode,
// normalizes transcript names, and selects the neoantigen-ready VCF for
// DNA tumor samples (status=1) based on params.neoantigen_input_source.
//
include { SALMON_QUANT        } from '../../../modules/local/salmon_quant/main'
include { QUANT_TSV_NORMALIZE } from '../../../modules/local/quant_tsv_normalize/main'

// Thin process to publish the neoantigen VCF + index to the output directory.
// publishDir cannot be set on a subworkflow directly, so a dedicated process is used.
process PUBLISH_NEOANTIGEN_VCF {
    tag "$meta.id"
    label 'process_single'

    publishDir "${params.outdir}/neoantigen/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${vcf}"), path("${tbi}")

    script:
    """
    """
}

workflow NEOANTIGEN_WORKFLOW {
    take:
    ch_fastq        // channel: [ [meta], [reads] ] — all samples
    ch_salmon_index // channel: path(index_dir)
    ch_consensus_vcf  // channel: [ [meta], vcf, tbi ] — all statuses
    ch_mutect2_vcf    // channel: [ [meta], vcf, tbi ] — all statuses

    main:
    versions = Channel.empty()

    // Filter to RNA tumor samples (status=2) for Salmon quantification
    ch_rna_fastq = ch_fastq.filter { meta, reads -> meta.status == 2 }

    // Run Salmon quasi-mapping quantification
    SALMON_QUANT(ch_rna_fastq, ch_salmon_index)

    // Normalize transcript names (strip pipe-delimited Gencode suffixes)
    QUANT_TSV_NORMALIZE(SALMON_QUANT.out.quant_sf)

    // Select VCF source based on neoantigen_input_source parameter
    ch_source_vcf = params.neoantigen_input_source == 'consensus' ? ch_consensus_vcf : ch_mutect2_vcf

    // Filter to DNA tumor samples (status=1) for neoantigen VCF output
    ch_neoantigen_vcf = ch_source_vcf.filter { meta, vcf, tbi -> meta.status == 1 }

    // Warn when no status=1 samples are present
    ch_neoantigen_vcf
        .ifEmpty { log.warn "NEOANTIGEN_WORKFLOW: No status=1 (DNA tumor) samples found. No neoantigen VCF will be published." }

    // Publish neoantigen VCF to ${outdir}/neoantigen/<sample_id>/
    PUBLISH_NEOANTIGEN_VCF(ch_neoantigen_vcf)

    // Collect software versions
    versions = versions.mix(SALMON_QUANT.out.versions)
    versions = versions.mix(QUANT_TSV_NORMALIZE.out.versions)

    emit:
    neoantigen_vcf = ch_neoantigen_vcf  // channel: [ [meta], vcf, tbi ] — status=1 only
    versions       = versions
}
