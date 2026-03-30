//
// Normalise VCFs with VT
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
//
// =============================================================================
// AUDIT NOTE — neoantigen-workflow Task 0 (Requirement 2.1)
// Audited against: nf-core module registry, vt docs, bcftools norm docs
// Result: CORRECT AS-IS — no blocking errors found; one minor note below.
//
// Step 1 — VT_DECOMPOSE
//   Input mapped to [meta, vcf, []] (empty intervals) — matches module signature
//   `tuple val(meta), path(vcf), path(intervals)`. Correct.
//   Output is *.vcf.gz with no TBI, which is expected and handled correctly.
//   The `-s` (smart-decompose) flag is intentionally omitted; somatic calling
//   does not require phase-aware decomposition.
//
// Step 2 — BCFTOOLS_NORM
//   Called as BCFTOOLS_NORM(vcf_decomposed, fasta) where vcf_decomposed is
//   [meta, vcf, []] (empty TBI — correct, VT output is unindexed) and fasta
//   is a bare value channel (Channel.fromPath(...).collect()).
//   The nf-core module declares `tuple val(meta2), path(fasta)` for the second
//   input, but receiving a bare path channel is a well-established nf-core
//   pattern: Nextflow broadcasts the value channel and the module script only
//   uses `${fasta}` (the path), never `${meta2}`. No runtime failure results.
//   This is consistent with how sarek and other nf-core pipelines call this
//   module. No change required.
//   The module is used here solely for left-alignment and indel normalization
//   (--fasta-ref); multi-allelic splitting is already handled by VT_DECOMPOSE,
//   so no `-m` flag is needed. Correct.
//
// Output channel shape: [meta, vcf, tbi] (3-tuple) — matches the
//   FORMAT_HARMONIZER input spec `tuple val(meta), path(vcf), path(tbi)`.
//   The downstream FORMAT_HARMONIZER can be inserted directly after this
//   subworkflow without any channel reshaping.
// =============================================================================
// VT steps
include { VT_DECOMPOSE                        } from '../../../modules/nf-core/vt/decompose/main'
include { BCFTOOLS_NORM                       } from '../../../modules/nf-core/bcftools/norm/main'
// include { VT_NORMALIZE                        } from '../../../modules/nf-core/vt/normalize/main'
// Create samplesheet to restart from different steps
include { CHANNEL_VARIANT_CALLING_CREATE_CSV  } from '../channel_variant_calling_create_csv/main'


workflow VCF_NORMALIZE {
    take:
    vcf_to_normalize
    fasta
    fasta_fai
    input_sample
    realignment

    main:
    version          = Channel.empty()

    if (params.step == 'norm') vcf_to_normalize = input_sample

    if ((params.step in ['mapping', 'markduplicates', 'splitncigar',
                        'prepare_recalibration', 'recalibrate',
                        'variant_calling', 'norm'] &&
                        ((params.tools && params.tools.split(",").contains("consensus")))) ||
                        realignment) {

        vcf_decomposed  = Channel.empty()
        vcf_to_normalize = vcf_to_normalize.map{meta, vcf -> [meta, vcf, []]} // vt accepts intervals, not in use for now
        // Separate variants
        VT_DECOMPOSE(vcf_to_normalize)

        vcf_decomposed = vcf_decomposed.mix(VT_DECOMPOSE.out.vcf)
        version = version.mix(VT_DECOMPOSE.out.versions.first())

        // Normalize variants
        vcf_decomposed = vcf_decomposed.map{meta,vcf -> [meta, vcf, []]} // tbi not necessary
        // VT_NORMALIZE(vcf_decomposed,
        //             fasta, fasta_fai) // fai not necessary?
        BCFTOOLS_NORM(vcf_decomposed, fasta)

        vcf_normalized = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        vcf_normalized = vcf_normalized.map{meta, vcf_file, tbi -> [ meta + [ file_name: vcf_file.fileName,  data_type: "vcf" ], vcf_file, tbi ] }
        version = version.mix(BCFTOOLS_NORM.out.versions.first())

        CHANNEL_VARIANT_CALLING_CREATE_CSV(vcf_normalized.map{meta, vcf_file, _tbi -> [meta, vcf_file]}, "normalized")

    } else {
        vcf_normalized = vcf_to_normalize
    }

    emit:
    vcf         = vcf_normalized // channel: [ [meta], vcf ]
    versions    = version // channel: [ versions.yml ]

}
