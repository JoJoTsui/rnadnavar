// BAM_REORDER_CONTIGS: Validate and, when safe, normalize an external BAM to
// the exact selected reference sequence dictionary.

include { CHECK_CONTIGS as CHECK_INPUT_DICTIONARY      } from '../../../modules/local/check_contigs/main'
include { CHECK_CONTIGS as CHECK_NORMALIZED_DICTIONARY } from '../../../modules/local/check_contigs/main'
include { PICARD_REORDERSAM                            } from '../../../modules/local/picard_reordersam/main'

workflow BAM_REORDER_CONTIGS {
    take:
    bam           // channel: [ meta, bam, bai ]
    fasta         // channel: path(fasta)
    fasta_fai     // channel: path(fai)
    dict          // channel: path(dict)
    policy        // value: normalize | strict

    main:
    versions = Channel.empty()

    CHECK_INPUT_DICTIONARY(
        bam,
        fasta,
        fasta_fai,
        dict,
        policy,
        'input'
    )
    versions = versions.mix(CHECK_INPUT_DICTIONARY.out.versions)

    bam_branched = CHECK_INPUT_DICTIONARY.out.alignment_with_status.branch { meta, bam_file, bai_file, needs_normalization ->
        needs_normalization: needs_normalization == "true"
            return [ meta, bam_file, bai_file ]
        compatible: needs_normalization == "false"
            return [ meta, bam_file, bai_file ]
    }

    PICARD_REORDERSAM(
        bam_branched.needs_normalization,
        fasta,
        fasta_fai,
        dict
    )
    versions = versions.mix(PICARD_REORDERSAM.out.versions)

    CHECK_NORMALIZED_DICTIONARY(
        PICARD_REORDERSAM.out.bam,
        fasta,
        fasta_fai,
        dict,
        'strict',
        'normalized'
    )
    versions = versions.mix(CHECK_NORMALIZED_DICTIONARY.out.versions)

    bam_out = Channel.empty().mix(
        CHECK_NORMALIZED_DICTIONARY.out.alignment_with_status.map { meta, bam_file, bai_file, ignored -> [ meta, bam_file, bai_file ] },
        bam_branched.compatible
    )
    audits = CHECK_INPUT_DICTIONARY.out.audit.mix(CHECK_NORMALIZED_DICTIONARY.out.audit)

    emit:
    bam      = bam_out       // channel: [ meta, bam, bai ]
    audit    = audits        // channel: [ meta, dictionary_audit.json ]
    versions = versions      // channel: [ versions.yml ]
}
