//
// BAM_REORDER_CONTIGS: Check and optionally reorder BAM contigs to match reference
//
// This subworkflow checks if BAM/CRAM files have contigs not present in the reference
// and reorders them using Picard ReorderSam if needed.
//

include { CHECK_CONTIGS      } from '../../../modules/local/check_contigs/main'
include { PICARD_REORDERSAM  } from '../../../modules/local/picard_reordersam/main'

workflow BAM_REORDER_CONTIGS {
    take:
    bam           // channel: [ meta, bam, bai ]
    fasta         // channel: path(fasta)
    fasta_fai     // channel: path(fai)
    dict          // channel: path(dict)

    main:
    versions = Channel.empty()

    // Check if BAM contigs match reference
    CHECK_CONTIGS(
        bam,
        fasta,
        fasta_fai,
        dict
    )
    versions = versions.mix(CHECK_CONTIGS.out.versions)

    // Branch based on whether reordering is needed
    bam_branched = CHECK_CONTIGS.out.bam_with_status.branch { meta, bam_file, bai_file, needs_reorder ->
        needs_reorder: needs_reorder == "true"
            return [ meta, bam_file, bai_file ]
        no_reorder: needs_reorder == "false"
            return [ meta, bam_file, bai_file ]
    }

    // Reorder BAMs that need it
    PICARD_REORDERSAM(
        bam_branched.needs_reorder,
        fasta,
        fasta_fai,
        dict
    )
    versions = versions.mix(PICARD_REORDERSAM.out.versions)

    // Combine reordered and non-reordered BAMs
    bam_out = Channel.empty().mix(
        PICARD_REORDERSAM.out.bam,
        bam_branched.no_reorder
    )

    emit:
    bam      = bam_out      // channel: [ meta, bam, bai ]
    versions = versions     // channel: [ versions.yml ]
}
