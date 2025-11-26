//
// VCF Filtering Subworkflow
//
include { VCF_FILTERING as VCF_FILTER } from '../../../modules/local/vcf_filtering/main'

workflow VCF_FILTERING {
    take:
    vcf_to_filter      // channel: [ [meta], vcf, tbi ]
    fasta              // channel: [ [meta], fasta ]
    input_sample
    realignment

    main:
    versions = Channel.empty()
    filtered_vcf = Channel.empty()
    filtered_vcf_stripped = Channel.empty()

    if ((params.step in ['mapping', 'markduplicates', 'splitncigar',
                        'prepare_recalibration', 'recalibrate', 'variant_calling', 'annotate',
                        'norm', 'consensus', 'filtering'] &&
                        ((params.tools && params.tools.split(",").contains("filtering")))) ||
                        realignment) {

        // Prepare whitelist and blacklist
        whitelist = params.whitelist ? file(params.whitelist) : []
        blacklist = params.blacklist ? file(params.blacklist) : []

        // Run filtering
        VCF_FILTER(vcf_to_filter, fasta, whitelist, blacklist)
        filtered_vcf          = VCF_FILTER.out.vcf
        filtered_vcf_stripped = VCF_FILTER.out.vcf_stripped
        versions              = versions.mix(VCF_FILTER.out.versions)
    }

    emit:
    vcf          = filtered_vcf           // channel: [ [meta], vcf, tbi ]
    vcf_stripped = filtered_vcf_stripped  // channel: [ [meta], vcf, tbi ]
    versions     = versions               // channel: [ versions.yml ]
}
