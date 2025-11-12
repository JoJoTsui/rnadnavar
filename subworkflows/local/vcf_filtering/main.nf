//
// VCF Filtering - parallel to MAF filtering
//
include { VCF_FILTERING as FILTERING_VCF } from '../../../modules/local/vcf_filtering/main'

workflow VCF_FILTERING {

    take:
    vcf_to_filter  // channel: [ [meta], vcf ]
    fasta
    input_sample
    realignment

    main:
    versions  = Channel.empty()
    vcf       = Channel.empty()
    
    if ((params.step in ['mapping', 'markduplicates', 'splitncigar',
                'prepare_recalibration', 'recalibrate', 'variant_calling', 'annotate',
                'norm', 'consensus', 'filtering'] &&
                ((params.tools && params.tools.split(",").contains("filtering")))) ||
                realignment) {

        if (params.step == 'filtering') vcf_to_filter = input_sample
        
        // BASIC VCF FILTERING
        FILTERING_VCF(vcf_to_filter, fasta)
        vcf      = FILTERING_VCF.out.vcf
        versions = versions.mix(FILTERING_VCF.out.versions)
    }

    emit:
    vcf        = vcf
    versions   = versions
}
