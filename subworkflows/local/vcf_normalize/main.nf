//
// Normalise VCFs with VT
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
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
