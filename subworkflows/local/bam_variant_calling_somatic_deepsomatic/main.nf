//
//
// DEEPSOMATIC: tumor-normal mode variant calling
//

include { DEEPSOMATIC                     } from '../../../modules/nf-core/deepsomatic/main'
include { GATK4_MERGEVCFS as MERGE_DEEPSOMATIC } from '../../../modules/nf-core/gatk4/mergevcfs/main'

workflow BAM_VARIANT_CALLING_SOMATIC_DEEPSOMATIC {
    take:
    input                     // channel: [ meta, [ input ], [ input_index ] ]
    fasta                     // channel: [ meta, fasta]
    fai                       // channel: [ meta, fai]
    dict                      // channel: [ meta, dict]
    intervals                 // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine input and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for DEEPSOMATIC module
        .map{ meta, input_list, input_index_list, intervls, num_intervals -> 
            [ meta + [ num_intervals:num_intervals ], input_list[0], input_index_list[0], input_list[1], input_index_list[1], intervls ] 
        }

    // Perform variant calling using DeepSomatic module
    DEEPSOMATIC( 
        input_intervals.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervls -> 
            [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai ] 
        },
        input_intervals.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervls -> 
            [ meta, intervls ] 
        },
        fasta,
        fai
    )

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = DEEPSOMATIC.out.vcf.branch{
        // Use meta.num_intervals to assess number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = DEEPSOMATIC.out.vcf_tbi.branch{
        // Use meta.num_intervals to assess number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals - merge VCFs
    vcf_to_merge = vcf_branch.intervals.map{ meta, vcf -> 
        [ groupKey(meta, meta.num_intervals), vcf ] 
    }.groupTuple()

    MERGE_DEEPSOMATIC(vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together and remove no longer necessary field: num_intervals
    vcf = Channel.empty().mix(MERGE_DEEPSOMATIC.out.vcf, vcf_branch.no_intervals).map{ meta, vcf ->
        [ meta - meta.subMap('num_intervals'), vcf ]
    }
    
    tbi = Channel.empty().mix(MERGE_DEEPSOMATIC.out.tbi, tbi_branch.no_intervals).map{ meta, tbi ->
        [ meta - meta.subMap('num_intervals'), tbi ]
    }

    versions = versions.mix(DEEPSOMATIC.out.versions)
    versions = versions.mix(MERGE_DEEPSOMATIC.out.versions)

    emit:
    vcf                                       // channel: [ meta, vcf ]
    tbi                                       // channel: [ meta, tbi ]
    
    vcf_filtered = vcf.map{ meta, vcf_file -> 
        [ meta.subMap('id', 'patient', 'status') + [ variantcaller:'deepsomatic' ], vcf_file ] 
    }                                         // channel: [ meta, vcf ] - standardized for downstream
    
    gvcf         = DEEPSOMATIC.out.gvcf.map{ meta, gvcf_file ->
        [ meta - meta.subMap('num_intervals'), gvcf_file ]
    }                                         // channel: [ meta, gvcf ]
    
    gvcf_tbi     = DEEPSOMATIC.out.gvcf_tbi.map{ meta, gvcf_tbi_file ->
        [ meta - meta.subMap('num_intervals'), gvcf_tbi_file ]
    }                                         // channel: [ meta, gvcf_tbi ]
    
    versions                                  // channel: [ versions.yml ]
}