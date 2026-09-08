// Exercise the real resource producer with a consumer per logical RNA sample.
include { PREPARE_GENOME } from '../../subworkflows/local/prepare_genome/main'

process USE_RESOURCES {
    input:
    val sample
    tuple val(index_meta), path(index)
    tuple val(splice_meta), path(splices)
    output:
    val sample
    script:
    """
    test -s ${splices}
    """
}

workflow {
    PREPARE_GENOME(
        Channel.empty(), Channel.value([file(params.fasta)]),
        Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty()
    )
    if (params.tools.contains('realignment')) {
        USE_RESOURCES(Channel.of('RT1', 'RT2'),
            PREPARE_GENOME.out.hisat2_index, PREPARE_GENOME.out.splicesites)
        USE_RESOURCES.out.toList().map { samples ->
            assert samples.sort() == ['RT1', 'RT2']: "Resource routing dropped samples: ${samples}"
            samples
        }.view { "RESOURCE_SAMPLES=${it}" }
    }
}
