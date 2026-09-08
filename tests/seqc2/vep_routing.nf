// Exercise both annotation consumers with the same collected cache resource.
include { VCF_ANNOTATE as FIRST_PASS; VCF_ANNOTATE as SECOND_PASS } from '../../subworkflows/local/vcf_annotate/main'

workflow {
    cache = Channel.fromPath(params.vep_cache).collect()
    fasta = Channel.value([[id: 'reference'], file(params.fasta)])
    vcfs = Channel.of('RT1', 'RT2').map { id ->
        [[id: id, patient: 'patient', status: 2, variantcaller: 'rescue'], file(params.input), []]
    }
    FIRST_PASS(vcfs, fasta, Channel.empty(), false, cache)
    SECOND_PASS(vcfs, fasta, Channel.empty(), params.test_realignment, cache)
    SECOND_PASS.out.vcf_ann.toList().view { "ANNOTATED_SAMPLES=${it.collect { row -> row[0].id }.sort()}" }
}
