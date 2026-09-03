//
// VCF Consensus Workflow - Standalone VCF processing
//
include { VCF_CONSENSUS                            } from '../../../modules/local/vcf_consensus/main'
include { VCF_RESCUE_WORKFLOW                      } from '../vcf_rescue_workflow/main'

workflow VCF_CONSENSUS_WORKFLOW {
    take:
    vcf_annotated      // channel: [ [meta], vcf, tbi ]
    input_sample
    realignment

    main:
    versions = Channel.empty()
    consensus_vcf = Channel.empty()
    consensus_vcf_rescue = Channel.empty()

    if ((params.step in ['mapping', 'markduplicates', 'splitncigar',
                        'prepare_recalibration', 'recalibrate', 'variant_calling', 'annotate',
                        'norm', 'consensus'] &&
                        ((params.tools && params.tools.split(",").contains("consensus")))) ||
                        realignment) {

        def active_callers = ['mutect2', 'strelka', 'deepsomatic']
        def expected_callers
        if (params.consensus_expected_callers) {
            expected_callers = params.consensus_expected_callers.split(',')
                .collect { it.trim().toLowerCase() }
                .findAll { it }
        } else {
            def configured = (realignment ||
                (params.step in ['consensus', 'annotate', 'filtering', 'rna_filtering'] &&
                 params.tools && params.tools.split(',').contains('realignment')))
                ? params.defaultvariantcallers
                : params.tools
            expected_callers = (configured ?: params.defaultvariantcallers ?: '')
                .split(',')
                .collect { it.trim().toLowerCase() }
                .findAll { it in active_callers }
        }
        if (!expected_callers || expected_callers.size() != expected_callers.unique().size()) {
            error "Consensus requires a non-empty, duplicate-free expected caller panel; got ${expected_callers}"
        }
        def unsupported_callers = expected_callers.findAll { !(it in active_callers) }
        if (unsupported_callers) {
            error "Unsupported callers in consensus expected panel: ${unsupported_callers}; active callers are ${active_callers}"
        }
        expected_callers = expected_callers.unique()
        def ncallers_expected = expected_callers.size()

        // Group to channel close, then validate the full configured panel.
        // A fixed-size groupKey can silently suppress incomplete groups.
        vcf_grouped = vcf_annotated
            .map { meta, vcf, tbi ->
                // Normalize tbi: some upstream code may wrap tbi in a singleton list
                def tbiFile = (tbi instanceof List && tbi.size()==1) ? tbi[0] : tbi
                // Reduce meta to essential fields and tag data_type
                def metaReduced = meta.subMap('id','patient','status') +
                    [data_type:'vcf', ncallers: ncallers_expected, expected_callers: expected_callers]
                [ metaReduced, vcf, tbiFile, (meta.variantcaller ?: 'unknown') ]
            }
            .map { metaReduced, vcf, tbiFile, variantcaller ->
                def key = metaReduced.subMap('id','patient','status') +
                    [ncallers: ncallers_expected, expected_callers: expected_callers]
                [ key, vcf, tbiFile, variantcaller.toString().toLowerCase() ]
            }
            .groupTuple() // [metaGrouped, [vcf...], [tbi...], [caller...]]
            .map { metaGrouped, vcfs, tbis, callers ->
                def metaMutable = metaGrouped.clone()
                def actual = callers.collect { it.toLowerCase() }
                def duplicates = actual.findAll { caller -> actual.count(caller) > 1 }.unique()
                def missing = expected_callers - actual
                def unexpected = actual - expected_callers
                if (duplicates || missing || unexpected || actual.size() != expected_callers.size()) {
                    error "Consensus caller panel mismatch id=${metaMutable.id} status=${metaMutable.status} " +
                        "expected=${expected_callers} actual=${actual} missing=${missing} " +
                        "unexpected=${unexpected} duplicates=${duplicates}"
                }
                metaMutable.ncallers = expected_callers.size()
                [ metaMutable, vcfs, tbis, actual, expected_callers ]
            }

        vcf_grouped.dump(tag:"vcf_grouped_for_consensus")
        
        // Run consensus
        VCF_CONSENSUS(vcf_grouped)
        consensus_vcf = VCF_CONSENSUS.out.vcf
        versions = versions.mix(VCF_CONSENSUS.out.versions)
        
        // RESCUE: Cross-modality variant rescue
        if (params.tools && params.tools.split(',').contains('rescue')) {
            // Separate DNA (status ≤ 1) and RNA (status = 2) consensus VCFs
            dna_consensus_vcf = consensus_vcf.filter { it[0].status <= 1 }
            rna_consensus_vcf = consensus_vcf.filter { it[0].status == 2 }
            
            // Separate DNA and RNA caller VCFs from grouped input
            dna_caller_vcfs = vcf_grouped
                .filter { it[0].status <= 1 }
                .flatMap { meta, vcfs, tbis, callers, expected ->
                    vcfs.indices.collect { i ->
                        [meta, vcfs[i], tbis[i], callers[i]]
                    }
                }
            
            rna_caller_vcfs = vcf_grouped
                .filter { it[0].status == 2 }
                .flatMap { meta, vcfs, tbis, callers, expected ->
                    vcfs.indices.collect { i ->
                        [meta, vcfs[i], tbis[i], callers[i]]
                    }
                }
            
            // Invoke VCF_RESCUE_WORKFLOW with consensus AND individual caller VCFs
            VCF_RESCUE_WORKFLOW(
                dna_consensus_vcf,
                rna_consensus_vcf,
                dna_caller_vcfs,
                rna_caller_vcfs
            )
            
            // Keep rescue results separate from consensus
            consensus_vcf_rescue = VCF_RESCUE_WORKFLOW.out.vcf
            versions = versions.mix(VCF_RESCUE_WORKFLOW.out.versions)
        }
    }

    emit:
    vcf        = consensus_vcf        // channel: [ [meta], vcf, tbi ]
    vcf_rescue = consensus_vcf_rescue // channel: [ [meta], vcf, tbi ] - rescued VCFs only
    versions   = versions             // channel: [ versions.yml ]
}
