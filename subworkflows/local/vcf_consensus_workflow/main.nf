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

        // Determine expected caller count for grouping
        def ncallers_expected
        
        if (realignment || (params.step in ['consensus', 'annotate','filtering', 'rna_filtering'] &&
                           params.tools && params.tools.split(',').contains("realignment"))) {
            // Realignment mode uses default callers
            def tools_list = params.defaultvariantcallers.split(',').toList()
            ncallers_expected = tools_list.unique().size()
        } else if (params.tools) {
            // Extract actual variant caller names from tools parameter
            def tools_list = params.tools.split(',').toList().findAll { it in ['sage', 'strelka', 'mutect2', 'deepsomatic'] }
            if (tools_list.size() > 0) {
                ncallers_expected = tools_list.unique().size()
            } else {
                // Tools specified but no caller names found - use dynamic grouping
                // This happens when starting from VCFs with --tools containing workflow steps
                ncallers_expected = 0
            }
        } else {
            // No tools specified - use dynamic grouping
            ncallers_expected = 0
        }

        // Align grouping pattern with MAF consensus: two maps then groupTuple
        vcf_grouped = vcf_annotated
            .map { meta, vcf, tbi ->
                // Normalize tbi: some upstream code may wrap tbi in a singleton list
                def tbiFile = (tbi instanceof List && tbi.size()==1) ? tbi[0] : tbi
                // Reduce meta to essential fields and tag data_type
                def metaReduced = meta.subMap('id','patient','status') + [data_type:'vcf', ncallers: ncallers_expected]
                [ metaReduced, vcf, tbiFile, (meta.variantcaller ?: 'unknown') ]
            }
            .map { metaReduced, vcf, tbiFile, variantcaller ->
                // Use groupKey only when we know the expected size (> 0)
                // Otherwise use regular grouping key for dynamic grouping
                if (ncallers_expected > 0) {
                    def key = groupKey(metaReduced.subMap('id','patient','status') + [ncallers: ncallers_expected], ncallers_expected)
                    [ key, vcf, tbiFile, variantcaller ]
                } else {
                    // Dynamic grouping - just use the metadata as key
                    [ metaReduced.subMap('id','patient','status') + [ncallers: 0], vcf, tbiFile, variantcaller ]
                }
            }
            .groupTuple() // [metaGrouped, [vcf...], [tbi...], [caller...]]
            .map { metaGrouped, vcfs, tbis, callers ->
                // Extract metadata - handle both GroupKey and regular Map
                def metaMutable
                if (metaGrouped instanceof nextflow.extension.GroupKey) {
                    // GroupKey from fixed-size grouping
                    metaMutable = metaGrouped.getGroupTarget() instanceof Map ? 
                                  metaGrouped.getGroupTarget().clone() : 
                                  [id: metaGrouped.id, patient: metaGrouped.patient, status: metaGrouped.status, ncallers: metaGrouped.ncallers]
                } else {
                    // Regular Map from dynamic grouping
                    metaMutable = metaGrouped instanceof Map ? metaGrouped.clone() : [:]
                }
                
                // Check for caller count mismatch
                def expected_ncallers = metaMutable.ncallers ?: 0
                def actual_ncallers = callers.size()
                
                if (expected_ncallers == 0) {
                    // Dynamic grouping mode - accept whatever callers are present
                    metaMutable.ncallers = actual_ncallers
                } else if (expected_ncallers != actual_ncallers) {
                    // Mismatch detected - warn but continue
                    metaMutable.ncallers_expected = expected_ncallers
                    metaMutable.ncallers = actual_ncallers
                    println "[CONSENSUS WARN] Caller count mismatch id=${metaMutable.id} expected=${expected_ncallers} actual=${actual_ncallers}" 
                } else {
                    // Expected matches actual - all good
                    metaMutable.ncallers = actual_ncallers
                }
                
                [ metaMutable, vcfs, tbis, callers ]
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
                .flatMap { meta, vcfs, tbis, callers ->
                    vcfs.indices.collect { i ->
                        [meta, vcfs[i], tbis[i], callers[i]]
                    }
                }
            
            rna_caller_vcfs = vcf_grouped
                .filter { it[0].status == 2 }
                .flatMap { meta, vcfs, tbis, callers ->
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
