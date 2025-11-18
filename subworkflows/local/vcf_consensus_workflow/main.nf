//
// VCF Consensus Workflow - Standalone VCF processing
//
include { VCF_CONSENSUS                            } from '../../../modules/local/vcf_consensus/main'
include { VCF_CONSENSUS as VCF_CONSENSUS_RESCUE    } from '../../../modules/local/vcf_consensus/main'
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

        // Determine tools list for grouping
        if (realignment || (params.step in ['consensus', 'annotate','filtering', 'rna_filtering'] && 
                           params.tools && params.tools.split(',').contains("realignment"))) {
            tools_list = params.defaultvariantcallers.split(',')
        } else {
            tools_list = params.tools.split(',').findAll { it in ['sage', 'strelka', 'mutect2', 'deepsomatic'] }
        }
        
        // Group VCFs by sample for consensus
        vcf_grouped = vcf_annotated.map { meta, vcf, tbi ->
                                def key = meta.subMap('id', 'patient', 'status')
                                [key, vcf, tbi, meta.variantcaller]
                            }.map { meta, vcf, tbi, variantcaller ->
                                def ncallers = tools_list.unique().size()
                                def key = groupKey(meta + [ncallers: ncallers], ncallers)
                                [key, vcf, tbi, variantcaller]
                            }.groupTuple()
        
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
            
            // Use rescue results as final consensus
            consensus_vcf_rescue = VCF_RESCUE_WORKFLOW.out.vcf
            versions = versions.mix(VCF_RESCUE_WORKFLOW.out.versions)
            
            // Optionally replace consensus with rescue results
            // For now, emit both - users can choose which to use downstream
            consensus_vcf = consensus_vcf.mix(consensus_vcf_rescue)
        }
    }

    emit:
    vcf      = consensus_vcf        // channel: [ [meta], vcf, tbi ]
    versions = versions             // channel: [ versions.yml ]
}
