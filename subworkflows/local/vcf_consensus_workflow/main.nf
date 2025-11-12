//
// VCF Consensus Workflow - Standalone VCF processing
//
include { VCF_CONSENSUS                            } from '../../../modules/local/vcf_consensus/main'
include { VCF_CONSENSUS as VCF_CONSENSUS_RESCUE    } from '../../../modules/local/vcf_consensus/main'

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
        
        // RESCUE: Cross DNA/RNA for second consensus
        if (params.tools && params.tools.split(',').contains('rescue')) {
            // Separate DNA and RNA
            vcf_dna = consensus_vcf.filter { it[0].status <= 1 }
            vcf_rna = consensus_vcf.filter { it[0].status == 2 }
            
            vcf_grouped_dna = vcf_grouped.filter { it[0].status <= 1 }
            vcf_grouped_rna = vcf_grouped.filter { it[0].status == 2 }
            
            // Cross DNA variants with RNA consensus
            dna_with_rna = vcf_grouped_dna
                .map { meta, vcfs, tbis, callers -> [meta.patient, meta, vcfs, tbis, callers] }
                .cross(vcf_rna.map { meta, vcf, tbi -> [meta.patient, meta, vcf, tbi] })
                .map { dna, rna ->
                    def meta = [:]
                    meta.patient = dna[0]
                    meta.dna_id = dna[1].id
                    meta.rna_id = rna[1].id
                    meta.status = dna[1].status
                    meta.id = "${meta.dna_id}_with_${meta.rna_id}".toString()
                    [meta, dna[2] + [rna[2]], dna[3] + [rna[3]], dna[4] + ['consensus']]
                }
            
            // Cross RNA variants with DNA consensus
            rna_with_dna = vcf_grouped_rna
                .map { meta, vcfs, tbis, callers -> [meta.patient, meta, vcfs, tbis, callers] }
                .cross(vcf_dna.map { meta, vcf, tbi -> [meta.patient, meta, vcf, tbi] })
                .map { rna, dna ->
                    def meta = [:]
                    meta.patient = rna[0]
                    meta.rna_id = rna[1].id
                    meta.dna_id = dna[1].id
                    meta.status = rna[1].status
                    meta.id = "${meta.rna_id}_with_${meta.dna_id}".toString()
                    [meta, rna[2] + [dna[2]], rna[3] + [dna[3]], rna[4] + ['consensus']]
                }
            
            // Run rescue consensus
            VCF_CONSENSUS_RESCUE(dna_with_rna.mix(rna_with_dna))
            consensus_vcf_rescue = VCF_CONSENSUS_RESCUE.out.vcf
            versions = versions.mix(VCF_CONSENSUS_RESCUE.out.versions)
            
            // Use rescue results as final consensus
            consensus_vcf = consensus_vcf_rescue
        }
    }

    emit:
    vcf      = consensus_vcf        // channel: [ [meta], vcf, tbi ]
    versions = versions             // channel: [ versions.yml ]
}
