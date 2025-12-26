//
// VCF Rescue Workflow - Cross-modality variant rescue
// Simplified version: Only uses consensus VCFs, not individual caller VCFs
//
include { VCF_RESCUE } from '../../../modules/local/vcf_rescue/main'

workflow VCF_RESCUE_WORKFLOW {
    take:
    dna_consensus_vcf      // channel: [ [meta], vcf, tbi ]
    rna_consensus_vcf      // channel: [ [meta], vcf, tbi ]
    dna_caller_vcfs        // channel: [ [meta], vcf, tbi, caller ]
    rna_caller_vcfs        // channel: [ [meta], vcf, tbi, caller ]
    
    main:
    versions = Channel.empty()
    rescued_vcf = Channel.empty()
    
    // Group DNA caller VCFs by patient
    dna_callers_grouped = dna_caller_vcfs
        .map { meta, vcf, tbi, caller -> 
            [meta.patient, vcf, tbi, caller] 
        }
        .groupTuple()
    
    // Group RNA caller VCFs by patient
    rna_callers_grouped = rna_caller_vcfs
        .map { meta, vcf, tbi, caller -> 
            [meta.patient, vcf, tbi, caller] 
        }
        .groupTuple()
    
    // Cross DNA and RNA consensus VCFs by patient, then join with caller VCFs
    // Refactored to use mix/groupTuple to handle empty caller channels gracefully
    
    // 1. Prepare Consensus Data
    ch_consensus_data = dna_consensus_vcf
        .map { meta, vcf, tbi -> [meta.patient, meta, vcf, tbi] }
        .cross(rna_consensus_vcf.map { meta, vcf, tbi -> 
            [meta.patient, meta, vcf, tbi] 
        })
        .map { dna, rna ->
            def meta = [:]
            meta.patient = dna[0]
            meta.dna_id = dna[1].id
            meta.rna_id = rna[1].id
            meta.id = "${meta.dna_id}_rescued_${meta.rna_id}"
            
            [meta.patient, [
                type: 'consensus',
                meta: meta,
                dna_cons_vcf: dna[2],
                dna_cons_tbi: dna[3],
                rna_cons_vcf: rna[2],
                rna_cons_tbi: rna[3]
            ]]
        }

    // 2. Prepare DNA Callers Data
    ch_dna_callers = dna_callers_grouped
        .map { patient, vcfs, tbis, callers ->
            [patient, [
                type: 'dna_callers',
                vcfs: vcfs,
                tbis: tbis,
                callers: callers
            ]]
        }

    // 3. Prepare RNA Callers Data
    ch_rna_callers = rna_callers_grouped
        .map { patient, vcfs, tbis, callers ->
            [patient, [
                type: 'rna_callers',
                vcfs: vcfs,
                tbis: tbis,
                callers: callers
            ]]
        }

    // 4. Mix, Group, and Combine
    rescue_input = ch_consensus_data
        .mix(ch_dna_callers)
        .mix(ch_rna_callers)
        .groupTuple()
        .flatMap { patient, data_list ->
            def dna_data = data_list.find { it.type == 'dna_callers' }
            def rna_data = data_list.find { it.type == 'rna_callers' }
            
            def consensus_items = data_list.findAll { it.type == 'consensus' }
            
            consensus_items.collect { consensus ->
                [
                    consensus.meta,
                    consensus.dna_cons_vcf,
                    consensus.dna_cons_tbi,
                    consensus.rna_cons_vcf,
                    consensus.rna_cons_tbi,
                    dna_data ? dna_data.vcfs : [],
                    dna_data ? dna_data.tbis : [],
                    dna_data ? dna_data.callers : [],
                    rna_data ? rna_data.vcfs : [],
                    rna_data ? rna_data.tbis : [],
                    rna_data ? rna_data.callers : []
                ]
            }
        }
    
    // Run rescue process with consensus AND individual caller VCFs
    VCF_RESCUE(rescue_input)
    rescued_vcf = VCF_RESCUE.out.vcf
    versions = versions.mix(VCF_RESCUE.out.versions)
    
    emit:
    vcf      = rescued_vcf
    versions = versions
}
