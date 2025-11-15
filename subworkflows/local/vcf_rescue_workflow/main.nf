//
// VCF Rescue Workflow - Cross-modality variant rescue
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
    
    // Cross DNA and RNA consensus VCFs by patient
    dna_rna_pairs = dna_consensus_vcf
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
            [meta, dna[2], dna[3], rna[2], rna[3]]
        }
    
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
    
    // Combine all inputs for rescue
    rescue_input = dna_rna_pairs
        .map { meta, dna_vcf, dna_tbi, rna_vcf, rna_tbi ->
            [meta.patient, meta, dna_vcf, dna_tbi, rna_vcf, rna_tbi]
        }
        .join(dna_callers_grouped)
        .join(rna_callers_grouped)
        .map { patient, meta, dna_vcf, dna_tbi, rna_vcf, rna_tbi,
               dna_vcfs, dna_tbis, dna_callers,
               rna_vcfs, rna_tbis, rna_callers ->
            [meta, dna_vcf, dna_tbi, rna_vcf, rna_tbi,
             dna_vcfs, dna_tbis, rna_vcfs, rna_tbis]
        }
    
    // Run rescue process
    VCF_RESCUE(rescue_input)
    rescued_vcf = VCF_RESCUE.out.vcf
    versions = versions.mix(VCF_RESCUE.out.versions)
    
    emit:
    vcf      = rescued_vcf
    versions = versions
}
