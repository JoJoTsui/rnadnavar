//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { DRAGMAP_ALIGN          } from '../../../modules/nf-core/dragmap/align/main'

workflow FASTQ_ALIGN {
    take:
    reads // channel: [mandatory] meta, reads
    index // channel: [mandatory] index
    fasta // channel: [mandatory] fasta
    sort  // boolean: [mandatory] true -> sort, false -> don't sort

    main:

    versions = Channel.empty()
    reports = Channel.empty()

    // Convert index and fasta to proper tuple format
    index_tuple = index.map { idx -> 
        def index_path = idx instanceof List ? idx[0] : idx
        [ [id:'index'], index_path ]
    }
    
    // Convert fasta to proper tuple format
    fasta_tuple = fasta.map { fa -> 
        if (fa instanceof List && fa.size() == 2) {
            // Already in [meta, path] format
            fa
        } else {
            // Single path, add meta
            def fasta_path = fa instanceof List ? fa[0] : fa
            [ [id:'fasta'], fasta_path ]
        }
    }
    
    // Only run the selected aligner based on params.aligner
    if (params.aligner == "bwa-mem") {
        BWAMEM1_MEM(reads, index_tuple, fasta_tuple, sort)
    }
    if (params.aligner == "bwa-mem2") {
        BWAMEM2_MEM(reads, index_tuple, fasta_tuple, sort)
    }
    if (params.aligner == "dragmap") {
        DRAGMAP_ALIGN(reads, index_tuple, fasta_tuple, sort)
    }

    // Get the bam files from the selected aligner
    bam = Channel.empty()
    if (params.aligner == "bwa-mem") {
        bam = bam.mix(BWAMEM1_MEM.out.bam)
        versions = versions.mix(BWAMEM1_MEM.out.versions)
    }
    if (params.aligner == "bwa-mem2") {
        bam = bam.mix(BWAMEM2_MEM.out.bam)
        versions = versions.mix(BWAMEM2_MEM.out.versions)
    }
    if (params.aligner == "dragmap") {
        bam = bam.mix(DRAGMAP_ALIGN.out.bam)
        reports = reports.mix(DRAGMAP_ALIGN.out.log)
        versions = versions.mix(DRAGMAP_ALIGN.out.versions)
    }

    emit:
    bam      // channel: [ [meta], bam ]
    reports
    versions // channel: [ versions.yml ]
}
