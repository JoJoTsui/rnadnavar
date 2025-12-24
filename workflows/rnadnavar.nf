/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                                                   } from '../modules/nf-core/multiqc'
include { samplesheetToList                                         } from 'plugin/nf-schema'
include { paramsSummaryMap                                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                    } from '../subworkflows/local/utils_nfcore_rnadnavar_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
// Build the genome index and other reference files
include { SAMPLESHEET_TO_CHANNEL                                    } from '../subworkflows/local/samplesheet_to_channel'
include { PREPARE_REFERENCE_AND_INTERVALS                           } from '../subworkflows/local/prepare_reference_and_intervals'
include { PREPARE_INTERVALS as PREPARE_INTERVALS_FOR_REALIGNMENT    } from '../subworkflows/local/prepare_intervals'
// Download annotation cache if needed
include { ENSEMBLVEP_DOWNLOAD                                       } from '../modules/nf-core/ensemblvep/download'
include { UNZIP as UNZIP_VEP_CACHE                                  } from '../modules/nf-core/unzip'
include { ANNOTATION_CACHE_INITIALISATION                           } from '../subworkflows/local/annotation_cache_initialisation'

// Alignment
include { BAM_ALIGN                                                 } from '../subworkflows/local/bam_align'

// Core subworkflows of the pipeline
include { BAM_VARIANT_CALLING_PRE_POST_PROCESSING as BAM_PROCESSING } from '../subworkflows/local/bam_variant_calling_pre_post_processing'

// Second run
include { BAM_EXTRACT_READS_HISAT2_ALIGN as PREPARE_REALIGNMENT     } from '../subworkflows/local/prepare_realignment'
include { BAM_VARIANT_CALLING_PRE_POST_PROCESSING as REALIGNMENT    } from '../subworkflows/local/bam_variant_calling_pre_post_processing'

// VCF-based realignment (new)
include { BAM_EXTRACT_READS_HISAT2_ALIGN_VCF as PREPARE_REALIGNMENT_VCF } from '../subworkflows/local/prepare_realignment_vcf'
include { RNA_REALIGNMENT_WORKFLOW                                  } from '../subworkflows/local/rna_realignment'
include { SECOND_RESCUE_WORKFLOW                                    } from '../subworkflows/local/second_rescue'

// Filter RNA
include { MAF_FILTERING_RNA                                         } from '../subworkflows/local/maf_rna_filtering'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNADNAVAR {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    // To gather all QC reports for MultiQC
    ch_multiqc_files = Channel.empty()
    multiqc_report = Channel.empty()
    reports = Channel.empty()
    // To gather used softwares versions for MultiQC
    versions = Channel.empty()

    // Set input, can either be from --input or from automatic retrieval in utils_nfcore_rnadnavar_pipeline
    if (params.input) {
        ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
    }
    else {
        ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
    }
    // Parse samplesheet
    SAMPLESHEET_TO_CHANNEL(ch_from_samplesheet)

    input_sample = SAMPLESHEET_TO_CHANNEL.out.input_sample


    // Initialise MULTIQC
    multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)

    // Download cache if needed
    if (params.download_cache) {
        ensemblvep_info = Channel.of([[id: "${params.vep_cache_version}_${params.vep_genome}"], params.vep_genome, params.vep_species, params.vep_cache_version])
        ENSEMBLVEP_DOWNLOAD(ensemblvep_info)
        vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect().map { meta, cache -> [cache] }
        versions = versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)
    }
    else if (params.vep_cache && params.vep_cache.endsWith(".zip")) {
        UNZIP_VEP_CACHE(Channel.fromPath(params.vep_cache).collect().map { it -> [[id: it[0].baseName], it] })
        vep_cache = UNZIP_VEP_CACHE.out.unzipped_archive.map { it[1] }
        versions = versions.mix(UNZIP_VEP_CACHE.out.versions)
    }
    else {
        // Assuming that if the cache is provided, the user has already downloaded it
        ANNOTATION_CACHE_INITIALISATION(
            (params.vep_cache && params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
            params.vep_cache,
            params.vep_species,
            params.vep_cache_version,
            params.vep_genome,
            params.vep_custom_args,
            "Please refer to https://nf-co.re/sarek/docs/usage/#how-to-customise-vep-annotation for more information.",
        )
        vep_cache = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache.map { it[1] }
    }

    // STEP 0: Build reference and indices if needed
    PREPARE_REFERENCE_AND_INTERVALS()
    versions = versions.mix(PREPARE_REFERENCE_AND_INTERVALS.out.versions)

    // Reference and intervals variables
    fasta = PREPARE_REFERENCE_AND_INTERVALS.out.fasta
    fasta_fai = PREPARE_REFERENCE_AND_INTERVALS.out.fasta_fai
    dict = PREPARE_REFERENCE_AND_INTERVALS.out.dict
    germline_resource = PREPARE_REFERENCE_AND_INTERVALS.out.germline_resource
    germline_resource_tbi = PREPARE_REFERENCE_AND_INTERVALS.out.germline_resource_tbi
    intervals = PREPARE_REFERENCE_AND_INTERVALS.out.intervals
    intervals_for_preprocessing = PREPARE_REFERENCE_AND_INTERVALS.out.intervals_for_preprocessing
    // specific for variant calling
    intervals_bed_combined = PREPARE_REFERENCE_AND_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi = PREPARE_REFERENCE_AND_INTERVALS.out.intervals_bed_gz_tbi
    intervals_bed_gz_tbi_combined = PREPARE_REFERENCE_AND_INTERVALS.out.intervals_bed_gz_tbi_combined
    dbsnp = PREPARE_REFERENCE_AND_INTERVALS.out.dbsnp
    dbsnp_tbi = PREPARE_REFERENCE_AND_INTERVALS.out.dbsnp_tbi
    pon = PREPARE_REFERENCE_AND_INTERVALS.out.pon
    pon_tbi = PREPARE_REFERENCE_AND_INTERVALS.out.pon_tbi
    known_sites_indels = PREPARE_REFERENCE_AND_INTERVALS.out.known_sites_indels
    known_sites_indels_tbi = PREPARE_REFERENCE_AND_INTERVALS.out.known_sites_indels_tbi
    known_sites_snps = PREPARE_REFERENCE_AND_INTERVALS.out.known_sites_snps
    known_sites_snps_tbi = PREPARE_REFERENCE_AND_INTERVALS.out.known_sites_snps_tbi

    intervals_and_num_intervals = intervals.map { interval, num_intervals ->
        if (num_intervals < 1) {
            [[], num_intervals]
        }
        else {
            [interval, num_intervals]
        }
    }
    // STEP 1: ALIGNMENT PREPROCESSING
    BAM_ALIGN(
        PREPARE_REFERENCE_AND_INTERVALS.out.bwa,
        PREPARE_REFERENCE_AND_INTERVALS.out.bwamem2,
        PREPARE_REFERENCE_AND_INTERVALS.out.dragmap,
        PREPARE_REFERENCE_AND_INTERVALS.out.star_index,
        PREPARE_REFERENCE_AND_INTERVALS.out.gtf,
        fasta,
        fasta_fai,
        input_sample,
    )

    reports = reports.mix(BAM_ALIGN.out.reports)
    versions = versions.mix(BAM_ALIGN.out.versions)
    // RNA editing parameters
    rediportal_vcf          = params.rediportal_vcf ? Channel.fromPath(params.rediportal_vcf).collect() : Channel.empty()
    rediportal_tbi          = params.rediportal_tbi ? Channel.fromPath(params.rediportal_tbi).collect() : Channel.empty()
    min_rna_support         = params.min_rna_support ?: 2
    enable_rna_annotation   = params.enable_rna_annotation ?: false

    // COSMIC/gnomAD annotation parameters
    cosmic_vcf                      = params.cosmic_database ? Channel.fromPath(params.cosmic_database).collect() : Channel.empty()
    cosmic_tbi                      = params.cosmic_database ? Channel.fromPath(params.cosmic_database + '.tbi').collect() : Channel.empty()
    gnomad_dir                      = params.gnomad_database ? Channel.fromPath(params.gnomad_database).collect() : Channel.empty()
    enable_cosmic_gnomad_annotation = params.enable_cosmic_gnomad_annotation ?: true
    cosmic_gnomad_verbose           = params.cosmic_gnomad_verbose ?: false

    // 5 MAIN STEPS: GATK PREPROCESING - VARIANT CALLING - NORMALIZATION - CONSENSUS - ANNOTATION
    BAM_PROCESSING(
        input_sample,
        BAM_ALIGN.out.bam_mapped,
        BAM_ALIGN.out.cram_mapped,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi,
        pon,
        pon_tbi,
        known_sites_indels,
        known_sites_indels_tbi,
        germline_resource,
        germline_resource_tbi,
        intervals,
        intervals_for_preprocessing,
        intervals_bed_gz_tbi,
        intervals_bed_combined,
        intervals_and_num_intervals,
        intervals_bed_gz_tbi_combined,
        vep_cache,
        null,
        null,
        false,
        params.no_intervals,
        rediportal_vcf,
        rediportal_tbi,
        min_rna_support,
        enable_rna_annotation,
        cosmic_vcf,
        cosmic_tbi,
        gnomad_dir,
        enable_cosmic_gnomad_annotation,
        cosmic_gnomad_verbose,
    )
    filtered_maf = BAM_PROCESSING.out.maf
    reports = reports.mix(BAM_PROCESSING.out.reports)
    versions = versions.mix(BAM_PROCESSING.out.versions)

    // === Store first-round outputs for second rescue ===
    // DNA consensus VCF (executed once only, status <= 1 for DNA normal+tumor)
    first_round_dna_consensus_vcf = BAM_PROCESSING.out.vcf.filter { it[0].status <= 1 }
    // First round rescued VCF (for comparison)
    first_round_rescued_vcf = BAM_PROCESSING.out.vcf_rescue

    // === Realignment branch ===
    realigned_rna_consensus_vcf = Channel.empty()
    realigned_filtered_rna_vcf = Channel.empty()
    second_rescued_vcf = Channel.empty()
    realigned_filtered_maf = Channel.empty()

    if (params.tools && params.tools.split(',').contains('realignment')) {
        // fastq will not be split when realignment
        params.split_fastq = 0
        // reset intervals to none (realignment files are small)
        PREPARE_INTERVALS_FOR_REALIGNMENT(fasta_fai, null, true)

        def mode = params.realignment_mode ?: 'maf'

        // === Step 1: Prepare realignment (extract reads + HISAT2) ===
        if (mode == 'vcf') {
            // VCF-based realignment: use RNA consensus VCF as input
            vcf_for_realignment = BAM_PROCESSING.out.vcf
                .filter { it[0].status == 2 }  // RNA samples only
                .map { meta, vcf, tbi -> [meta, vcf, tbi] }

            PREPARE_REALIGNMENT_VCF(
                input_sample,
                vcf_for_realignment,
                BAM_PROCESSING.out.cram_variant_calling,
                fasta,
                fasta_fai,
                dict,
                PREPARE_REFERENCE_AND_INTERVALS.out.hisat2_index,
                PREPARE_REFERENCE_AND_INTERVALS.out.splicesites,
                BAM_PROCESSING.out.dna_consensus_maf,
                BAM_PROCESSING.out.dna_varcall_mafs,
            )
            versions = versions.mix(PREPARE_REALIGNMENT_VCF.out.versions)
            realigned_bam = PREPARE_REALIGNMENT_VCF.out.bam_mapped
        } else {
            // MAF-based realignment (existing)
            PREPARE_REALIGNMENT(
                input_sample,
                filtered_maf,
                BAM_PROCESSING.out.cram_variant_calling,
                fasta,
                fasta_fai,
                dict,
                PREPARE_REFERENCE_AND_INTERVALS.out.hisat2_index,
                PREPARE_REFERENCE_AND_INTERVALS.out.splicesites,
                BAM_PROCESSING.out.dna_consensus_maf,
                BAM_PROCESSING.out.dna_varcall_mafs,
            )
            versions = versions.mix(PREPARE_REALIGNMENT.out.versions)
            realigned_bam = PREPARE_REALIGNMENT.out.bam_mapped
        }

        // === Step 2: RNA variant calling on realigned BAM ===
        if (mode == 'vcf') {
            // VCF mode: use dedicated RNA realignment workflow
            RNA_REALIGNMENT_WORKFLOW(
                input_sample,
                realigned_bam,
                fasta,
                fasta_fai,
                dict,
                dbsnp,
                dbsnp_tbi,
                pon,
                pon_tbi,
                known_sites_indels,
                known_sites_indels_tbi,
                germline_resource,
                germline_resource_tbi,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_for_preprocessing,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed_gz_tbi,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed_combined,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_and_num_intervals,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed_gz_tbi_combined,
                true  // realignment = true
            )
            versions = versions.mix(RNA_REALIGNMENT_WORKFLOW.out.versions)
            reports = reports.mix(RNA_REALIGNMENT_WORKFLOW.out.reports)

            // Capture realigned RNA outputs
            realigned_rna_consensus_vcf = RNA_REALIGNMENT_WORKFLOW.out.rna_consensus_vcf
            realigned_filtered_rna_vcf = RNA_REALIGNMENT_WORKFLOW.out.filtered_rna_vcf

            // === Step 3: Second rescue (DNA + realigned RNA) ===
            if (params.tools && params.tools.split(',').contains('rescue')) {
                SECOND_RESCUE_WORKFLOW(
                    first_round_dna_consensus_vcf,
                    realigned_rna_consensus_vcf,
                    Channel.empty(),  // DNA caller VCFs (consensus-only for now)
                    Channel.empty(),  // RNA caller VCFs (consensus-only for now)
                    fasta,
                    rediportal_vcf,
                    rediportal_tbi,
                    min_rna_support,
                    enable_rna_annotation,
                    cosmic_vcf,
                    cosmic_tbi,
                    gnomad_dir,
                    enable_cosmic_gnomad_annotation,
                    cosmic_gnomad_verbose
                )
                versions = versions.mix(SECOND_RESCUE_WORKFLOW.out.versions)
                second_rescued_vcf = SECOND_RESCUE_WORKFLOW.out.second_rescued_vcf
            }
        } else {
            // MAF mode: use existing realignment workflow
            REALIGNMENT(
                Channel.empty(),
                realigned_bam,
                Channel.empty(),
                fasta,
                fasta_fai,
                dict,
                dbsnp,
                dbsnp_tbi,
                pon,
                pon_tbi,
                known_sites_indels,
                known_sites_indels_tbi,
                germline_resource,
                germline_resource_tbi,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_for_preprocessing,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed_gz_tbi,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed_combined,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_and_num_intervals,
                PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed_gz_tbi_combined,
                vep_cache,
                PREPARE_REALIGNMENT.out.dna_consensus_maf,
                PREPARE_REALIGNMENT.out.dna_varcall_mafs,
                true,
                true,
                rediportal_vcf,
                rediportal_tbi,
                min_rna_support,
                enable_rna_annotation,
                cosmic_vcf,
                cosmic_tbi,
                gnomad_dir,
                enable_cosmic_gnomad_annotation,
                cosmic_gnomad_verbose,
            )
            reports = reports.mix(REALIGNMENT.out.reports)
            versions = versions.mix(REALIGNMENT.out.versions)
            realigned_filtered_maf = REALIGNMENT.out.maf
        }
    }
    filtered_maf.dump(tag: "filtered_maf")
    realigned_filtered_maf.dump(tag: "realigned_filtered_maf")
    // MAF filtering only needed for MAF mode (VCF mode is VCF-only workflow)
    if (params.realignment_mode == 'maf') {
        MAF_FILTERING_RNA(
            filtered_maf.branch {
                dna: it[0].status < 2
                rna: it[0].status >= 2
            }.rna,
            realigned_filtered_maf.branch {
                dna: it[0].status < 2
                rna: it[0].status >= 2
            }.rna,
            fasta,
            fasta_fai,
            input_sample,
        )
        versions = versions.mix(MAF_FILTERING_RNA.out.versions)
    }
    //
    // REPORTING
    //
    version_yaml = Channel.empty()
    if (!(params.skip_tools && params.skip_tools.split(',').contains('versions'))) {
        version_yaml = softwareVersionsToYAML(versions).collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_rnadnavar_software_mqc_versions.yml', sort: true, newLine: true)
    }

    if (!(params.skip_tools && params.skip_tools.split(',').contains('multiqc'))) {
        ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(version_yaml)
        ch_multiqc_files = ch_multiqc_files.mix(reports)
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        MULTIQC(
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            [],
        )
        multiqc_report = MULTIQC.out.report.toList()
    }

    emit:
    multiqc_report = multiqc_report // channel: /path/to/multiqc_report.html
    versions       = versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        def InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        def Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        def BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    def String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    }
    else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}
