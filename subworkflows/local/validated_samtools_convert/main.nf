//
// VALIDATED_SAMTOOLS_CONVERT: Samtools convert with command validation
//

include { COMMAND_VALIDATION } from '../command_validation/main'
include { SAMTOOLS_CONVERT_ENHANCED } from '../../../modules/local/samtools_convert_enhanced/main'

workflow VALIDATED_SAMTOOLS_CONVERT {
    take:
        input_channel   // [meta, input_file, index_file]
        fasta          // reference fasta file
        fasta_fai      // reference fasta index file
        
    main:
        versions = Channel.empty()
        
        // Validate samtools convert commands before execution
        validated_commands = input_channel
            .map { meta, input_file, index_file ->
                // Determine file type and output format
                def file_type = input_file.getExtension()
                def output_fmt = ""
                def output_ext = ""
                
                if (file_type == "bam") {
                    output_fmt = "CRAM"
                    output_ext = "cram"
                } else if (file_type == "cram") {
                    output_fmt = "BAM"
                    output_ext = "bam"
                } else if (file_type == "sam") {
                    output_fmt = "BAM"
                    output_ext = "bam"
                } else {
                    throw new IllegalArgumentException("Unsupported input file type: ${file_type} for sample ${meta.id}")
                }
                
                // Construct command parts for validation
                def command_parts = [
                    executable: "samtools",
                    args: [
                        "view",
                        "-o", "${meta.id}.${output_ext}",
                        input_file.toString()
                    ],
                    input_files: [input_file],
                    output_files: ["${meta.id}.${output_ext}"],
                    resources: [
                        memory: "4.GB",
                        cpus: 2,
                        time: "2.h"
                    ]
                ]
                
                // Add reference files if provided
                if (fasta && fasta_fai) {
                    command_parts.args.addAll(1, ["--reference", fasta.toString()])
                    command_parts.input_files.addAll([fasta, fasta_fai])
                }
                
                // Add index file if provided
                if (index_file && index_file.toString() != "null") {
                    command_parts.input_files.add(index_file)
                }
                
                return [meta, input_file, index_file, command_parts]
            }
        
        // Perform command validation
        command_validation_results = validated_commands
            .map { meta, input_file, index_file, command_parts ->
                // Create a channel for command validation
                return [
                    process_name: "SAMTOOLS_CONVERT_${meta.id}",
                    command_parts: command_parts,
                    meta: meta,
                    input_file: input_file,
                    index_file: index_file
                ]
            }
        
        // Validate each command
        validated_results = command_validation_results
            .map { validation_data ->
                try {
                    // Perform inline command validation
                    def process_name = validation_data.process_name
                    def cmd_parts = validation_data.command_parts
                    
                    log.info "VALIDATED_SAMTOOLS_CONVERT: Validating command for ${process_name}"
                    
                    // Validate command structure
                    if (!cmd_parts || !(cmd_parts instanceof Map)) {
                        throw new IllegalArgumentException("Invalid command parts for ${process_name}")
                    }
                    
                    // Validate executable
                    if (!cmd_parts.executable || cmd_parts.executable != "samtools") {
                        throw new IllegalArgumentException("Invalid executable for ${process_name}: ${cmd_parts.executable}")
                    }
                    
                    // Validate arguments
                    if (!cmd_parts.args || !(cmd_parts.args instanceof List)) {
                        throw new IllegalArgumentException("Invalid arguments for ${process_name}")
                    }
                    
                    // Validate input files exist
                    cmd_parts.input_files.each { file ->
                        if (!file.exists()) {
                            throw new IllegalArgumentException("Input file does not exist for ${process_name}: ${file}")
                        }
                        if (!file.canRead()) {
                            throw new IllegalArgumentException("Input file is not readable for ${process_name}: ${file}")
                        }
                    }
                    
                    // Validate resource requirements
                    if (cmd_parts.resources) {
                        if (cmd_parts.resources.memory && !cmd_parts.resources.memory.matches(/^\d+(\.\d+)?\.(GB|MB|KB)$/)) {
                            throw new IllegalArgumentException("Invalid memory format for ${process_name}: ${cmd_parts.resources.memory}")
                        }
                        if (cmd_parts.resources.cpus && cmd_parts.resources.cpus < 1) {
                            throw new IllegalArgumentException("Invalid CPU count for ${process_name}: ${cmd_parts.resources.cpus}")
                        }
                    }
                    
                    // Construct full command for logging
                    def full_command = "${cmd_parts.executable} ${cmd_parts.args.join(' ')}"
                    
                    log.info "VALIDATED_SAMTOOLS_CONVERT: Command validation passed for ${process_name}"
                    log.info "  - Full command: ${full_command}"
                    log.info "  - Input files: ${cmd_parts.input_files.size()}"
                    log.info "  - Output files: ${cmd_parts.output_files.size()}"
                    
                    // Return validated data
                    return [
                        validation_data.meta,
                        validation_data.input_file,
                        validation_data.index_file,
                        cmd_parts + [
                            validated: true,
                            full_command: full_command,
                            validation_timestamp: new Date().toString()
                        ]
                    ]
                    
                } catch (Exception e) {
                    log.error "VALIDATED_SAMTOOLS_CONVERT: Command validation failed for ${validation_data.process_name}: ${e.getMessage()}"
                    throw new RuntimeException("Command validation failed for ${validation_data.process_name}: ${e.getMessage()}", e)
                }
            }
        
        // Extract validated inputs for samtools convert
        validated_inputs = validated_results.map { meta, input_file, index_file, validated_cmd ->
            log.info "VALIDATED_SAMTOOLS_CONVERT: Proceeding with validated conversion for ${meta.id}"
            return [meta, input_file, index_file]
        }
        
        // Create reference file channels
        fasta_channel = Channel.of([[id: "fasta"], fasta])
        fasta_fai_channel = Channel.of([[id: "fasta_fai"], fasta_fai])
        
        // Execute samtools convert with validated commands
        SAMTOOLS_CONVERT_ENHANCED(
            validated_inputs,
            fasta_channel,
            fasta_fai_channel
        )
        
        versions = versions.mix(SAMTOOLS_CONVERT_ENHANCED.out.versions)
        
        // Log successful completion
        completed_conversions = SAMTOOLS_CONVERT_ENHANCED.out.bam
            .mix(SAMTOOLS_CONVERT_ENHANCED.out.cram)
            .map { meta, output_file ->
                log.info "VALIDATED_SAMTOOLS_CONVERT: Successfully completed conversion for ${meta.id}: ${output_file.getSimpleName()}"
                return [meta, output_file]
            }
    
    emit:
        bam = SAMTOOLS_CONVERT_ENHANCED.out.bam
        cram = SAMTOOLS_CONVERT_ENHANCED.out.cram
        sam = SAMTOOLS_CONVERT_ENHANCED.out.sam
        bai = SAMTOOLS_CONVERT_ENHANCED.out.bai
        csi = SAMTOOLS_CONVERT_ENHANCED.out.csi
        crai = SAMTOOLS_CONVERT_ENHANCED.out.crai
        versions = versions
}