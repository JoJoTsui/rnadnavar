//
// VALIDATED_SAMTOOLS_CONVERT: Samtools convert with command validation
//

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
                    input_files: [],
                    output_files: ["${meta.id}.${output_ext}"],
                    resources: [
                        memory: "4.GB",
                        cpus: 2,
                        time: "2.h"
                    ]
                ]
                
                // Add input file (only if it's a proper file object)
                if (input_file) {
                    command_parts.input_files.add(input_file)
                }
                
                // Add reference files if provided
                if (fasta && fasta_fai) {
                    command_parts.args.addAll(1, ["--reference", fasta.toString()])
                    // Only add to input_files if they're proper file objects
                    if (fasta.hasProperty('exists')) {
                        command_parts.input_files.add(fasta)
                    }
                    if (fasta_fai.hasProperty('exists')) {
                        command_parts.input_files.add(fasta_fai)
                    }
                }
                
                // Add index file if provided
                if (index_file && index_file.toString() != "null" && index_file.hasProperty('exists')) {
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
                    
                    if (params.debug_verbose) { log.info "VALIDATED_SAMTOOLS_CONVERT: Validating command for ${process_name}" }
                    
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
                    
                    // Validate input files exist (only for proper file objects)
                    if (cmd_parts.input_files && cmd_parts.input_files.size() > 0) {
                        cmd_parts.input_files.each { file_obj ->
                            if (file_obj == null) {
                                return // Skip null files
                            }
                            
                            try {
                                // Only validate if it's a file object with exists() method
                                if (file_obj.hasProperty('exists') && file_obj.respondsTo('exists')) {
                                    if (!file_obj.exists()) {
                                        throw new IllegalArgumentException("Input file does not exist for ${process_name}: ${file_obj}")
                                    }
                                    if (file_obj.respondsTo('canRead') && !file_obj.canRead()) {
                                        throw new IllegalArgumentException("Input file is not readable for ${process_name}: ${file_obj}")
                                    }
                                    log.debug "VALIDATED_SAMTOOLS_CONVERT: File validation passed for ${file_obj}"
                                } else {
                                    log.debug "VALIDATED_SAMTOOLS_CONVERT: Skipping validation for non-file object: ${file_obj.getClass().simpleName}"
                                }
                            } catch (Exception e) {
                                log.warn "VALIDATED_SAMTOOLS_CONVERT: File validation error for ${file_obj}: ${e.getMessage()}"
                            }
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
                    
                    if (params.debug_verbose) {
                        log.info "VALIDATED_SAMTOOLS_CONVERT: Command validation passed for ${process_name}"
                        log.info "  - Full command: ${full_command}"
                        log.info "  - Input files: ${cmd_parts.input_files.size()}"
                        log.info "  - Output files: ${cmd_parts.output_files.size()}"
                    }
                    
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
            if (params.debug_verbose) { log.info "VALIDATED_SAMTOOLS_CONVERT: Proceeding with validated conversion for ${meta.id}" }
            return [meta, input_file, index_file]
        }
        
        // Execute samtools convert with validated commands
        // Ensure both fasta and fasta_fai are proper tuples as expected by the module
        SAMTOOLS_CONVERT_ENHANCED(
            validated_inputs,
            fasta,  // Already a tuple [meta, file] from PREPARE_GENOME
            fasta_fai.map{ fai -> [[id: "fasta_fai"], fai] }  // Convert file to tuple
        )
        
        versions = versions.mix(SAMTOOLS_CONVERT_ENHANCED.out.versions)
        
        // Log successful completion (only when debug_verbose enabled)
        completed_conversions = SAMTOOLS_CONVERT_ENHANCED.out.bam
            .mix(SAMTOOLS_CONVERT_ENHANCED.out.cram)
            .map { meta, output_file ->
                if (params.debug_verbose) { log.info "VALIDATED_SAMTOOLS_CONVERT: Successfully completed conversion for ${meta.id}: ${output_file.getSimpleName()}" }
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