//
// COMMAND_VALIDATION: Validate command syntax and parameters before execution
//

workflow COMMAND_VALIDATION {
    take:
        process_name    // String: Name of the process for validation context
        command_parts   // Map: Command components to validate
        
    main:
        versions = Channel.empty()
        
        // Validate command structure and syntax
        validated_command = Channel.of(command_parts)
            .map { cmd_parts ->
                try {
                    if (params.debug_verbose) { log.info "COMMAND_VALIDATION: Validating command for ${process_name}" }
                    
                    // Validate command structure
                    if (!cmd_parts) {
                        throw new IllegalArgumentException("Command parts are null or empty for ${process_name}")
                    }
                    
                    if (!(cmd_parts instanceof Map)) {
                        throw new IllegalArgumentException("Command parts must be a Map for ${process_name}")
                    }
                    
                    // Validate required command components
                    def requiredFields = ['executable', 'args']
                    def missingFields = requiredFields.findAll { field ->
                        !cmd_parts.containsKey(field) || cmd_parts[field] == null
                    }
                    
                    if (missingFields.size() > 0) {
                        throw new IllegalArgumentException("Missing required command fields for ${process_name}: ${missingFields}")
                    }
                    
                    // Validate executable
                    def executable = cmd_parts.executable
                    if (!(executable instanceof String) || executable.trim().isEmpty()) {
                        throw new IllegalArgumentException("Executable must be a non-empty string for ${process_name}")
                    }
                    
                    // Check for potentially dangerous commands
                    def dangerousCommands = ['rm', 'rmdir', 'del', 'format', 'fdisk', 'mkfs']
                    if (dangerousCommands.any { dangerous -> executable.toLowerCase().contains(dangerous) }) {
                        log.warn "POTENTIALLY DANGEROUS COMMAND detected in ${process_name}: ${executable}"
                    }
                    
                    // Validate arguments
                    def args = cmd_parts.args
                    if (args != null) {
                        if (args instanceof List) {
                            // Validate each argument in the list
                            args.eachWithIndex { arg, idx ->
                                if (arg == null) {
                                    throw new IllegalArgumentException("Argument ${idx} is null in ${process_name}")
                                }
                                if (!(arg instanceof String) && !(arg instanceof Number) && !(arg instanceof Boolean)) {
                                    log.warn "Argument ${idx} in ${process_name} is not a basic type: ${arg.getClass().getSimpleName()}"
                                }
                            }
                        } else if (args instanceof String) {
                            // Single string argument - validate it's not empty
                            if (args.trim().isEmpty()) {
                                log.warn "Empty string argument detected in ${process_name}"
                            }
                        } else {
                            throw new IllegalArgumentException("Arguments must be a List or String for ${process_name}")
                        }
                    }
                    
                    // Validate input files if specified
                    if (cmd_parts.input_files) {
                        validateInputFiles(cmd_parts.input_files, process_name)
                    }
                    
                    // Validate output files if specified
                    if (cmd_parts.output_files) {
                        validateOutputFiles(cmd_parts.output_files, process_name)
                    }
                    
                    // Validate resource requirements if specified
                    if (cmd_parts.resources) {
                        validateResourceRequirements(cmd_parts.resources, process_name)
                    }
                    
                    // Construct and validate the full command
                    def full_command = constructFullCommand(executable, args)
                    validateCommandSyntax(full_command, process_name)
                    
                    // Log successful validation
                    if (params.debug_verbose) {
                        log.info "COMMAND_VALIDATION: Successfully validated command for ${process_name}"
                        log.info "  - Executable: ${executable}"
                        log.info "  - Arguments: ${args}"
                        log.info "  - Full command: ${full_command}"
                    }
                    
                    // Return validated command parts with additional metadata
                    return cmd_parts + [
                        validated: true,
                        full_command: full_command,
                        validation_timestamp: new Date().toString()
                    ]
                    
                } catch (Exception e) {
                    log.error "COMMAND_VALIDATION: Validation failed for ${process_name}: ${e.getMessage()}"
                    throw new RuntimeException("Command validation failed for ${process_name}: ${e.getMessage()}", e)
                }
            }
    
    emit:
        validated_command = validated_command
        versions = versions
}

// Helper function to validate input files
def validateInputFiles(inputFiles, processName) {
    if (inputFiles instanceof List) {
        inputFiles.eachWithIndex { file, idx ->
            validateSingleFile(file, "input file ${idx}", processName, true)
        }
    } else if (inputFiles instanceof String || inputFiles instanceof java.io.File) {
        validateSingleFile(inputFiles, "input file", processName, true)
    } else {
        throw new IllegalArgumentException("Input files must be a List, String, or File for ${processName}")
    }
}

// Helper function to validate output files
def validateOutputFiles(outputFiles, processName) {
    if (outputFiles instanceof List) {
        outputFiles.eachWithIndex { file, idx ->
            validateSingleFile(file, "output file ${idx}", processName, false)
        }
    } else if (outputFiles instanceof String || outputFiles instanceof java.io.File) {
        validateSingleFile(outputFiles, "output file", processName, false)
    } else {
        throw new IllegalArgumentException("Output files must be a List, String, or File for ${processName}")
    }
}

// Helper function to validate a single file
def validateSingleFile(file, fileType, processName, mustExist) {
    if (file == null) {
        throw new IllegalArgumentException("${fileType} is null for ${processName}")
    }
    
    def filePath
    if (file instanceof String) {
        filePath = file
    } else if (file instanceof java.io.File) {
        filePath = file.getAbsolutePath()
    } else {
        throw new IllegalArgumentException("${fileType} must be String or File for ${processName}")
    }
    
    // Validate file path format
    if (filePath.trim().isEmpty()) {
        throw new IllegalArgumentException("${fileType} path is empty for ${processName}")
    }
    
    // Check for potentially dangerous paths
    def dangerousPaths = ['/etc/', '/bin/', '/sbin/', '/usr/bin/', '/usr/sbin/', 'C:\\Windows\\', 'C:\\Program Files\\']
    if (dangerousPaths.any { dangerous -> filePath.startsWith(dangerous) }) {
        log.warn "POTENTIALLY DANGEROUS PATH detected for ${fileType} in ${processName}: ${filePath}"
    }
    
    // Check file existence if required
    if (mustExist) {
        def fileObj = new java.io.File(filePath)
        if (!fileObj.exists()) {
            throw new IllegalArgumentException("${fileType} does not exist for ${processName}: ${filePath}")
        }
        if (!fileObj.canRead()) {
            throw new IllegalArgumentException("${fileType} is not readable for ${processName}: ${filePath}")
        }
        if (fileObj.length() == 0) {
            log.warn "${fileType} is empty for ${processName}: ${filePath}"
        }
    } else {
        // For output files, check if parent directory exists and is writable
        def fileObj = new java.io.File(filePath)
        def parentDir = fileObj.getParentFile()
        if (parentDir != null && !parentDir.exists()) {
            log.warn "Parent directory does not exist for ${fileType} in ${processName}: ${parentDir.getAbsolutePath()}"
        } else if (parentDir != null && !parentDir.canWrite()) {
            throw new IllegalArgumentException("Parent directory is not writable for ${fileType} in ${processName}: ${parentDir.getAbsolutePath()}")
        }
    }
    
    log.debug "Validated ${fileType} for ${processName}: ${filePath}"
}

// Helper function to validate resource requirements
def validateResourceRequirements(resources, processName) {
    if (!(resources instanceof Map)) {
        throw new IllegalArgumentException("Resources must be a Map for ${processName}")
    }
    
    // Validate memory requirements
    if (resources.memory) {
        def memory = resources.memory
        if (!(memory instanceof String)) {
            throw new IllegalArgumentException("Memory requirement must be a String for ${processName}")
        }
        
        // Check memory format (e.g., "4.GB", "1024.MB")
        if (!memory.matches(/^\d+(\.\d+)?\.(GB|MB|KB)$/)) {
            throw new IllegalArgumentException("Invalid memory format for ${processName}: ${memory}. Expected format: number.UNIT (e.g., 4.GB)")
        }
        
        // Extract numeric value and check if reasonable
        def memoryValue = memory.replaceAll(/\.(GB|MB|KB)$/, '').toDouble()
        def unit = memory.replaceAll(/^\d+(\.\d+)?\./, '')
        
        def memoryInMB = 0
        switch (unit) {
            case 'GB':
                memoryInMB = memoryValue * 1024
                break
            case 'MB':
                memoryInMB = memoryValue
                break
            case 'KB':
                memoryInMB = memoryValue / 1024
                break
        }
        
        if (memoryInMB < 100) {
            log.warn "Very low memory requirement for ${processName}: ${memory}"
        } else if (memoryInMB > 64000) {
            log.warn "Very high memory requirement for ${processName}: ${memory}"
        }
    }
    
    // Validate CPU requirements
    if (resources.cpus) {
        def cpus = resources.cpus
        if (!(cpus instanceof Integer) && !(cpus instanceof String)) {
            throw new IllegalArgumentException("CPU requirement must be Integer or String for ${processName}")
        }
        
        def cpuCount = 0
        if (cpus instanceof String) {
            try {
                cpuCount = cpus.toInteger()
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Invalid CPU format for ${processName}: ${cpus}")
            }
        } else {
            cpuCount = cpus
        }
        
        if (cpuCount < 1) {
            throw new IllegalArgumentException("CPU count must be at least 1 for ${processName}: ${cpuCount}")
        } else if (cpuCount > 64) {
            log.warn "Very high CPU requirement for ${processName}: ${cpuCount}"
        }
    }
    
    // Validate time requirements
    if (resources.time) {
        def time = resources.time
        if (!(time instanceof String)) {
            throw new IllegalArgumentException("Time requirement must be a String for ${processName}")
        }
        
        // Check time format (e.g., "1.h", "30.m", "3600.s")
        if (!time.matches(/^\d+(\.\d+)?\.(h|m|s)$/)) {
            throw new IllegalArgumentException("Invalid time format for ${processName}: ${time}. Expected format: number.UNIT (e.g., 1.h)")
        }
    }
    
    log.debug "Validated resource requirements for ${processName}: ${resources}"
}

// Helper function to construct full command
def constructFullCommand(executable, args) {
    def command = executable
    
    if (args != null) {
        if (args instanceof List) {
            // Join list arguments with spaces
            def argString = args.collect { arg ->
                // Quote arguments that contain spaces
                if (arg.toString().contains(' ')) {
                    return "\"${arg}\""
                } else {
                    return arg.toString()
                }
            }.join(' ')
            command = "${executable} ${argString}"
        } else if (args instanceof String) {
            command = "${executable} ${args}"
        }
    }
    
    return command
}

// Helper function to validate command syntax
def validateCommandSyntax(fullCommand, processName) {
    if (fullCommand == null || fullCommand.trim().isEmpty()) {
        throw new IllegalArgumentException("Full command is empty for ${processName}")
    }
    
    // Check for potentially dangerous command patterns
    def dangerousPatterns = [
        /rm\s+-rf\s+\//, // rm -rf /
        />\s*\/dev\/null\s*2>&1/, // Redirecting to /dev/null (might hide errors)
        /;\s*rm\s+/, // Command chaining with rm
        /\|\s*sh/, // Piping to shell
        /\$\(.*\)/, // Command substitution
        /`.*`/ // Backtick command substitution
    ]
    
    dangerousPatterns.each { pattern ->
        if (fullCommand =~ pattern) {
            log.warn "POTENTIALLY DANGEROUS COMMAND PATTERN detected in ${processName}: ${pattern}"
        }
    }
    
    // Check for basic syntax issues
    def openQuotes = fullCommand.count('"') % 2
    if (openQuotes != 0) {
        throw new IllegalArgumentException("Unmatched quotes in command for ${processName}: ${fullCommand}")
    }
    
    def openSingleQuotes = fullCommand.count("'") % 2
    if (openSingleQuotes != 0) {
        throw new IllegalArgumentException("Unmatched single quotes in command for ${processName}: ${fullCommand}")
    }
    
    // Check for excessive length
    if (fullCommand.length() > 8192) {
        log.warn "Very long command detected for ${processName}: ${fullCommand.length()} characters"
    }
    
    log.debug "Command syntax validation passed for ${processName}: ${fullCommand}"
}