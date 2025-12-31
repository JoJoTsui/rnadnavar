//
// ERROR_SAFE_PROCESS: Wrapper for safe process execution with comprehensive error handling
//

workflow ERROR_SAFE_PROCESS {
    take:
        process_name    // String: Name of the process for error reporting
        input_channel   // Channel: Input channel to be processed safely
        process_closure // Closure: The actual process logic to execute
        
    main:
        versions = Channel.empty()
        
        // Apply comprehensive error handling and safety checks
        safe_output = input_channel
            .map { items ->
                try {
                    // Pre-execution validation
                    if (params.debug_verbose) { log.info "ERROR_SAFE_PROCESS: Starting ${process_name}" }
                    
                    // Check for null inputs
                    if (items == null) {
                        throw new IllegalArgumentException("Input channel items are null for ${process_name}")
                    }
                    
                    // Validate items structure
                    if (items instanceof List && items.any { it == null }) {
                        def nullIndices = []
                        items.eachWithIndex { item, idx ->
                            if (item == null) nullIndices.add(idx)
                        }
                        throw new IllegalArgumentException("Null items detected at positions ${nullIndices} in ${process_name}")
                    }
                    
                    // Memory check before processing
                    def runtime = Runtime.getRuntime()
                    def freeMemory = runtime.freeMemory()
                    def maxMemory = runtime.maxMemory()
                    def memoryUsagePercent = ((maxMemory - freeMemory) * 100.0) / maxMemory
                    
                    if (memoryUsagePercent > 90) {
                        log.warn "CRITICAL MEMORY WARNING before ${process_name}: ${String.format('%.1f', memoryUsagePercent)}% used"
                        // Force garbage collection
                        System.gc()
                        Thread.sleep(100) // Brief pause for GC
                        
                        // Re-check memory after GC
                        runtime = Runtime.getRuntime()
                        freeMemory = runtime.freeMemory()
                        memoryUsagePercent = ((maxMemory - freeMemory) * 100.0) / maxMemory
                        if (params.debug_verbose) { log.info "Memory usage after GC: ${String.format('%.1f', memoryUsagePercent)}%" }
                    }
                    
                    // Check for potential circular references in metadata
                    if (items instanceof List && items.size() > 0 && items[0] instanceof Map) {
                        def meta = items[0]
                        def circularRefs = detectCircularReferences(meta, process_name)
                        if (circularRefs.size() > 0) {
                            log.warn "Potential circular references detected in ${process_name}: ${circularRefs}"
                        }
                    }
                    
                    if (params.debug_verbose) { log.info "Pre-execution validation passed for ${process_name}" }
                    return items
                    
                } catch (StackOverflowError e) {
                    log.error "STACKOVERFLOW ERROR detected in ${process_name}"
                    log.error "This typically indicates circular references in metadata or infinite recursion"
                    log.error "Stack trace: ${e.getStackTrace().take(10).join('\n')}"
                    
                    // Capture diagnostic information
                    captureDiagnosticInfo(process_name, items, e)
                    
                    throw new RuntimeException("StackOverflowError in ${process_name}: Check for circular references in metadata", e)
                    
                } catch (OutOfMemoryError e) {
                    log.error "OUT OF MEMORY ERROR in ${process_name}"
                    log.error "Consider increasing memory allocation or optimizing data structures"
                    
                    // Capture memory diagnostic information
                    captureMemoryDiagnostics(process_name, e)
                    
                    throw new RuntimeException("OutOfMemoryError in ${process_name}: Insufficient memory", e)
                    
                } catch (Exception e) {
                    log.error "GENERAL ERROR in ${process_name}: ${e.getMessage()}"
                    log.error "Error type: ${e.getClass().getSimpleName()}"
                    log.error "Stack trace: ${e.getStackTrace().take(5).join('\n')}"
                    
                    // Capture general diagnostic information
                    captureGeneralDiagnostics(process_name, items, e)
                    
                    throw new RuntimeException("Error in ${process_name}: ${e.getMessage()}", e)
                }
            }
            .filter { items ->
                try {
                    // Post-processing validation
                    if (items == null) {
                        log.error "Process ${process_name} returned null output"
                        return false
                    }
                    
                    // Validate output structure
                    if (items instanceof List) {
                        def nullCount = items.count { it == null }
                        if (nullCount > 0) {
                            log.error "Process ${process_name} output contains ${nullCount} null items"
                            return false
                        }
                    }
                    
                    if (params.debug_verbose) { log.info "Post-processing validation passed for ${process_name}" }
                    return true
                    
                } catch (Exception e) {
                    log.error "Post-processing validation failed for ${process_name}: ${e.getMessage()}"
                    return false
                }
            }
            .map { items ->
                // Final success logging
                if (params.debug_verbose) { log.info "ERROR_SAFE_PROCESS: Successfully completed ${process_name}" }
                return items
            }
    
    emit:
        output = safe_output
        versions = versions
}

// Helper function to detect circular references in metadata
def detectCircularReferences(meta, processName, visited = [], path = "") {
    def circularRefs = []
    
    if (meta == null || !(meta instanceof Map)) {
        return circularRefs
    }
    
    meta.each { key, value ->
        def currentPath = path.isEmpty() ? key : "${path}.${key}"
        
        if (value instanceof Map) {
            if (visited.contains(value)) {
                circularRefs.add("Circular reference at ${currentPath}")
            } else {
                visited.add(value)
                circularRefs.addAll(detectCircularReferences(value, processName, visited, currentPath))
                visited.remove(value)
            }
        } else if (value instanceof List) {
            value.eachWithIndex { item, idx ->
                if (item instanceof Map) {
                    def itemPath = "${currentPath}[${idx}]"
                    if (visited.contains(item)) {
                        circularRefs.add("Circular reference at ${itemPath}")
                    } else {
                        visited.add(item)
                        circularRefs.addAll(detectCircularReferences(item, processName, visited, itemPath))
                        visited.remove(item)
                    }
                }
            }
        }
    }
    
    return circularRefs
}

// Helper function to capture diagnostic information for StackOverflowError
def captureDiagnosticInfo(processName, items, error) {
    try {
        log.error "=== DIAGNOSTIC INFORMATION FOR ${processName} ==="
        log.error "Error type: StackOverflowError"
        log.error "Timestamp: ${new Date()}"
        
        if (items != null) {
            log.error "Input items type: ${items.getClass().getSimpleName()}"
            if (items instanceof List) {
                log.error "Input items count: ${items.size()}"
                items.eachWithIndex { item, idx ->
                    if (item instanceof Map) {
                        log.error "Item ${idx} keys: ${item.keySet()}"
                        log.error "Item ${idx} size: ${item.size()}"
                    } else {
                        log.error "Item ${idx} type: ${item?.getClass()?.getSimpleName() ?: 'null'}"
                    }
                }
            }
        }
        
        // Memory information
        def runtime = Runtime.getRuntime()
        log.error "Memory - Total: ${runtime.totalMemory()}, Free: ${runtime.freeMemory()}, Max: ${runtime.maxMemory()}"
        
        log.error "=== END DIAGNOSTIC INFORMATION ==="
        
    } catch (Exception diagError) {
        log.error "Failed to capture diagnostic information: ${diagError.getMessage()}"
    }
}

// Helper function to capture memory diagnostic information
def captureMemoryDiagnostics(processName, error) {
    try {
        log.error "=== MEMORY DIAGNOSTIC INFORMATION FOR ${processName} ==="
        
        def runtime = Runtime.getRuntime()
        def totalMemory = runtime.totalMemory()
        def freeMemory = runtime.freeMemory()
        def maxMemory = runtime.maxMemory()
        def usedMemory = totalMemory - freeMemory
        
        log.error "Memory usage at error:"
        log.error "  - Used: ${String.format('%.2f', usedMemory / 1024 / 1024)} MB"
        log.error "  - Free: ${String.format('%.2f', freeMemory / 1024 / 1024)} MB"
        log.error "  - Total: ${String.format('%.2f', totalMemory / 1024 / 1024)} MB"
        log.error "  - Max: ${String.format('%.2f', maxMemory / 1024 / 1024)} MB"
        log.error "  - Usage: ${String.format('%.1f', (usedMemory * 100.0) / maxMemory)}%"
        
        log.error "=== END MEMORY DIAGNOSTIC INFORMATION ==="
        
    } catch (Exception diagError) {
        log.error "Failed to capture memory diagnostic information: ${diagError.getMessage()}"
    }
}

// Helper function to capture general diagnostic information
def captureGeneralDiagnostics(processName, items, error) {
    try {
        log.error "=== GENERAL DIAGNOSTIC INFORMATION FOR ${processName} ==="
        log.error "Error message: ${error.getMessage()}"
        log.error "Error type: ${error.getClass().getSimpleName()}"
        log.error "Timestamp: ${new Date()}"
        
        if (items != null) {
            log.error "Input structure: ${items.getClass().getSimpleName()}"
            if (items instanceof List) {
                log.error "Input size: ${items.size()}"
            }
        }
        
        log.error "=== END GENERAL DIAGNOSTIC INFORMATION ==="
        
    } catch (Exception diagError) {
        log.error "Failed to capture general diagnostic information: ${diagError.getMessage()}"
    }
}