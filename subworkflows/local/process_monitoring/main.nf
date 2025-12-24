//
// PROCESS_MONITORING: Comprehensive logging and monitoring for channel operations
//

workflow PROCESS_MONITORING {
    take:
        process_name    // String: Name of the process being monitored
        input_channel   // Channel: Any channel to be monitored
        
    main:
        versions = Channel.empty()
        
        // Monitor and log channel operations with comprehensive state inspection
        monitored = input_channel
            .map { items ->
                // Log process start with detailed channel state
                log.info "=== PROCESS_MONITORING: ${process_name} ==="
                log.info "Channel items count: ${items.size()}"
                
                // Inspect and log metadata structure if present
                if (items.size() > 0 && items[0] instanceof Map) {
                    def meta = items[0]
                    log.info "Meta structure for ${process_name}:"
                    log.info "  - ID: ${meta.id ?: 'MISSING'}"
                    log.info "  - Patient: ${meta.patient ?: 'MISSING'}"
                    log.info "  - Sample: ${meta.sample ?: 'MISSING'}"
                    log.info "  - Status: ${meta.status != null ? meta.status : 'MISSING'}"
                    log.info "  - Data type: ${meta.data_type ?: 'UNKNOWN'}"
                    
                    // Check for potential circular references or problematic fields
                    def problematicFields = []
                    meta.each { key, value ->
                        if (value instanceof java.io.File || value instanceof java.nio.file.Path) {
                            problematicFields.add("${key}: File/Path reference detected")
                        } else if (value instanceof Map && value.containsKey(key)) {
                            problematicFields.add("${key}: Potential circular reference")
                        } else if (value == null) {
                            problematicFields.add("${key}: null value")
                        }
                    }
                    
                    if (problematicFields.size() > 0) {
                        log.warn "Potential issues detected in ${process_name} metadata:"
                        problematicFields.each { issue ->
                            log.warn "  - ${issue}"
                        }
                    }
                }
                
                // Log file information if present
                for (int i = 1; i < items.size(); i++) {
                    def item = items[i]
                    if (item instanceof java.io.File || item instanceof java.nio.file.Path) {
                        // Handle both File and Path objects
                        def fileName = item instanceof java.nio.file.Path ? item.getFileName().toString() : item.name
                        def fileExists = item instanceof java.nio.file.Path ? java.nio.file.Files.exists(item) : item.exists()
                        def fileSize = fileExists ? (item instanceof java.nio.file.Path ? java.nio.file.Files.size(item) : item.length()) : 0
                        log.info "File ${i}: ${fileName} (${fileExists ? 'EXISTS' : 'MISSING'}, ${fileExists ? "${fileSize} bytes" : 'N/A'})"
                    } else if (item instanceof List) {
                        log.info "File list ${i}: ${item.size()} files"
                        item.eachWithIndex { subItem, idx ->
                            if (subItem instanceof java.io.File || subItem instanceof java.nio.file.Path) {
                                def fileName = subItem instanceof java.nio.file.Path ? subItem.getFileName().toString() : subItem.name
                                def fileExists = subItem instanceof java.nio.file.Path ? java.nio.file.Files.exists(subItem) : subItem.exists()
                                log.info "  - File ${idx}: ${fileName} (${fileExists ? 'EXISTS' : 'MISSING'})"
                            }
                        }
                    } else {
                        log.info "Item ${i}: ${item?.getClass()?.simpleName ?: 'null'} - ${item}"
                    }
                }
                
                // Monitor memory usage
                def runtime = Runtime.getRuntime()
                def totalMemory = runtime.totalMemory()
                def freeMemory = runtime.freeMemory()
                def usedMemory = totalMemory - freeMemory
                def maxMemory = runtime.maxMemory()
                
                log.info "Memory usage for ${process_name}:"
                log.info "  - Used: ${String.format('%.2f', usedMemory / 1024 / 1024)} MB"
                log.info "  - Free: ${String.format('%.2f', freeMemory / 1024 / 1024)} MB"
                log.info "  - Total: ${String.format('%.2f', totalMemory / 1024 / 1024)} MB"
                log.info "  - Max: ${String.format('%.2f', maxMemory / 1024 / 1024)} MB"
                log.info "  - Usage: ${String.format('%.1f', (usedMemory * 100.0) / maxMemory)}%"
                
                // Check for memory pressure
                def memoryUsagePercent = (usedMemory * 100.0) / maxMemory
                if (memoryUsagePercent > 80) {
                    log.warn "HIGH MEMORY USAGE WARNING in ${process_name}: ${String.format('%.1f', memoryUsagePercent)}%"
                } else if (memoryUsagePercent > 60) {
                    log.info "Moderate memory usage in ${process_name}: ${String.format('%.1f', memoryUsagePercent)}%"
                }
                
                log.info "=== END MONITORING: ${process_name} ==="
                return items
            }
            .filter { items ->
                // Validate that no items are null
                def nullItems = []
                items.eachWithIndex { item, idx ->
                    if (item == null) {
                        nullItems.add(idx)
                    }
                }
                
                if (nullItems.size() > 0) {
                    log.error "NULL ITEMS DETECTED in ${process_name} at positions: ${nullItems}"
                    log.error "This may indicate upstream processing issues"
                    return false
                }
                
                // Validate metadata structure if present
                if (items.size() > 0 && items[0] instanceof Map) {
                    def meta = items[0]
                    def requiredFields = ['id', 'patient']
                    def missingFields = requiredFields.findAll { field ->
                        !meta.containsKey(field) || meta[field] == null
                    }
                    
                    if (missingFields.size() > 0) {
                        log.error "MISSING REQUIRED FIELDS in ${process_name}: ${missingFields}"
                        return false
                    }
                }
                
                log.info "Validation passed for ${process_name}"
                return true
            }
            .map { items ->
                // Final logging before passing items downstream
                if (items.size() > 0 && items[0] instanceof Map) {
                    log.info "Successfully processed ${process_name} for sample: ${items[0].id}"
                } else {
                    log.info "Successfully processed ${process_name}"
                }
                return items
            }
    
    emit:
        output = monitored
        versions = versions
}