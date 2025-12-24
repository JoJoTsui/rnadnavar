/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Channel Inspector - Debug tool for analyzing Nextflow channels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNEL INSPECTION UTILITIES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Inspect channel contents and structure
workflow INSPECT_CHANNEL {
    take:
        input_channel
        channel_name
        
    main:
        inspected = input_channel
            .map { items ->
                def inspection = [
                    channel_name: channel_name,
                    timestamp: new Date().toString(),
                    item_count: items instanceof List ? items.size() : 1,
                    item_types: items instanceof List ? 
                        items.collect { it?.getClass()?.getSimpleName() ?: 'null' } : 
                        [items?.getClass()?.getSimpleName() ?: 'null'],
                    meta_structure: items instanceof List && items.size() > 0 && items[0] instanceof Map ? 
                        analyzeMetaStructure(items[0]) : null,
                    file_info: items instanceof List ? 
                        items.findAll { it instanceof File }.collect { analyzeFileInfo(it) } : 
                        items instanceof File ? [analyzeFileInfo(items)] : []
                ]
                
                log.info "Channel Inspection: ${channel_name}"
                log.info "  Items: ${inspection.item_count}"
                log.info "  Types: ${inspection.item_types.join(', ')}"
                if (inspection.meta_structure) {
                    log.info "  Meta fields: ${inspection.meta_structure.fields.join(', ')}"
                    log.info "  Meta depth: ${inspection.meta_structure.max_depth}"
                }
                if (inspection.file_info.size() > 0) {
                    inspection.file_info.each { fileInfo ->
                        log.info "  File: ${fileInfo.name} (${fileInfo.size} bytes, exists: ${fileInfo.exists})"
                    }
                }
                
                return [inspection, items]
            }
    
    emit:
        inspected = inspected.map { inspection, items -> items }
        report = inspected.map { inspection, items -> inspection }
}

// Analyze metadata structure for circular references
def analyzeMetaStructure(meta, visited = [], depth = 0) {
    if (depth > 10) {
        return [
            fields: ['CIRCULAR_REFERENCE_DETECTED'],
            max_depth: depth,
            circular: true
        ]
    }
    
    if (meta == null || !(meta instanceof Map)) {
        return [
            fields: [],
            max_depth: depth,
            circular: false
        ]
    }
    
    def objectId = System.identityHashCode(meta)
    if (visited.contains(objectId)) {
        return [
            fields: ['CIRCULAR_REFERENCE_DETECTED'],
            max_depth: depth,
            circular: true
        ]
    }
    
    visited = visited + [objectId]
    
    def fields = []
    def maxDepth = depth
    def hasCircular = false
    
    meta.each { key, value ->
        fields << key
        
        if (value instanceof Map) {
            def subAnalysis = analyzeMetaStructure(value, visited, depth + 1)
            maxDepth = Math.max(maxDepth, subAnalysis.max_depth)
            hasCircular = hasCircular || subAnalysis.circular
        } else if (value instanceof List) {
            value.each { item ->
                if (item instanceof Map) {
                    def subAnalysis = analyzeMetaStructure(item, visited, depth + 1)
                    maxDepth = Math.max(maxDepth, subAnalysis.max_depth)
                    hasCircular = hasCircular || subAnalysis.circular
                }
            }
        }
    }
    
    return [
        fields: fields,
        max_depth: maxDepth,
        circular: hasCircular
    ]
}

// Analyze file information
def analyzeFileInfo(file) {
    if (file == null || !(file instanceof File)) {
        return [
            name: 'null',
            exists: false,
            size: 0,
            readable: false,
            writable: false
        ]
    }
    
    return [
        name: file.getName(),
        path: file.getAbsolutePath(),
        exists: file.exists(),
        size: file.exists() ? file.length() : 0,
        readable: file.exists() ? file.canRead() : false,
        writable: file.exists() ? file.canWrite() : false,
        directory: file.exists() ? file.isDirectory() : false
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNEL STATE MONITORING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MONITOR_CHANNEL_STATE {
    take:
        input_channel
        monitor_name
        
    main:
        // Count items passing through
        item_counter = Channel.value(0)
        
        monitored = input_channel
            .map { items ->
                def timestamp = System.currentTimeMillis()
                def itemId = System.identityHashCode(items)
                
                log.info "Channel Monitor [${monitor_name}]: Item ${itemId} at ${new Date(timestamp)}"
                
                // Check for potential issues
                if (items == null) {
                    log.warn "Channel Monitor [${monitor_name}]: NULL item detected"
                } else if (items instanceof List && items.any { it == null }) {
                    log.warn "Channel Monitor [${monitor_name}]: List contains NULL elements"
                } else if (items instanceof List && items.size() == 0) {
                    log.warn "Channel Monitor [${monitor_name}]: Empty list detected"
                }
                
                return [timestamp, itemId, items]
            }
            .map { timestamp, itemId, items -> items }
    
    emit:
        monitored = monitored
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNEL JOIN DEBUGGING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEBUG_CHANNEL_JOIN {
    take:
        left_channel
        right_channel
        join_key_extractor
        
    main:
        // Extract and log join keys
        left_keyed = left_channel
            .map { items ->
                def key = join_key_extractor(items)
                log.info "Left channel key: ${key} for items: ${items}"
                return [key, items]
            }
        
        right_keyed = right_channel
            .map { items ->
                def key = join_key_extractor(items)
                log.info "Right channel key: ${key} for items: ${items}"
                return [key, items]
            }
        
        // Collect all keys for analysis
        left_keys = left_keyed.map { key, items -> key }.collect()
        right_keys = right_keyed.map { key, items -> key }.collect()
        
        // Analyze join compatibility
        join_analysis = left_keys.combine(right_keys)
            .map { leftKeys, rightKeys ->
                def commonKeys = leftKeys.intersect(rightKeys)
                def leftOnlyKeys = leftKeys - rightKeys
                def rightOnlyKeys = rightKeys - leftKeys
                
                log.info "Join Analysis:"
                log.info "  Left keys: ${leftKeys}"
                log.info "  Right keys: ${rightKeys}"
                log.info "  Common keys: ${commonKeys}"
                log.info "  Left-only keys: ${leftOnlyKeys}"
                log.info "  Right-only keys: ${rightOnlyKeys}"
                log.info "  Expected join results: ${commonKeys.size()}"
                
                return [
                    left_keys: leftKeys,
                    right_keys: rightKeys,
                    common_keys: commonKeys,
                    left_only: leftOnlyKeys,
                    right_only: rightOnlyKeys
                ]
            }
        
        // Perform the actual join
        joined = left_keyed.join(right_keyed, by: 0, failOnMismatch: false)
            .map { key, leftItems, rightItems ->
                log.info "Successful join for key: ${key}"
                return [leftItems, rightItems]
            }
    
    emit:
        joined = joined
        analysis = join_analysis
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MEMORY USAGE MONITORING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MONITOR_MEMORY_USAGE {
    take:
        input_channel
        process_name
        
    main:
        monitored = input_channel
            .map { items ->
                def runtime = Runtime.getRuntime()
                def totalMemory = runtime.totalMemory()
                def freeMemory = runtime.freeMemory()
                def usedMemory = totalMemory - freeMemory
                def maxMemory = runtime.maxMemory()
                
                def memoryInfo = [
                    process: process_name,
                    timestamp: new Date().toString(),
                    used_mb: Math.round(usedMemory / 1024 / 1024),
                    free_mb: Math.round(freeMemory / 1024 / 1024),
                    total_mb: Math.round(totalMemory / 1024 / 1024),
                    max_mb: Math.round(maxMemory / 1024 / 1024),
                    usage_percent: Math.round((usedMemory / maxMemory) * 100)
                ]
                
                log.info "Memory Usage [${process_name}]: ${memoryInfo.used_mb}MB used (${memoryInfo.usage_percent}% of max)"
                
                if (memoryInfo.usage_percent > 80) {
                    log.warn "High memory usage detected: ${memoryInfo.usage_percent}%"
                }
                
                return [memoryInfo, items]
            }
            .map { memoryInfo, items -> items }
    
    emit:
        monitored = monitored
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPREHENSIVE CHANNEL DEBUGGER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEBUG_CHANNEL_COMPREHENSIVE {
    take:
        input_channel
        debug_name
        
    main:
        // Apply all debugging tools
        INSPECT_CHANNEL(input_channel, debug_name)
        MONITOR_CHANNEL_STATE(INSPECT_CHANNEL.out.inspected, debug_name)
        MONITOR_MEMORY_USAGE(MONITOR_CHANNEL_STATE.out.monitored, debug_name)
        
        // Collect comprehensive debug report
        debug_report = INSPECT_CHANNEL.out.report
            .collect()
            .map { reports ->
                def summary = [
                    debug_name: debug_name,
                    total_items: reports.size(),
                    timestamp: new Date().toString(),
                    circular_references: reports.any { it.meta_structure?.circular },
                    file_issues: reports.any { report ->
                        report.file_info.any { !it.exists || !it.readable }
                    },
                    type_consistency: reports.collect { it.item_types }.unique().size() == 1
                ]
                
                log.info "Debug Summary [${debug_name}]:"
                log.info "  Total items processed: ${summary.total_items}"
                log.info "  Circular references detected: ${summary.circular_references}"
                log.info "  File issues detected: ${summary.file_issues}"
                log.info "  Type consistency: ${summary.type_consistency}"
                
                return summary
            }
    
    emit:
        debugged = MONITOR_MEMORY_USAGE.out.monitored
        report = debug_report
}