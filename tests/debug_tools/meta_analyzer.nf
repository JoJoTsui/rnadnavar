/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Meta Data Structure Analyzer - Debug tool for analyzing metadata structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    META DATA ANALYSIS UTILITIES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Deep analysis of metadata structure
workflow ANALYZE_META_STRUCTURE {
    take:
        input_channel
        analysis_name
        
    main:
        analyzed = input_channel
            .map { items ->
                def meta = items instanceof List && items.size() > 0 ? items[0] : items
                
                if (!(meta instanceof Map)) {
                    log.warn "Meta Analysis [${analysis_name}]: No metadata found"
                    return [items, null]
                }
                
                def analysis = performDeepMetaAnalysis(meta, analysis_name)
                
                // Log analysis results
                log.info "Meta Analysis [${analysis_name}]:"
                log.info "  Fields: ${analysis.field_count}"
                log.info "  Max depth: ${analysis.max_depth}"
                log.info "  Circular refs: ${analysis.has_circular_refs}"
                log.info "  File objects: ${analysis.file_object_count}"
                log.info "  Memory footprint: ~${analysis.estimated_size_bytes} bytes"
                
                if (analysis.issues.size() > 0) {
                    log.warn "Meta Analysis Issues:"
                    analysis.issues.each { issue ->
                        log.warn "  - ${issue}"
                    }
                }
                
                return [items, analysis]
            }
    
    emit:
        analyzed = analyzed.map { items, analysis -> items }
        reports = analyzed.map { items, analysis -> analysis }.filter { it != null }
}

// Perform deep analysis of metadata
def performDeepMetaAnalysis(meta, analysisName, visited = [], path = "", depth = 0) {
    def analysis = [
        analysis_name: analysisName,
        field_count: 0,
        max_depth: depth,
        has_circular_refs: false,
        file_object_count: 0,
        estimated_size_bytes: 0,
        issues: [],
        field_types: [:],
        nested_structures: []
    ]
    
    if (depth > 20) {
        analysis.issues << "Maximum analysis depth exceeded (possible infinite recursion)"
        analysis.has_circular_refs = true
        return analysis
    }
    
    if (meta == null) {
        analysis.issues << "Null metadata at path: ${path}"
        return analysis
    }
    
    if (!(meta instanceof Map)) {
        analysis.issues << "Non-map metadata at path: ${path}"
        return analysis
    }
    
    def objectId = System.identityHashCode(meta)
    if (visited.contains(objectId)) {
        analysis.issues << "Circular reference detected at path: ${path}"
        analysis.has_circular_refs = true
        return analysis
    }
    
    visited = visited + [objectId]
    
    meta.each { key, value ->
        analysis.field_count++
        def currentPath = path.isEmpty() ? key : "${path}.${key}"
        
        // Analyze field type
        def valueType = value?.getClass()?.getSimpleName() ?: 'null'
        analysis.field_types[valueType] = (analysis.field_types[valueType] ?: 0) + 1
        
        // Estimate size
        analysis.estimated_size_bytes += estimateObjectSize(value)
        
        // Check for file objects
        if (value instanceof File) {
            analysis.file_object_count++
            analysis.issues << "File object found at path: ${currentPath}"
        }
        
        // Recursive analysis for nested structures
        if (value instanceof Map) {
            def nestedAnalysis = performDeepMetaAnalysis(value, analysisName, visited, currentPath, depth + 1)
            analysis.max_depth = Math.max(analysis.max_depth, nestedAnalysis.max_depth)
            analysis.has_circular_refs = analysis.has_circular_refs || nestedAnalysis.has_circular_refs
            analysis.file_object_count += nestedAnalysis.file_object_count
            analysis.estimated_size_bytes += nestedAnalysis.estimated_size_bytes
            analysis.issues.addAll(nestedAnalysis.issues)
            analysis.nested_structures << currentPath
        } else if (value instanceof List) {
            value.eachWithIndex { item, index ->
                if (item instanceof Map) {
                    def itemPath = "${currentPath}[${index}]"
                    def nestedAnalysis = performDeepMetaAnalysis(item, analysisName, visited, itemPath, depth + 1)
                    analysis.max_depth = Math.max(analysis.max_depth, nestedAnalysis.max_depth)
                    analysis.has_circular_refs = analysis.has_circular_refs || nestedAnalysis.has_circular_refs
                    analysis.file_object_count += nestedAnalysis.file_object_count
                    analysis.estimated_size_bytes += nestedAnalysis.estimated_size_bytes
                    analysis.issues.addAll(nestedAnalysis.issues)
                }
            }
        }
    }
    
    return analysis
}

// Estimate object size in bytes
def estimateObjectSize(obj) {
    if (obj == null) return 8
    if (obj instanceof String) return obj.length() * 2 + 24
    if (obj instanceof Integer) return 24
    if (obj instanceof Long) return 24
    if (obj instanceof Boolean) return 16
    if (obj instanceof File) return 48 + (obj.getAbsolutePath()?.length() ?: 0) * 2
    if (obj instanceof List) return 24 + obj.size() * 8
    if (obj instanceof Map) return 32 + obj.size() * 16
    return 32 // Default estimate
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    META DATA VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VALIDATE_META_REQUIREMENTS {
    take:
        input_channel
        required_fields
        validation_name
        
    main:
        validated = input_channel
            .map { items ->
                def meta = items instanceof List && items.size() > 0 ? items[0] : items
                
                if (!(meta instanceof Map)) {
                    log.error "Meta Validation [${validation_name}]: No metadata map found"
                    return [items, false, ["No metadata map found"]]
                }
                
                def validation = validateMetaRequirements(meta, required_fields, validation_name)
                
                if (!validation.valid) {
                    log.error "Meta Validation [${validation_name}]: FAILED"
                    validation.errors.each { error ->
                        log.error "  - ${error}"
                    }
                } else {
                    log.info "Meta Validation [${validation_name}]: PASSED"
                }
                
                return [items, validation.valid, validation.errors]
            }
    
    emit:
        validated = validated.map { items, valid, errors -> items }
        results = validated.map { items, valid, errors -> [valid: valid, errors: errors] }
}

// Validate metadata against requirements
def validateMetaRequirements(meta, requiredFields, validationName) {
    def validation = [
        valid: true,
        errors: []
    ]
    
    // Check required fields
    requiredFields.each { field ->
        if (!meta.containsKey(field)) {
            validation.valid = false
            validation.errors << "Missing required field: ${field}"
        } else if (meta[field] == null) {
            validation.valid = false
            validation.errors << "Required field is null: ${field}"
        } else if (meta[field] instanceof String && meta[field].trim().isEmpty()) {
            validation.valid = false
            validation.errors << "Required field is empty: ${field}"
        }
    }
    
    // Check for common issues
    if (meta.containsKey('patient') && meta.patient == null) {
        validation.valid = false
        validation.errors << "Patient ID is null (required for channel joins)"
    }
    
    if (meta.containsKey('status') && !(meta.status instanceof Integer)) {
        validation.valid = false
        validation.errors << "Status field must be integer (0=DNA normal, 1=DNA tumor, 2=RNA)"
    }
    
    // Check for file objects in metadata (potential circular reference cause)
    meta.each { key, value ->
        if (value instanceof File) {
            validation.errors << "File object in metadata field '${key}' may cause circular references"
        }
    }
    
    return validation
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    META DATA SANITIZATION TESTING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TEST_META_SANITIZATION {
    take:
        input_channel
        test_name
        
    main:
        // Test before sanitization
        ANALYZE_META_STRUCTURE(input_channel, "${test_name}_before")
        
        // Apply sanitization
        sanitized = ANALYZE_META_STRUCTURE.out.analyzed
            .map { items ->
                def meta = items instanceof List && items.size() > 0 ? items[0] : items
                
                if (!(meta instanceof Map)) {
                    return items
                }
                
                // Create sanitized metadata
                def sanitized_meta = sanitizeMetadata(meta)
                
                // Replace metadata in items
                if (items instanceof List) {
                    def newItems = [sanitized_meta] + items.drop(1)
                    return newItems
                } else {
                    return sanitized_meta
                }
            }
        
        // Test after sanitization
        ANALYZE_META_STRUCTURE(sanitized, "${test_name}_after")
        
        // Compare before and after
        comparison = ANALYZE_META_STRUCTURE.out.reports
            .collect()
            .map { reports ->
                def beforeReports = reports.findAll { it.analysis_name.endsWith('_before') }
                def afterReports = reports.findAll { it.analysis_name.endsWith('_after') }
                
                def comparison = [
                    test_name: test_name,
                    before_issues: beforeReports.sum { it.issues.size() } ?: 0,
                    after_issues: afterReports.sum { it.issues.size() } ?: 0,
                    before_circular_refs: beforeReports.any { it.has_circular_refs },
                    after_circular_refs: afterReports.any { it.has_circular_refs },
                    before_file_objects: beforeReports.sum { it.file_object_count } ?: 0,
                    after_file_objects: afterReports.sum { it.file_object_count } ?: 0
                ]
                
                log.info "Sanitization Test [${test_name}]:"
                log.info "  Issues: ${comparison.before_issues} -> ${comparison.after_issues}"
                log.info "  Circular refs: ${comparison.before_circular_refs} -> ${comparison.after_circular_refs}"
                log.info "  File objects: ${comparison.before_file_objects} -> ${comparison.after_file_objects}"
                
                def improvement = comparison.before_issues > comparison.after_issues ||
                                 (comparison.before_circular_refs && !comparison.after_circular_refs) ||
                                 comparison.before_file_objects > comparison.after_file_objects
                
                log.info "  Sanitization effective: ${improvement}"
                
                return comparison
            }
    
    emit:
        sanitized = ANALYZE_META_STRUCTURE.out.analyzed
        comparison = comparison
}

// Sanitize metadata by removing problematic elements
def sanitizeMetadata(meta) {
    if (!(meta instanceof Map)) {
        return meta
    }
    
    def sanitized = [:]
    
    // Essential fields to preserve
    def essentialFields = ['id', 'patient', 'sample', 'status', 'single_end', 'data_type']
    
    essentialFields.each { field ->
        if (meta.containsKey(field) && meta[field] != null) {
            // Only copy primitive values, not objects
            def value = meta[field]
            if (value instanceof String || value instanceof Number || value instanceof Boolean) {
                sanitized[field] = value
            }
        }
    }
    
    return sanitized
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPREHENSIVE META ANALYZER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COMPREHENSIVE_META_ANALYSIS {
    take:
        input_channel
        analysis_name
        
    main:
        // Required fields for VCF realignment workflow
        def required_fields = ['id', 'patient', 'sample', 'status']
        
        // Run all analyses
        ANALYZE_META_STRUCTURE(input_channel, analysis_name)
        VALIDATE_META_REQUIREMENTS(ANALYZE_META_STRUCTURE.out.analyzed, required_fields, analysis_name)
        TEST_META_SANITIZATION(VALIDATE_META_REQUIREMENTS.out.validated, analysis_name)
        
        // Generate comprehensive report
        comprehensive_report = ANALYZE_META_STRUCTURE.out.reports
            .combine(VALIDATE_META_REQUIREMENTS.out.results)
            .combine(TEST_META_SANITIZATION.out.comparison)
            .map { structure_report, validation_result, sanitization_comparison ->
                def report = [
                    analysis_name: analysis_name,
                    timestamp: new Date().toString(),
                    structure_analysis: structure_report,
                    validation_result: validation_result,
                    sanitization_test: sanitization_comparison,
                    overall_health: assessOverallHealth(structure_report, validation_result, sanitization_comparison)
                ]
                
                log.info "Comprehensive Meta Analysis [${analysis_name}]:"
                log.info "  Overall health: ${report.overall_health}"
                
                return report
            }
    
    emit:
        analyzed = TEST_META_SANITIZATION.out.sanitized
        report = comprehensive_report
}

// Assess overall metadata health
def assessOverallHealth(structureReport, validationResult, sanitizationComparison) {
    def score = 100
    
    // Deduct points for issues
    if (structureReport.has_circular_refs) score -= 30
    if (structureReport.file_object_count > 0) score -= 20
    if (!validationResult.valid) score -= 25
    if (structureReport.max_depth > 5) score -= 10
    if (structureReport.issues.size() > 0) score -= (structureReport.issues.size() * 5)
    
    // Bonus points for good sanitization
    if (sanitizationComparison.after_issues < sanitizationComparison.before_issues) score += 10
    if (!sanitizationComparison.after_circular_refs && sanitizationComparison.before_circular_refs) score += 15
    
    score = Math.max(0, Math.min(100, score))
    
    if (score >= 90) return "EXCELLENT"
    if (score >= 75) return "GOOD"
    if (score >= 60) return "FAIR"
    if (score >= 40) return "POOR"
    return "CRITICAL"
}