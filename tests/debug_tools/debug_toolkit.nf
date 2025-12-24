/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Comprehensive Debug Toolkit for VCF Realignment Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DEBUG_CHANNEL_COMPREHENSIVE } from './channel_inspector'
include { COMPREHENSIVE_META_ANALYSIS } from './meta_analyzer'
include { COMPREHENSIVE_PERFORMANCE_PROFILE } from './performance_profiler'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INTEGRATED DEBUGGING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEBUG_VCF_REALIGNMENT_WORKFLOW {
    take:
        vcf_channel
        cram_channel
        debug_session_name
        
    main:
        // Debug VCF channel
        DEBUG_CHANNEL_COMPREHENSIVE(vcf_channel, "${debug_session_name}_vcf")
        COMPREHENSIVE_META_ANALYSIS(DEBUG_CHANNEL_COMPREHENSIVE.out.debugged, "${debug_session_name}_vcf_meta")
        COMPREHENSIVE_PERFORMANCE_PROFILE(COMPREHENSIVE_META_ANALYSIS.out.analyzed, "${debug_session_name}_vcf_perf")
        
        // Debug CRAM channel
        DEBUG_CHANNEL_COMPREHENSIVE(cram_channel, "${debug_session_name}_cram")
        COMPREHENSIVE_META_ANALYSIS(DEBUG_CHANNEL_COMPREHENSIVE.out.debugged, "${debug_session_name}_cram_meta")
        COMPREHENSIVE_PERFORMANCE_PROFILE(COMPREHENSIVE_META_ANALYSIS.out.analyzed, "${debug_session_name}_cram_perf")
        
        // Generate integrated debug report
        integrated_report = COMPREHENSIVE_META_ANALYSIS.out.report
            .combine(COMPREHENSIVE_PERFORMANCE_PROFILE.out.report)
            .collect()
            .map { reports ->
                generateIntegratedDebugReport(reports, debug_session_name)
            }
    
    emit:
        debugged_vcf = COMPREHENSIVE_PERFORMANCE_PROFILE.out.profiled
        debugged_cram = COMPREHENSIVE_PERFORMANCE_PROFILE.out.profiled
        debug_report = integrated_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STACK OVERFLOW DEBUGGING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEBUG_STACKOVERFLOW_ISSUES {
    take:
        input_channel
        debug_name
        
    main:
        // Comprehensive analysis focused on StackOverflow causes
        DEBUG_CHANNEL_COMPREHENSIVE(input_channel, "${debug_name}_stackoverflow")
        COMPREHENSIVE_META_ANALYSIS(DEBUG_CHANNEL_COMPREHENSIVE.out.debugged, "${debug_name}_meta_stackoverflow")
        
        // Specific StackOverflow analysis
        stackoverflow_analysis = COMPREHENSIVE_META_ANALYSIS.out.report
            .map { meta_report ->
                def analysis = [
                    debug_name: debug_name,
                    timestamp: new Date().toString(),
                    circular_references_detected: meta_report.structure_analysis.has_circular_refs,
                    file_objects_in_meta: meta_report.structure_analysis.file_object_count > 0,
                    excessive_nesting: meta_report.structure_analysis.max_depth > 10,
                    large_metadata_size: meta_report.structure_analysis.estimated_size_bytes > 1_000_000,
                    stackoverflow_risk: assessStackOverflowRisk(meta_report),
                    recommendations: generateStackOverflowRecommendations(meta_report)
                ]
                
                log.info "StackOverflow Analysis [${debug_name}]:"
                log.info "  Circular references: ${analysis.circular_references_detected}"
                log.info "  File objects in meta: ${analysis.file_objects_in_meta}"
                log.info "  Excessive nesting: ${analysis.excessive_nesting}"
                log.info "  Large metadata: ${analysis.large_metadata_size}"
                log.info "  StackOverflow risk: ${analysis.stackoverflow_risk}"
                
                if (analysis.stackoverflow_risk != "LOW") {
                    log.warn "StackOverflow risk detected: ${analysis.stackoverflow_risk}"
                    analysis.recommendations.each { rec ->
                        log.warn "  Recommendation: ${rec}"
                    }
                }
                
                return analysis
            }
    
    emit:
        analyzed = COMPREHENSIVE_META_ANALYSIS.out.analyzed
        stackoverflow_analysis = stackoverflow_analysis
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNEL JOIN DEBUGGING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEBUG_CHANNEL_JOIN_ISSUES {
    take:
        left_channel
        right_channel
        debug_name
        
    main:
        // Debug both channels individually
        DEBUG_CHANNEL_COMPREHENSIVE(left_channel, "${debug_name}_left")
        DEBUG_CHANNEL_COMPREHENSIVE(right_channel, "${debug_name}_right")
        
        // Analyze join compatibility
        join_analysis = DEBUG_CHANNEL_COMPREHENSIVE.out.debugged
            .combine(DEBUG_CHANNEL_COMPREHENSIVE.out.debugged)
            .map { left_items, right_items ->
                analyzeJoinCompatibility(left_items, right_items, debug_name)
            }
        
        // Test actual join operation
        test_join = left_channel
            .map { items -> 
                def meta = items instanceof List ? items[0] : items
                def patient_id = meta instanceof Map ? meta.patient : null
                log.info "Left channel patient ID: ${patient_id}"
                return [patient_id, items]
            }
            .join(
                right_channel.map { items ->
                    def meta = items instanceof List ? items[0] : items
                    def patient_id = meta instanceof Map ? meta.patient : null
                    log.info "Right channel patient ID: ${patient_id}"
                    return [patient_id, items]
                },
                by: 0,
                failOnMismatch: false
            )
            .map { patient_id, left_items, right_items ->
                log.info "Successful join for patient: ${patient_id}"
                return [left_items, right_items]
            }
    
    emit:
        join_result = test_join
        join_analysis = join_analysis
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UTILITY FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def generateIntegratedDebugReport(reports, sessionName) {
    def report = [
        session_name: sessionName,
        timestamp: new Date().toString(),
        total_reports: reports.size(),
        meta_analyses: reports.findAll { it.containsKey('structure_analysis') },
        performance_profiles: reports.findAll { it.containsKey('workflow_performance') },
        overall_health: "UNKNOWN"
    ]
    
    // Calculate overall health
    def meta_health_scores = report.meta_analyses.collect { 
        it.structure_analysis.overall_health 
    }
    def perf_health_scores = report.performance_profiles.collect { 
        it.overall_performance 
    }
    
    def all_scores = meta_health_scores + perf_health_scores
    def critical_count = all_scores.count { it == "CRITICAL" }
    def poor_count = all_scores.count { it == "POOR" }
    def fair_count = all_scores.count { it == "FAIR" }
    def good_count = all_scores.count { it == "GOOD" }
    def excellent_count = all_scores.count { it == "EXCELLENT" }
    
    if (critical_count > 0) {
        report.overall_health = "CRITICAL"
    } else if (poor_count > all_scores.size() / 2) {
        report.overall_health = "POOR"
    } else if (fair_count > all_scores.size() / 2) {
        report.overall_health = "FAIR"
    } else if (good_count > all_scores.size() / 2) {
        report.overall_health = "GOOD"
    } else if (excellent_count > 0) {
        report.overall_health = "EXCELLENT"
    }
    
    log.info "Integrated Debug Report [${sessionName}]:"
    log.info "  Total reports: ${report.total_reports}"
    log.info "  Meta analyses: ${report.meta_analyses.size()}"
    log.info "  Performance profiles: ${report.performance_profiles.size()}"
    log.info "  Overall health: ${report.overall_health}"
    
    return report
}

def assessStackOverflowRisk(metaReport) {
    def risk_score = 0
    
    if (metaReport.structure_analysis.has_circular_refs) risk_score += 40
    if (metaReport.structure_analysis.file_object_count > 0) risk_score += 30
    if (metaReport.structure_analysis.max_depth > 10) risk_score += 20
    if (metaReport.structure_analysis.estimated_size_bytes > 1_000_000) risk_score += 10
    
    if (risk_score >= 70) return "CRITICAL"
    if (risk_score >= 40) return "HIGH"
    if (risk_score >= 20) return "MEDIUM"
    return "LOW"
}

def generateStackOverflowRecommendations(metaReport) {
    def recommendations = []
    
    if (metaReport.structure_analysis.has_circular_refs) {
        recommendations << "Remove circular references from metadata structures"
        recommendations << "Use SANITIZE_CHANNELS workflow to clean metadata"
    }
    
    if (metaReport.structure_analysis.file_object_count > 0) {
        recommendations << "Remove File objects from metadata to prevent circular references"
        recommendations << "Store file paths as strings instead of File objects"
    }
    
    if (metaReport.structure_analysis.max_depth > 10) {
        recommendations << "Reduce metadata nesting depth"
        recommendations << "Flatten complex metadata structures"
    }
    
    if (metaReport.structure_analysis.estimated_size_bytes > 1_000_000) {
        recommendations << "Reduce metadata size by removing unnecessary fields"
        recommendations << "Consider using references instead of embedding large data"
    }
    
    if (recommendations.isEmpty()) {
        recommendations << "Metadata structure appears healthy for StackOverflow prevention"
    }
    
    return recommendations
}

def analyzeJoinCompatibility(leftItems, rightItems, debugName) {
    def leftMeta = leftItems instanceof List ? leftItems[0] : leftItems
    def rightMeta = rightItems instanceof List ? rightItems[0] : rightItems
    
    def analysis = [
        debug_name: debugName,
        left_has_meta: leftMeta instanceof Map,
        right_has_meta: rightMeta instanceof Map,
        left_patient_id: leftMeta instanceof Map ? leftMeta.patient : null,
        right_patient_id: rightMeta instanceof Map ? rightMeta.patient : null,
        patient_ids_match: false,
        join_compatible: false
    ]
    
    if (analysis.left_has_meta && analysis.right_has_meta) {
        analysis.patient_ids_match = analysis.left_patient_id == analysis.right_patient_id
        analysis.join_compatible = analysis.patient_ids_match && 
                                  analysis.left_patient_id != null && 
                                  analysis.right_patient_id != null
    }
    
    log.info "Join Compatibility Analysis [${debugName}]:"
    log.info "  Left patient ID: ${analysis.left_patient_id}"
    log.info "  Right patient ID: ${analysis.right_patient_id}"
    log.info "  Patient IDs match: ${analysis.patient_ids_match}"
    log.info "  Join compatible: ${analysis.join_compatible}"
    
    return analysis
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN DEBUG TOOLKIT WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEBUG_TOOLKIT_COMPREHENSIVE {
    take:
        vcf_channel
        cram_channel
        toolkit_session_name
        
    main:
        // Run comprehensive debugging
        DEBUG_VCF_REALIGNMENT_WORKFLOW(vcf_channel, cram_channel, toolkit_session_name)
        
        // Run specific issue debugging
        DEBUG_STACKOVERFLOW_ISSUES(vcf_channel, "${toolkit_session_name}_vcf_stackoverflow")
        DEBUG_STACKOVERFLOW_ISSUES(cram_channel, "${toolkit_session_name}_cram_stackoverflow")
        DEBUG_CHANNEL_JOIN_ISSUES(vcf_channel, cram_channel, "${toolkit_session_name}_join")
        
        // Generate final toolkit report
        toolkit_report = DEBUG_VCF_REALIGNMENT_WORKFLOW.out.debug_report
            .combine(DEBUG_STACKOVERFLOW_ISSUES.out.stackoverflow_analysis)
            .combine(DEBUG_CHANNEL_JOIN_ISSUES.out.join_analysis)
            .collect()
            .map { reports ->
                def final_report = [
                    toolkit_session: toolkit_session_name,
                    timestamp: new Date().toString(),
                    comprehensive_debug: reports[0],
                    stackoverflow_analyses: reports.findAll { it.containsKey('stackoverflow_risk') },
                    join_analyses: reports.findAll { it.containsKey('join_compatible') },
                    summary: generateToolkitSummary(reports)
                ]
                
                log.info "Debug Toolkit Summary [${toolkit_session_name}]:"
                log.info "  ${final_report.summary}"
                
                return final_report
            }
    
    emit:
        debugged_vcf = DEBUG_VCF_REALIGNMENT_WORKFLOW.out.debugged_vcf
        debugged_cram = DEBUG_VCF_REALIGNMENT_WORKFLOW.out.debugged_cram
        toolkit_report = toolkit_report
}

def generateToolkitSummary(reports) {
    def issues = []
    def recommendations = []
    
    // Analyze all reports for issues
    reports.each { report ->
        if (report.containsKey('stackoverflow_risk') && report.stackoverflow_risk != "LOW") {
            issues << "StackOverflow risk: ${report.stackoverflow_risk}"
            recommendations.addAll(report.recommendations)
        }
        
        if (report.containsKey('join_compatible') && !report.join_compatible) {
            issues << "Channel join compatibility issue"
            recommendations << "Check patient ID consistency between channels"
        }
        
        if (report.containsKey('overall_health') && report.overall_health in ["CRITICAL", "POOR"]) {
            issues << "Overall health: ${report.overall_health}"
        }
    }
    
    def summary = "Debug analysis complete. "
    if (issues.isEmpty()) {
        summary += "No critical issues detected."
    } else {
        summary += "${issues.size()} issue(s) found: ${issues.join(', ')}"
    }
    
    return summary
}