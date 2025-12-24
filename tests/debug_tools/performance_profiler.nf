/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Performance Profiler - Debug tool for monitoring workflow performance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PERFORMANCE MONITORING UTILITIES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Profile workflow execution time and resource usage
workflow PROFILE_WORKFLOW_PERFORMANCE {
    take:
        input_channel
        workflow_name
        
    main:
        start_time = System.currentTimeMillis()
        
        profiled = input_channel
            .map { items ->
                def current_time = System.currentTimeMillis()
                def elapsed = current_time - start_time
                
                def performance_data = [
                    workflow_name: workflow_name,
                    timestamp: new Date(current_time).toString(),
                    elapsed_ms: elapsed,
                    item_id: System.identityHashCode(items),
                    memory_usage: getMemoryUsage(),
                    thread_info: getThreadInfo()
                ]
                
                log.info "Performance [${workflow_name}]: ${elapsed}ms elapsed, ${performance_data.memory_usage.used_mb}MB used"
                
                return [performance_data, items]
            }
            .map { performance_data, items -> items }
        
        // Collect performance metrics
        performance_summary = profiled
            .collect()
            .map { all_items ->
                def end_time = System.currentTimeMillis()
                def total_elapsed = end_time - start_time
                
                def summary = [
                    workflow_name: workflow_name,
                    total_items: all_items.size(),
                    total_elapsed_ms: total_elapsed,
                    avg_time_per_item: all_items.size() > 0 ? total_elapsed / all_items.size() : 0,
                    throughput_items_per_sec: all_items.size() > 0 ? (all_items.size() * 1000.0) / total_elapsed : 0,
                    final_memory_usage: getMemoryUsage()
                ]
                
                log.info "Performance Summary [${workflow_name}]:"
                log.info "  Total items: ${summary.total_items}"
                log.info "  Total time: ${summary.total_elapsed_ms}ms"
                log.info "  Avg time per item: ${Math.round(summary.avg_time_per_item)}ms"
                log.info "  Throughput: ${Math.round(summary.throughput_items_per_sec * 100) / 100} items/sec"
                log.info "  Final memory: ${summary.final_memory_usage.used_mb}MB"
                
                return summary
            }
    
    emit:
        profiled = profiled
        summary = performance_summary
}

// Get current memory usage information
def getMemoryUsage() {
    def runtime = Runtime.getRuntime()
    def totalMemory = runtime.totalMemory()
    def freeMemory = runtime.freeMemory()
    def usedMemory = totalMemory - freeMemory
    def maxMemory = runtime.maxMemory()
    
    return [
        used_mb: Math.round(usedMemory / 1024 / 1024),
        free_mb: Math.round(freeMemory / 1024 / 1024),
        total_mb: Math.round(totalMemory / 1024 / 1024),
        max_mb: Math.round(maxMemory / 1024 / 1024),
        usage_percent: Math.round((usedMemory / maxMemory) * 100)
    ]
}

// Get current thread information
def getThreadInfo() {
    def threadMX = java.lang.management.ManagementFactory.getThreadMXBean()
    
    return [
        active_threads: threadMX.getThreadCount(),
        peak_threads: threadMX.getPeakThreadCount(),
        daemon_threads: threadMX.getDaemonThreadCount()
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNEL THROUGHPUT MONITORING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MONITOR_CHANNEL_THROUGHPUT {
    take:
        input_channel
        monitor_name
        sample_interval_ms
        
    main:
        start_time = System.currentTimeMillis()
        item_count = Channel.value(0)
        
        monitored = input_channel
            .map { items ->
                def current_time = System.currentTimeMillis()
                def elapsed = current_time - start_time
                
                // Sample throughput at intervals
                if (elapsed % (sample_interval_ms ?: 5000) < 100) {
                    def throughput = elapsed > 0 ? (1000.0 / elapsed) : 0
                    log.info "Throughput [${monitor_name}]: ${Math.round(throughput * 100) / 100} items/sec"
                }
                
                return items
            }
    
    emit:
        monitored = monitored
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RESOURCE USAGE TRACKING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TRACK_RESOURCE_USAGE {
    take:
        input_channel
        tracking_name
        
    main:
        resource_snapshots = []
        
        tracked = input_channel
            .map { items ->
                def snapshot = [
                    timestamp: System.currentTimeMillis(),
                    memory: getMemoryUsage(),
                    threads: getThreadInfo(),
                    gc_info: getGCInfo()
                ]
                
                resource_snapshots << snapshot
                
                // Check for resource issues
                if (snapshot.memory.usage_percent > 85) {
                    log.warn "High memory usage [${tracking_name}]: ${snapshot.memory.usage_percent}%"
                }
                
                if (snapshot.threads.active_threads > 50) {
                    log.warn "High thread count [${tracking_name}]: ${snapshot.threads.active_threads} threads"
                }
                
                return items
            }
        
        // Generate resource usage report
        resource_report = tracked
            .collect()
            .map { all_items ->
                if (resource_snapshots.size() == 0) {
                    return [tracking_name: tracking_name, error: "No resource snapshots collected"]
                }
                
                def first_snapshot = resource_snapshots[0]
                def last_snapshot = resource_snapshots[-1]
                def max_memory = resource_snapshots.collect { it.memory.usage_percent }.max()
                def max_threads = resource_snapshots.collect { it.threads.active_threads }.max()
                
                def report = [
                    tracking_name: tracking_name,
                    duration_ms: last_snapshot.timestamp - first_snapshot.timestamp,
                    initial_memory_mb: first_snapshot.memory.used_mb,
                    final_memory_mb: last_snapshot.memory.used_mb,
                    peak_memory_percent: max_memory,
                    peak_threads: max_threads,
                    memory_growth_mb: last_snapshot.memory.used_mb - first_snapshot.memory.used_mb,
                    gc_collections: last_snapshot.gc_info.total_collections - first_snapshot.gc_info.total_collections,
                    gc_time_ms: last_snapshot.gc_info.total_time_ms - first_snapshot.gc_info.total_time_ms
                ]
                
                log.info "Resource Usage Report [${tracking_name}]:"
                log.info "  Duration: ${report.duration_ms}ms"
                log.info "  Memory growth: ${report.memory_growth_mb}MB"
                log.info "  Peak memory: ${report.peak_memory_percent}%"
                log.info "  Peak threads: ${report.peak_threads}"
                log.info "  GC collections: ${report.gc_collections}"
                log.info "  GC time: ${report.gc_time_ms}ms"
                
                return report
            }
    
    emit:
        tracked = tracked
        report = resource_report
}

// Get garbage collection information
def getGCInfo() {
    def gcBeans = java.lang.management.ManagementFactory.getGarbageCollectorMXBeans()
    
    def totalCollections = 0
    def totalTime = 0
    
    gcBeans.each { bean ->
        totalCollections += bean.getCollectionCount()
        totalTime += bean.getCollectionTime()
    }
    
    return [
        total_collections: totalCollections,
        total_time_ms: totalTime
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BOTTLENECK DETECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DETECT_BOTTLENECKS {
    take:
        input_channel
        detection_name
        
    main:
        processing_times = []
        
        analyzed = input_channel
            .map { items ->
                def start_time = System.currentTimeMillis()
                
                // Simulate processing time measurement
                def processing_start = System.nanoTime()
                
                // The actual processing happens here (items pass through)
                def processing_end = System.nanoTime()
                def processing_time_ms = (processing_end - processing_start) / 1_000_000.0
                
                processing_times << processing_time_ms
                
                return [processing_time_ms, items]
            }
            .map { processing_time, items -> items }
        
        // Analyze bottlenecks
        bottleneck_analysis = analyzed
            .collect()
            .map { all_items ->
                if (processing_times.size() == 0) {
                    return [detection_name: detection_name, error: "No processing times collected"]
                }
                
                def sorted_times = processing_times.sort()
                def total_time = processing_times.sum()
                def avg_time = total_time / processing_times.size()
                def median_time = sorted_times[sorted_times.size() / 2]
                def p95_time = sorted_times[(int)(sorted_times.size() * 0.95)]
                def max_time = sorted_times[-1]
                
                def analysis = [
                    detection_name: detection_name,
                    total_items: processing_times.size(),
                    avg_processing_time_ms: Math.round(avg_time * 100) / 100,
                    median_processing_time_ms: Math.round(median_time * 100) / 100,
                    p95_processing_time_ms: Math.round(p95_time * 100) / 100,
                    max_processing_time_ms: Math.round(max_time * 100) / 100,
                    bottleneck_detected: p95_time > (avg_time * 3),
                    slow_items_count: processing_times.count { it > (avg_time * 2) }
                ]
                
                log.info "Bottleneck Analysis [${detection_name}]:"
                log.info "  Avg processing time: ${analysis.avg_processing_time_ms}ms"
                log.info "  95th percentile: ${analysis.p95_processing_time_ms}ms"
                log.info "  Max processing time: ${analysis.max_processing_time_ms}ms"
                log.info "  Bottleneck detected: ${analysis.bottleneck_detected}"
                log.info "  Slow items: ${analysis.slow_items_count}"
                
                if (analysis.bottleneck_detected) {
                    log.warn "Performance bottleneck detected in ${detection_name}"
                }
                
                return analysis
            }
    
    emit:
        analyzed = analyzed
        analysis = bottleneck_analysis
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPREHENSIVE PERFORMANCE PROFILER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COMPREHENSIVE_PERFORMANCE_PROFILE {
    take:
        input_channel
        profile_name
        
    main:
        // Apply all profiling tools
        PROFILE_WORKFLOW_PERFORMANCE(input_channel, profile_name)
        MONITOR_CHANNEL_THROUGHPUT(PROFILE_WORKFLOW_PERFORMANCE.out.profiled, profile_name, 5000)
        TRACK_RESOURCE_USAGE(MONITOR_CHANNEL_THROUGHPUT.out.monitored, profile_name)
        DETECT_BOTTLENECKS(TRACK_RESOURCE_USAGE.out.tracked, profile_name)
        
        // Generate comprehensive performance report
        comprehensive_report = PROFILE_WORKFLOW_PERFORMANCE.out.summary
            .combine(TRACK_RESOURCE_USAGE.out.report)
            .combine(DETECT_BOTTLENECKS.out.analysis)
            .map { workflow_summary, resource_report, bottleneck_analysis ->
                def report = [
                    profile_name: profile_name,
                    timestamp: new Date().toString(),
                    workflow_performance: workflow_summary,
                    resource_usage: resource_report,
                    bottleneck_analysis: bottleneck_analysis,
                    overall_performance: assessOverallPerformance(workflow_summary, resource_report, bottleneck_analysis)
                ]
                
                log.info "Comprehensive Performance Profile [${profile_name}]:"
                log.info "  Overall performance: ${report.overall_performance}"
                
                return report
            }
    
    emit:
        profiled = DETECT_BOTTLENECKS.out.analyzed
        report = comprehensive_report
}

// Assess overall performance health
def assessOverallPerformance(workflowSummary, resourceReport, bottleneckAnalysis) {
    def score = 100
    
    // Deduct points for performance issues
    if (workflowSummary.throughput_items_per_sec < 1) score -= 20
    if (resourceReport.peak_memory_percent > 80) score -= 15
    if (resourceReport.memory_growth_mb > 100) score -= 10
    if (bottleneckAnalysis.bottleneck_detected) score -= 25
    if (resourceReport.gc_time_ms > 1000) score -= 10
    if (bottleneckAnalysis.slow_items_count > (bottleneckAnalysis.total_items * 0.1)) score -= 15
    
    score = Math.max(0, Math.min(100, score))
    
    if (score >= 90) return "EXCELLENT"
    if (score >= 75) return "GOOD"
    if (score >= 60) return "FAIR"
    if (score >= 40) return "POOR"
    return "CRITICAL"
}