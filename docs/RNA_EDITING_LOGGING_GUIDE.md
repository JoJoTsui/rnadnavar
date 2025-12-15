# RNA Editing Annotation Logging and Monitoring Guide

This guide describes the comprehensive logging and monitoring integration implemented for the RNA editing annotation pipeline in nf-core/rnadnavar.

## Overview

The RNA editing annotation module now includes comprehensive logging and monitoring capabilities designed for operational use and integration with Nextflow's built-in logging and reporting system.

## Features

### 1. Comprehensive Process Logging
- **Structured logging format** compatible with Nextflow reporting
- **Performance metrics collection** (processing time, memory usage, file sizes)
- **Processing statistics reporting** (annotation rates, evidence classification)
- **Diagnostic information** for troubleshooting annotation failures
- **Configurable log levels** for different operational environments

### 2. Operational Monitoring
- **Resource usage tracking** (memory, CPU, disk space)
- **Processing rate monitoring** (variants per second)
- **Error and warning categorization**
- **Success/failure metrics**
- **Performance benchmarking**

### 3. Nextflow Integration
- **Automatic log file publishing** to output directories
- **Metrics export** in JSON format for monitoring systems
- **Environment variable configuration** for different deployment scenarios
- **Graceful error handling** with detailed troubleshooting information

## Configuration Parameters

### Basic Configuration
```nextflow
params {
    // Enable RNA editing annotation
    enable_rna_annotation = true
    rediportal_vcf = "/path/to/rediportal.vcf.gz"
    min_rna_support = 2
}
```

### Logging Configuration
```nextflow
params {
    // Logging and monitoring options
    rna_annotation_log_level = 'INFO'           // DEBUG, INFO, WARNING, ERROR
    rna_annotation_verbose = false              // Enable verbose logging
    rna_annotation_enable_metrics = true       // Enable metrics collection
    rna_annotation_enable_performance = true   // Enable performance monitoring
    rna_annotation_save_reports = false        // Save detailed reports
}
```

### Environment Variables
The logging system can also be configured using environment variables:
- `RNA_EDITING_LOG_LEVEL`: Set logging level (DEBUG, INFO, WARNING, ERROR)
- `RNA_EDITING_LOG_DIR`: Directory for log files
- `RNA_EDITING_ENABLE_METRICS`: Enable metrics collection (true/false)
- `RNA_EDITING_ENABLE_PERFORMANCE`: Enable performance monitoring (true/false)

## Output Files

### Log Files
The RNA editing annotation process generates several types of log files:

1. **Process Log** (`rna_editing_process_<sample>_<timestamp>.log`)
   - Detailed process execution log
   - System resource information
   - Input/output file validation
   - Error messages and troubleshooting information

2. **Metrics Log** (`rna_editing_metrics_<sample>_<timestamp>.log`)
   - Performance metrics
   - Resource usage statistics
   - Processing rates and timing

3. **Statistics Log** (`rna_editing_stats_<sample>_<timestamp>.log`)
   - Processing statistics
   - Annotation rates
   - Evidence classification results

### Metrics Files
1. **JSON Metrics Report** (`rna_editing_metrics_<timestamp>.json`)
   - Structured metrics data
   - Performance indicators
   - Operational status
   - Integration-ready format for monitoring systems

### Published Locations
```
outdir/
├── annotation/rna_editing/<sample>/     # Main output files
│   ├── <sample>.rescue.rna_annotated.vcf.gz
│   └── <sample>.rescue.rna_annotated.vcf.gz.tbi
├── logs/rna_editing/<sample>/           # Log files
│   ├── rna_editing_process.log
│   ├── rna_editing_metrics.json
│   └── rna_editing_statistics.log
└── reports/rna_editing/<sample>/        # Detailed reports (optional)
    └── comprehensive_metrics.json
```

## Monitoring Integration

### Key Performance Indicators (KPIs)
The logging system tracks several KPIs for operational monitoring:

- **Processing Rate**: Variants processed per second
- **Annotation Rate**: Percentage of variants successfully annotated
- **Success Rate**: Percentage of successful pipeline executions
- **Resource Efficiency**: Memory and CPU utilization
- **Error Rate**: Number of errors per execution

### Metrics Collection
```json
{
  "pipeline_info": {
    "name": "rna_editing_annotation",
    "version": "1.0.0",
    "timestamp": 1765776464,
    "total_duration_seconds": 120.5
  },
  "performance_metrics": {
    "processing_rate_variants_per_second": 850,
    "total_variants": 102400,
    "annotated_variants": 15680,
    "annotation_rate_percent": 15.3
  },
  "quality_metrics": {
    "success": true,
    "error_count": 0,
    "warning_count": 2
  }
}
```

## Troubleshooting

### Common Issues and Solutions

#### 1. High Memory Usage
**Symptoms**: Process killed with exit code 137
**Solution**: Increase memory allocation
```nextflow
params.rna_annotation_memory = '8.GB'
```

#### 2. Slow Processing
**Symptoms**: Processing rate < 100 variants/second
**Solutions**:
- Increase CPU allocation: `params.rna_annotation_cpus = 4`
- Check disk I/O performance
- Verify REDIportal database accessibility

#### 3. Missing Log Files
**Symptoms**: Log files not published to output directory
**Solution**: Check publishing configuration
```nextflow
params.rna_annotation_save_reports = true
```

### Error Codes
- **Exit Code 2**: File not found or invalid arguments
- **Exit Code 13**: Permission denied
- **Exit Code 124**: Process timeout
- **Exit Code 127**: Required tools not found
- **Exit Code 137**: Out of memory (SIGKILL)

## Best Practices

### 1. Production Deployment
```nextflow
params {
    rna_annotation_log_level = 'INFO'
    rna_annotation_enable_metrics = true
    rna_annotation_enable_performance = true
    rna_annotation_save_reports = true
}
```

### 2. Development/Debugging
```nextflow
params {
    rna_annotation_log_level = 'DEBUG'
    rna_annotation_verbose = true
    rna_annotation_enable_metrics = true
    rna_annotation_enable_performance = true
}
```

### 3. Resource Optimization
- Monitor processing rates and adjust CPU/memory allocation accordingly
- Use metrics data to optimize resource allocation for different file sizes
- Set appropriate timeouts based on expected processing times

### 4. Monitoring Integration
- Export metrics to monitoring systems (Prometheus, Grafana, etc.)
- Set up alerts for error rates and performance degradation
- Use log aggregation systems for centralized monitoring

## Integration with Monitoring Systems

### Prometheus Integration
The JSON metrics files can be easily integrated with Prometheus using file-based service discovery or custom exporters.

### Grafana Dashboards
Key metrics to monitor:
- Processing rate trends
- Error rate over time
- Resource utilization
- Success/failure ratios

### Log Aggregation
The structured log format is compatible with:
- ELK Stack (Elasticsearch, Logstash, Kibana)
- Fluentd/Fluent Bit
- Splunk
- CloudWatch Logs

## Requirements Satisfied

This implementation satisfies the following requirements:
- **8.1**: Comprehensive logging to RNA editing annotation module
- **8.2**: Processing statistics reporting (annotation rates, evidence classification)
- **8.3**: Performance metrics logging (processing time, memory usage, file sizes)
- **8.4**: Integration with Nextflow's built-in logging and reporting system
- **8.5**: Diagnostic information for troubleshooting annotation failures
- **8.5**: Log levels and output formatting for operational monitoring

## Support

For issues related to RNA editing annotation logging and monitoring:
1. Check the process log files for detailed error information
2. Review the metrics JSON files for performance indicators
3. Verify configuration parameters are set correctly
4. Ensure all required tools are available in the environment
5. Check resource allocation (memory, CPU, disk space)