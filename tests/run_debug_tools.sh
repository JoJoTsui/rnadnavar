#!/bin/bash

#
# Debug Tools Runner for VCF Realignment Optimization
# Provides comprehensive debugging and diagnostic capabilities
#

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DEBUG_CONFIG="$SCRIPT_DIR/config/property_test.config"
RESULTS_DIR="$SCRIPT_DIR/results/debug_tools"
LOG_FILE="$RESULTS_DIR/debug_tools_run.log"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_FILE"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOG_FILE"
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1" | tee -a "$LOG_FILE"
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1" | tee -a "$LOG_FILE"
}

info() {
    echo -e "${CYAN}[INFO]${NC} $1" | tee -a "$LOG_FILE"
}

# Setup
setup_debug_environment() {
    log "Setting up debug tools environment..."
    
    # Create results directory
    mkdir -p "$RESULTS_DIR"
    
    # Initialize log file
    echo "Debug Tools Run - $(date)" > "$LOG_FILE"
    echo "========================================" >> "$LOG_FILE"
    
    # Check dependencies
    if ! command -v nextflow &> /dev/null; then
        error "Nextflow is not installed or not in PATH"
        exit 1
    fi
    
    log "Debug environment setup complete"
}

# Channel Inspector
run_channel_inspector() {
    local test_name="$1"
    local vcf_input="$2"
    local cram_input="$3"
    
    log "Running channel inspector for: $test_name"
    
    local output_dir="$RESULTS_DIR/channel_inspector_$test_name"
    mkdir -p "$output_dir"
    
    # Create a simple test workflow to inspect channels
    cat > "$output_dir/test_inspector.nf" << EOF
#!/usr/bin/env nextflow

include { DEBUG_CHANNEL_COMPREHENSIVE } from '${SCRIPT_DIR}/debug_tools/channel_inspector'

workflow {
    // Create test channels
    vcf_ch = Channel.fromPath("$vcf_input")
        .map { vcf -> 
            def meta = [id: "test_sample", patient: "test_patient", status: 2]
            return [meta, vcf, file("\${vcf}.tbi")]
        }
    
    cram_ch = Channel.fromPath("$cram_input")
        .map { cram ->
            def meta = [id: "test_sample", patient: "test_patient", status: 2]
            return [meta, cram, file("\${cram}.crai")]
        }
    
    // Run channel inspection
    DEBUG_CHANNEL_COMPREHENSIVE(vcf_ch, "vcf_inspection")
    DEBUG_CHANNEL_COMPREHENSIVE(cram_ch, "cram_inspection")
    
    // Output results
    DEBUG_CHANNEL_COMPREHENSIVE.out.report.view { report ->
        println "Channel Inspection Report: \$report"
    }
}
EOF

    if nextflow run "$output_dir/test_inspector.nf" \
        -config "$DEBUG_CONFIG" \
        -work-dir "$output_dir/work" \
        >> "$LOG_FILE" 2>&1; then
        success "Channel inspector completed for $test_name"
        return 0
    else
        error "Channel inspector failed for $test_name"
        return 1
    fi
}

# Meta Analyzer
run_meta_analyzer() {
    local test_name="$1"
    local input_file="$2"
    
    log "Running meta analyzer for: $test_name"
    
    local output_dir="$RESULTS_DIR/meta_analyzer_$test_name"
    mkdir -p "$output_dir"
    
    # Create a test workflow for meta analysis
    cat > "$output_dir/test_meta_analyzer.nf" << EOF
#!/usr/bin/env nextflow

include { COMPREHENSIVE_META_ANALYSIS } from '${SCRIPT_DIR}/debug_tools/meta_analyzer'

workflow {
    // Create test channel with complex metadata
    test_ch = Channel.fromPath("$input_file")
        .map { file ->
            def complex_meta = [
                id: "test_sample_\${System.currentTimeMillis()}",
                patient: "patient_001",
                sample: "sample_001",
                status: 2,
                single_end: false,
                data_type: "vcf",
                nested_data: [
                    processing: [
                        timestamp: new Date().toString(),
                        version: "1.0.0"
                    ],
                    files: [
                        input: file.toString(),
                        output: "processed_\${file.name}"
                    ]
                ],
                // Potential circular reference test
                self_ref: null
            ]
            
            // Add self-reference to test circular detection
            complex_meta.self_ref = complex_meta
            
            return [complex_meta, file]
        }
    
    // Run comprehensive meta analysis
    COMPREHENSIVE_META_ANALYSIS(test_ch, "meta_analysis_$test_name")
    
    // Output results
    COMPREHENSIVE_META_ANALYSIS.out.report.view { report ->
        println "Meta Analysis Report: \$report"
    }
}
EOF

    if nextflow run "$output_dir/test_meta_analyzer.nf" \
        -config "$DEBUG_CONFIG" \
        -work-dir "$output_dir/work" \
        >> "$LOG_FILE" 2>&1; then
        success "Meta analyzer completed for $test_name"
        return 0
    else
        error "Meta analyzer failed for $test_name"
        return 1
    fi
}

# Performance Profiler
run_performance_profiler() {
    local test_name="$1"
    local input_file="$2"
    local iterations="${3:-50}"
    
    log "Running performance profiler for: $test_name ($iterations iterations)"
    
    local output_dir="$RESULTS_DIR/performance_profiler_$test_name"
    mkdir -p "$output_dir"
    
    # Create a test workflow for performance profiling
    cat > "$output_dir/test_performance_profiler.nf" << EOF
#!/usr/bin/env nextflow

include { COMPREHENSIVE_PERFORMANCE_PROFILE } from '${SCRIPT_DIR}/debug_tools/performance_profiler'

workflow {
    // Create test channel with multiple items
    test_ch = Channel.from(1..$iterations)
        .map { i ->
            def meta = [
                id: "test_sample_\$i",
                patient: "patient_\${i % 5}",
                status: i % 3,
                iteration: i
            ]
            return [meta, file("$input_file")]
        }
    
    // Run comprehensive performance profiling
    COMPREHENSIVE_PERFORMANCE_PROFILE(test_ch, "performance_profile_$test_name")
    
    // Output results
    COMPREHENSIVE_PERFORMANCE_PROFILE.out.report.view { report ->
        println "Performance Profile Report: \$report"
    }
}
EOF

    if nextflow run "$output_dir/test_performance_profiler.nf" \
        -config "$DEBUG_CONFIG" \
        -work-dir "$output_dir/work" \
        >> "$LOG_FILE" 2>&1; then
        success "Performance profiler completed for $test_name"
        return 0
    else
        error "Performance profiler failed for $test_name"
        return 1
    fi
}

# Comprehensive Debug Toolkit
run_debug_toolkit() {
    local test_name="$1"
    local vcf_input="$2"
    local cram_input="$3"
    
    log "Running comprehensive debug toolkit for: $test_name"
    
    local output_dir="$RESULTS_DIR/debug_toolkit_$test_name"
    mkdir -p "$output_dir"
    
    # Create a comprehensive test workflow
    cat > "$output_dir/test_debug_toolkit.nf" << EOF
#!/usr/bin/env nextflow

include { DEBUG_TOOLKIT_COMPREHENSIVE } from '${SCRIPT_DIR}/debug_tools/debug_toolkit'

workflow {
    // Create VCF channel
    vcf_ch = Channel.fromPath("$vcf_input")
        .map { vcf -> 
            def meta = [
                id: "test_vcf_sample",
                patient: "patient_001",
                sample: "vcf_sample_001",
                status: 2,
                single_end: false,
                data_type: "vcf"
            ]
            return [meta, vcf, file("\${vcf}.tbi")]
        }
    
    // Create CRAM channel
    cram_ch = Channel.fromPath("$cram_input")
        .map { cram ->
            def meta = [
                id: "test_cram_sample",
                patient: "patient_001",
                sample: "cram_sample_001", 
                status: 2,
                single_end: false,
                data_type: "cram"
            ]
            return [meta, cram, file("\${cram}.crai")]
        }
    
    // Run comprehensive debug toolkit
    DEBUG_TOOLKIT_COMPREHENSIVE(vcf_ch, cram_ch, "toolkit_$test_name")
    
    // Output results
    DEBUG_TOOLKIT_COMPREHENSIVE.out.toolkit_report.view { report ->
        println "Debug Toolkit Report: \$report"
    }
}
EOF

    if nextflow run "$output_dir/test_debug_toolkit.nf" \
        -config "$DEBUG_CONFIG" \
        -work-dir "$output_dir/work" \
        >> "$LOG_FILE" 2>&1; then
        success "Debug toolkit completed for $test_name"
        return 0
    else
        error "Debug toolkit failed for $test_name"
        return 1
    fi
}

# Generate debug report
generate_debug_report() {
    log "Generating comprehensive debug report..."
    
    local report_file="$RESULTS_DIR/debug_tools_report.md"
    
    cat > "$report_file" << EOF
# Debug Tools Report

**Generated:** $(date)
**Debug Configuration:** $DEBUG_CONFIG
**Results Directory:** $RESULTS_DIR

## Debug Tools Overview

This report contains the results of running comprehensive debugging and diagnostic tools for the VCF Realignment Optimization feature.

### Available Debug Tools

1. **Channel Inspector** - Analyzes channel contents, structure, and flow
2. **Meta Analyzer** - Deep analysis of metadata structures for circular references
3. **Performance Profiler** - Monitors resource usage, throughput, and bottlenecks
4. **Debug Toolkit** - Comprehensive debugging combining all tools

### Tool Capabilities

#### Channel Inspector
- Channel content inspection
- Channel state monitoring
- Channel join debugging
- Memory usage monitoring

#### Meta Analyzer
- Deep metadata structure analysis
- Circular reference detection
- File object detection
- Metadata validation
- Sanitization testing

#### Performance Profiler
- Workflow execution profiling
- Channel throughput monitoring
- Resource usage tracking
- Bottleneck detection

#### Debug Toolkit
- Integrated debugging workflow
- StackOverflow issue analysis
- Channel join issue debugging
- Comprehensive reporting

### Usage Examples

\`\`\`bash
# Run channel inspector
./run_debug_tools.sh channel_inspector test_name /path/to/vcf /path/to/cram

# Run meta analyzer
./run_debug_tools.sh meta_analyzer test_name /path/to/input

# Run performance profiler
./run_debug_tools.sh performance_profiler test_name /path/to/input 100

# Run comprehensive toolkit
./run_debug_tools.sh debug_toolkit test_name /path/to/vcf /path/to/cram

# Run all tools
./run_debug_tools.sh all
\`\`\`

### Requirements Addressed

- **Requirements 6.5:** Debugging and diagnostic tools
- **Requirements 1.1-1.5:** StackOverflow error analysis and prevention
- **Requirements 3.1-3.5:** Channel validation and monitoring
- **Requirements 5.1-5.5:** Error handling and command validation

### Results

See individual tool logs in the results directory for detailed information.

EOF

    log "Debug report generated: $report_file"
}

# Run all debug tools
run_all_debug_tools() {
    log "Running all debug tools with test data..."
    
    # Create test data if it doesn't exist
    local test_vcf="$RESULTS_DIR/test_data/test.vcf"
    local test_cram="$RESULTS_DIR/test_data/test.cram"
    
    mkdir -p "$RESULTS_DIR/test_data"
    
    # Create minimal test VCF
    if [[ ! -f "$test_vcf" ]]; then
        cat > "$test_vcf" << EOF
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	1000	.	A	T	60	PASS	.	GT	0/1
chr1	2000	.	G	C	60	PASS	.	GT	0/1
EOF
        touch "$test_vcf.tbi"
    fi
    
    # Create minimal test CRAM
    if [[ ! -f "$test_cram" ]]; then
        echo "mock_cram_data" > "$test_cram"
        touch "$test_cram.crai"
    fi
    
    local exit_code=0
    
    # Run all tools
    if ! run_channel_inspector "all_test" "$test_vcf" "$test_cram"; then
        exit_code=1
    fi
    
    if ! run_meta_analyzer "all_test" "$test_vcf"; then
        exit_code=1
    fi
    
    if ! run_performance_profiler "all_test" "$test_vcf" 25; then
        exit_code=1
    fi
    
    if ! run_debug_toolkit "all_test" "$test_vcf" "$test_cram"; then
        exit_code=1
    fi
    
    return $exit_code
}

# Main execution
main() {
    local command="${1:-help}"
    
    case "$command" in
        "channel_inspector")
            setup_debug_environment
            run_channel_inspector "${2:-test}" "${3:-}" "${4:-}"
            ;;
        "meta_analyzer")
            setup_debug_environment
            run_meta_analyzer "${2:-test}" "${3:-}"
            ;;
        "performance_profiler")
            setup_debug_environment
            run_performance_profiler "${2:-test}" "${3:-}" "${4:-50}"
            ;;
        "debug_toolkit")
            setup_debug_environment
            run_debug_toolkit "${2:-test}" "${3:-}" "${4:-}"
            ;;
        "all")
            setup_debug_environment
            if run_all_debug_tools; then
                success "All debug tools completed successfully"
            else
                error "Some debug tools failed"
                exit 1
            fi
            ;;
        "report")
            generate_debug_report
            ;;
        "help"|*)
            echo "Debug Tools Runner for VCF Realignment Optimization"
            echo ""
            echo "Usage: $0 <command> [arguments...]"
            echo ""
            echo "Commands:"
            echo "  channel_inspector <name> <vcf_path> <cram_path>  - Run channel inspector"
            echo "  meta_analyzer <name> <input_path>               - Run meta analyzer"
            echo "  performance_profiler <name> <input_path> [iter] - Run performance profiler"
            echo "  debug_toolkit <name> <vcf_path> <cram_path>     - Run comprehensive toolkit"
            echo "  all                                              - Run all debug tools"
            echo "  report                                           - Generate debug report"
            echo "  help                                             - Show this help"
            echo ""
            echo "Examples:"
            echo "  $0 all                                           # Run all tools with test data"
            echo "  $0 channel_inspector test /path/to/test.vcf /path/to/test.cram"
            echo "  $0 meta_analyzer test /path/to/test.vcf"
            echo "  $0 performance_profiler test /path/to/test.vcf 100"
            echo "  $0 debug_toolkit test /path/to/test.vcf /path/to/test.cram"
            echo ""
            exit 0
            ;;
    esac
    
    # Generate report after any successful run
    if [[ $? -eq 0 && "$command" != "help" && "$command" != "report" ]]; then
        generate_debug_report
    fi
}

# Execute main function with all arguments
main "$@"