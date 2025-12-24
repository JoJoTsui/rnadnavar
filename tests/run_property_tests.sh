#!/bin/bash

#
# Property-Based Test Runner for VCF Realignment Optimization
# Runs all property tests with minimum 100 iterations each
#

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
TEST_CONFIG="$SCRIPT_DIR/config/property_test.config"
RESULTS_DIR="$SCRIPT_DIR/results/property_tests"
LOG_FILE="$RESULTS_DIR/property_test_run.log"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

# Setup
setup_test_environment() {
    log "Setting up property test environment..."
    
    # Create results directory
    mkdir -p "$RESULTS_DIR"
    
    # Initialize log file
    echo "Property-Based Test Run - $(date)" > "$LOG_FILE"
    echo "========================================" >> "$LOG_FILE"
    
    # Check dependencies
    if ! command -v nextflow &> /dev/null; then
        error "Nextflow is not installed or not in PATH"
        exit 1
    fi
    
    if ! command -v nf-test &> /dev/null; then
        error "nf-test is not installed or not in PATH"
        exit 1
    fi
    
    log "Environment setup complete"
}

# Run individual property test
run_property_test() {
    local test_file="$1"
    local test_name="$(basename "$test_file" .nf.test)"
    
    log "Running property test: $test_name"
    
    local test_output_dir="$RESULTS_DIR/$test_name"
    mkdir -p "$test_output_dir"
    
    # Run the test with nf-test
    if nf-test test "$test_file" \
        --config "$TEST_CONFIG" \
        --outdir "$test_output_dir" \
        --verbose \
        >> "$LOG_FILE" 2>&1; then
        success "Property test $test_name completed successfully"
        return 0
    else
        error "Property test $test_name failed"
        return 1
    fi
}

# Run all property tests
run_all_property_tests() {
    log "Starting property-based test suite..."
    
    local test_files=(
        "$SCRIPT_DIR/property_tests/stackoverflow_prevention.nf.test"
        "$SCRIPT_DIR/property_tests/null_condition_handling.nf.test"
        "$SCRIPT_DIR/property_tests/channel_join_integrity.nf.test"
        "$SCRIPT_DIR/property_tests/file_path_validation.nf.test"
        "$SCRIPT_DIR/property_tests/command_construction.nf.test"
    )
    
    local total_tests=${#test_files[@]}
    local passed_tests=0
    local failed_tests=0
    
    log "Found $total_tests property test files"
    
    for test_file in "${test_files[@]}"; do
        if [[ -f "$test_file" ]]; then
            if run_property_test "$test_file"; then
                ((passed_tests++))
            else
                ((failed_tests++))
            fi
        else
            warning "Test file not found: $test_file"
            ((failed_tests++))
        fi
    done
    
    # Summary
    echo ""
    log "Property Test Suite Summary:"
    log "  Total tests: $total_tests"
    log "  Passed: $passed_tests"
    log "  Failed: $failed_tests"
    log "  Success rate: $(( passed_tests * 100 / total_tests ))%"
    
    if [[ $failed_tests -eq 0 ]]; then
        success "All property tests passed!"
        return 0
    else
        error "$failed_tests property test(s) failed"
        return 1
    fi
}

# Run comprehensive property test suite
run_comprehensive_suite() {
    log "Running comprehensive property test suite with framework..."
    
    local suite_output_dir="$RESULTS_DIR/comprehensive_suite"
    mkdir -p "$suite_output_dir"
    
    # Run the main property test suite workflow
    if nf-test test "$SCRIPT_DIR/property_test_framework.nf" \
        --config "$TEST_CONFIG" \
        --outdir "$suite_output_dir" \
        --verbose \
        --entry PROPERTY_TEST_SUITE \
        >> "$LOG_FILE" 2>&1; then
        success "Comprehensive property test suite completed successfully"
        return 0
    else
        error "Comprehensive property test suite failed"
        return 1
    fi
}

# Generate test report
generate_report() {
    log "Generating property test report..."
    
    local report_file="$RESULTS_DIR/property_test_report.md"
    
    cat > "$report_file" << EOF
# Property-Based Test Report

**Generated:** $(date)
**Test Configuration:** $TEST_CONFIG
**Results Directory:** $RESULTS_DIR

## Test Summary

This report contains the results of property-based testing for the VCF Realignment Optimization feature.

### Properties Tested

1. **StackOverflow Prevention** - Validates that channel operations complete without StackOverflowError
2. **Null Condition Handling** - Validates graceful handling of null task.ext.when conditions
3. **Channel Join Integrity** - Validates patient ID matching and metadata preservation
4. **File Path Validation** - Validates file existence and accessibility checks
5. **Command Construction Validation** - Validates command syntax and parameter correctness

### Test Configuration

- **Iterations per property:** 100 minimum
- **Test data seed:** 42 (reproducible)
- **Maximum test samples:** 10
- **Maximum test patients:** 5
- **Timeout:** 30 minutes per test

### Results

See individual test logs in the results directory for detailed information.

### Requirements Validation

- **Requirements 1.1, 1.2, 1.3:** StackOverflow Prevention
- **Requirements 1.5:** Null Condition Handling
- **Requirements 3.1, 3.2:** Channel Join Integrity
- **Requirements 3.3:** File Path Validation
- **Requirements 5.5:** Command Construction Validation

EOF

    log "Report generated: $report_file"
}

# Main execution
main() {
    log "Starting VCF Realignment Optimization Property-Based Test Suite"
    
    setup_test_environment
    
    local exit_code=0
    
    # Run individual property tests
    if ! run_all_property_tests; then
        exit_code=1
    fi
    
    # Run comprehensive suite
    if ! run_comprehensive_suite; then
        exit_code=1
    fi
    
    # Generate report
    generate_report
    
    if [[ $exit_code -eq 0 ]]; then
        success "Property-based test suite completed successfully"
    else
        error "Property-based test suite completed with failures"
    fi
    
    log "Test results available in: $RESULTS_DIR"
    log "Test log available in: $LOG_FILE"
    
    exit $exit_code
}

# Handle command line arguments
case "${1:-}" in
    "setup")
        setup_test_environment
        ;;
    "run")
        main
        ;;
    "report")
        generate_report
        ;;
    *)
        echo "Usage: $0 {setup|run|report}"
        echo "  setup  - Setup test environment only"
        echo "  run    - Run complete property test suite"
        echo "  report - Generate test report only"
        exit 1
        ;;
esac