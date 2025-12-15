#!/bin/bash

# Final Integration and Compatibility Test Suite Runner
# This script runs comprehensive final integration tests for RNA editing annotation

set -euo pipefail

echo "=== Final Integration and Compatibility Test Suite ==="
echo "Starting comprehensive final integration testing..."

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test results tracking
PASSED_TESTS=0
FAILED_TESTS=0
TOTAL_TESTS=0
TEST_OUTPUT_DIR="test_results_$(date +%Y%m%d_%H%M%S)"

# Create test output directory
mkdir -p "$TEST_OUTPUT_DIR"

# Function to run a test and track results
run_test() {
    local test_name="$1"
    local test_file="$2"
    local test_category="${3:-integration}"
    
    echo -e "\n${BLUE}=== $test_category Test ===${NC}"
    echo -e "${YELLOW}Running: $test_name${NC}"
    echo "Test file: $test_file"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    local log_file="$TEST_OUTPUT_DIR/$(basename "$test_file" .nf.test)_$(date +%H%M%S).log"
    
    if nf-test test "$test_file" --verbose > "$log_file" 2>&1; then
        echo -e "${GREEN}✓ PASSED: $test_name${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        echo "Log: $log_file"
    else
        echo -e "${RED}✗ FAILED: $test_name${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        echo "Error log: $log_file"
        echo -e "${RED}Last 20 lines of error log:${NC}"
        tail -20 "$log_file" || true
    fi
}

# Function to run downstream compatibility validation
run_downstream_validation() {
    local output_dir="$1"
    local test_name="Downstream Tool Compatibility"
    
    echo -e "\n${BLUE}=== Downstream Compatibility Test ===${NC}"
    echo -e "${YELLOW}Running: $test_name${NC}"
    echo "Output directory: $output_dir"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    local log_file="$TEST_OUTPUT_DIR/downstream_compatibility_$(date +%H%M%S).log"
    
    if ./tests/validate_downstream_compatibility.sh "$output_dir" > "$log_file" 2>&1; then
        echo -e "${GREEN}✓ PASSED: $test_name${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        echo -e "${RED}✗ FAILED: $test_name${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        echo "Error log: $log_file"
        echo -e "${RED}Last 10 lines of error log:${NC}"
        tail -10 "$log_file" || true
    fi
}

# Function to validate test data exists
validate_test_data() {
    echo -e "\n${BLUE}=== Validating Test Data ===${NC}"
    
    local required_files=(
        "tests/data/test_rescue.vcf.gz"
        "tests/data/test_rediportal.vcf.gz"
        "tests/data/test_large_rescue.vcf.gz"
        "tests/data/test_malformed.vcf.gz"
        "tests/csv/3.0/fastq_pair.csv"
    )
    
    local missing_files=()
    
    for file in "${required_files[@]}"; do
        if [ ! -f "$file" ]; then
            missing_files+=("$file")
        else
            echo -e "${GREEN}✓ Found: $file${NC}"
        fi
    done
    
    if [ ${#missing_files[@]} -gt 0 ]; then
        echo -e "${RED}Missing required test files:${NC}"
        for file in "${missing_files[@]}"; do
            echo -e "${RED}  - $file${NC}"
        done
        return 1
    fi
    
    echo -e "${GREEN}All required test data files found${NC}"
    return 0
}

# Function to check required tools
check_required_tools() {
    echo -e "\n${BLUE}=== Checking Required Tools ===${NC}"
    
    local required_tools=("nf-test" "nextflow")
    local optional_tools=("bcftools" "tabix" "vcf-validator")
    
    for tool in "${required_tools[@]}"; do
        if command -v "$tool" &> /dev/null; then
            echo -e "${GREEN}✓ Found required tool: $tool${NC}"
        else
            echo -e "${RED}✗ Missing required tool: $tool${NC}"
            return 1
        fi
    done
    
    for tool in "${optional_tools[@]}"; do
        if command -v "$tool" &> /dev/null; then
            echo -e "${GREEN}✓ Found optional tool: $tool${NC}"
        else
            echo -e "${YELLOW}⚠ Optional tool not found: $tool${NC}"
        fi
    done
    
    return 0
}

# Function to generate test report
generate_test_report() {
    local report_file="$TEST_OUTPUT_DIR/final_integration_test_report.md"
    
    cat > "$report_file" << EOF
# Final Integration Test Report

**Date:** $(date)
**Test Suite:** RNA Editing Nextflow Integration - Final Integration and Compatibility

## Summary

- **Total Tests:** $TOTAL_TESTS
- **Passed:** $PASSED_TESTS
- **Failed:** $FAILED_TESTS
- **Success Rate:** $(( PASSED_TESTS * 100 / TOTAL_TESTS ))%

## Test Categories

### 1. Final Integration Tests
- Complete workflow with production-like data
- Backward compatibility validation
- Existing workflow configuration compatibility
- Downstream tool compatibility validation
- Error handling with malformed inputs
- Resource management and scalability
- Configuration parameter validation

### 2. Production Compatibility Tests
- VCF format validation
- Annotation field preservation
- Large dataset simulation
- Multi-sample processing
- Workflow reporting and monitoring
- Configuration profiles

### 3. Downstream Tool Compatibility
- bcftools compatibility
- tabix compatibility
- VCF format compliance
- RNA editing annotation preservation

## Requirements Validation

This test suite validates the following requirements:

- **Requirement 2.4:** Workflow integration preservation and backward compatibility
- **Requirement 3.5:** VCF rescue filtering compatibility with annotated VCF
- **Requirement 6.5:** Configuration parameter integration and validation
- **Requirement 7.5:** Error handling and resource management

## Test Results

EOF

    if [ $FAILED_TESTS -eq 0 ]; then
        cat >> "$report_file" << EOF
### ✅ All Tests Passed

The RNA editing annotation integration is ready for production use. All compatibility tests passed successfully.

EOF
    else
        cat >> "$report_file" << EOF
### ❌ Some Tests Failed

$FAILED_TESTS out of $TOTAL_TESTS tests failed. Please review the individual test logs for details.

**Failed Test Logs:**
EOF
        find "$TEST_OUTPUT_DIR" -name "*.log" -exec grep -l "FAILED\|ERROR" {} \; | while read -r log_file; do
            echo "- $(basename "$log_file")" >> "$report_file"
        done
    fi
    
    echo -e "\n${BLUE}Test report generated: $report_file${NC}"
}

# Main execution
main() {
    echo "Final Integration Test Suite for RNA Editing Annotation"
    echo "Test output directory: $TEST_OUTPUT_DIR"
    
    # Pre-flight checks
    if ! check_required_tools; then
        echo -e "${RED}Required tools missing. Please install missing tools and try again.${NC}"
        exit 1
    fi
    
    if ! validate_test_data; then
        echo -e "${RED}Required test data missing. Please ensure test data files exist.${NC}"
        exit 1
    fi
    
    echo -e "\n${BLUE}=== Phase 1: Final Integration Tests ===${NC}"
    
    # Run final integration tests
    run_test "Final Integration - Complete Workflow" "tests/rna_editing_final_integration.nf.test" "Final Integration"
    
    echo -e "\n${BLUE}=== Phase 2: Production Compatibility Tests ===${NC}"
    
    # Run production compatibility tests
    run_test "Production Compatibility - VCF Format and Downstream Tools" "tests/rna_editing_production_compatibility.nf.test" "Production Compatibility"
    
    echo -e "\n${BLUE}=== Phase 3: Comprehensive Validation ===${NC}"
    
    # Run existing comprehensive tests to ensure nothing broke
    run_test "Comprehensive RNA Editing Validation" "tests/rna_editing_comprehensive.nf.test" "Comprehensive"
    
    # Run integration tests to ensure backward compatibility
    run_test "RNA Editing Integration Validation" "tests/rna_editing_integration.nf.test" "Integration"
    
    echo -e "\n${BLUE}=== Phase 4: Error Handling Validation ===${NC}"
    
    # Run error handling tests
    run_test "Error Handling Validation" "tests/rna_editing_error_handling.nf.test" "Error Handling"
    
    echo -e "\n${BLUE}=== Phase 5: Downstream Tool Compatibility ===${NC}"
    
    # Look for workflow output directories to test downstream compatibility
    local output_dirs
    output_dirs=$(find . -name "output*" -type d 2>/dev/null | head -3 || true)
    
    if [ -n "$output_dirs" ]; then
        while IFS= read -r output_dir; do
            if [ -d "$output_dir/variant_calling" ] || [ -d "$output_dir" ]; then
                run_downstream_validation "$output_dir"
                break
            fi
        done <<< "$output_dirs"
    else
        echo -e "${YELLOW}No workflow output directories found for downstream compatibility testing${NC}"
        echo -e "${YELLOW}This is expected if this is the first test run${NC}"
    fi
    
    # Generate comprehensive test report
    generate_test_report
    
    # Final summary
    echo -e "\n${BLUE}=== Final Integration Test Suite Summary ===${NC}"
    echo "Total tests run: $TOTAL_TESTS"
    echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
    echo -e "Failed: ${RED}$FAILED_TESTS${NC}"
    
    if [ $FAILED_TESTS -eq 0 ]; then
        echo -e "\n${GREEN}🎉 All final integration tests passed!${NC}"
        echo -e "${GREEN}RNA editing annotation is ready for production deployment.${NC}"
        echo -e "\n${BLUE}✅ Requirements Validated:${NC}"
        echo -e "  - Requirement 2.4: Workflow integration and backward compatibility"
        echo -e "  - Requirement 3.5: VCF filtering compatibility with annotations"
        echo -e "  - Requirement 6.5: Configuration parameter integration"
        echo -e "  - Requirement 7.5: Error handling and resource management"
        exit 0
    else
        echo -e "\n${RED}❌ Some final integration tests failed.${NC}"
        echo -e "${RED}Please review the test logs in $TEST_OUTPUT_DIR${NC}"
        echo -e "\n${YELLOW}Common issues to check:${NC}"
        echo -e "  - Test data file availability"
        echo -e "  - Resource allocation (memory/CPU limits)"
        echo -e "  - Network connectivity for remote test data"
        echo -e "  - Tool dependencies (bcftools, tabix, etc.)"
        exit 1
    fi
}

# Cleanup function
cleanup() {
    echo -e "\n${YELLOW}Cleaning up temporary files...${NC}"
    # Add any cleanup logic here if needed
}

# Set up cleanup trap
trap cleanup EXIT

# Run main function
main "$@"