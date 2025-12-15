#!/bin/bash

# RNA Editing Annotation Test Suite Runner
# This script runs all RNA editing related tests in the correct order

set -euo pipefail

echo "=== RNA Editing Annotation Test Suite ==="
echo "Starting comprehensive test execution..."

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test results tracking
PASSED_TESTS=0
FAILED_TESTS=0
TOTAL_TESTS=0

# Function to run a test and track results
run_test() {
    local test_name="$1"
    local test_file="$2"
    
    echo -e "\n${YELLOW}Running: $test_name${NC}"
    echo "Test file: $test_file"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    if nf-test test "$test_file" --verbose; then
        echo -e "${GREEN}✓ PASSED: $test_name${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        echo -e "${RED}✗ FAILED: $test_name${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
}

# Function to run performance tests with timing
run_performance_test() {
    local test_name="$1"
    local test_file="$2"
    
    echo -e "\n${YELLOW}Running Performance Test: $test_name${NC}"
    echo "Test file: $test_file"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    start_time=$(date +%s)
    if nf-test test "$test_file" --verbose; then
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        echo -e "${GREEN}✓ PASSED: $test_name (${duration}s)${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        echo -e "${RED}✗ FAILED: $test_name (${duration}s)${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
}

echo -e "\n${YELLOW}=== Phase 1: Unit Tests ===${NC}"

# Test RNA editing annotation module
run_test "RNA Editing Annotation Module" "modules/local/rna_editing_annotation/tests/main.nf.test"

echo -e "\n${YELLOW}=== Phase 2: Integration Tests ===${NC}"

# Test VCF rescue post-processing subworkflow
run_test "VCF Rescue Post-Processing Subworkflow" "subworkflows/local/vcf_rescue_post_processing/tests/main.nf.test"

# Test full pipeline integration
run_test "RNA Editing Pipeline Integration" "tests/rna_editing_integration.nf.test"

echo -e "\n${YELLOW}=== Phase 3: Error Handling Tests ===${NC}"

# Test error handling scenarios
run_test "RNA Editing Error Handling" "tests/rna_editing_error_handling.nf.test"

echo -e "\n${YELLOW}=== Phase 4: Performance Tests ===${NC}"

# Test performance with different VCF sizes
run_performance_test "RNA Editing Performance" "tests/rna_editing_performance.nf.test"

echo -e "\n${YELLOW}=== Phase 5: Comprehensive Tests ===${NC}"

# Run comprehensive validation tests
run_test "RNA Editing Comprehensive Validation" "tests/rna_editing_comprehensive.nf.test"

# Summary
echo -e "\n${YELLOW}=== Test Suite Summary ===${NC}"
echo "Total tests run: $TOTAL_TESTS"
echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
echo -e "Failed: ${RED}$FAILED_TESTS${NC}"

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "\n${GREEN}🎉 All tests passed! RNA editing annotation is ready for production.${NC}"
    exit 0
else
    echo -e "\n${RED}❌ Some tests failed. Please review the failures above.${NC}"
    exit 1
fi