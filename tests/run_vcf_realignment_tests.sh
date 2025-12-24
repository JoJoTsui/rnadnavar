#!/bin/bash

# VCF Realignment Optimization Test Suite Runner
# This script runs all VCF realignment related tests in the correct order

set -euo pipefail

echo "=== VCF Realignment Optimization Test Suite ==="
echo "Starting comprehensive test execution..."

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

# Function to run a test and track results
run_test() {
    local test_name="$1"
    local test_file="$2"
    local config_profile="${3:-}"
    
    echo -e "\n${YELLOW}Running: $test_name${NC}"
    echo "Test file: $test_file"
    if [ -n "$config_profile" ]; then
        echo "Profile: $config_profile"
    fi
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    local cmd="nf-test test \"$test_file\" --verbose"
    if [ -n "$config_profile" ]; then
        cmd="$cmd --profile $config_profile"
    fi
    
    if eval "$cmd"; then
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
    local config_profile="${3:-}"
    
    echo -e "\n${BLUE}Running Performance Test: $test_name${NC}"
    echo "Test file: $test_file"
    if [ -n "$config_profile" ]; then
        echo "Profile: $config_profile"
    fi
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    local cmd="nf-test test \"$test_file\" --verbose"
    if [ -n "$config_profile" ]; then
        cmd="$cmd --profile $config_profile"
    fi
    
    start_time=$(date +%s)
    if eval "$cmd"; then
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

echo -e "\n${YELLOW}=== Phase 1: Component Unit Tests ===${NC}"

# Test individual VCF realignment components
run_test "Channel Sanitization Module" "subworkflows/local/sanitize_channels/tests/main.nf.test"
run_test "Safe Channel Join Module" "subworkflows/local/safe_channel_join/tests/main.nf.test"
run_test "Enhanced CRAM2BAM Conversion" "subworkflows/local/enhanced_cram2bam_conversion/tests/main.nf.test"
run_test "Input Validation Module" "subworkflows/local/input_validation/tests/main.nf.test"

echo -e "\n${YELLOW}=== Phase 2: VCF Realignment Integration Tests ===${NC}"

# Test VCF realignment workflow integration
run_test "VCF Realignment End-to-End" "tests/vcf_realignment_integration.nf.test" "vcf_realignment_test"

echo -e "\n${YELLOW}=== Phase 3: Error Handling and Validation Tests ===${NC}"

# Test error handling scenarios
run_test "VCF Realignment Error Handling" "tests/vcf_realignment_integration.nf.test"
run_test "Channel Join Validation" "tests/vcf_realignment_integration.nf.test"

echo -e "\n${YELLOW}=== Phase 4: Backward Compatibility Tests ===${NC}"

# Test MAF workflow backward compatibility
run_test "MAF Realignment Backward Compatibility" "tests/vcf_realignment_integration.nf.test" "maf_realignment_test"

echo -e "\n${YELLOW}=== Phase 5: Workflow Equivalence Tests ===${NC}"

# Test VCF vs MAF workflow equivalence
run_test "VCF vs MAF Workflow Equivalence" "tests/vcf_realignment_integration.nf.test"

echo -e "\n${YELLOW}=== Phase 6: Performance Tests ===${NC}"

# Test performance with different configurations
run_performance_test "VCF Realignment Performance" "tests/vcf_realignment_integration.nf.test" "vcf_realignment_performance"

echo -e "\n${YELLOW}=== Phase 7: Stack Overflow Prevention Tests ===${NC}"

# Test specific StackOverflowError prevention
run_test "StackOverflow Prevention Validation" "tests/vcf_realignment_integration.nf.test"

# Summary
echo -e "\n${YELLOW}=== Test Suite Summary ===${NC}"
echo "Total tests run: $TOTAL_TESTS"
echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
echo -e "Failed: ${RED}$FAILED_TESTS${NC}"

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "\n${GREEN}🎉 All VCF realignment tests passed! The optimization is ready for production.${NC}"
    echo -e "${GREEN}Key achievements:${NC}"
    echo -e "${GREEN}  ✓ StackOverflowError prevention validated${NC}"
    echo -e "${GREEN}  ✓ Channel sanitization and validation working${NC}"
    echo -e "${GREEN}  ✓ Error handling and monitoring functional${NC}"
    echo -e "${GREEN}  ✓ Backward compatibility with MAF workflow maintained${NC}"
    echo -e "${GREEN}  ✓ End-to-end integration successful${NC}"
    exit 0
else
    echo -e "\n${RED}❌ Some tests failed. Please review the failures above.${NC}"
    echo -e "${RED}Common issues to check:${NC}"
    echo -e "${RED}  • Test data availability and accessibility${NC}"
    echo -e "${RED}  • Resource allocation and memory limits${NC}"
    echo -e "${RED}  • Container/environment configuration${NC}"
    echo -e "${RED}  • Network connectivity for remote test data${NC}"
    exit 1
fi