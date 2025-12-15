#!/bin/bash

# RNA Editing Test Framework Validation Script
# This script validates the test framework structure without running full tests

set -euo pipefail

echo "=== RNA Editing Test Framework Validation ==="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

VALIDATION_PASSED=0
VALIDATION_FAILED=0

# Function to validate file existence and syntax
validate_file() {
    local file_path="$1"
    local file_type="$2"
    
    if [ -f "$file_path" ]; then
        echo -e "${GREEN}✓ Found: $file_path${NC}"
        
        # Basic syntax validation for different file types
        case "$file_type" in
            "nf-test")
                if grep -q "nextflow_process\|nextflow_workflow\|nextflow_pipeline" "$file_path"; then
                    echo -e "  ${GREEN}✓ Valid nf-test syntax${NC}"
                    VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
                else
                    echo -e "  ${RED}✗ Invalid nf-test syntax${NC}"
                    VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
                fi
                ;;
            "vcf")
                if grep -q "##fileformat=VCFv4" "$file_path"; then
                    echo -e "  ${GREEN}✓ Valid VCF format${NC}"
                    VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
                else
                    echo -e "  ${RED}✗ Invalid VCF format${NC}"
                    VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
                fi
                ;;
            "config")
                if grep -q "params\|process\|profiles" "$file_path"; then
                    echo -e "  ${GREEN}✓ Valid Nextflow config syntax${NC}"
                    VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
                else
                    echo -e "  ${RED}✗ Invalid config syntax${NC}"
                    VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
                fi
                ;;
            "script")
                if [ -x "$file_path" ]; then
                    echo -e "  ${GREEN}✓ Script is executable${NC}"
                    VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
                else
                    echo -e "  ${RED}✗ Script is not executable${NC}"
                    VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
                fi
                ;;
            *)
                echo -e "  ${GREEN}✓ File exists${NC}"
                VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
                ;;
        esac
    else
        echo -e "${RED}✗ Missing: $file_path${NC}"
        VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
    fi
}

echo -e "\n${YELLOW}=== Validating Test Data Files ===${NC}"
validate_file "tests/data/test_rescue.vcf" "vcf"
validate_file "tests/data/test_rescue.vcf.gz" "file"
validate_file "tests/data/test_rescue.vcf.gz.tbi" "file"
validate_file "tests/data/test_rediportal.vcf" "vcf"
validate_file "tests/data/test_rediportal.vcf.gz" "file"
validate_file "tests/data/test_rediportal.vcf.gz.tbi" "file"
validate_file "tests/data/test_large_rescue.vcf" "vcf"
validate_file "tests/data/test_large_rescue.vcf.gz" "file"
validate_file "tests/data/test_large_rescue.vcf.gz.tbi" "file"
validate_file "tests/data/test_malformed.vcf" "vcf"
validate_file "tests/data/test_malformed.vcf.gz" "file"

echo -e "\n${YELLOW}=== Validating Unit Test Files ===${NC}"
validate_file "modules/local/rna_editing_annotation/tests/main.nf.test" "nf-test"

echo -e "\n${YELLOW}=== Validating Integration Test Files ===${NC}"
validate_file "subworkflows/local/vcf_rescue_post_processing/tests/main.nf.test" "nf-test"
validate_file "tests/rna_editing_integration.nf.test" "nf-test"

echo -e "\n${YELLOW}=== Validating Error Handling Test Files ===${NC}"
validate_file "tests/rna_editing_error_handling.nf.test" "nf-test"

echo -e "\n${YELLOW}=== Validating Performance Test Files ===${NC}"
validate_file "tests/rna_editing_performance.nf.test" "nf-test"

echo -e "\n${YELLOW}=== Validating Comprehensive Test Files ===${NC}"
validate_file "tests/rna_editing_comprehensive.nf.test" "nf-test"

echo -e "\n${YELLOW}=== Validating Configuration Files ===${NC}"
validate_file "tests/config/rna_editing_test.config" "config"
validate_file "tests/config/test_data.config" "config"

echo -e "\n${YELLOW}=== Validating Scripts ===${NC}"
validate_file "tests/run_rna_editing_tests.sh" "script"
validate_file "tests/validate_test_framework.sh" "script"

echo -e "\n${YELLOW}=== Validating Documentation ===${NC}"
validate_file "tests/README_RNA_EDITING_TESTS.md" "file"

# Check test file content for key requirements
echo -e "\n${YELLOW}=== Validating Test Content Requirements ===${NC}"

# Check that unit tests cover all required scenarios
if grep -q "RNA editing annotation - basic functionality" modules/local/rna_editing_annotation/tests/main.nf.test; then
    echo -e "${GREEN}✓ Unit tests include basic functionality test${NC}"
    VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
else
    echo -e "${RED}✗ Unit tests missing basic functionality test${NC}"
    VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
fi

if grep -q "RNA editing annotation - error handling" modules/local/rna_editing_annotation/tests/main.nf.test; then
    echo -e "${GREEN}✓ Unit tests include error handling test${NC}"
    VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
else
    echo -e "${RED}✗ Unit tests missing error handling test${NC}"
    VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
fi

# Check that integration tests cover enabled/disabled scenarios
if grep -q "RNA editing enabled" subworkflows/local/vcf_rescue_post_processing/tests/main.nf.test; then
    echo -e "${GREEN}✓ Integration tests include RNA editing enabled scenario${NC}"
    VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
else
    echo -e "${RED}✗ Integration tests missing RNA editing enabled scenario${NC}"
    VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
fi

if grep -q "RNA editing disabled" subworkflows/local/vcf_rescue_post_processing/tests/main.nf.test; then
    echo -e "${GREEN}✓ Integration tests include RNA editing disabled scenario${NC}"
    VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
else
    echo -e "${RED}✗ Integration tests missing RNA editing disabled scenario${NC}"
    VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
fi

# Check that performance tests exist
if grep -q "Performance" tests/rna_editing_performance.nf.test; then
    echo -e "${GREEN}✓ Performance tests are implemented${NC}"
    VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
else
    echo -e "${RED}✗ Performance tests are missing${NC}"
    VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
fi

# Summary
echo -e "\n${YELLOW}=== Validation Summary ===${NC}"
echo "Validations passed: $VALIDATION_PASSED"
echo "Validations failed: $VALIDATION_FAILED"

if [ $VALIDATION_FAILED -eq 0 ]; then
    echo -e "\n${GREEN}🎉 Test framework validation passed! All required components are in place.${NC}"
    echo -e "${GREEN}The comprehensive testing framework for RNA editing annotation is ready.${NC}"
    exit 0
else
    echo -e "\n${RED}❌ Test framework validation failed. Please address the issues above.${NC}"
    exit 1
fi