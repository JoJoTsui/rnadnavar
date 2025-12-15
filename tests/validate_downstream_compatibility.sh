#!/bin/bash

# Downstream Tool Compatibility Validation Script
# This script validates that RNA editing annotated VCF files are compatible with common downstream tools

set -euo pipefail

echo "=== Downstream Tool Compatibility Validation ==="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test results tracking
PASSED_TESTS=0
FAILED_TESTS=0
TOTAL_TESTS=0

# Function to run a compatibility test
run_compatibility_test() {
    local test_name="$1"
    local test_command="$2"
    local vcf_file="$3"
    
    echo -e "\n${YELLOW}Testing: $test_name${NC}"
    echo "Command: $test_command"
    echo "VCF file: $vcf_file"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    if [ ! -f "$vcf_file" ]; then
        echo -e "${RED}✗ FAILED: VCF file not found: $vcf_file${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    if eval "$test_command" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASSED: $test_name${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    else
        echo -e "${RED}✗ FAILED: $test_name${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
}

# Function to check if a tool is available
check_tool() {
    local tool="$1"
    if command -v "$tool" &> /dev/null; then
        echo -e "${GREEN}✓ Found: $tool${NC}"
        return 0
    else
        echo -e "${YELLOW}⚠ Not found: $tool (skipping related tests)${NC}"
        return 1
    fi
}

# Function to find VCF files in output directory
find_vcf_files() {
    local output_dir="$1"
    find "$output_dir" -name "*.vcf.gz" -o -name "*.vcf" 2>/dev/null || true
}

# Main validation function
validate_downstream_compatibility() {
    local output_dir="$1"
    
    echo "Searching for VCF files in: $output_dir"
    
    # Find all VCF files
    local vcf_files
    vcf_files=$(find_vcf_files "$output_dir")
    
    if [ -z "$vcf_files" ]; then
        echo -e "${RED}No VCF files found in output directory${NC}"
        return 1
    fi
    
    echo "Found VCF files:"
    echo "$vcf_files"
    
    # Test each VCF file with available tools
    while IFS= read -r vcf_file; do
        echo -e "\n${YELLOW}=== Testing VCF file: $(basename "$vcf_file") ===${NC}"
        
        # Test 1: bcftools compatibility
        if check_tool "bcftools"; then
            run_compatibility_test "bcftools view" "bcftools view -h '$vcf_file'" "$vcf_file"
            run_compatibility_test "bcftools stats" "bcftools stats '$vcf_file'" "$vcf_file"
            run_compatibility_test "bcftools query" "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' '$vcf_file' | head -5" "$vcf_file"
        fi
        
        # Test 2: tabix compatibility (for compressed files)
        if [[ "$vcf_file" == *.gz ]] && check_tool "tabix"; then
            run_compatibility_test "tabix query" "tabix -h '$vcf_file' chr1:1-1000000 || tabix -h '$vcf_file' 1:1-1000000" "$vcf_file"
        fi
        
        # Test 3: VCF format validation with vcf-validator (if available)
        if check_tool "vcf-validator"; then
            run_compatibility_test "vcf-validator" "vcf-validator '$vcf_file'" "$vcf_file"
        fi
        
        # Test 4: Basic file format checks
        run_compatibility_test "VCF header check" "head -100 '$vcf_file' | grep -q '^##fileformat=VCF' || zhead -100 '$vcf_file' | grep -q '^##fileformat=VCF'" "$vcf_file"
        
        # Test 5: Check for RNA editing annotations (if this is an annotated file)
        if [[ "$(basename "$vcf_file")" == *"rescue"* ]] || [[ "$(basename "$vcf_file")" == *"rna"* ]]; then
            echo -e "\n${YELLOW}Checking for RNA editing annotations...${NC}"
            
            # Check for REDI_* INFO fields in header
            if [[ "$vcf_file" == *.gz ]]; then
                if zgrep -q "##INFO=.*REDI_" "$vcf_file" 2>/dev/null; then
                    echo -e "${GREEN}✓ Found RNA editing INFO fields in header${NC}"
                    PASSED_TESTS=$((PASSED_TESTS + 1))
                else
                    echo -e "${YELLOW}⚠ No RNA editing INFO fields found (may be expected)${NC}"
                fi
                TOTAL_TESTS=$((TOTAL_TESTS + 1))
                
                # Check for RNAedit FILTER
                if zgrep -q "##FILTER=.*RNAedit" "$vcf_file" 2>/dev/null; then
                    echo -e "${GREEN}✓ Found RNAedit FILTER in header${NC}"
                    PASSED_TESTS=$((PASSED_TESTS + 1))
                else
                    echo -e "${YELLOW}⚠ No RNAedit FILTER found (may be expected)${NC}"
                fi
                TOTAL_TESTS=$((TOTAL_TESTS + 1))
            else
                if grep -q "##INFO=.*REDI_" "$vcf_file" 2>/dev/null; then
                    echo -e "${GREEN}✓ Found RNA editing INFO fields in header${NC}"
                    PASSED_TESTS=$((PASSED_TESTS + 1))
                else
                    echo -e "${YELLOW}⚠ No RNA editing INFO fields found (may be expected)${NC}"
                fi
                TOTAL_TESTS=$((TOTAL_TESTS + 1))
                
                if grep -q "##FILTER=.*RNAedit" "$vcf_file" 2>/dev/null; then
                    echo -e "${GREEN}✓ Found RNAedit FILTER in header${NC}"
                    PASSED_TESTS=$((PASSED_TESTS + 1))
                else
                    echo -e "${YELLOW}⚠ No RNAedit FILTER found (may be expected)${NC}"
                fi
                TOTAL_TESTS=$((TOTAL_TESTS + 1))
            fi
        fi
        
        # Test 6: Check file integrity
        if [[ "$vcf_file" == *.gz ]]; then
            run_compatibility_test "gzip integrity" "gzip -t '$vcf_file'" "$vcf_file"
        fi
        
    done <<< "$vcf_files"
}

# Function to test with sample data
test_with_sample_data() {
    echo -e "\n${YELLOW}=== Testing with sample data ===${NC}"
    
    local test_data_dir="tests/data"
    
    if [ -d "$test_data_dir" ]; then
        local sample_vcfs
        sample_vcfs=$(find "$test_data_dir" -name "*.vcf.gz" -o -name "*.vcf" 2>/dev/null || true)
        
        if [ -n "$sample_vcfs" ]; then
            echo "Testing sample VCF files for baseline compatibility:"
            while IFS= read -r vcf_file; do
                echo -e "\n${YELLOW}Testing sample file: $(basename "$vcf_file")${NC}"
                
                if check_tool "bcftools"; then
                    run_compatibility_test "Sample bcftools view" "bcftools view -h '$vcf_file'" "$vcf_file"
                fi
                
                run_compatibility_test "Sample VCF header" "head -10 '$vcf_file' | grep -q '^##fileformat=VCF' || zhead -10 '$vcf_file' | grep -q '^##fileformat=VCF'" "$vcf_file"
                
            done <<< "$sample_vcfs"
        fi
    fi
}

# Main execution
main() {
    local output_dir="${1:-}"
    
    if [ -z "$output_dir" ]; then
        echo "Usage: $0 <output_directory>"
        echo "Example: $0 results/variant_calling"
        exit 1
    fi
    
    if [ ! -d "$output_dir" ]; then
        echo -e "${RED}Output directory not found: $output_dir${NC}"
        exit 1
    fi
    
    echo "Validating downstream tool compatibility for RNA editing annotated VCF files"
    echo "Output directory: $output_dir"
    
    # Test with sample data first (baseline)
    test_with_sample_data
    
    # Test with workflow output
    validate_downstream_compatibility "$output_dir"
    
    # Summary
    echo -e "\n${YELLOW}=== Compatibility Test Summary ===${NC}"
    echo "Total tests run: $TOTAL_TESTS"
    echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
    echo -e "Failed: ${RED}$FAILED_TESTS${NC}"
    
    if [ $FAILED_TESTS -eq 0 ]; then
        echo -e "\n${GREEN}🎉 All compatibility tests passed! VCF files are compatible with downstream tools.${NC}"
        exit 0
    else
        echo -e "\n${RED}❌ Some compatibility tests failed. Please review the failures above.${NC}"
        exit 1
    fi
}

# Helper function for zhead (head for gzipped files)
zhead() {
    local lines="$1"
    local file="$2"
    zcat "$file" | head -n "$lines"
}

# Run main function with all arguments
main "$@"