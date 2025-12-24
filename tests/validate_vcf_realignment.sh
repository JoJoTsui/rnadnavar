#!/bin/bash
# VCF Realignment Validation Script
#
# This script performs automated validation checks for the VCF-based RNA realignment feature.
# It verifies that all required modules, subworkflows, and configurations are in place.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

echo "=========================================="
echo "VCF Realignment Validation"
echo "=========================================="
echo ""

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Track test results
PASS=0
FAIL=0

check_pass() {
    echo -e "${GREEN}   ✓${NC} $1"
    ((PASS++))
}

check_fail() {
    echo -e "${RED}   ✗${NC} $1"
    ((FAIL++))
}

check_warn() {
    echo -e "${YELLOW}   ⚠${NC} $1"
}

echo "1. Checking VCF2BED module..."
if [ -f "modules/local/vcf2bed/main.nf" ]; then
    check_pass "VCF2BED module exists"
else
    check_fail "VCF2BED module missing"
fi

echo ""
echo "2. Checking subworkflows..."
for wf in prepare_realignment_vcf rna_realignment second_rescue; do
    if [ -f "subworkflows/local/${wf}/main.nf" ]; then
        check_pass "${wf} exists"
    else
        check_fail "${wf} missing"
    fi
done

echo ""
echo "3. Checking VCF2BED module config..."
if [ -f "conf/modules/prepare_realignment/vcf2bed.config" ]; then
    check_pass "vcf2bed.config exists"
else
    check_fail "vcf2bed.config missing"
fi

echo ""
echo "4. Checking nf-core module configs..."
for cfg in picard/filtersamreads.config samtools/convert.config; do
    if [ -f "conf/modules/nf-core/${cfg}" ]; then
        check_pass "${cfg} exists"
    else
        check_warn "${cfg} not found (may be optional)"
    fi
done

echo ""
echo "5. Checking nextflow.config configuration..."
if grep -q "realignment_mode = 'vcf'" nextflow.config; then
    check_pass "VCF mode is set as default"
else
    check_fail "VCF mode not set as default in nextflow.config"
fi

echo ""
echo "6. Checking workflow imports..."
imports=(
    "PREPARE_REALIGNMENT_VCF"
    "RNA_REALIGNMENT_WORKFLOW"
    "SECOND_RESCUE_WORKFLOW"
)

for import_name in "${imports[@]}"; do
    if grep -q "include.*${import_name}" workflows/rnadnavar.nf; then
        check_pass "${import_name} import found"
    else
        check_fail "${import_name} import missing"
    fi
done

echo ""
echo "7. Checking nextflow_schema.json..."
if [ -f "nextflow_schema.json" ]; then
    if grep -q '"realignment_mode"' nextflow_schema.json; then
        check_pass "realignment_mode in schema"
    else
        check_fail "realignment_mode not in schema"
    fi

    if grep -q '"enum".*\[.*"maf".*"vcf"' nextflow_schema.json; then
        check_pass "enum values (maf, vcf) defined"
    else
        check_fail "enum values not properly defined"
    fi
else
    check_fail "nextflow_schema.json not found"
fi

echo ""
echo "8. Checking realignment branch in main workflow..."
if grep -q "mode == 'vcf'" workflows/rnadnavar.nf; then
    check_pass "VCF mode branch logic found"
else
    check_fail "VCF mode branch logic missing"
fi

if grep -q "mode == 'maf'" workflows/rnadnavar.nf; then
    check_pass "MAF mode branch logic found"
else
    check_fail "MAF mode branch logic missing"
fi

echo ""
echo "9. Checking MAF_FILTERING_RNA conditional..."
if grep -q "if (params.realignment_mode == 'maf')" workflows/rnadnavar.nf; then
    check_pass "MAF_FILTERING_RNA properly wrapped for MAF mode only"
else
    check_fail "MAF_FILTERING_RNA conditional not found"
fi

echo ""
echo "10. Checking Python rescue script..."
if [ -f "bin/run_rescue_vcf.py" ]; then
    check_pass "run_rescue_vcf.py exists"
else
    check_fail "run_rescue_vcf.py missing"
fi

echo ""
echo "=========================================="
echo "Validation Summary"
echo "=========================================="
echo -e "${GREEN}Passed:${NC} $PASS"
echo -e "${RED}Failed:${NC}  $FAIL"
echo ""

if [ $FAIL -eq 0 ]; then
    echo -e "${GREEN}All validation checks passed!${NC}"
    exit 0
else
    echo -e "${RED}Validation failed with $FAIL error(s).${NC}"
    exit 1
fi
