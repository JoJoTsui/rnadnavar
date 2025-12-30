#!/bin/bash
# Verify output directory structure and metadata correctness

OUTDIR="/t9k/mnt/hdd/work/Vax/sequencing/aim_exp/rdv_test/COO8801.subset" # Adjust if your test output goes elsewhere

echo "Verifying output structure in $OUTDIR..."

# 1. Check MAF Results Directory
if [ -d "$OUTDIR/maf_results" ]; then
    echo "[PASS] maf_results directory exists."
    
    if [ -d "$OUTDIR/maf_results/vcf2maf" ]; then
        echo "[PASS] maf_results/vcf2maf exists."
    else
        echo "[FAIL] maf_results/vcf2maf missing."
    fi

    if [ -d "$OUTDIR/maf_results/consensus" ]; then
        echo "[PASS] maf_results/consensus exists."
    else
        echo "[FAIL] maf_results/consensus missing."
    fi

    if [ -d "$OUTDIR/maf_results/filtered" ]; then
        echo "[PASS] maf_results/filtered exists."
    else
        echo "[FAIL] maf_results/filtered missing."
    fi
else
    echo "[FAIL] maf_results directory missing."
fi

# 2. Check VCF Realignment Directory
if [ -d "$OUTDIR/vcf_realignment" ]; then
    echo "[PASS] vcf_realignment directory exists."
    
    # Check for GATK outputs in the correct location
    if [ -d "$OUTDIR/vcf_realignment/preprocessing/markduplicates" ]; then
        echo "[PASS] vcf_realignment/preprocessing/markduplicates exists."
    else
        echo "[FAIL] vcf_realignment/preprocessing/markduplicates missing."
    fi

    if [ -d "$OUTDIR/vcf_realignment/preprocessing/splitncigarreads" ]; then
        echo "[PASS] vcf_realignment/preprocessing/splitncigarreads exists."
    else
        echo "[FAIL] vcf_realignment/preprocessing/splitncigarreads missing."
    fi
    
    # Check for variant calling outputs
    if [ -d "$OUTDIR/vcf_realignment/variant_calling" ]; then
        echo "[PASS] vcf_realignment/variant_calling exists."
    else
        echo "[FAIL] vcf_realignment/variant_calling missing."
    fi

else
    echo "[FAIL] vcf_realignment directory missing."
fi

# 3. Check Metadata Suffix (_realign)
echo "Checking for _realign suffix in filenames..."
FOUND_REALIGN=$(find "$OUTDIR" -name "*_realign*" | head -n 1)

if [ -n "$FOUND_REALIGN" ]; then
    echo "[PASS] Found files with _realign suffix: $FOUND_REALIGN"
else
    echo "[WARN] No files with _realign suffix found. Check if realignment ran successfully."
fi

# 4. Check VCF Realignment Subdirectories
echo "\nChecking vcf_realignment subdirectories..."
for subdir in variant_calling normalized consensus filtered; do
    if [ -d "$OUTDIR/vcf_realignment/$subdir" ]; then
        echo "[PASS] vcf_realignment/$subdir exists."
        # Check for _realign suffix in files within this directory
        REALIGN_FILES=$(find "$OUTDIR/vcf_realignment/$subdir" -type f -name "*_realign*" | wc -l)
        if [ "$REALIGN_FILES" -gt 0 ]; then
            echo "  [PASS] Found $REALIGN_FILES files with _realign suffix."
        else
            echo "  [WARN] No files with _realign suffix found in $subdir."
        fi
    else
        echo "[FAIL] vcf_realignment/$subdir missing."
    fi
done

# 5. Check Second Rescue Outputs
echo "\nChecking for second rescue outputs..."
if [ -d "$OUTDIR/rescue" ]; then
    FIRST_RESCUE=$(find "$OUTDIR/rescue" -name "*rescued_RNA_TUMOR_vs*" -type f | head -n 1)
    SECOND_RESCUE=$(find "$OUTDIR/rescue" -name "*rescued_RNA_TUMOR_realign_vs*" -type f | head -n 1)
    
    if [ -n "$FIRST_RESCUE" ]; then
        echo "[PASS] Found first-round rescue: $(basename "$FIRST_RESCUE")"
    fi
    
    if [ -n "$SECOND_RESCUE" ]; then
        echo "[PASS] Found second-round rescue (with realigned RNA): $(basename "$SECOND_RESCUE")"
    else
        echo "[WARN] Second-round rescue with _realign suffix not found."
    fi
fi

echo "\nVerification complete."
