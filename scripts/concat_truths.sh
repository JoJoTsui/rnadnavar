#!/usr/bin/env bash

# IO
WD="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/benchmark/seqc2/somatic/release/latest"
TRUTH_SNV="${WD}/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
TRUTH_INDEL="${WD}/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz"
TRUTH_SEQC2="${WD}/high-confidence_sSNV+INDEL_in_HC_regions_v1.2.1.vcf.gz"

# concat
bcftools concat \
    --allow-overlaps \
    --remove-duplicates \
    "${TRUTH_SNV}" "${TRUTH_INDEL}" \
    | bcftools sort -Oz -o "${TRUTH_SEQC2}"

# Index the new merged truth set
tabix -p vcf "${TRUTH_SEQC2}"
